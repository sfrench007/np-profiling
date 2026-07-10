# =============================================================================
# 05.fuzzyRF_OvO_roc_duck.R
#
# Purpose:  Train one-vs-one (OvO) Random Forest classifiers for each
#           annotation level (drug target, ATC levels 1-4) using the
#           annotated training data prepared by script 04.  Saves a
#           per-level RF model object to results/ as an Rds file.
#
# Inputs:
#   data/output/full_db_training_scale_dba_remap.Rdata  (script 04 output)
#   legends/drugbank_curated_annotations.xlsx            (curated target remaps)
#   functions/functions.R                                (shared helper functions)
#
# Outputs:
#   results/RF_<level>.Rds   one file per annotation level (target, level_1..4)
#   data/output/RF_training_annotations.Rds   compiled class-level summary
#
# CLI args:
#   --wd <path>   set the working directory (defaults to current directory)
#
# Notes:
#   - Parallel processing uses makePSOCKcluster(); adjust the core count to
#     match your hardware (default 50 — reduce on workstations).
#   - Run from the project root so relative paths resolve correctly.
# =============================================================================

# Standard option declarations (consistent with scripts 01-04)
options(stringsAsFactors=FALSE) # Otherwise we need to force them as strings repeatedly
options(echo=FALSE)             # Rscript needs options(echo=TRUE) to make an output file
options(cli.progress_show_after=0)
options(cli.condition="always")
options(warn=(-1))              # Suppress NA-coercion warnings
arguments <- commandArgs(trailingOnly=FALSE) # Capture all arguments for later inspection

library(cli)          # Progress bars with ETA
library(stringi)
library(stringr)
library(progress)
library(reshape2)
library(ggplot2)
library(openxlsx)
library(caret)
library(gbm)
library(pROC)
library(GGally)
library(doParallel)
library(MLmetrics)
library(Metrics)
library(RecordLinkage)

# Header line
cli_h1("Training Random Forest classifiers (OvO)")

# Set working directory from --wd argument if provided
if(any(arguments=="--wd")) {
    arg_ind <- which(arguments=="--wd")
    trailing_arg <- arguments[arg_ind[1]+1]
    if(length(arg_ind)>=1) {
        invisible(tryCatch(setwd(trailing_arg),error=function(e) e, finally=function(x) setwd(trailing_arg)))
    }
    cli_alert_success(" Set working directory to {trailing_arg}")
}

# Multiclass ROC/logLoss summary function adapted from:
# https://github.com/topepo/caret/issues/107
multiClassSummaryGit <- function (data, lev = NULL, model = NULL) {
  #Load Libraries
  require(Metrics)
  require(caret)
  
  #Check data 
  if (!all(levels(data[, "pred"]) == levels(data[, "obs"]))) 
    stop("levels of observed and predicted data do not match")
  
  #Calculate custom one-vs-all stats for each class
  prob_stats <- lapply(levels(data[, "pred"]), function(class){
    
    #Grab one-vs-all data for the class
    pred <- ifelse(data[, "pred"] == class, 1, 0)
    obs  <- ifelse(data[,  "obs"] == class, 1, 0)
    prob <- data[,class]
    
    #Calculate one-vs-all AUC and logLoss and return
    cap_prob <- pmin(pmax(prob, .000001), .999999)
    prob_stats <- c(auc(obs, prob), logLoss(obs, cap_prob))
    names(prob_stats) <- c('ROC', 'logLoss')
    return(prob_stats) 
  })
  prob_stats <- do.call(rbind, prob_stats)
  rownames(prob_stats) <- paste('Class:', levels(data[, "pred"]))
  
  #Calculate confusion matrix-based statistics
  CM <- confusionMatrix(data[, "pred"], data[, "obs"])
  
  #Aggregate and average class-wise stats
  #Todo: add weights
  class_stats <- cbind(CM$byClass, prob_stats)
  class_stats <- colMeans(class_stats)
  
  #Aggregate overall stats
  overall_stats <- c(CM$overall)
  
  #Combine overall with class-wise stats and remove some stats we don't want 
  stats <- c(overall_stats, class_stats)
  stats <- stats[! names(stats) %in% c('AccuracyNull', 
                                       'Prevalence', 'Detection Prevalence')]
  
  #Clean names and return
  names(stats) <- gsub('[[:blank:]]+', '_', names(stats))
  return(stats)
  
}

# Create a classlist and training set (same nrows) for model training
createTrainTest <- function(sourceList,mainClasses,trainSplit=0.8,level="target") {
  # Create a large list of everything first
  classes_train <- c(); classes_test <- c()
  for(ii in 1:length(mainClasses)) {
    cc_rows_t <- sourceList[[mainClasses[ii]]]$index
    # working_t <- cc[cc_rows_t,]

    working_t <- full_db_training[[level]]$median[cc_rows_t,]
    train_t <- sample(1:nrow(working_t),round(nrow(working_t)*trainSplit,0))
    test_t <-  which(1:nrow(working_t) %in% train_t==FALSE)
    
    if(ii==1) {
      trainOutput <- working_t[train_t,]
      testOutput <- working_t[test_t,]
    } else {
      trainOutput <- rbind(trainOutput,working_t[train_t,])
      testOutput <- rbind(testOutput,working_t[test_t,])
    }
    classes_train <- c(classes_train,rep(names(mainClasses)[ii],nrow(working_t[train_t,])))
    classes_test <- c(classes_test,rep(names(mainClasses)[ii],nrow(working_t[test_t,])))
  }
  outputList <- list()
  outputList$training <- trainOutput
  outputList$classes_training <- classes_train
  outputList$test <- testOutput
  outputList$classes_test <- classes_test
  return(outputList)
}

# Subset lists from the main cp list
concList <- function(meta_in,data_in,conc) {
  conc_temp <- which(meta_in[,"Image_Metadata_CPD_MMOL_CONC"]==conc)
  outmatrix <- cbind(
    meta_in[conc_temp,],
    data_in[conc_temp,]
  )
  return(outmatrix)
}

# Identify the top classes based on frequency
gatherTargetInfo <- function(an_c,dataIn=cc_mesh,simval=0.8,smilesim=NA) {
  print(paste0("Compiling ",colnames(dataIn)[an_c],"..."))
  # Need to break up the ;s to see what individual classes are present (some are multiple, ; separated)
  total_annots <- sapply(dataIn[,an_c],function(x) strsplit(x,";")[[1]])
  
  # First compile all indvs
  ta <- c()
  for(ii in 1:length(total_annots)) ta <- c(ta,total_annots[[ii]])
  # Remove nonsense
  if(any(is.na(ta))) ta <- ta[-which(is.na(ta))]
  if(any(ta=="",na.rm=TRUE)) ta <- ta[-which(ta=="")]
  
  # Combine and count targets that are pretty much the same (70% similar by NAME otherwise)
  #   ta_sim <- substr(ta,1,nchar(ta)-round(nchar(ta)*(1-simval),0))
  #   unq_ta_sim <- sort(unique(ta_sim))
  
  # Force uppercase to ensure no case issues, for fuzzy search purposes
  patterns <- unique(ta)
  #   pmatch <- toupper(unq_ta_sim)
  
  # Compile as list for ease
  output_list <- list()
  pb_t <- txtProgressBar(min=1,max=length(patterns),style=3)
  for(ii in 1:length(patterns)) {
    # Get the matches from the OG strings - if any of the subitems is a match, log it
    total_annots_ind <- which(sapply(total_annots,function(x) any(x==patterns[ii])))
    if(!is.na(smilesim)) total_annots_ind <- total_annots_ind[-which(cc_mesh[total_annots_ind,15]<smilesim)]
    unq_name <- strsplit(names(total_annots_ind),";")
    names(total_annots_ind) <- NULL
    output_list[[patterns[ii]]]$name <- patterns[ii]
    output_list[[patterns[ii]]]$counts <- length(total_annots_ind)
    output_list[[patterns[ii]]]$index <- total_annots_ind
    output_list[[patterns[ii]]]$pubchem_cid <- cc_mesh[total_annots_ind,1]
    output_list[[patterns[ii]]]$compounds <- cc_mesh[total_annots_ind,2]
    output_list[[patterns[ii]]]$simscore <- cc_mesh[total_annots_ind,15]
    output_list[[patterns[ii]]]$full_list <- cc_mesh[total_annots_ind,an_c]
    output_list[[patterns[ii]]]$smiles <- cc_mesh[total_annots_ind,14]
    setTxtProgressBar(pb_t,ii)
  }
  return(output_list)
}

# Identify the top classes based on frequency

# Identify the top classes based on frequency
refineTarget <- function(dataIn=full_db_training,level="target",simval=0.8) {
    print(paste0("Compiling ",level,"..."))
    # Need to break up the ;s to see what individual classes are present (some are multiple, ; separated)
    total_annots <- dataIn[[level]]$db
    ta <<- total_annots[,2]

    # Combine and count targets that are pretty much the same (70% similar by NAME otherwise)
    if(level!="target") {
      xtemp <- list()
      for(iii in 1:length(unique(ta))) {
          xsim <- levenshteinSim(unique(ta),unique(ta)[iii])
          xtemp[[unique(ta)[iii]]] <- unique(ta)[which(xsim>=simval)]
          if(iii==1) corout <- xsim
          else corout <- cbind(corout,xsim)
      }
      colnames(corout) <- unique(ta);rownames(corout) <- unique(ta)

      # Need to make all the entries the same that have similar xtemp annotations
      ta2 <- ta
      for(iii in 1:length(xtemp)) {
          ta2[which(ta2 %in% xtemp[[iii]])] <- names(xtemp)[iii]
      }
      patterns <- unique(ta2)
    } else {
      ta2 <- c()
      for(iii in 1:length(ta)) {
        t0 <- which(remaps[,"Target_Name"]==ta[iii])
        if(length(t0)>0) ta2[iii] <- remaps[t0,"Functional_Class"][1]
        else ta2[iii] <- ta[iii]
      }
      patterns <- unique(ta2)
    }

    # Keep the conversion key as a global variable
    conversion_key <<- cbind(ta,ta2); colnames(conversion_key) <- c("original_target","matched_target")
    # Compile as list for ease
    output_list <- list()
    pb_t <- txtProgressBar(min=1,max=length(patterns),style=3)
    for(ii in 1:length(patterns)) {
        # Get the matches from the OG strings - if any of the subitems is a match, log it
        # total_annots_ind <- total_annots[which(ta2==patterns[ii]),1]
        total_annots_ind <- which(ta2==patterns[ii])
        # names(total_annots_ind) <- NULL
        output_list[[patterns[ii]]]$name <- patterns[ii]
        output_list[[patterns[ii]]]$counts <- length(total_annots_ind)
        output_list[[patterns[ii]]]$index <- total_annots_ind
        output_list[[patterns[ii]]]$index_key <- cbind(total_annots[which(ta2==patterns[ii]),1],total_annots_ind)
        colnames(output_list[[patterns[ii]]]$index_key) <- c("mesh_row","annot_ind")
    }
    return(output_list)
}


expandedSummary <- function(data,lev=NULL,model=NULL) {
    a1 <- twoClassSummary(data,lev,model)
    b1 <- defaultSummary(data,lev,model)
    c1 <- prSummary(data,lev,model)
    out <- c(a1,b1,c1)
    return(out)
}

load("data/output/full_db_training_scale_dba_remap.Rdata")
remaps <- openxlsx::read.xlsx("legends/drugbank_curated_annotations.xlsx","targets-final",colNames=TRUE)

# Make lists summarizing each class
sv <- 1
target_list <- refineTarget(dataIn=full_db_training,level="target",simval=sv) 
level_1_list <- refineTarget(dataIn=full_db_training,level="level_1",simval=sv)  
level_2_list <- refineTarget(dataIn=full_db_training,level="level_2",simval=sv) 
level_3_list <- refineTarget(dataIn=full_db_training,level="level_3",simval=sv) 
level_4_list <- refineTarget(dataIn=full_db_training,level="level_4",simval=sv) 

final_training <- list()
final_training[["level_1"]] <- level_1_list
final_training[["level_2"]] <- level_2_list
final_training[["level_3"]] <- level_3_list
final_training[["level_4"]] <- level_4_list
final_training[["target"]] <- target_list
# final_training[["mesh"]] <- mesh_list
saveRDS(final_training,"data/output/RF_training_annotations.Rds",version=2)

# Identify which members of each list have >=x counts for RF training
ncounts <- 20
target_mainClasses <- which(sapply(target_list,function(x) x$counts)>=ncounts)
level_1_mainClasses <- which(sapply(level_1_list,function(x) x$counts)>=ncounts)
level_2_mainClasses <- which(sapply(level_2_list,function(x) x$counts)>=ncounts)
level_3_mainClasses <- which(sapply(level_3_list,function(x) x$counts)>=ncounts)
level_4_mainClasses <- which(sapply(level_4_list,function(x) x$counts)>=ncounts)

# Try random forest stuff - first set a seed, same for all would be ideal so set it up top here
seed <- sample(1:100000,1)
print(paste0("Trying seed ",seed)); set.seed(seed)

# Populate the list starting with seed
source("functions/functions.R")
RF <- list()
tag <<- generateTag(hqp="Anonymous")
RF[["seed"]] <- seed

# Protect RF
RF2 <- RF

# Set up parallel processing
cl <- makePSOCKcluster(50, outfile="")
registerDoParallel(cl)
# stopCluster(cl)

# filename <- "targets" # Can be one of:  targets, pathway, level_1, level_2, level_3, level_4
for(bb in 1:5) {
    RF <- RF2
    filename <- c("target","level_1","level_2","level_3","level_4","mesh")[bb]
    print(paste0("Training on ",filename,"..."))

    # Split training and test sets, and create labels
    if(filename=="mesh") {
        if(length(mesh_mainClasses)==0) break()
        targets_sets <- createTrainTest(mesh_list,mesh_mainClasses,level="mesh")
    } else if(filename=="target") {
        if(length(target_mainClasses)==0) break()
        # targets_sets <- createTrainTest(target_list,newTargs)
        targets_sets <- createTrainTest(target_list,target_mainClasses,level="target")
    } else if(filename=="level_1") {
        if(length(level_1_mainClasses)==0) break()
        targets_sets <- createTrainTest(level_1_list,level_1_mainClasses,level="level_1")
    } else if(filename=="level_2") {
        if(length(level_2_mainClasses)==0) break()
        targets_sets <- createTrainTest(level_2_list,level_2_mainClasses,level="level_2")
    } else if(filename=="level_3") {
        if(length(level_3_mainClasses)==0) break()
        targets_sets <- createTrainTest(level_3_list,level_3_mainClasses,level="level_3")
    } else if(filename=="level_4") {
        if(length(level_4_mainClasses)==0) break()
        targets_sets <- createTrainTest(level_4_list,level_4_mainClasses,level="level_4")
    }
  
    targets_sets$seed <- seed # Inject seed
    RF[["train_test"]] <- targets_sets
    # save(targets_sets,file=paste0("D:/rf_",filename,"_train_test.Rdata"),version=2) # save for a later date

    # Define the dataset and annotations
    Class <- as.factor(make.names(targets_sets$classes_training))
    ProcessedData <- targets_sets$training
  
    # Define trainControl values using the custom summary function from GitHub (above)
    ctrl <- trainControl(method = "repeatedcv",
                        number = 10,
                        repeats = 2,
                        search = "random",
                        summaryFunction = multiClassSummaryGit,
                        classProbs = TRUE,
                        savePredictions = TRUE,
                        returnResamp = "all",
                        verboseIter = TRUE,
                        allowParallel=TRUE)

    mtryValues <- c(2,4,8,10,20,40,80,100)
#   mtryValues <- c(2,4,8,10,12,14,16,18,20,15,30,40,60,80,100)
  
  # OvO first (fastest)
  # Train random forest model, low ntree value given the limited training samples
    rf_model <- train(x = ProcessedData,y = Class,
                        method = "rf",
                        ntree = 100,
                        tuneGrid = data.frame(mtry = mtryValues),
                        importance = TRUE,
                        metric = "ROC",
                        trControl = ctrl)

    # Save the model
    RF[["ovo"]]$model <- rf_model

    # Create a confusion matrix for the model (rf_cm)
    rfp <- merge(rf_model$pred,rf_model$bestTune)
    rf_cm <- confusionMatrix(rf_model, norm = "none")
    RF[["ovo"]]$cm <- rf_cm

    # Set up the roc QC
    rf_response <- rf_model$pred$obs
    rf_predictor <- rf_model$pred[,3:(ncol(rf_model$pred)-3)] 
    rf_roc_ovo <- multiclass.roc(response=rf_response,predictor=rf_predictor,levels=rev(levels(rf_response)))
    RF[["ovo"]]$roc <- rf_roc_ovo

    # Save it to have it, OvR takes a while, this lets us use OvO right away
    RF$tag <- tag
    saveRDS(RF,file=paste0("results/RF_",filename,".Rds"),version=2)

    print("...Done!")
}

stopCluster(cl)
