# =============================================================================
# 07.refine_predictions_OvO.R
# Cell Painting Pipeline — Refined OvO RF Training
# =============================================================================
#
# PURPOSE
#   Loads the initial OvO predictions from script 06 and removes low-quality
#   or biologically uninformative annotation classes. Then retrains One-vs-One
#   Random Forest models on the refined class list, producing updated model
#   objects for downstream re-prediction.
#
# WORKFLOW
#   1. Load training data, OvO predictions, and annotation lists
#   2. Filter annotation classes by AUROC quality (good_auc)
#   3. Identify classes meeting the minimum sample count threshold
#   4. Retrain OvO RF models on the refined class set (all five levels)
#   5. Save per-level RF model objects to results/
#
# INPUTS  (all relative to --wd or the working directory)
#   data/output/full_db_training_scale_dba_remap.Rdata  — annotated training set
#   results/ovo_preds_fractions_auto_roc.Rds            — OvO predictions (script 06)
#   results/RF_training_annotations.Rds                 — class annotation lists
#   functions/functions.R                               — shared helper functions
#
# OUTPUTS
#   RF_training_annotations_cutoffRefined_roc.Rds       — refined annotation list
#   RF_<level>_cutoffRefined_roc.Rds   one file per level (target, level_1..4)
#
# COMMAND-LINE FLAGS
#   --wd <path>   Set working directory
#
# DEPENDENCIES
#   Script 01 (database), Script 04 (training annotations),
#   Script 05 (initial RF models), Script 06 (OvO predictions)
# =============================================================================

# Initial option declarations
options(stringsAsFactors=FALSE) # Otherwise we need to force them as strings repeatedly
options(echo=FALSE)             # Rscript needs options(echo=TRUE) to make an output file
options(cli.progress_show_after=0) # For progress bars
options(cli.condition="always")    # For progress bars
options(warn=(-1))              # Ignore NA warnings
arguments <- commandArgs(trailingOnly=FALSE) # Will capture all arguments, so can search these later

library(caret)         # RF training, confusionMatrix, trainControl
library(pROC)          # multiclass.roc()
library(MLmetrics)     # logLoss(), auc() — used in multiClassSummaryGit
library(Metrics)       # auc() fallback
library(doParallel)    # makePSOCKcluster, registerDoParallel
library(RecordLinkage) # levenshteinSim() — used in refineTarget
library(cli)           # CLI progress bars and alerts
library(DBI)           # Database interface
library(duckdb)        # Embedded analytical database

# let's set --wd as the working directory, defaulting at the directory that it's run in
if(any(arguments=="--wd")) {
    arg_ind <- which(arguments=="--wd")
    trailing_arg <- arguments[arg_ind[1]+1]
    if(length(arg_ind)>=1) {
        invisible(tryCatch(setwd(trailing_arg),error=function(e) e, finally=function(x) setwd(trailing_arg)))
    }
    cli_alert_success(" Set working directory to {trailing_arg}")
}

# Custom multiclass summary function for caret trainControl
# Source: https://github.com/topepo/caret/issues/107
multiClassSummaryGit <- function(data, lev=NULL, model=NULL) {
    require(Metrics)
    require(caret)

    if (!all(levels(data[, "pred"]) == levels(data[, "obs"])))
        stop("levels of observed and predicted data do not match")

    # Calculate custom one-vs-all stats for each class
    prob_stats <- lapply(levels(data[, "pred"]), function(class) {
        pred <- ifelse(data[, "pred"] == class, 1, 0)
        obs  <- ifelse(data[,  "obs"] == class, 1, 0)
        prob <- data[, class]
        cap_prob   <- pmin(pmax(prob, .000001), .999999)
        prob_stats <- c(auc(obs, prob), logLoss(obs, cap_prob))
        names(prob_stats) <- c("ROC", "logLoss")
        return(prob_stats)
    })
    prob_stats <- do.call(rbind, prob_stats)
    rownames(prob_stats) <- paste("Class:", levels(data[, "pred"]))

    # Calculate confusion matrix-based statistics
    CM          <- confusionMatrix(data[, "pred"], data[, "obs"])
    class_stats <- colMeans(cbind(CM$byClass, prob_stats))
    overall_stats <- c(CM$overall)

    stats <- c(overall_stats, class_stats)
    stats <- stats[!names(stats) %in% c("AccuracyNull", "Prevalence", "Detection Prevalence")]
    names(stats) <- gsub("[[:blank:]]+", "_", names(stats))
    return(stats)
}

# Create a classlist and training set (same nrows) for model training
createTrainTest <- function(sourceList, mainClasses, trainSplit=0.8, level="target") {
    classes_train <- c(); classes_test <- c()
    for(ii in 1:length(mainClasses)) {
        cc_rows_t <- sourceList[[mainClasses[ii]]]$index
        working_t <- full_db_training[[level]]$median[cc_rows_t,]
        train_t   <- sample(1:nrow(working_t), round(nrow(working_t)*trainSplit, 0))
        test_t    <- which(1:nrow(working_t) %in% train_t == FALSE)
        if(ii==1) {
            trainOutput <- working_t[train_t,]
            testOutput  <- working_t[test_t,]
        } else {
            trainOutput <- rbind(trainOutput, working_t[train_t,])
            testOutput  <- rbind(testOutput,  working_t[test_t,])
        }
        classes_train <- c(classes_train, rep(names(mainClasses)[ii], nrow(working_t[train_t,])))
        classes_test  <- c(classes_test,  rep(names(mainClasses)[ii], nrow(working_t[test_t,])))
    }
    outputList <- list(
        training         = trainOutput,
        classes_training = classes_train,
        test             = testOutput,
        classes_test     = classes_test
    )
    return(outputList)
}

# Identify unique annotation classes and group similar targets by Levenshtein similarity
refineTarget <- function(dataIn=full_db_training, level="target", simval=0.8) {
    cli_alert_info(" Compiling {level}...")
    total_annots <- dataIn[[level]]$db
    ta <- total_annots[,2]

    # Group entries with similar names (simval threshold)
    xtemp <- list()
    for(iii in 1:length(unique(ta))) {
        xsim <- levenshteinSim(unique(ta), unique(ta)[iii])
        xtemp[[unique(ta)[iii]]] <- unique(ta)[which(xsim >= simval)]
        if(iii==1) corout <- xsim
        else corout <- cbind(corout, xsim)
    }
    colnames(corout) <- unique(ta); rownames(corout) <- unique(ta)

    ta2 <- ta
    for(iii in 1:length(xtemp)) {
        ta2[which(ta2 %in% xtemp[[iii]])] <- names(xtemp)[iii]
    }
    patterns <- unique(ta2)

    output_list <- list()
    for(ii in 1:length(patterns)) {
        total_annots_ind <- which(ta2 == patterns[ii])
        output_list[[patterns[ii]]]$name      <- patterns[ii]
        output_list[[patterns[ii]]]$counts    <- length(total_annots_ind)
        output_list[[patterns[ii]]]$index     <- total_annots_ind
        output_list[[patterns[ii]]]$index_key <- cbind(
            total_annots[which(ta2 == patterns[ii]), 1],
            total_annots_ind
        )
        colnames(output_list[[patterns[ii]]]$index_key) <- c("mesh_row", "annot_ind")
    }
    return(output_list)
}

# Tuning parameters — adjust to match available hardware and desired stringency
n_cores  <- 50  # Number of parallel workers for caret
n_counts <- 25  # Minimum class sample count for inclusion in RF training

# ── Data loading ──────────────────────────────────────────────────────────────
source("functions/functions.R")
load("data/output/full_db_training_scale_dba_remap.Rdata")
preds_ovo  <- readRDS("results/ovo_preds_fractions_auto_roc.Rds")
annots_all <- readRDS("results/RF_training_annotations.Rds")

# ── Filter annotation classes ─────────────────────────────────────────────────
# Remove classes with known biological ambiguity or poor specificity
good_auc <- names(preds_ovo$target)[7:length(names(preds_ovo$target))]
bad_pattern <- "nterleukin|ucin|yc.proto|lactamase|eptidoglycan|uter.membrane|rothrombin"
if(any(grepl(bad_pattern, good_auc))) {
    good_auc <- good_auc[-which(grepl(bad_pattern, tolower(good_auc)))]
}

# Apply the filter to the annotations list
annots <- make.names(names(annots_all$target))
annots_all$target <- annots_all$target[which(annots %in% good_auc)]

target_list  <- annots_all$target
level_1_list <- annots_all$level_1
level_2_list <- annots_all$level_2
level_3_list <- annots_all$level_3
level_4_list <- annots_all$level_4

saveRDS(annots_all, "results/RF_training_annotations_cutoffRefined_roc.Rds")

# ── Identify classes meeting the minimum sample count ─────────────────────────
target_mainClasses  <- which(sapply(target_list,  function(x) x$counts) >= n_counts)
level_1_mainClasses <- which(sapply(level_1_list, function(x) x$counts) >= n_counts)
level_2_mainClasses <- which(sapply(level_2_list, function(x) x$counts) >= n_counts)
level_3_mainClasses <- which(sapply(level_3_list, function(x) x$counts) >= n_counts)
level_4_mainClasses <- which(sapply(level_4_list, function(x) x$counts) >= n_counts)

# ── RF training ───────────────────────────────────────────────────────────────
seed <- sample(1:100000, 1)
cli_alert_info(" Using seed {seed}"); set.seed(seed)

source("functions/functions.R")
RF      <- list()
RF$tag  <- generateTag()
RF$seed <- seed

cl <- makePSOCKcluster(n_cores, outfile="")
registerDoParallel(cl)

for(bb in 1:5) {
    filename <- c("target", "level_1", "level_2", "level_3", "level_4")[bb]
    cli_alert_info(" Training on {filename}...")

    # Split training and test sets, and create labels
    if(filename=="target") {
        if(length(target_mainClasses)==0) next
        targets_sets <- createTrainTest(target_list,  target_mainClasses,  level="target")
    } else if(filename=="level_1") {
        if(length(level_1_mainClasses)==0) next
        targets_sets <- createTrainTest(level_1_list, level_1_mainClasses, level="level_1")
    } else if(filename=="level_2") {
        if(length(level_2_mainClasses)==0) next
        targets_sets <- createTrainTest(level_2_list, level_2_mainClasses, level="level_2")
    } else if(filename=="level_3") {
        if(length(level_3_mainClasses)==0) next
        targets_sets <- createTrainTest(level_3_list, level_3_mainClasses, level="level_3")
    } else if(filename=="level_4") {
        if(length(level_4_mainClasses)==0) next
        targets_sets <- createTrainTest(level_4_list, level_4_mainClasses, level="level_4")
    }

    targets_sets$seed  <- seed
    RF[["train_test"]] <- targets_sets

    Class         <- as.factor(make.names(targets_sets$classes_training))
    ProcessedData <- targets_sets$training

    ctrl <- trainControl(
        method          = "repeatedcv",
        number          = 10,
        repeats         = 10,
        search          = "random",
        summaryFunction = multiClassSummaryGit,
        classProbs      = TRUE,
        savePredictions = TRUE,
        returnResamp    = "all",
        verboseIter     = TRUE,
        allowParallel   = TRUE
    )

    mtryValues <- c(2, 4, 8, 10, 20, 40, 80, 100)

    # Train OvO Random Forest model
    rf_model <- train(
        x         = ProcessedData,
        y         = Class,
        method    = "rf",
        ntree     = 100,
        tuneGrid  = data.frame(mtry = mtryValues),
        importance = TRUE,
        metric    = "ROC",
        trControl = ctrl
    )

    RF[["ovo"]]$model <- rf_model

    # Confusion matrix
    rfp   <- merge(rf_model$pred, rf_model$bestTune)
    rf_cm <- confusionMatrix(rf_model, norm="none")
    RF[["ovo"]]$cm <- rf_cm

    # ROC
    rf_response  <- rf_model$pred$obs
    rf_predictor <- rf_model$pred[, 3:(ncol(rf_model$pred)-3)]
    RF[["ovo"]]$roc <- multiclass.roc(
        response  = rf_response,
        predictor = rf_predictor,
        levels    = rev(levels(rf_response))
    )

    saveRDS(RF, file=paste0("results/RF_", filename, "_cutoffRefined_roc.Rds"), version=2)
    cli_alert_success(" {filename} done")
}

stopCluster(cl)
