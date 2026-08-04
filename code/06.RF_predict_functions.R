# =============================================================================
# 06.RF_predict_functions.R
# Cell Painting Pipeline — RF Prediction & Analysis
# =============================================================================
#
# PURPOSE
#   Runs trained Random Forest models (OvO and OvR strategies) on the
#   "fractions" (test/new compound) subset of the Cell Painting DuckDB
#   database, then filters and ranks predictions by AUROC and probability
#   score to produce final hit lists.
#
# WORKFLOW
#   1. Load training and fractions data from cellpainting.duckdb
#   2. predict_ovo()  — apply One-vs-One RF models; save raw predictions
#   3. predict_ovr()  — apply One-vs-Rest RF models; save raw predictions
#   4. analyze_ovo()  — filter OvO hits by AUROC/accuracy; rank and export
#   5. analyze_ovr()  — filter OvR hits by AUROC/Kappa; rank and export
#
# INPUTS  (all relative to --wd or the working directory)
#   data/db/cellpainting.duckdb    — source database (read-only)
#   results/RF_<level>.Rds         — trained RF model objects (from script 05)
#
# OUTPUTS
#   results/RF_Preds_fractions_complete_ovo_auto_roc.Rds
#   results/ovo_preds_fractions_auto_roc.Rds
#
# COMMAND-LINE FLAGS
#   --wd <path>   Set working directory
#   --force       Re-run even if the processing log shows script 06 has run
#
# DEPENDENCIES
#   Script 01 (database), Script 04 (training annotations), Script 05 (RF models)
# =============================================================================

# Initial option declarations
options(stringsAsFactors=FALSE) # Otherwise we need to force them as strings repeatedly
options(echo=FALSE)             # Rscript needs options(echo=TRUE) to make an output file
options(cli.progress_show_after=0) # For progress bars
options(cli.condition="always")    # For progress bars
options(warn=(-1))              # Ignore NA warnings
arguments <- commandArgs(trailingOnly=FALSE) # Will capture all arguments, so can search these later

# Initial library loading
library(caret)      # RF training, confusionMatrix, thresholder
library(pROC)       # roc(), multiclass.roc()
library(MLmetrics)  # logLoss(), auc() — used in multiClassSummaryGit
library(Metrics)    # auc() fallback
library(duckdb)     # Embedded analytical database
library(cli)        # CLI progress bars and alerts
library(DBI)        # Database interface

# Get the peak of the density distribution
denPeak <- function(dataIn,adjval=0.7) {
    denT <- density(dataIn,adjust=adjval,na.rm=TRUE)
    denP <- denT$x[which(denT$y==max(denT$y,na.rm=TRUE))]
    return(denP)
}

# First try the OvO models because it's very fast
predict_ovo <- function(ovo_dir="results",outfile="RF_Preds_fractions_complete_ovo.Rds") {
    # Make a master output file for predictions
    Preds <- list()
    Preds[["fractions"]]$data <- fractions$features_1pass_3
    Preds[["fractions"]]$metadata <- fractions$metadata
    
    ovo_nums <- 0
    if(any(grepl("level_|target",list.files(ovo_dir,pattern=".Rds")))) {
        ovo_nums <- ovo_nums + length(which(grepl("level_|target",list.files(ovo_dir,pattern=".Rds"))))
    }

    for(bb in 1:ovo_nums) {
        filename <- c("target","level_1","level_2","level_3","level_4","mesh")[bb]
        print(paste0("Predicting on ",filename))

        # filename <- "level_4" # Can be one of:  targets, pathway, level_1, level_2, level_3, level_4
        RF <- readRDS(paste0(ovo_dir,"/RF_",filename,".Rds"))

        # Predict on previous test set 20% pulled out
        rf_model <- RF$ovo$model
        cc_t_cmpd <- training$metadata[match(rownames(RF$train_test$test),rownames(training$features_1pass_3)),8]
        cc_t_mahal <- training$mahalanobis[match(rownames(RF$train_test$test),rownames(training$features_1pass_3))]
        cc_t_names <- RF$train_test$classes_test
        pred_ovo_t_names <- predict(rf_model,RF$train_test$test)
        pred_ovo_t <- predict(rf_model,RF$train_test$test,type="prob")
        top_prob_t <- apply(pred_ovo_t,1,max)
        output_ovo_t <- cbind(
            cc_t_cmpd,
            cc_t_mahal,
            RF$train_test$classes_test,
            pred_ovo_t_names,
            top_prob_t,
            pred_ovo_t
        )
        Preds[[filename]]$ovo$test <- output_ovo_t

        # Predict fractions on current model/level
        ft0 <<- colnames(RF$train_test$test)
        pred_ovo_f_names <- predict(rf_model,fractions$features_1pass_3)
        pred_ovo_f_mahal <- fractions$mahalanobis
        pred_ovo_f <- predict(rf_model,fractions$features_1pass_3,type="prob")
        top_prob_f <- apply(pred_ovo_f,1,max)
        output_ovo_f <- cbind(
            pred_ovo_f_names,
            pred_ovo_f_mahal,
            top_prob_f,
            pred_ovo_f
        )
        Preds[[filename]]$ovo$fractions <- output_ovo_f

        # OvO holdout set
        holdout <- Preds[[filename]]$ovo$test
        rownames(holdout) <- NULL
        colnames(holdout)[1:5] <- c("compound","pheno.activity","actual.class","pred.class","pred.score")
        holdout <- holdout[order(holdout[,4],decreasing=TRUE),]

        # Create a big confusion matrix for the ovo classes
        y_actual <- as.factor(make.names(holdout[,3]))
        y_predict <- as.factor(holdout[,4])
        ovo_cm <- confusionMatrix(y_actual,y_predict)

        # OvO predictions
        # rnames1 <- rownames(Preds[[filename]]$ovo$fractions)
        rnames <- apply(fractions$metadata,1,function(x) paste(x[3],x[2],paste0("R",x[11]),sep="-"))
        predictions <- cbind(rnames,Preds[[filename]]$ovo$fractions)
        rownames(predictions) <- NULL
        colnames(predictions)[1:4] <- c("compound","pheno.activity","pred.class","pred.score")
        # predictions <- predictions[order(predictions[,4],decreasing=TRUE),]

        Preds[[filename]]$ovo_holdout <- holdout
        Preds[[filename]]$ovo_predictions <- predictions
        Preds[[filename]]$ovo_cm <- ovo_cm
    }

    # Save the predictions
    saveRDS(Preds,file=outfile,version=2) 
}
 
# Then move onto the OvR models if we have them calculated
predict_ovr <- function(ovr_dir="results",outfile="RF_Preds_fractions_complete_ovr.Rds") {
    Preds <- list()
    Preds[["fractions"]]$data <- fractions$features_1pass_3
    Preds[["fractions"]]$metadata <- fractions$metadata

    # Fix the NAs
    norows <- which(apply(fractions$features_1pass_3,1,function(x) length(which(is.na(x))))>100)
    if(length(norows)>0) {
        fractions$features_1pass_3 <- fractions$features_1pass_3[-norows,]
        fractions$metadata <- fractions$metadata[-norows,]
    }
    nocols <- unique(which(is.na(fractions$features_1pass_3),arr.ind=TRUE)[,2])
    if(length(nocols)>0) {
        tvals <- apply(fractions$features_1pass_3[,nocols],2,function(x) median(x,na.rm=TRUE))
        for(kk in 1:length(nocols)) {
            fractions$features_1pass_3[which(is.na(fractions$features_1pass_3[,nocols[kk]])),nocols[kk]] <- tvals[kk]
        }
    }
    
    ovr_nums <- 0
    if(any(grepl("level_|target",list.files(ovr_dir,pattern=".Rds")))) {
        ovr_nums <- ovr_nums + length(which(grepl("level_|target",list.files(ovr_dir,pattern=".Rds"))))
    }

    for(bb in 1:ovr_nums) {
        filename <- c("target","level_1","level_2","level_3","level_4","mesh")[bb]
        print(paste0("Predicting on ",filename))

        # filename <- "level_4" # Can be one of:  targets, pathway, level_1, level_2, level_3, level_4
        RF <- readRDS(paste0(ovr_dir,"/RF_",filename,".Rds"))

        # Next move onto the OvR - predict on each level/filename across each model
        ovr_holdout <- list(); ovr_cm <- list(); hname_all <- c()
        cli_progress_bar(paste0("OvR predictions [" , filename, "]"), total=length(RF$ovr), clear=FALSE)
        for(k in 1:length(RF$ovr)) {
            cli_progress_update()
            # Reinitialize values just to be safe
            rf_model <- RF$ovr[[k]]$model
            xNames <- rf_model$finalModel$xNames
            cc_t <- RF$train_test$test
            cc_t_mahal <- training$mahalanobis[match(rownames(RF$train_test$test),rownames(training$features_1pass_3))]
            cc_t_names <- RF$train_test$classes_test
            wname <- names(RF$ovr)[k]

            # For test set 
            pred_ovr_t_names <- predict(rf_model,cc_t)
            pred_ovr_t <- predict(rf_model,cc_t,type="prob")
            fixednames <- make.names(tolower(RF$train_test$classes_test))
            fixednames[which(fixednames %in% wname==FALSE)] <- "other"
            output_ovr_t <- cbind(
                fixednames,
                cc_t_mahal,
                pred_ovr_t_names,
                pred_ovr_t
            )
            Preds[[filename]]$ovr[[wname]]$test <- output_ovr_t

            # And fractions
            pred_ovr_f_names <- predict(rf_model,fractions$features_1pass_3)
            pred_ovr_f_mahal <- fractions$mahalanobis
            pred_ovr_f <- predict(rf_model,fractions$features_1pass_3,type="prob")
            output_ovr_f <- cbind(
                pred_ovr_f_names,
                pred_ovr_f_mahal,                
                pred_ovr_f
            )
            Preds[[filename]]$ovr[[wname]]$fractions <- output_ovr_f

            # OvR holdout
            hname_all[k] <- names(RF$ovr)[k]
            ovr_holdout[[k]] <- Preds[[filename]]$ovr[[k]][[1]]
            colnames(ovr_holdout[[k]]) <- c("compound","pheno.activity","prediction","active","inactive")

            # Create a confusion matrix for each ovr class
            y_actual <- sapply(ovr_holdout[[k]][,1],function(xx) ifelse(xx=="other","Inactive","Active"))
            names(y_actual) <- NULL
            y_predict <- ovr_holdout[[k]][,3]
            ovr_cm[[k]] <- confusionMatrix(as.factor(y_actual), as.factor(y_predict))
        }
        cli_progress_done()
        names(ovr_cm) <- hname_all

        # OvR predictions
        rnames <- apply(fractions$metadata,1,function(x) paste(x[3],x[2],paste0("R",x[11]),sep="-"))
        for(i in 1:length(Preds[[filename]]$ovr)) {
            rnames <- apply(fractions$metadata,1,function(x) paste(x[3],x[2],paste0("R",x[11]),sep="-"))
            if(i==1) {
                ovr_preds <- cbind(rnames,Preds[[filename]]$ovr[[i]][[2]][,2:3])
            } else {
                ovr_preds <- cbind(ovr_preds,Preds[[filename]]$ovr[[i]][[2]][,3])
            }
        }
        colnames(ovr_preds) <- c("Compound","Activity",names(Preds[[filename]]$ovr))

        top_p <- apply(ovr_preds[,3:ncol(ovr_preds)],1,function(x) max(as.numeric(x)))
        top_p2 <- sapply(1:length(top_p),function(x) {
            qq <- colnames(ovr_preds)[which(ovr_preds[x,]==top_p[x])]
            ifelse(length(qq)>1,paste(qq,collapse=" | "),qq)
        })

        ovr_preds <- cbind(ovr_preds[,1:2],top_p,top_p2,ovr_preds[,3:ncol(ovr_preds)])

        # ovr_holdout <- ovr_holdout[order(rownames(Preds[[filename]]$ovr[[i]][[2]])),]
        # ovr_preds <- ovr_preds[order(rownames(Preds[[filename]]$ovr[[i]][[2]])),]
        ovr_preds <- ovr_preds[order(ovr_preds[,3],decreasing=TRUE),]
        colnames(ovr_preds)[1:4] <- c("compound","pheno.activity","pred.score","top.pred")

        Preds[[filename]]$ovr_holdout <- ovr_holdout
        Preds[[filename]]$ovr_predictions <- ovr_preds
        Preds[[filename]]$ovr_cm <- ovr_cm
    }

    # Save the predictions
    saveRDS(Preds,file=outfile,version=2) 
}

# OvO results analysis
analyze_ovo <- function(fileIn="RF_Preds_fractions_complete_ovo.Rds",fc_mult=1,ovo_dir="results",outfile="ovo_preds_fractions.Rds",useAcc=TRUE,QC=FALSE,iqr_cutoff=0) {
    Preds <- readRDS(paste0(ovo_dir,"/",fileIn))
    # OvO tests
    ovo_final <- list()
    for(q in 2:length(Preds)) {
        # > names(Preds)
        # [1] "fractions" "target"    "level_4"   "level_3"   "level_2"   "level_1"
        ovo_out <- Preds[[q]]$ovo$fractions
        ovo_out[,4:ncol(ovo_out)] <- sweep(ovo_out[,4:ncol(ovo_out)],2,apply(ovo_out[,4:ncol(ovo_out)],2,median))
        ovo_out <- cbind(apply(ovo_out[,4:ncol(ovo_out)],1,max),ovo_out)

        # Fix rownames 
        rnames <- Preds[[q]][[3]][,1]
        ovo_final[[names(Preds)[q]]] <- cbind(rnames,ovo_out)
    }

    # Check for fold changes vs the norm for each column of interest
    # fc_cutoff <- 0.5
    mahal <- Preds[[2]]$ovo$fraction[,2]
    ovo_preds <- list()
    for(q in 1:length(ovo_final)) {
        roc_cutoff <- c(0.65,0.65,0.65,0.60,0.55)[q]
        # Load the model
        RF <- readRDS(paste0(ovo_dir,"/RF_",names(Preds)[q+1],".Rds"))

        # Get AUROCs from OvO
        allrocs <- sapply(names(RF$ovo[[3]][[5]]),function(x) strsplit(x,"/")[[1]][1])
        unq_t <- unique(allrocs)
        all_rocs <- c()
        cli_progress_bar(paste0("Computing AUROCs [" , names(Preds)[q+1], "]"), total=length(unq_t), clear=FALSE)
        for(v in 1:length(unq_t)) {
            roc_ind <- which(allrocs==unq_t[v])
            roc_means <- sapply(RF$ovo[[3]][[5]][roc_ind],function(x) roc(x[[1]][[13]],x[[1]][[12]],quiet=TRUE)[[9]][[1]])
            all_rocs[v] <-  round(mean(roc_means,na.rm=TRUE),3)
            cli_progress_update()
        }
        cli_progress_done()
        names(all_rocs) <- unq_t

        # Get the acceptable accuracy classes (adjusted accuracy from CM)
        working_acc <- Preds[[q+1]]$ovo_cm[[4]][,11]
        if(length(working_acc)==(length(all_rocs)+1)) working_acc <- working_acc[-length(working_acc)]
        if(any(is.na(working_acc))) working_acc[which(is.na(working_acc))] <- (Preds[[q+1]]$ovo_cm[[4]][which(is.na(working_acc)),2]-0.9)*10

        # Get the good ROC targets and accuracy targets
        if(!useAcc) {
            goodrocs <- which(all_rocs>roc_cutoff | grepl("Actin",names(all_rocs)))
            names(goodrocs) <- names(all_rocs)[goodrocs]
            if(any(grepl(badClasses,names(goodrocs)))) {
                goodrocs <- goodrocs[-which(grepl(badClasses,tolower(names(goodrocs))))]
            }
        } else {
            goodrocs <- sort(unique(c(which(all_rocs>roc_cutoff & working_acc>roc_cutoff),which(grepl("Actin",names(all_rocs))))))
            names(goodrocs) <- names(all_rocs)[goodrocs]
            if(any(grepl(badClasses,names(goodrocs)))) {
                goodrocs <- goodrocs[-which(grepl(badClasses,tolower(names(goodrocs))))]
            }
        }
        
        ovo_final_t <- ovo_final[[q]][,c(1:3,which(colnames(ovo_final[[q]]) %in% names(goodrocs)))]
        rownames(ovo_final_t) <- gsub("[.]","-",rownames(ovo_out))
        colnames(ovo_final_t)[1:3] <- c("Code","Prediction","Top.FC")
        
        # Must scrub out the carpenter controls *before* the quantile calculations otherwise they'll skew
        # Split the names column by "-" and then search for something oerlapping with a well string
        wells384 <- paste0(rep(LETTERS[1:16],24),sprintf("%02d",sort(rep(1:24,16))))
        allwells <- sapply(ovo_final_t[,1],function(x) {
            split_temp <- strsplit(x,"-")[[1]]
            if(any(split_temp %in% wells384)) {
                split_temp[which(split_temp %in% wells384)][1]
            } else {
                NA
            }
        })
        rT_ovo <- which(as.numeric(substr(allwells,2,3)) %in% c(1,2,3,4,21,22,23,24))
        ovo_final_t <- ovo_final_t[-rT_ovo,]

        # Go through the list with the good ROC values and update with fold-change
        # QC <- FALSE
        actives <- c(); preds <- c()
        nvals <- 0

        if(!is.na(names(ovo_final)[q])) {
            if(names(ovo_final)[q]=="target") {
                targ_classes <<- c(); targ_iqr <<- c(); targ_max.iqr <<- c()
            }
        }
        
        if(QC!=FALSE) {
            dev.new()
            # par(mfrow=c(10,12))
            par(mfrow=c(8,10))
            par(mar=c(2,1,2,1))
        }
        for(i in 4:ncol(ovo_final_t)) {
            ovo_final_t[,i] <- ovo_final_t[,i]-min(ovo_final_t[,i])
            ovo_temp <- ovo_final_t[,i]

            # Let's try using Cohen's Kappa coefficient to find the proper thresholds to use here from the CM
            # Update:  For OvO this doesn't work as well as it does with OvR, the multiclass gives just a single Kappa
            # Solution: Stick with quantile for now, but need to figure this out

            # The quantile cutoff is a fine start, but some classes have a lot of hits that are missed this way          
            iqr <- quantile(ovo_temp,probs=c(0,0.1,0.5,0.9,1))[4]-quantile(ovo_temp,probs=c(0,0.1,0.5,0.9,1))[2]
            if(!is.na(names(ovo_final)[q])) {
                if(names(ovo_final)[q]=="target") {
                    targ_classes <<- c(targ_classes,colnames(ovo_final_t)[i])
                    targ_iqr <<- c(targ_iqr,iqr)
                    targ_max.iqr <<- c(targ_max.iqr,(max(ovo_temp)-iqr))                
                }
            }
        
            if(iqr>0.10) {
                fc_cutoff <- quantile(ovo_temp,probs=0.96)
            } else if(iqr>=0.05) {
                fc_cutoff <- quantile(ovo_temp,probs=0.97)
            } else if(iqr>=0.02) {
                fc_cutoff <- quantile(ovo_temp,probs=0.98)
            } else if(iqr>0) {
                fc_cutoff <- quantile(ovo_temp,probs=0.99)
            } else {
                fc_cutoff <- quantile(ovo_temp,probs=0.995)
            }
            fc_cutoff <- fc_cutoff*fc_mult
            print(paste0(colnames(ovo_final_t)[i]," | iqr=",iqr," | max-iqr=",max(ovo_temp)-iqr," | fc_cutoff=",fc_cutoff))
            if(iqr>=iqr_cutoff) {
                act_temp <- which(ovo_temp>fc_cutoff)
            } else {
                act_temp <- integer(0)
            }
            
            pre_temp <- rep(colnames(ovo_final_t)[i],length(act_temp))
            actives <- c(actives,act_temp)
            preds <- c(preds,pre_temp)
            nval <- length(act_temp)
            nvals <- nvals+nval

            # QC portion if needed
            if(QC!=FALSE) {
                # print(paste0(substr(colnames(ovo_final_t)[i],1,15)," | iqr=",iqr," | shapiro=",round(shapiro.test(sample(ovo_temp,5000,replace=TRUE))[[1]],3)," | ",round(shapiro.test(sample(ovo_temp,5000,replace=TRUE))[[1]]*iqr,3)," | n=",nval))
                plot(ovo_temp[sample(1:length(ovo_temp),length(ovo_temp))],pch=19,col=rgb(0,0,0,0.2),ylab="Probability")#,ylim=c(0,30))
                abline(h=fc_cutoff,col="blue",lwd=3)
                fdata <- density(ovo_temp[sample(1:length(ovo_temp),length(ovo_temp))],adjust=2)
                plot(y=fdata$x,x=fdata$y,type="l")
                polygon(y=fdata$x,x=fdata$y,fill="grey")
                text(x=1000,y=.020,labels=paste0("n=",nval),col="red")
                # hist(ovo_temp[sample(1:length(ovo_temp),length(ovo_temp))],breaks=100,col="grey",border=NA)
                # text(x=0,y=100,labels=i,col="red")
                # abline(v=fc_cutoff,col=rgb(0,0,1,0.5),lwd=3)                
            }
        }

        # Get the max FC score
        ovo_final_t[,3] <- apply(ovo_final_t[,4:ncol(ovo_final_t)],1,function(x) max(as.numeric(x),na.rm=TRUE))

        # Subset actives are above the new cutoff
        ovo_final2 <- ovo_final_t[actives,]
        ovo_final2[,2] <- preds
        mahal <- mahal[-rT_ovo]; mahal <- mahal[actives]

        # Add in the ovo ROC scores
        rocval <- c()
        for(s in 1:nrow(ovo_final2)) {
            nval <- as.numeric(all_rocs[which(toupper(names(all_rocs))==toupper(ovo_final2[s,2]))])
            if(length(nval)>0) rocval[s] <- nval
            else rocval[s] <- NA
        }

        # Recompile everything together
        ovo_final3 <- cbind(ovo_final2[,1:3],ovo_final_t[actives,2],rocval,mahal,ovo_final2[,4:ncol(ovo_final2)])
        colnames(ovo_final3)[4:6] <- c("OG.Pred","AUROC","Pheno.Act")
        rownames(ovo_final3) <- NULL
        
        # Order the various hits by name, and then by probability score
        finalOrder_t <- ovo_final3[order(as.numeric(ovo_final3[,3]),decreasing=TRUE),]
        finalOrder <- finalOrder_t[order(finalOrder_t[,1]),]

        # And a final filter if desired
        # finalOrder <- finalOrder[-which(finalOrder[,4]<0.04 | finalOrder[,3]<0.07),]

        # Save the predictions    
        ovo_preds[[names(Preds)[q+1]]] <- finalOrder
        # ovo_preds[[names(Preds)[q+1]]] <- ovo_final3[finalOrder,]
        # write.table(ovo_final3,"ovo_final_level2.tsv",sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
    }

    # View(ovo_preds[[1]])
    saveRDS(ovo_preds,file=outfile,version=2)
}

# OvR results analysis
analyze_ovr <- function(fileIn="RF_Preds_fractions_complete_ovr.Rds",fc_mult=1,ovr_dir="results",outfile="results/ovr_preds_fractions.Rds",useAcc=FALSE,useKappas=TRUE,QC=FALSE) {
    Preds <- readRDS(paste0(ovr_dir,"/",fileIn))

    # OvR tests
    mahal <- Preds[[2]]$ovr[[1]]$fraction[,2]
    ovr_final <- list()
    for(q in 2:length(Preds)) {
        # > names(Preds)
        # [1] "fractions" "target"    "level_4"   "level_3"   "level_2"   "level_1"
        ovr_out <- sapply(Preds[[q]]$ovr,function(x) x$fractions[,3])
        rownames(ovr_out) <- rownames(Preds[[q]]$ovr[[1]]$fractions)
        # ovr_out[,3:ncol(ovr_out)] <- sweep(ovr_out[,3:ncol(ovr_out)],2,apply(ovr_out[,3:ncol(ovr_out)],2,median))
        pred <- apply(ovr_out,1,function(x) colnames(ovr_out)[which(x==max(x))[1]])
        ovr_out <- cbind(apply(ovr_out,1,max),ovr_out)

        # Fix rownames 
        rnames <- Preds[[q]][[3]][,1]
        ovr_final[[names(Preds)[q]]] <- matrix(data=cbind(rnames,pred,ovr_out),ncol=ncol(cbind(rnames,pred,ovr_out)),dimnames=list(c(),c("code","prediction","max.pred",colnames(ovr_out)[2:ncol(ovr_out)])))
    }

    # Check for fold changes vs the norm for each column of interest
    # fc_cutoff <- 0.5
    ovr_preds <- list()
    for(q in 1:length(ovr_final)) {
        roc_cutoff <- c(0.65,0.65,0.65,0.60,0.55)[q]
        # Load the model
        RF <- readRDS(paste0(ovr_dir,"/RF_",names(Preds)[q+1],".Rds"))

        # Get AUROCs from OvR
        all_rocs <-  round(as.numeric(sapply(RF$ovr,function(x) x$roc$auc)),3)
        names(all_rocs) <- names(RF$ovr)

        # Get accuracies
        working_acc <- sapply(Preds[[q+1]]$ovr_cm,function(xx) xx[[4]][11])
        if(length(working_acc)==(length(all_rocs)+1)) working_acc <- working_acc[-length(working_acc)]
        if(any(is.na(working_acc))) working_acc[which(is.na(working_acc))] <- ((sapply(Preds[[q+1]]$ovr_cm,function(xx) xx[[4]][2])-0.9)*10)[which(is.na(working_acc))]

        # Get the good ROC and accuracy targets
        if(!useAcc) {
            goodrocs <- which(all_rocs>roc_cutoff)
        } else {
            goodrocs <- which(all_rocs>roc_cutoff & working_acc>roc_cutoff)
        }
        # Get default kappas just for QC purposes
        kappas <- sapply(Preds[[q+1]]$ovr_cm,function(xx) xx[[3]][2])[goodrocs]

        ovr_final_t <- ovr_final[[q]][,c(1:3,which(colnames(ovr_final[[q]]) %in% names(goodrocs)))]

        rownames(ovr_final_t) <- gsub("[.]","-",rownames(ovr_out))
        colnames(ovr_final_t)[1:3] <- c("code","prediction","top.score")

        # Must scrub out the carpenter controls *before* the quantile calculations otherwise they'll skew
        wells384 <- paste0(rep(LETTERS[1:16],24),sprintf("%02d",sort(rep(1:24,16))))
        allwells <- sapply(ovr_final_t[,1],function(x) {
            split_temp <- strsplit(x,"-")[[1]]
            if(any(split_temp %in% wells384)) {
                split_temp[which(split_temp %in% wells384)][1]
            } else {
                NA
            }
        })

        rT_ovr <- which(as.numeric(substr(allwells,2,3)) %in% c(1,2,3,4,21,22,23,24))
        ovr_final_t <- ovr_final_t[-rT_ovr,]

        # Get the max FC score
        ovr_final_t[,3] <- apply(ovr_final_t[,4:ncol(ovr_final_t)],1,function(x) max(as.numeric(x),na.rm=TRUE))

        # Go through the list with the good ROC values and update with fold-change
        # QC <- FALSE
        actives <- c(); preds <- c()
        nvals <- 0
        if(QC!=FALSE) {
            dev.new()
            par(mfrow=c(9,10))
            par(mar=c(2,1,2,1))
        }
        for(i in 4:ncol(ovr_final_t)) {
            # ovr_final_t[,i] <- as.numeric(ovr_final_t[,i])-min(as.numeric(ovr_final_t[,i]))
            ovr_temp <- as.numeric(ovr_final_t[,i])
            # fc_cutoff <- quantile(ovr_temp,probs=0.99)

            # Try Kappa coefficients here for each class's CM as thresholds
            if(useKappas==TRUE) {
                # Get the proper model
                ovr_model <- RF[[3]][[which(names(RF[[3]])==colnames(ovr_final_t)[i])]][[1]]
                ovr_model_kappa <- thresholder(ovr_model,threshold=seq(quantile(ovr_temp,probs=0.98),max(ovr_temp),by=0.01),statistics=c("Kappa"),final=FALSE)      
                fc_cutoff <- ovr_model_kappa[which(ovr_model_kappa[,3]==max(ovr_model_kappa[,3]))[1],2]
            } else {
                iqr <- quantile(ovr_temp,probs=c(0,0.1,0.5,0.9,1))[4]-quantile(ovr_temp,probs=c(0,0.1,0.5,0.9,1))[2]
                # Same as above, set some thresholds here
                if(iqr>0.10) {
                    fc_cutoff <- quantile(ovr_temp,probs=0.955)
                } else if(iqr>0.05) {
                    fc_cutoff <- quantile(ovr_temp,probs=0.96)
                } else if(iqr>0.02) {
                    fc_cutoff <- quantile(ovr_temp,probs=0.97)
                } else if(iqr>0) {
                    fc_cutoff <- quantile(ovr_temp,probs=0.99)
                } else {
                    fc_cutoff <- quantile(ovr_temp,probs=0.995)
                }
            }

            fc_cutoff <- fc_cutoff*fc_mult
            act_temp <- which(ovr_temp>fc_cutoff)
            pre_temp <- rep(colnames(ovr_final_t)[i],length(act_temp))
            actives <- c(actives,act_temp)
            preds <- c(preds,pre_temp)
            ovr_final_t[,i] <- ovr_temp
            nval <- length(which(ovr_temp>fc_cutoff))
            nvals <- nvals+nval

            # QC portion if needed
            if(QC!=FALSE) {
                if(useKappas==FALSE) {
                    print(paste0(substr(colnames(ovr_final_t)[i],1,15)," | iqr=",iqr," | shapiro=",round(shapiro.test(sample(ovr_temp,5000,replace=TRUE))[[1]],3)," | ",round(shapiro.test(sample(ovr_temp,5000,replace=TRUE))[[1]]*iqr,3)," | n=",nval))
                } else {
                    print(paste0(substr(colnames(ovr_final_t)[i],1,15)," | kap=",round(fc_cutoff,4)," | shapiro=",round(shapiro.test(sample(ovr_temp,5000,replace=TRUE))[[1]],3)," | ",round(shapiro.test(sample(ovr_temp,5000,replace=TRUE))[[1]]*fc_cutoff,3)," | n=",nval))
                }
                plot(ovr_temp[sample(1:length(ovr_temp),length(ovr_temp))],pch=19,col=rgb(0,0,0,0.2),main=colnames(ovr_final_t)[i])#,ylim=c(0,30))
                abline(h=fc_cutoff,col="blue",lwd=3)
                text(x=1000,y=.020,labels=paste0("n=",nval),col="red")
            }
        }

        # Subset actives are above the new cutoff
        ovr_final2 <- ovr_final_t[actives,]
        ovr_final2[,2] <- preds
        mahal <- mahal[-rT_ovr]
        mahal <- mahal[actives]

        # Add in the ovr ROC scores
        rocval <- c()
        for(s in 1:nrow(ovr_final2)) {
            nval <- as.numeric(all_rocs[which(toupper(names(all_rocs))==toupper(ovr_final2[s,2]))])
            if(length(nval)>0) rocval[s] <- nval
            else rocval[s] <- NA
        }

        # Recompile everything together
        ovr_final3 <- cbind(ovr_final2[,1:3],ovr_final_t[actives,2],rocval,mahal,ovr_final2[,4:ncol(ovr_final2)])
        colnames(ovr_final3)[4:6] <- c("OG.Pred","AUROC","Pheno.Act")
        rownames(ovr_final3) <- NULL
        
        # Order the various hits by name, and then by probability score
        finalOrder_t <- ovr_final3[order(as.numeric(ovr_final3[,3]),decreasing=TRUE),]
        finalOrder <- finalOrder_t[order(finalOrder_t[,1]),]

        # And a final filter if desired
        # finalOrder <- finalOrder[-which(finalOrder[,4]<0.04 | finalOrder[,3]<0.07),]

        # Save the predictions    
        ovr_preds[[names(Preds)[q+1]]] <- finalOrder
        # ovr_preds[[names(Preds)[q+1]]] <- ovr_final3[finalOrder,]
        # write.table(ovr_final3,"ovr_final_level2.tsv",sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
    }

    # View(ovr_preds[[1]])
    saveRDS(ovr_preds,file=outfile,version=2)
}

# Make a function that subsets the same rows for each list table
subsetRows <- function(listIn,rows=which(cp_compiled$sample_type=="training")) {
    listOut <- listIn
    for(iii in 1:length(listOut)) {
        if(!is.null(nrow(listOut[[iii]]))) {
            listOut[[iii]] <- listOut[[iii]][rows,]
        } else {
            listOut[[iii]] <- listOut[[iii]][rows]
        }
    }
    return(listOut)
}

# ===================================================================================================================================

# Load core functions from functions.R to save space here
source("functions/functions.R")

# Connect to the Cell Painting database with DuckDB
db_path <- "data/db/cellpainting.duckdb"
con <- dbConnect(duckdb(), dbdir=db_path, read_only=TRUE)
avail_tbls <- dbListTables(con)

tbls_to_load <- c("metadata", "mesh", "db", "features_1pass_3", "sample_type", "phenotypic_active", "mahalanobis")
cli_progress_bar("Loading database tables", total=length(tbls_to_load), clear=FALSE)
cp_compiled <- list()
for (tbl in tbls_to_load) {
    if (tbl %in% avail_tbls) cp_compiled[[tbl]] <- as.matrix(dbReadTable(con, tbl))
    cli_progress_update()
}
dbDisconnect(con, shutdown=TRUE)

training <- subsetRows(cp_compiled, rows=which(cp_compiled$sample_type=="training" & cp_compiled$phenotypic_active[,2]=="active"))
training <- global_rownames(training,newnames=training$metadata[,1])
fractions <- subsetRows(cp_compiled, rows=which(cp_compiled$sample_type!="training" & cp_compiled$phenotypic_active[,2]=="active"))
fractions <- global_rownames(fractions,newnames=fractions$metadata[,1])

cli_alert_success(" Data loaded!")

odir <- "results/"
badClasses <- "mucin|interleukin|myc|enicillin|actamase|rothrombin" # Manually remove certain classes

# Set the rows as extract_fraction because CP names aren't set up yet for fractions
rownames(fractions$features_1pass_3) <- fractions$metadata[,"Image_Metadata_SOURCE_COMPOUND_NAME"]

# Fix the NAs
norows <- which(apply(fractions$features_1pass_3,1,function(x) length(which(is.na(x))))>100)
if(length(norows)>0) {
    fractions$features_1pass_3 <- fractions$features_1pass_3[-norows,]
    fractions$metadata <- fractions$metadata[-norows,]
}
nocols <- unique(which(is.na(fractions$features_1pass_3),arr.ind=TRUE)[,2])
if(length(nocols)>0) {
    tvals <- apply(fractions$features_1pass_3[,nocols],2,function(x) median(x,na.rm=TRUE))
    for(kk in 1:length(nocols)) {
        fractions$features_1pass_3[which(is.na(fractions$features_1pass_3[,nocols[kk]])),nocols[kk]] <- tvals[kk]
    }
}

filein <- "RF_Preds_fractions_complete_ovo_auto_roc.Rds"
fileout <- "ovo_preds_fractions_auto_roc.Rds"

# let's set --wd as the working directory, defaulting at the directory that it's run in
if(any(arguments=="--wd")) {
    arg_ind <- which(arguments=="--wd")
    trailing_arg <- arguments[arg_ind[1]+1]
    if(length(arg_ind)>=1) {
        invisible(tryCatch(setwd(trailing_arg),error=function(e) e, finally=function(x) setwd(trailing_arg)))
    }
    cli_alert_success(" Set working directory to {trailing_arg}")
}

# Check if script 06 has already been applied to this database.
# Pass --force on the command line to bypass and re-run anyway.
already_processed <- FALSE
if ("_processing_log" %in% avail_tbls) {
    proc_log <- dbReadTable(con, "_processing_log")
    if (any(proc_log$script == "06.RF_predict_functions")) {
        already_processed <- TRUE
        prev_ts <- tail(proc_log$timestamp[proc_log$script == "06.RF_predict_functions"], 1)
        cli_alert_warning(" Database already processed by script 06 on {prev_ts}")
    }
}
if (already_processed && !any(arguments == "--force")) {
    cli_alert_info(" Skipping. Use --force to re-run script 06 on this database.")
    quit(save="no", status=0)
}

# Predict using the above functions
predict_ovo(ovo_dir=odir,outfile=paste0(odir,filein))

# Analyze using the above functions
analyze_ovo(fileIn=paste0(odir,filein),fc_mult=0.9,ovo_dir=odir,outfile=paste0(odir,fileout),useAcc=FALSE,QC=TRUE,iqr_cutoff=0.02)

# Log this run so re-runs are skipped unless --force is passed
con_log <- dbConnect(duckdb(), dbdir=db_path, read_only=FALSE)
dbWriteTable(con_log, "_processing_log",
             data.frame(
                 script    = "06.RF_predict_functions",
                 timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                 row.names = NULL
             ),
             append=TRUE, overwrite=FALSE)
dbDisconnect(con_log, shutdown=TRUE)
cli_alert_success(" Processing log updated")
