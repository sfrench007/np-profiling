# =============================================================================
# functions/functions.R
#
# Shared helper functions sourced by the main pipeline scripts (01-05).
# Source this file at the top of each script with: source("functions/functions.R")
#
# DEPENDENCY NOTE:
#   Lines ~95-100 load a Carpenter controls Excel file (c_controls.xlsx) from
#   the legends/ directory.  Ensure this file is present before sourcing.
#
# Functions defined in this file:
#
#   --- Metadata injection (fill missing well/compound info from plate legends) ---
#   acJInjectMetadata()          Repair metadata for Anticancer library (J-format plates)
#   ac2InjectMetadata()          Repair metadata for Anticancer library (custom sourced)
#   apinInjectMetadata()         Repair metadata for APIN inhibitor plates
#   bioactivesInjectMetadata()   Repair metadata for Bioactives2 library plates
#   caithnessInjectMetadata()    Repair metadata for Caithness library plates
#   fractionInjectMetadata()     Repair metadata for NP fraction plates
#   pkInjectMetadata()           Repair metadata for pharmacokinetics library plates
#   psychoactivesInjectMetadata() Repair metadata for psychoactives library plates
#   trkbInjectMetadata()         Repair metadata for TrkB library plates
#   unknownsInjectMetadata()     Repair metadata for unknown/uncharacterized plates
#
#   --- Data processing and normalization ---
#   compressData()    Collapse 9 per-well images into a single well-level summary
#   dataScrub()       Center data by subtracting column medians (MAD-scaled)
#   fixNA_zeroSum()   Impute NAs with column medians; optionally remove zero-variance columns
#   rescale()         Plate-level normalization (density peak, IQM, or mean method)
#   subsetAndCombine() Subset and bind cell painting data by type/concentration
#   stackData()       Stack a list of data frames into a single matrix
#
#   --- Quality control and feature selection ---
#   getInactive()     Compute Mahalanobis distances; flag inactive wells
#   getRedundant()    Identify highly correlated (redundant) feature columns
#   plateScan()       Summarize per-plate SD and mean for QC
#   removeRows()      Remove a set of row indices from every element of a cp list
#
#   --- Classification and prediction ---
#   getPredictions()  Predict compound class using a trained model
#   getMatches()      Find nearest tSNE neighbours for unknown samples
#   nTSNE_preds()     Run tSNE and generate neighbour-based predictions
#   makeNNlist()      Build a nearest-neighbour list from a correlation matrix
#   phenomap()        Map rownames to phenotypic annotation
#
#   --- Annotation and labelling ---
#   generateTag()     Create a provenance tag (user, timestamp, version)
#   injectTag()       Attach a provenance tag to a cp list
#   labelScrub()      Clean up class label strings for model compatibility
#   nameScrub()       Simplify compound name strings
#   reverseMatch()    Map predicted labels back to original annotation strings
#   global_rownames() Assign CP_Index values as rownames across all list elements
#
#   --- Utilities ---
#   euclidean()       Euclidean distance between two vectors
#   iqm()             Interquartile mean of a vector
#   generateWells()   Generate a well-plate layout from compound list
#   make_cormat()     Compute and plot a feature correlation matrix
#   makeXYZ()         Export a feature matrix to XYZ coordinate format
#   get_lower_tri()   Extract lower triangle of a correlation matrix
#   get_upper_tri()   Extract upper triangle of a correlation matrix
#   reorder_cormat()  Reorder a correlation matrix by hierarchical clustering
# =============================================================================

# Repair missing metadata for anticancer library CP outputs (J-format plate layout)
acJInjectMetadata <- function(metawithNAs) {
    newsmiles <- c(); newcpdnames <- c()
    for(kk in 1:nrow(metawithNAs)) {
        wrow <- which(acJ_legend[,"J_WELL"]==metawithNAs[kk,"Image_Metadata_CPD_WELL_POSITION"])
        if(length(wrow)>=1) {
            wrow <- wrow[1]
            newsmiles[kk] <- acJ_legend[wrow,"SMILES"]
            newcpdnames[kk] <- acJ_legend[wrow,"Compound"]
        } else {
            newsmiles[kk] <- NA
            newcpdnames[kk] <- NA            
        }
    }

    newcpdnames[which(is.na(newcpdnames))] <- "DMSO"
    metanoNAs <- metawithNAs
    metanoNAs[,3] <- newsmiles
    metanoNAs[,7] <- newcpdnames
    return(metanoNAs)
}

# The AC library that was custom sourced - NEED TO ADD IN SMILES CODES
ac2InjectMetadata <- function(metawithNAs) {
    # if(nchar(dirList[i])==52 & substr(dirList[i],4,8)!="Robot") {
    #     activeLegend <- pa500_1
    # } else if(nchar(dirList[i])!=52 & substr(dirList[i],4,8)!="Robot") {
    #     activeLegend <- pa5_10
    # } else {
        activeLegend <- ac_legend
    # }

    for(kk in 1:nrow(activeLegend)) {
        wellTemp <- which(metawithNAs[,"Image_Metadata_CPD_WELL_POSITION"]==activeLegend[kk,1])
        metawithNAs[wellTemp,"Image_Metadata_CPD_MMOL_CONC"] <- activeLegend[kk,7]   # Concentration
        metawithNAs[wellTemp,"Image_Metadata_SOURCE_COMPOUND_NAME"] <- activeLegend[kk,4]   # Treatment name
    }

    metanoNAs <- metawithNAs
    return(metanoNAs)
}

# The APIN inhibitors that were custom sourced
apinInjectMetadata <- function(metawithNAs) {
    # if(nchar(dirList[i])==52 & substr(dirList[i],4,8)!="Robot") {
    #     activeLegend <- pa500_1
    # } else if(nchar(dirList[i])!=52 & substr(dirList[i],4,8)!="Robot") {
    #     activeLegend <- pa5_10
    # } else {
        activeLegend <- apin_legend
    # }

    for(kk in 1:nrow(activeLegend)) {
        wellTemp <- which(metawithNAs[,"Image_Metadata_CPD_WELL_POSITION"]==activeLegend[kk,1])
        metawithNAs[wellTemp,"Image_Metadata_CPD_MMOL_CONC"] <- activeLegend[kk,7]   # Concentration
        metawithNAs[wellTemp,"Image_Metadata_SOURCE_COMPOUND_NAME"] <- activeLegend[kk,4]   # Treatment name
    }

    metanoNAs <- metawithNAs
    return(metanoNAs)
}

# Quick little function to repair the missing metadata from bioactives2 library CP outputs
# Carpenter controls
library(xlsx)
c_controls <- read.xlsx("legends/c_controls.xlsx",1,startRow=4,header=TRUE)

ac_controls <- unique(c_controls[,6])[-c(1,2,8)]
jc_controls <- c_controls[which(c_controls[,6] %in% ac_controls),]
bioactivesInjectMetadata <- function(metawithNAs) {
    metawithNAs[,"Image_Metadata_CPD_MMOL_CONC"] <- working_conc   # Concentration
    
    for(kk in 1:nrow(metawithNAs)) {
        wellTemp <- which(bioa2_legend[,21]==metawithNAs[kk,"Image_Metadata_CPD_WELL_POSITION"] & bioa2_legend[,20]==metawithNAs[kk,"Image_Metadata_PlateID"])
        # In case it's not there
        if(length(wellTemp)==0) {
            # It might be a carpenter control...
            if(any(jc_controls[,1]==metawithNAs[kk,"Image_Metadata_CPD_WELL_POSITION"])) {
                cc_work <- jc_controls[which(jc_controls[,1]==metawithNAs[kk,"Image_Metadata_CPD_WELL_POSITION"]),]
                metawithNAs[kk,"Image_Metadata_SOURCE_COMPOUND_NAME"] <- cc_work[6]
                metawithNAs[kk,"Image_Metadata_CPD_MMOL_CONC"] <- cc_work[7]
            } else {
                # Or else otherwise it's DMSO
                metawithNAs[kk,"Image_Metadata_SOURCE_COMPOUND_NAME"] <- "DMSO"
                metawithNAs[kk,"Image_Metadata_CPD_MMOL_CONC"] <- "0.5%"
                metawithNAs[kk,"Image_Metadata_CPD_SMILES"] <- "CS(=O)C"
            }
        } else {
            metawithNAs[kk,"Image_Metadata_CPD_SMILES"] <- bioa2_legend[wellTemp,3][1]
        }
    }

    metanoNAs <- metawithNAs
    return(metanoNAs)
}

# Quick little funciton to repair the missing metadata from caithness library CP outputs
caithnessInjectMetadata <- function(metawithNAs) {
    # if(nchar(dirList[i])==52 & substr(dirList[i],4,8)!="Robot") {
    #     activeLegend <- pa500_1
    # } else if(nchar(dirList[i])!=52 & substr(dirList[i],4,8)!="Robot") {
    #     activeLegend <- pa5_10
    # } else {
        activeLegend <- cn_legend
    # }

    for(kk in 1:nrow(activeLegend)) {
        wellTemp <- which(metawithNAs[,"Image_Metadata_CPD_WELL_POSITION"]==activeLegend[kk,1])
        metawithNAs[wellTemp,"Image_Metadata_CPD_MMOL_CONC"] <- activeLegend[kk,7]   # Concentration
        metawithNAs[wellTemp,"Image_Metadata_SOURCE_COMPOUND_NAME"] <- activeLegend[kk,4]   # Treatment name
    }

    metanoNAs <- metawithNAs
    return(metanoNAs)
}


# Merge 9 images of data here, otherwise it will be insane
compressData <- function(dataIn,tWells=d1[,"Image_Metadata_CPD_WELL_POSITION"],message="the",type=NA,all=FALSE) {
    print(paste0("Combining ",message," dataset replicates..."))
    # There unfortunately isn't always 9 images....so just take what we have for each well
    # for(iii in seq(1,nrow(dataIn),by=9)) {
    tunq <- unique(tWells)
    if(length(tunq)!=0) {
        pb <- txtProgressBar(min=1,max=length(tunq),style=3)
    }
    wellstemp <<- c()
    for(iii in 1:length(tunq)) {
        awells <- which(tWells==tunq[iii])
        wellstemp[iii] <<- awells[1]
        if(all==FALSE) {
            if(is.na(type)) {
                lineOut <- apply(dataIn[awells,],2,function(x) {
                    mt <- mean(x,na.rm=TRUE)
                    if(is.na(mt)) {
                        if(any(!is.na(x))) {
                            mt <- x[which(!is.na(x))[1]]
                        } else {
                            mt <- NA
                        }
                    }
                    mt
                })
            } else {
                lineOut <- apply(dataIn[awells,],2,function(x) sum(x,na.rm=TRUE))
            }
        } else {
            lineOut <- dataIn[awells,]
        }

        if(iii==1) dataOut <- lineOut
        else dataOut <- rbind(dataOut,lineOut)
        if(length(tunq)!=0) {
            setTxtProgressBar(pb,iii)
        }
    }
    return(dataOut)
}

# Scrub out column means if desired (centers data)
dataScrub <- function(dataSetIn,rn) {
	# Scrub the colMeans to better normalize the dataset
	oo <- dataSetIn
	rm <- apply(oo,2,function(x) median(x,na.rm=TRUE))
	oo_sweep <- sweep(oo,2,rm)
	sx <- apply(oo_sweep,2,function(x) mad(x,na.rm=TRUE))
	oo_sweep <- sweep(oo_sweep,2,sx,"/")
	rownames(oo_sweep) <- rn
	if(any(is.na(oo_sweep))) {
		temp <- which(is.na(oo_sweep),arr.ind=TRUE)
		colmeans_t <- apply(oo_sweep,2,function(xx) mean(xx,na.rm=TRUE))
        for(ii in 1:nrow(temp)) {
			cmTemp <- colmeans_t[temp[ii,2]]
			oo_sweep[temp[ii,1],temp[ii,2]] <- cmTemp
        }
	}
	return(oo_sweep)
}

# Quick Euclidean distance code
euclidean <- function(a,b) {
    ed <- sqrt(sum((a-b)^2))
    return(ed)
}

# Check for zerosum columns and remove them
fixNA_zeroSum <- function(dataIn=cp_features,removeZeroVar=TRUE,repairRowNA=FALSE) {
    medTemp <- apply(dataIn,2,function(x) median(x,na.rm=TRUE))
    dataOut <- dataIn
    # Remove the columns (not the treatments) that have NAs
    dataNAs <<- as.numeric(na.omit(as.numeric(names(which(summary(as.factor((which(is.na(dataIn),arr.ind=TRUE)[,2])))>100)))))
    if(length(dataNAs)>0) dataOut <- dataOut[,-dataNAs]
    zsRows <<- unique(which(is.na(dataOut),arr.ind=TRUE)[,1])
    if(length(zsRows)>0 & repairRowNA==FALSE) {
        dataOut <- dataOut[-zsRows,]
    } else if(length(zsRows)>0 & repairRowNA==TRUE) {
        zsCols <- which(is.na(dataOut[zsRows,]),arr.ind=TRUE)[,2]
        dataOut[zsRows,zsCols] <- medTemp[zsCols]
    }

    # Remove zerosums
    zerosum <- which(apply(dataOut,2,sd)==0)
    if(length(zerosum>0)) dataOut <- dataOut[,-zerosum]

    if(removeZeroVar==TRUE) {
        zv <- apply(dataOut,2,function(x) sd(x,na.rm=TRUE))
        if(any(zv==0)) {
            zv_rm <- which(zv==0)
            dataOut <- dataOut[,-zv_rm]
        }
    }

    return(dataOut)
}

# Quick little function to repair the missing metadata from present fraction legends
fractionInjectMetadata <- function(metawithNAs) {
    # These are the metadata columns:
    # Image_Metadata_CPD_MMOL_CONC
    # Image_Metadata_CPD_PLATE_MAP_NAME
    # Image_Metadata_CPD_SMILES
    # Image_Metadata_CPD_WELL_POSITION
    # Image_Metadata_MAC_ID
    # Image_Metadata_PlateID
    # Image_Metadata_SOURCE_NAME
    # Image_Metadata_Site

    # > head(legend_working)
    #   Well Row Column      Code Fraction   Treatment Concentration Units
    # 1  A01   A     01      DMSO                 DMSO           0.5     %
    # 2  A02   A     02      DMSO                 DMSO           0.5     %
    # 3  A03   A     03 KCCC-1071        C KCCC-1071_C           0.5     %
    # 4  A04   A     04 KCCC-1088        C KCCC-1088_C           0.5     %
    # 5  A05   A     05 KCCC-1071        1 KCCC-1071_1           0.5     %
    # 6  A06   A     06 KCCC-1088        1 KCCC-1088_1           0.5     %

    if(any(ls()=="legend_working")) {
        for(kk in 1:nrow(legend_working)) {
            wellTemp <- which(metawithNAs[,"Image_Metadata_CPD_WELL_POSITION"]==legend_working[kk,1])
            metawithNAs[wellTemp,"Image_Metadata_CPD_MMOL_CONC"] <- legend_working[kk,7]   # Concentration
            metawithNAs[wellTemp,"Image_Metadata_CPD_PLATE_MAP_NAME"] <- legend_working[kk,4]   # Plate name
            metawithNAs[wellTemp,"Image_Metadata_PlateID"] <- legend_working[kk,4]   # Plate name
            metawithNAs[wellTemp,"Image_Metadata_SOURCE_COMPOUND_NAME"] <- legend_working[kk,6]   # Treatment name
        }
    } else {
        newname <- paste(metawithNAs[,"Image_Metadata_PlateID"],metawithNAs[,"Image_Metadata_CPD_WELL_POSITION"],sep="-")
        metawithNAs[,"Image_Metadata_CPD_PLATE_MAP_NAME"] <- newname   # Plate name
        metawithNAs[,"Image_Metadata_SOURCE_COMPOUND_NAME"] <- newname   # Treatment name        
    }
    metanoNAs <- metawithNAs
    return(metanoNAs)
}

generateTag <- function(hqp="Anonymous") {
    # Generate tag IDs
    tag_id <- list()
    stime <- Sys.time()
    tag_id$hash <- cli::hash_sha1(stime)
    tag_id$date <- stime
    tag_id$hqp <- hqp
    return(tag_id)
}

# Quick QC function to make a mock plate
generateWells <- function(compounds,platerows) {
	# 1-2, 23:24 are DMSO
	# 5-6 are crude
	# 7-22 are F1-F8 (n=2; F1,F1,F2,F2,...,F8,F8)
	useWells <- cp_fractions$metadata[which(substr(cp_fractions$metadata[,5],1,1) %in% platerows & substr(cp_fractions$metadata[,5],2,3) %in% sprintf("%02d",c(1:2,5:24))),5]
	cmpds <- c()
	for(ii in 1:length(compounds)) {
		cmpds <- c(cmpds,rep("DMSO",2),rep(paste0(compounds[ii],"_C"),2),rep(paste0(compounds[ii],"_F1"),2),rep(paste0(compounds[ii],"_F2"),2),rep(paste0(compounds[ii],"_F3"),2),rep(paste0(compounds[ii],"_F4"),2),rep(paste0(compounds[ii],"_F5"),2),rep(paste0(compounds[ii],"_F6"),2),rep(paste0(compounds[ii],"_F7"),2),rep(paste0(compounds[ii],"_F8"),2),rep("DMSO",2))
	}
	return(cbind(cmpds,useWells))
}

# PCA test to capture inactive molecules via Malanobis distance
getInactive <- function(dataIn=cc,sd_cutoff=2,plotQC=FALSE,quantiles=0.25) {
    prtemp <<- princomp(dataIn,cor=TRUE)
    proportions <- prtemp$sdev^2/sum(prtemp$sdev^2) 
    csum_proportions <- cumsum(proportions)
    # Get how many columns to capture 90% of variances
    v90 <- which(csum_proportions>=0.90)[1]
    vals_t0 <- predict(prtemp)[,1:v90]
    vals <<- mahalanobis(vals_t0,apply(vals_t0,2,mean),cov(vals_t0),tol=1e-50)

    vq <- quantile(vals,probs=seq(0,1,quantiles))
    midpoint <- iqm(vals[which(vals>vq[2] & vals<vq[4])])
    vsd <- sd(vals[which(vals>vq[2] & vals<vq[4])])
    high_c <- midpoint+(vsd*sd_cutoff)
    inactive_temp <- which(vals<high_c) 

    # Plot index plot for QC and check variances if desired
    if(plotQC==TRUE) {
        par(mfrow=c(1,2)); par(mar=c(4,4,1,1))
        plot(prtemp)
        plot(vals)
        abline(h=high_c,col="red");abline(h=midpoint,col="green")
    }
    
    return(inactive_temp)
}

# Subset lower triangle of correlation matrix
get_lower_tri <- function(cormat) {
	cormat[upper.tri(cormat)] <- NA
	return(cormat)
}

# Function to get the "priority" actives using euclidean distances
getMatches <- function(tsneIn=cp_tsne$Y,labels=rownames(cc_training),whichUnknowns=fractions,cutoff=2) {
	keepers <- c()
	matches <- c()
	distances <- c()
	fract <- c()
	rncc3 <- make.unique(labels)

	# Compare everything to generate a massive list of predictions based on cutoff
    counter <- 1
    for(iii in 1:length(whichUnknowns)) {
        yy <- tsneIn[whichUnknowns[iii],]
        edist <- apply(tsneIn[-whichUnknowns,],1,function(x) euclidean(x,yy))

        if(any(edist<=cutoff)) {
            keepers <- whichUnknowns[iii]
            matches <- rncc3[-whichUnknowns][which(edist<=cutoff)]
            distances <- edist[which(edist<=cutoff)]
            fract <- rncc3[whichUnknowns][iii]
            if(counter==1) {
            	thelist <- cbind(keepers,matches,distances,fract)
                counter <- counter+1
            } else {
            	thelist <- rbind(thelist,cbind(keepers,matches,distances,fract))
            }
        }
    }
	
	# Generate a list of matches based on the euclidean distances
    # if(any(substr(fract,nchar(fract)-1,nchar(fract)-1)==".")) {
    #     fixthese <- which(substr(fract,nchar(fract)-1,nchar(fract)-1)==".")
    #     fract[fixthese] <- substr(fract[fixthese],1,nchar(fract[fixthese])-2)
    # }
	thelist <- thelist[order(thelist[,3],decreasing=FALSE),]
    checkbad <- substr(thelist[,2],nchar(thelist[,2])-3,nchar(thelist[,2])-3) 
    if(any(checkbad=="_" | str_detect(toupper(thelist[,2]),"KCCC") | str_detect(toupper(thelist[,2]),"KCHC"))) {
        thelist <- thelist[-which(checkbad=="_" | str_detect(toupper(thelist[,2]),"KCCC") | str_detect(toupper(thelist[,2]),"KCHC")),]
    }
    return(thelist)
}

# Subset upper triangle of correlation matrix
get_upper_tri <- function(cormat) {
	cormat[lower.tri(cormat)]<- NA
	return(cormat)
}

# Dendrogram-based predictions.  Might be best just to look at the heatmap, but how to subset a smaller area?
t_cormat <- NA
cormat_reorder <- NA
getPredictions <- function(drug,processedSet=cp_scale,cutoff=0.6,plotTree=FALSE,plotHeat=FALSE) {
	# This will need to fully process each time from an Rdata set...
	# Also use fuzzy searching with stringr/stringi for the drug
	# Prepare the correlation matrix if it doesn't exist, and make.unique the names
	if(length(t_cormat)<=1) {
		rownames(processedSet) <- make.unique(toupper(rownames(processedSet)))
		t_cormat <<- cor(na.omit(t(processedSet)))
		t_dd <<- as.dist((1-t_cormat)/2)
		t_hc <<- hclust(t_dd)
		cormat_reorder <<- t_cormat[t_hc$order,t_hc$order]
	}

	# Cluster the dendrogram by the cutoff
	avg_dend_obj <- as.dendrogram(t_hc)
	t_groups <- cutree(avg_dend_obj,h=cutoff)
	if(plotTree==TRUE) {
		avg_col_dend <- color_branches(avg_dend_obj,h=cutoff)
		plot(avg_col_dend)
		abline(h=cutoff,lwd=2,lty=2,col=rgb(1,0,0,0.5))
	}

	# Then, search for the drug (fuzzy) in the cormat, get the highest correlations, and the dendrogram matches
	findDrug_t_hc <- which(stri_detect_fixed(t_hc$labels,pattern=toupper(drug)))
	drugpreds <- list()
	for(ii in 1:length(findDrug_t_hc)) {
		tgtemp <- t_groups[findDrug_t_hc][ii]
		tmatches <- names(which(t_groups==tgtemp))
		drugpreds[[t_hc$labels[findDrug_t_hc[ii]]]]$correlations <- tmatches
		if(plotHeat==TRUE) {
			cc_t <- cormat_reorder[which(rownames(cormat_reorder) %in% tmatches),]
			tcols <- colorRampPalette(c(rgb(0,0,0),rgb(0.0,0.0,0.3),rgb(0.9,0.9,0)))
			dev.new()
			heatmap.2(as.matrix(t(cc_t)),trace="none",Rowv=FALSE,Colv=FALSE,dendrogram="none",labCol="",col=tcols(200),keysize=0.5,margins=c(2,13))
		}
	}
	return(drugpreds)
}

# Identify redundant features based on Pearson correlation cutoffs
getRedundant <- function(dataIn=cc,r_cutoff=0.99) {
    cormat_temp <- cor(na.omit(dataIn))
    redundant_temp <- which(cormat_temp>r_cutoff & cormat_temp<1,arr.ind=TRUE)
    redundant_out <- c()
    for(iii in unique(redundant_temp[,1])) {
        moreF <- which(cormat_temp[,iii]>r_cutoff & cormat_temp[,iii]<1)
        redundant_out <- c(redundant_out,moreF)
    }
    redundant_out <- sort(unique(redundant_out))
    return(redundant_out)
}

# Change rownames of list items to match a common key
global_rownames <- function(listIn,listItems=1:length(listIn),newnames=listIn$metadata[,1]) {
    qt <- sapply(listIn[listItems],function(x) {
        if(length(nrow(x))==0) {
            length(x)
        } else {
            nrow(x)
        }
    })
    if(any(qt!=length(newnames) & names(listIn)!="tag")) {
        bads <- which(qt!=length(newnames))
        cat("\nRow length mismatch! Check list items",bads,"and retry\n")
        return(NULL)
    } else {
        listOut <- listIn; listItems2 <- listItems
        if(any(names(listIn)=="tag")) {
            listItems2 <- listItems[-which(names(listIn)=="tag")]
        }
        for(iii in 1:length(listItems2)) {
            if(length(nrow(listOut[[listItems2[iii]]]))==0) {
                names(listOut[[listItems2[iii]]]) <- newnames
            } else {
                rownames(listOut[[listItems2[iii]]]) <- newnames
            }
        }
    }
    return(listOut)
}

# Add in a tag
injectTag <- function(listIn,hqp="Anonymous") {
    tempList <- listIn
    listOut <- list()
    listOut$tag <- generateTag(hqp=hqp)
    for(iii in 1:length(tempList)) {
        listOut[[iii+1]] <- tempList[[iii]]
    }
    names(listOut)[2:length(listOut)] <- names(listIn)
    return(listOut)
}

# Interquartile mean
iqm <- function(dataIn) {
    qT <- quantile(dataIn,na.rm=TRUE)
    iqmOut <- mean(dataIn[which(dataIn>qT[2] & dataIn<qT[4])],na.rm=TRUE)
    return(iqmOut)
}

# Remove periods and numbers after calling make.unique
labelScrub <- function(labelsIn) {
    labelsOut <- sapply(labelsIn,function(x) {
        tval <- strsplit(x,",")[[1]]
        toString(tval[length(tval)])
    })
    return(labelsOut)
}

# Create a correlation matrix based on a data frame (correlates columns)
make_cormat <- function(dataIn,fontSize=8,colscale=c(rgb(0,0,0),rgb(0.0,0.0,0.3),rgb(0.9,0.9,0)),square=FALSE) {
	# Make and melt correlation matrix
	cormat <<- cor(na.omit(dataIn))
	melted_cormat <- melt(reorder_cormat(cormat))
	melted_cormat <- na.omit(melted_cormat)

	# Default color scale (black=low, blue=mid, yellow=high)
	lowCol <- colscale[1]
	midCol <- rgb(0.0,0.0,0.3)
	highCol <- rgb(0.9,0.9,0)

	vJust <- 0.5
	hJust <- 1

	if(square==TRUE) {
		ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
		geom_tile() +
		# geom_tile(color = rgb(0,0,0,0.3)) +
		scale_fill_gradient2(low = lowCol, mid = midCol, high = highCol,limit = c(-1,1), name="Pearson\nCorrelation") +
		theme_minimal() +
		theme(axis.ticks = element_blank(), axis.text.y = element_blank()) +
		theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
		theme(axis.text.x = element_text(angle = 90, vjust = vJust, size = fontSize, hjust = hJust)) +	theme(axis.text.y = element_text(angle = 0, vjust = vJust, size = fontSize, hjust = hJust)) +
		coord_fixed()
	} else {
		ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
		geom_tile() +
		# geom_tile(color = rgb(0,0,0,0.3)) +
		scale_fill_gradient2(low = lowCol, mid = midCol, high = highCol,limit = c(-1,1), name="Pearson\nCorrelation") +
		theme_minimal() +
		theme(axis.ticks = element_blank(), axis.text.y = element_blank()) +
		theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
		theme(axis.text.x = element_text(angle = 90, vjust = vJust, size = fontSize, hjust = hJust)) +
		theme(axis.text.y = element_text(angle = 0, vjust = vJust, size = fontSize, hjust = hJust))
	}
}

# Network list generation
makeNNlist <- function(corMatrixIn,cutoff=0.8) {
    cnlist <- colnames(corMatrixIn)
    tcounter <- 1
    for(iii in 1:ncol(corMatrixIn)) {
        node1 <- cnlist[iii]
        n2temp <- which(as.numeric(corMatrixIn[,iii])>cutoff)
        if(length(n2temp)>0) {
            node2 <- rownames(corMatrixIn)[n2temp]
            dists <- corMatrixIn[n2temp,iii]
            tlineout <- cbind(node1,node2,dists)
            if(tcounter==1) {
                nodeList <- tlineout
                tcounter <- tcounter + 1
            } else {
                nodeList <- rbind(nodeList,tlineout)
            }
        }
    }
    colnames(nodeList) <- c("Node1","Node2","Distance")
    rownames(nodeList) <- NULL
    nodeList <- nodeList[-which(nodeList[,"Distance"]==1),]
    return(nodeList)
}

# Make pymol-compatible XYZ files for raytraced rendering
makeXYZ <- function(dataIn,classes,filename="temp.xyz",randomize=FALSE) {
	library(PeriodicTable)
	data(periodicTable)
	con1 <- file(filename, "w", encoding = "UTF-8")
	# Write the number of rows in the file first as line 1, then a blank line
	cat(nrow(dataIn), file=con1)
	cat("\n\n", file=con1)

	# Calculate the 'elements' used
	useEl <- periodicTable[3:nrow(periodicTable),2]
	useEl <- useEl[-which(useEl %in% c("C","Fe","N","O"))]  # For the existing .pse file, these elements have special maps
	if(randomize==TRUE) useEl <- sample(useEl,length(useEl))
    unqC <- unique(classes)
	element <- c()
	for(ii in 1:length(classes)) {
		element[ii] <- useEl[which(unqC==classes[ii])]
	}

	# Write the data, starting with the 'element'
	for(ii in 1:nrow(dataIn)) {
		cat(paste(element[ii],dataIn[ii,1],dataIn[ii,2],dataIn[ii,3],sep="\t"),file=con1)
		cat("\n", file=con1)
	}
	close(con1)
}

# Remove the 'make.unique" .X suffix, or a leading " "
nameScrub <- function(nameIn) {
	nameOut <- nameIn
	if(any(substr(nameIn,nchar(nameIn)-1,nchar(nameIn)-1)==".")) {
		fixEm <- which(substr(nameIn,nchar(nameIn)-1,nchar(nameIn)-1)==".")
		nameOut[fixEm] <- substr(nameIn[fixEm],1,nchar(nameIn[fixEm])-2)
	}
	if(any(substr(nameOut,1,1)==" ")) {
		fixEm <- which(substr(nameOut,1,1)==" ")
		nameOut[fixEm] <- substr(nameOut[fixEm],2,nchar(nameOut[fixEm]))
	}
	return(nameOut)
}

# Get phenotypes from names
phenomap <- function(rownamesIn) {
	startX <- which(colnames(legend)=="X1")
	rownamesIn <- nameScrub(rownamesIn)
	everythingT <- matrix(data=NA,ncol=5,nrow=length(rownamesIn))

	for(ii in 1:length(rownamesIn)) {
		nameT <- which(legend[,"Synonyms"]==rownamesIn[ii])[1]
		# If it finds a match, get the bioactivity, and the NUMBER of activities it has - return EVERYTHING
		colnames(everythingT) <- c("Name","Main_MOA","MOA_Color","Total_MOA","Complete_MOA")
		name_temp <- rownamesIn[ii]

		# If it exists in the pharmakon...
		if(!is.na(nameT)) {
			mainMOA_temp <- toString(legend[nameT,"X1"])
			totalMOA_temp <- 8-length(which(legend[nameT,startX:(startX+7)]==""))
			completeMOA_temp <- toString(legend[nameT,"Bioactivity"])
		} else {
		# Otherwise tap the AC legend
			nameT <- which(ac_legend[,"Compound"]==rownamesIn[ii])[1]
			mainMOA_temp <- "anticancer"
			totalMOA_temp <- 1
			completeMOA_temp <- toString(ac_legend[nameT,"Pathway"])
		}
		everythingT[ii,] <- c(name_temp,mainMOA_temp,NA,totalMOA_temp,completeMOA_temp)
	}
	npheno <- length(unique(everythingT[,2]))
	ncolors <- rainbow(npheno)
	randcol <- sample(ncolors,npheno)	
	unqT <- unique(everythingT[,2])
	newcols <- everythingT[,3]
	for(ii in 1:length(unqT)) {
		newcols[which(everythingT[,2]==unqT[ii])] <- randcol[ii]
	}
	everythingT[,3] <- newcols
	return(everythingT)
}

# Quick little funciton to repair the missing metadata from pharmakon library CP outputs 
pkInjectMetadata <- function(metawithNAs,pk_plate) {
    activeLegend <- pk_legend[which(pk_legend[,"Plate"]==toString(pk_plate)),]

    for(kk in 1:nrow(activeLegend)) {
        wellTemp <- which(metawithNAs[,"Image_Metadata_CPD_WELL_POSITION"]==activeLegend[kk,"X384.well"])
        metawithNAs[wellTemp,"Image_Metadata_MAC_ID"] <- activeLegend[kk,"Batch.Compound.Batch.ID"]   # MAC-ID
        metawithNAs[wellTemp,"Image_Metadata_CPD_SMILES"] <- activeLegend[kk,"SMILES"]   # SMILES
        metawithNAs[wellTemp,"Image_Metadata_CPD_MMOL_CONC"] <- 5   # Concentration
        metawithNAs[wellTemp,"Image_Metadata_CPD_PLATE_MAP_NAME"] <- activeLegend[kk,"UID.384"]   # Plate name
        metawithNAs[wellTemp,"Image_Metadata_PlateID"] <- toString(pk_plate)   # Plate name
        metawithNAs[wellTemp,"Image_Metadata_SOURCE_COMPOUND_NAME"] <- activeLegend[kk,"Compound"]   # Treatment name
    }
    metanoNAs <- metawithNAs
    return(metanoNAs)
}

# Quick little funciton to repair the missing metadata from psychoactive library CP outputs
psychoactivesInjectMetadata <- function(metawithNAs) {
    # if(nchar(dirList[i])==52 & substr(dirList[i],4,8)!="Robot") {
    #     activeLegend <- pa500_1
    # } else if(nchar(dirList[i])!=52 & substr(dirList[i],4,8)!="Robot") {
    #     activeLegend <- pa5_10
    # } else {
        activeLegend <- pa_legend
    # }

    for(kk in 1:nrow(activeLegend)) {
        wellTemp <- which(metawithNAs[,"Image_Metadata_CPD_WELL_POSITION"]==activeLegend[kk,1])
        metawithNAs[wellTemp,"Image_Metadata_CPD_MMOL_CONC"] <- activeLegend[kk,7]   # Concentration
        metawithNAs[wellTemp,"Image_Metadata_SOURCE_COMPOUND_NAME"] <- activeLegend[kk,4]   # Treatment name
    }

    metanoNAs <- metawithNAs
    return(metanoNAs)
}

# Function to randomize a vector to simplify the sampling syntax
randomize <- function(vectorIn) {
    vectorTemp <- sample(1:length(vectorIn),length(vectorIn))
    vectorOut <- vectorIn[vectorTemp]
    names(vectorOut) <- vectorTemp
    return(vectorOut)
}

# Function to update lists by removing the same number of members from each matrix/vector
# (Note: Does not work with sub-lists yet)
removeRows <- function(listIn,removeThese) {
    listOut <- listIn
    for(iii in 1:length(listIn)) {
        if(names(listIn)[iii]=="tag") {
            listOut[[iii]] <- listIn[[iii]]
        } else {
            if(!is.null(nrow(listIn[[iii]]))) {
                listOut[[iii]] <- listIn[[iii]][-removeThese,]
            } else {
                listOut[[iii]] <- listIn[[iii]][-removeThese]
            }
        }
    }
    return(listOut)
}

reorder_cormat <- function(cormat) {
	# Use correlation between variables as distance
	dd <- as.dist((1-cormat)/2)
	hc <- hclust(dd)
	cormat <- cormat[hc$order, hc$order]
	return(cormat)
}

# Make a custom rescaling method to test the 0-1 scaling
ccon <- c_controls[which(c_controls[,4]=="DMSO" | substr(c_controls[,1],2,3) %in% c("23","24")),1]
rescale <- function(dataIn,cp_features2=cp_features,adjval=0.5,meth="iqm") {
    if(any(ls()=="cp_meta")) {
        mval <- cp_meta[which(rownames(cp_features2) %in% rownames(dataIn)),5]
    } else {
        mval <- cp_working$metadata[which(rownames(cp_features2) %in% rownames(dataIn)),5]
    }
    ewells <- which(mval %in% ccon)
    if(meth=="median") {
        mval2 <- apply(dataIn[ewells,],2,function(xx) median(xx,na.rm=TRUE))
    } else if(meth=="iqm") {
        mval2 <- apply(dataIn[ewells,],2,function(xx) iqm(xx))
    } else if(meth=="mean") {
        mval2 <- apply(dataIn[ewells,],2,function(xx) mean(xx,na.rm=TRUE))
    }
    
    if(meth=="density") {
        dataOut <- apply(dataIn,2,function(xx) {
            xxx <- as.numeric(xx)
            # Try dividing by the peak of the density distribution
            dtemp <- density(xxx,na.rm=TRUE,adjust=adjval)
            as.numeric(xxx/(dtemp$x[which(dtemp$y==max(dtemp$y))][[1]]))
        })
    } else {
        dataOut <- sapply(1:ncol(dataIn),function(xx) {
            as.numeric(dataIn[,xx]/mval2[[xx]])
        })
    }
    return(dataOut)
}

# Needs to be updated for pharmakon and not just AC library
reverseMatch <- function(newLabsIn) {
	# Match to the original label before scrubbing
	oglabel <- cp_pharmakon$metadata[which(labelScrub(cp_pharmakon$metadata[,8])==newLabsIn)[1],8]
	acmatch <- which(ac_legend[,"Compound"]==oglabel)
	labInfo <- ac_legend[acmatch,]
	return(labInfo)
}

stackData <- function(dataIn) {
    for(ii in 1:nrow(dataIn)) {
        wout <- t(as.matrix(dataIn[ii,]))
        if(ii==1) dataOut <- wout
        else dataOut <- rbind(dataOut,wout)
    }
    return(dataOut)
}

subsetAndCombine <- function(inputSet,dtype,pkplate=NA,cell_cutoff=NA,ptAll=FALSE) {
    # As a test, start with median only
    meta_ind <- which(substr(colnames(inputSet),1,10)=="Image_Meta")
    counts_ind <- which(substr(colnames(inputSet),1,5)=="Count")
    med_ind <- which(substr(colnames(inputSet),1,6)=="Median")
    # mean_ind <- which(substr(colnames(inputSet),1,4)=="Mean")
    text_ind <- which(substr(colnames(inputSet),1,7)=="Texture")
    sd_ind <- which(substr(colnames(inputSet),1,5)=="StDev")

    if(dtype=="PK") {
        if(is.na(pkplate)) {
            d1_meta <- inputSet[,meta_ind]
        } else {
            d1_meta <- pkInjectMetadata(inputSet[,meta_ind],pk_plate=pkplate)
        }
    } else if(dtype=="AC") {
        d1_meta <- acInjectMetadata(inputSet[,meta_ind])
    } else if(dtype=="ACJ") {
        d1_meta <- acJInjectMetadata(inputSet[,meta_ind])
    } else if(dtype=="FR") {
        d1_meta <- fractionInjectMetadata(inputSet[,meta_ind])
    } else if(dtype=="PA") {
        d1_meta <- psychoactivesInjectMetadata(inputSet[,meta_ind])
    } else if(dtype=="BA") {
        d1_meta <- bioactivesInjectMetadata(inputSet[,meta_ind])
    } else if(dtype=="UN") {
        d1_meta <- unknownsInjectMetadata(inputSet[,meta_ind])
    } else if(dtype=="CA") {
        d1_meta <- caithnessInjectMetadata(inputSet[,meta_ind])
    } else if(dtype=="TR") {
        d1_meta <- trkbInjectMetadata(inputSet[,meta_ind])
    } else if(dtype=="AC2") {
        d1_meta <- ac2InjectMetadata(inputSet[,meta_ind])
    } else if(dtype=="AP") {
        d1_meta <- apinInjectMetadata(inputSet[,meta_ind])
    }
    
    d1_counts <- inputSet[,counts_ind]
    d1_median <- inputSet[,med_ind]
    # d1_mean <- inputSet[,mean_ind]
    d1_text <- inputSet[,text_ind]
    d1_sdT <- inputSet[,sd_ind]

    # Scrub out any wells that have *too many* cells
    if(!is.na(cell_cutoff)) {
        cellCountsOnly <- as.numeric(d1_counts[,"Count_Cells"])
        if(any(cellCountsOnly>=cell_cutoff)) {
            # First check that we're not eliminating *all* of the images for a treatment
            scrubThese <- which(cellCountsOnly>=cell_cutoff)
            # Get the wells that are now MISSING, and remove them from the purge
            if(length(scrubThese)>0) {
                d1_meta <- d1_meta[-scrubThese,]
                d1_counts <- d1_counts[-scrubThese,]
                d1_median <- d1_median[-scrubThese,]
                # d1_mean <- d1_mean[-scrubThese,]
                d1_text <- d1_text[-scrubThese,]
                d1_sdT <- d1_sdT[-scrubThese,]   

                # Combine the technical replicates for each of the treatements
                d1_counts_c <- compressData(d1_counts,message="counts",type="counts",all=ptAll,tWells=d1[-scrubThese,"Image_Metadata_CPD_WELL_POSITION"])
                d1_med_c <- compressData(d1_median,message="median",all=ptAll,tWells=d1[-scrubThese,"Image_Metadata_CPD_WELL_POSITION"])
                # d1_mea_c <- compressData(d1_mean,message="mean",all=ptAll,tWells=d1[-scrubThese,"Image_Metadata_CPD_WELL_POSITION"])
                d1_text_c <- compressData(d1_text,message="texture",all=ptAll,tWells=d1[-scrubThese,"Image_Metadata_CPD_WELL_POSITION"])
                d1_sd_c <- compressData(d1_sdT,message="sd",all=ptAll,tWells=d1[-scrubThese,"Image_Metadata_CPD_WELL_POSITION"])
            }
        } else {
            d1_counts_c <- compressData(d1_counts,message="counts",type="counts",all=ptAll)
            d1_med_c <- compressData(d1_median,message="median",all=ptAll)
            # d1_mea_c <- compressData(d1_mean,message="mean",all=ptAll)
            d1_text_c <- compressData(d1_text,message="texture",all=ptAll)
            d1_sd_c <- compressData(d1_sdT,message="sd",all=ptAll)        
        }

    } else {
        # Combine the technical replicates for each of the treatements
        d1_counts_c <- compressData(d1_counts,message="counts",type="counts",all=ptAll)
        d1_med_c <- compressData(d1_median,message="median",all=ptAll)
        # d1_mea_c <- compressData(d1_mean,message="mean",all=ptAll)
        d1_text_c <- compressData(d1_text,message="texture",all=ptAll)
        d1_sd_c <- compressData(d1_sdT,message="sd",all=ptAll)
    }

    if(ptAll==FALSE) {
        d1_meta_c <- d1_meta[wellstemp,]
    } else {
        d1_meta_c <- d1_meta
    }

    # Combine measurements
    d1_out <- list()
    d1_out$metadata <- d1_meta_c
    d1_out$counts <- d1_counts_c
    d1_out$median <- d1_med_c
    # d1_out$mean <- d1_mea_c
    d1_out$texture <- d1_text_c
    d1_out$sd <- d1_sd_c

    return(d1_out)
}

# The TrkB inhibitors that were custom sourced
trkbInjectMetadata <- function(metawithNAs) {
    # if(nchar(dirList[i])==52 & substr(dirList[i],4,8)!="Robot") {
    #     activeLegend <- pa500_1
    # } else if(nchar(dirList[i])!=52 & substr(dirList[i],4,8)!="Robot") {
    #     activeLegend <- pa5_10
    # } else {
        activeLegend <- trkb_legend
    # }

    for(kk in 1:nrow(activeLegend)) {
        wellTemp <- which(metawithNAs[,"Image_Metadata_CPD_WELL_POSITION"]==activeLegend[kk,1])
        metawithNAs[wellTemp,"Image_Metadata_CPD_MMOL_CONC"] <- activeLegend[kk,7]   # Concentration
        metawithNAs[wellTemp,"Image_Metadata_SOURCE_COMPOUND_NAME"] <- activeLegend[kk,4]   # Treatment name
    }

    metanoNAs <- metawithNAs
    return(metanoNAs)
}

# Quick little function to repair the missing metadata from present fraction legends
unknownsInjectMetadata <- function(metawithNAs) {
    # These are the metadata columns:
    # Image_Metadata_CPD_MMOL_CONC
    # Image_Metadata_CPD_PLATE_MAP_NAME
    # Image_Metadata_CPD_SMILES
    # Image_Metadata_CPD_WELL_POSITION
    # Image_Metadata_MAC_ID
    # Image_Metadata_PlateID
    # Image_Metadata_SOURCE_NAME
    # Image_Metadata_Site

    # > head(legend_working)
    #   Well Row Column      Code Fraction   Treatment Concentration Units
    # 1  A01   A     01      DMSO                 DMSO           0.5     %
    # 2  A02   A     02      DMSO                 DMSO           0.5     %
    # 3  A03   A     03 KCCC-1071        C KCCC-1071_C           0.5     %
    # 4  A04   A     04 KCCC-1088        C KCCC-1088_C           0.5     %
    # 5  A05   A     05 KCCC-1071        1 KCCC-1071_1           0.5     %
    # 6  A06   A     06 KCCC-1088        1 KCCC-1088_1           0.5     %
    # if(any(ls()=="legend_working")) {
        for(kk in 1:nrow(legend_working)) {
            wellTemp <- which(metawithNAs[,"Image_Metadata_CPD_WELL_POSITION"]==legend_working[kk,1])
            metawithNAs[wellTemp,"Image_Metadata_CPD_MMOL_CONC"] <- legend_working[kk,7]   # Concentration
            metawithNAs[wellTemp,"Image_Metadata_CPD_PLATE_MAP_NAME"] <- legend_working[kk,4]   # Plate name
            metawithNAs[wellTemp,"Image_Metadata_PlateID"] <- legend_working[kk,4]   # Plate name
            metawithNAs[wellTemp,"Image_Metadata_SOURCE_COMPOUND_NAME"] <- legend_working[kk,6]   # Treatment name
        }
    # } else {
    #         newname <- paste(metawithNAs[,"Image_Metadata_PlateID"],metawithNAs[,"Image_Metadata_CPD_WELL_POSITION"],sep="-")
    #         metawithNAs[,"Image_Metadata_CPD_PLATE_MAP_NAME"] <- newname   # Plate name
    #         metawithNAs[,"Image_Metadata_SOURCE_COMPOUND_NAME"] <- newname   # Treatment name        
    # }
    metanoNAs <- metawithNAs
    return(metanoNAs)
}


