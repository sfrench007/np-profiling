# Initial option declarations
options(stringsAsFactors=FALSE) # Otherwise we need to force them as strings repeatedly
options(echo=FALSE) # Rscript needs options(echo=TRUE) to make an output file
options(cli.progress_show_after=0) # For progress bars
options(cli.condition="always") # For progress bars
options(warn=(-1)) # Ignore NA warnings
arguments <- commandArgs(trailingOnly=FALSE) # Will capture all arguments, so can search these later

library(cli)
library(DBI)    # Database interface
library(duckdb) # Embedded analytical database
library(openxlsx)

# Header line
cli_h1("Processing compounds and building reference set")

# Scan plates for SD and mean
plateScan <- function(dataIn=cp_features_1pass,plateNames=allplates) {
    unique_plates <- unique(plateNames)
    cli_progress_bar("Scanning plates", total=length(unique_plates), clear=FALSE)
    for(i in 1:length(unique_plates)) {
        cli_progress_update()
        wt <- which(plateNames==unique_plates[i])
        wt_out <- sum(apply(dataIn[wt,],2,function(x) sd(x,na.rm=TRUE)),na.rm=TRUE)
        wt_out2 <- sum(apply(dataIn[wt,],2,function(x) mean(x,na.rm=TRUE)),na.rm=TRUE)
        eo <- c(unique_plates[i],wt_out,wt_out2)
        if(i==1) sd_out <- eo
        else sd_out <- rbind(sd_out,eo)
    }
    cli_progress_done()
    colnames(sd_out) <- c("plate","sd","mean")
    return(sd_out)
}

# Get the peak of the density distribution
denPeak <- function(dataIn,adjval=0.7) {
    denT <- density(dataIn,adjust=adjval,na.rm=TRUE)
    denP <- denT$x[which(denT$y==max(denT$y,na.rm=TRUE))]
    return(denP)
}

# Try to pick the more targeted annotation
nameScan <- function(setIn) {
    setIn <- na.omit(setIn)
    # Check for non-C,D starting letters
    firstLetter <- substr(setIn,1,1)
    # Get length of code
    justCodes <- sapply(setIn,function(x) strsplit(x," - ")[[1]][1])
    codeLength <- nchar(justCodes)
    # Check if it's SARS
    isSARS <- substr(justCodes,1,4)
    # Check if the last 2 characters are letters
    lastTwo <- substr(justCodes,nchar(justCodes)-1,nchar(justCodes))
    # Choose the code ending in 2 letters if available
    if(any(is.na(as.numeric(lastTwo)))) {
        useCode <- setIn[which(is.na(as.numeric(lastTwo)) & !is.na(lastTwo))]
    # Otherwise choose a non-SARS,C,D alpha first letter
    } else if(any(firstLetter!="C" & firstLetter!="D" & isSARS!="SARS")) {
        useCode <- setIn[which(firstLetter!="C" & firstLetter!="D" & isSARS!="SARS")[1]]
    # Otherwise use the longest "D" annotation
    } else if(any(firstLetter=="D")) {
        useCode <- names(which.max(codeLength[which(firstLetter=="D")]))
    # Otherwise use the longest C annotation
    } else if(any(firstLetter=="C")) {
        useCode <- names(which.max(codeLength[which(firstLetter=="C")]))
    # Otherwise just use the last one
    } else {
        useCode <- setIn[length(setIn)]
    }
    return(useCode)
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

# let's set --wd as the working directory, defaulting at the directory that it's run in
if(any(arguments=="--wd")) {
    arg_ind <- which(arguments=="--wd")
    trailing_arg <- arguments[arg_ind[1]+1]
    if(length(arg_ind)>=1) {
        invisible(tryCatch(setwd(trailing_arg),error=function(e) e, finally=function(x) setwd(trailing_arg)))
    }
    cli_alert_success(" Set working directory to {trailing_arg}")
}

# Load functions from functions.R to save space here
source("functions/functions.R")

required_dirs <- c(
    "data",
    "data/db",
    "data/output",
    "logs",
    "logs/images",
    "functions",
    "drugbank"
)
for (d in required_dirs) {
    if (!dir.exists(d)) {
        cli_alert_warning("Directory {d} does not exist.")
    }
}

# Load in the compiled data from DuckDB, applying the deadwells_filesize quality
# filter written by script 02 (equivalent to the old noDBA Rds).
db_path <- "data/db/cellpainting.duckdb"
con <- dbConnect(duckdb(), dbdir=db_path, read_only=TRUE)
avail_tbls <- dbListTables(con)

# Check if script 03 has already been applied to this database.
# Pass --force on the command line to bypass and re-run anyway.
already_processed <- FALSE
if ("_processing_log" %in% avail_tbls) {
    proc_log <- dbReadTable(con, "_processing_log")
    if (any(proc_log$script == "03.create_cc_comp_fractions")) {
        already_processed <- TRUE
        prev_ts <- tail(proc_log$timestamp[proc_log$script == "03.create_cc_comp_fractions"], 1)
        cli_alert_warning(" Database already processed by script 03 on {prev_ts}")
    }
}
if (already_processed && !any(arguments == "--force")) {
    dbDisconnect(con, shutdown=TRUE)
    cli_alert_info(" Skipping. Use --force to re-run script 03 on this database.")
    quit(save="no", status=0)
}

tbls_to_load <- c("metadata", "mesh", "db", "counts", "sd", "url", "method",
                  "features_raw", "features_1pass_1", "features_1pass_2", "features_1pass_3")
cli_progress_bar("Loading database tables", total=length(tbls_to_load), clear=FALSE)
cp_compiled <- list()
for (tbl in tbls_to_load) {
    if (tbl %in% avail_tbls) cp_compiled[[tbl]] <- as.matrix(dbReadTable(con, tbl))
    cli_progress_update()
}
if ("sample_type" %in% avail_tbls) {
    st <- dbReadTable(con, "sample_type")
    cp_compiled[["type"]] <- st$sample_type
}
if ("_tag" %in% avail_tbls) {
    tag_df <- dbReadTable(con, "_tag")
    cp_compiled[["tag"]] <- setNames(tag_df$value, tag_df$key)
}
if ("deadwells_filesize" %in% avail_tbls) {
    cp_compiled[["deadwells_filesize"]] <- as.matrix(dbReadTable(con, "deadwells_filesize"))
}
dbDisconnect(con, shutdown=TRUE)
cli_alert_success(" All data loaded, proceeding to compile")

# Repair the CP_Index for some of the tables
cp_compiled$sd[,1] <- cp_compiled$metadata[,1]
cp_compiled$url[,1] <- cp_compiled$metadata[,1]
cp_compiled$method[,1] <- cp_compiled$metadata[,1]

# Convert everything to numeric where appropriate and bind CP_Index to names/rownames
for(i in seq_along(cp_compiled)) {
    nm  <- names(cp_compiled)[i]
    el  <- cp_compiled[[i]]
    if(nm == "type") {
        # Named vector: assign CP_Index as element names
        names(cp_compiled[[i]]) <- cp_compiled$metadata[,1]
    } else if(is.matrix(el) && nrow(el) == nrow(cp_compiled$metadata)) {
        rownames(cp_compiled[[i]]) <- cp_compiled$metadata[,1]
        # Tables that keep all columns intact (no CP_Index strip, no numeric coercion)
        keep_as_is <- c("metadata","mesh","db","url","method","features_raw","deadwells_filesize")
        if(!nm %in% keep_as_is) {
            cp_compiled[[i]] <- cp_compiled[[i]][,-1]
            cp_compiled[[i]] <- apply(cp_compiled[[i]],2,as.numeric)
        }
    }
    # tag (short named vector), deadwells_filesize (different row count),
    # and any other non-conforming elements are silently skipped
}

# Everything is set up, now check for bad plates (outrageous SD) that throw things off
allplates <- cp_compiled$metadata[,"Image_Metadata_ArbPlate"]
colnum <- 1110 # Arbitrary feature, but specific to the cutoff chosen
cli_alert_info(" Performing data QC and redundancy checks")

# Find out which have iqms < 1 (they shouldn't be < 1)...something is up
unq_nct <- unique(cp_compiled$metadata[,12])
crange <- 0.1
badones <- c()
for(i in 1:length(unq_nct)) {
    wplate <- which(cp_compiled$metadata[,12]==unq_nct[i])
    if(length(wplate)>3) {
        wval <- iqm(cp_compiled$features_1pass_3[wplate,colnum])
        if(wval<(-crange) | wval>(crange)) {
            badones <- c(badones,i)
        }
    }
}
wonkysamples <- unique(badones)

# Reconstruct good dataset
temp <- fixNA_zeroSum(dataIn=cp_compiled$features_1pass_3,removeZeroVar=TRUE); rm(temp)
redundant <- getRedundant(dataIn=cp_compiled$features_1pass_3,r_cutoff=0.95)

# Archive metadata for rows that will be removed, before removing them
rows_to_remove <- if(length(badones)==0) unique(zsRows) else unique(c(wonkysamples,zsRows))
removed_rows_metadata <- cp_compiled$metadata[rows_to_remove, , drop=FALSE]
removed_rows_metadata <- cbind(
    removed_rows_metadata,
    removal_reason = ifelse(
        rows_to_remove %in% wonkysamples & rows_to_remove %in% zsRows,
        "wonky_and_zsRows",
        ifelse(rows_to_remove %in% wonkysamples, "wonky_plate", "zsRows")
    )
)

# Remove redundant features
if(length(badones)==0) {
    cp_purged <- removeRows(cp_compiled,unique(zsRows))
} else {
    cp_purged <- removeRows(cp_compiled,unique(c(wonkysamples,zsRows)))
}
# features_raw has the same rows removed as all other tables, but no columns are removed

# Check for features that are obscenely high, or have very very low SD, but are not technically NZV
removeThese <- unique(c(
    which(apply(cp_compiled$features_1pass_3,2,max)>130),
    which(apply(cp_compiled$features_1pass_3,2,sd)<0.001)
))

cp_purged$features_1pass_1 <- cp_purged$features_1pass_1[,-unique(c(redundant,removeThese))]
cp_purged$features_1pass_2 <- cp_purged$features_1pass_2[,-unique(c(redundant,removeThese))]
cp_purged$features_1pass_3 <- cp_purged$features_1pass_3[,-unique(c(redundant,removeThese))]
# Final NA/inf removal
finalNA <- unique(which(is.na(cp_purged$features_1pass_3),arr.ind=TRUE)[,2])
if(length(finalNA)>0) {
    cp_purged$features_1pass_1 <- cp_purged$features_1pass_1[,-finalNA]
    cp_purged$features_1pass_2 <- cp_purged$features_1pass_2[,-finalNA]
    cp_purged$features_1pass_3 <- cp_purged$features_1pass_3[,-finalNA]
}
# Manual checking here, check the loadings and make sure there are no features pulling things in crazy directions
getInactive(dataIn=scale(cp_purged$features_1pass_3),sd_cutoff=3,plot=FALSE)
wonkyFeatures <- which(prtemp$loadings[,2]<(-0.1))
if(length(wonkyFeatures)!=0) {
    cp_purged$features_1pass_1 <- cp_purged$features_1pass_1[,-wonkyFeatures]
    cp_purged$features_1pass_2 <- cp_purged$features_1pass_2[,-wonkyFeatures]
    cp_purged$features_1pass_3 <- cp_purged$features_1pass_3[,-wonkyFeatures]
}
# Log the removed features
cp_purged$removed_features <- colnames(cp_purged$features_1pass_3)[unique(sort(c(redundant,removeThese,finalNA,wonkyFeatures)))]
# Archive metadata
cp_purged$removed_samples <- removed_rows_metadata
cli_alert_success(" Operations complete, writing updated data to tables")

# Check where activities lie, we want controls to be active, so pick a cutoff that makes sense
options(echo=TRUE)
png(filename="logs/images/phenotypic_activity_mahalanobis.png",res=600,width=10,height=7,bg="transparent",units="in")
inactive <- getInactive(dataIn=scale(cp_purged$features_1pass_3),sd_cutoff=3,plot=TRUE,quantiles=0.3)
dev.off()
options(echo=FALSE)
cp_purged$mahalanobis <- vals
cp_purged$phenotypic_active <- rep("active",length(cp_purged$mahalanobis))
cp_purged$phenotypic_active[inactive] <- "inactive"
names(cp_purged$phenotypic_active) <- rownames(cp_purged$metadata)

# If any updated annotations exist, add them in here
remaps <- openxlsx::read.xlsx("legends/drugbank_curated_annotations.xlsx",sheet="targets-new",colNames=TRUE)
cli_progress_bar("Remapping targets", total=nrow(cp_purged$db), clear=FALSE)
functional_class <- c()
for(i in 1:nrow(cp_purged$db)) {
    targs <- strsplit(cp_purged$db[i,"name"],";")[[1]]
    targs <- gsub("; ",";",targs)
    if(is.na(targs[1])) {
        functional_class[i] <- NA
    } else if(targs[1]=="") {
        functional_class[i] <- NA
    } else {
        targs_remapped <- c()
        for(j in 1:length(targs)) {
            if(targs[j] %in% remaps$Target_Name) {
                targs_remapped[j] <- remaps$Functional_Class[which(remaps$Target_Name==targs[j])]
            }
        }
        functional_class[i] <-paste(targs_remapped,collapse=";")
    }
    cli_progress_update()
}
cp_purged$db <- cbind(cp_purged$db,functional_class)

# Export master metadata flat file (optional, triggered by --masterlist flag)
if(any(arguments=="--masterlist")) {
    options(echo=TRUE)
    write.table(cbind(cp_purged$metadata,cp_purged$mesh,cp_purged$db,cp_purged$counts,cp_purged$url,cp_purged$method,cp_purged$mahalanobis,cp_purged$phenotypic_active,cp_purged$deadwells_filesize),"data/output/masterdata.tsv",row.names=FALSE,col.names=TRUE,sep="\t")
    options(echo=FALSE)
    cli_alert_success(" Master metadata log stored at [data/output/masterdata.tsv]")
}
cli_alert_success(" All data compiled and quality-checked, ready for downstream analysis")

# Write cp_purged back to DuckDB 
cli_alert_info(" Writing data to database")
options(echo=TRUE)
con_write <- dbConnect(duckdb(), dbdir=db_path, read_only=FALSE)

# Matrix-backed tables: overwrite with purged versions
matrix_tbls <- c("metadata","mesh","db","counts","sd","url","method",
                 "features_raw","features_1pass_1","features_1pass_2","features_1pass_3")
cli_progress_bar("Writing database tables", total=length(matrix_tbls), clear=FALSE)
for (tbl in matrix_tbls) {
    if (!is.null(cp_purged[[tbl]]))
        dbWriteTable(con_write, tbl, as.data.frame(cp_purged[[tbl]]), overwrite=TRUE)
    cli_progress_update()
}
cli_progress_done()

# sample_type vector → sample_type table
if (!is.null(cp_purged$type))
    dbWriteTable(con_write, "sample_type",
                 data.frame(sample_type=cp_purged$type), overwrite=TRUE)

# Named tag vector → _tag table (key / value)
if (!is.null(cp_purged$tag))
    dbWriteTable(con_write, "_tag",
                 data.frame(key=names(cp_purged$tag), value=cp_purged$tag,
                            row.names=NULL), overwrite=TRUE)

# Mahalanobis distances
if (!is.null(cp_purged$mahalanobis))
    dbWriteTable(con_write, "mahalanobis",
                 data.frame(CP_Index=rownames(cp_purged$metadata),
                            mahalanobis=as.numeric(cp_purged$mahalanobis),
                            row.names=NULL), overwrite=TRUE)

# Phenotypic activity calls
if (!is.null(cp_purged$phenotypic_active))
    dbWriteTable(con_write, "phenotypic_active",
                 data.frame(CP_Index=names(cp_purged$phenotypic_active),
                            phenotypic_active=cp_purged$phenotypic_active,
                            row.names=NULL), overwrite=TRUE)

# Removed features archive
if (!is.null(cp_purged$removed_features))
    dbWriteTable(con_write, "removed_features",
                 data.frame(feature=cp_purged$removed_features), overwrite=TRUE)

# Removed samples metadata archive
if (!is.null(cp_purged$removed_samples))
    dbWriteTable(con_write, "removed_samples",
                 as.data.frame(cp_purged$removed_samples), overwrite=TRUE)

# Log this run so re-runs are skipped unless --force is passed
dbWriteTable(con_write, "_processing_log",
             data.frame(
                 script    = "03.create_cc_comp_fractions",
                 timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                 row.names = NULL
             ),
             append=TRUE, overwrite=FALSE)

dbDisconnect(con_write, shutdown=TRUE)
options(echo=FALSE)
cli_alert_success(" Database updated!")

