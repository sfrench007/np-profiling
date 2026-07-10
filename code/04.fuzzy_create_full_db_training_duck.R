# Initial option declarations
options(stringsAsFactors=FALSE) # Otherwise we need to force them as strings repeatedly
options(echo=FALSE) # Rscript needs options(echo=TRUE) to make an output file
options(cli.progress_show_after=0) # For progress bars
options(cli.condition="always") # For progress bars
options(warn=(-1)) # Ignore NA warnings
arguments <- commandArgs(trailingOnly=FALSE) # Will capture all arguments, so can search these later

library(stringi)
library(stringr)
library(duckdb)
library(cli)
library(DBI)

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

# let's set --wd as the working directory, defaulting at the directory that it's run in
if(any(arguments=="--wd")) {
    arg_ind <- which(arguments=="--wd")
    trailing_arg <- arguments[arg_ind[1]+1]
    if(length(arg_ind)>=1) {
        invisible(tryCatch(setwd(trailing_arg),error=function(e) e, finally=function(x) setwd(trailing_arg)))
    }
    cli_alert_success(" Set working directory to {trailing_arg}")
}

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

source("functions/functions.R")

# And drugbank rip
db_all <- readRDS("drugbank/db_all.Rds")

# Subset the PubChem CIDs to connect annotation databases
db_cid <- db_all$drugs$external_identifiers[which(db_all$drugs$external_identifiers[,1]=="PubChem Compound"),]

# Connect to the Cell Painting database with DuckDB
db_path <- "data/db/cellpainting.duckdb"
con <- dbConnect(duckdb(), dbdir=db_path, read_only=TRUE)
avail_tbls <- dbListTables(con)

# Check if script 04 has already been applied to this database.
# Pass --force on the command line to bypass and re-run anyway.
already_processed <- FALSE
if ("_processing_log" %in% avail_tbls) {
    proc_log <- dbReadTable(con, "_processing_log")
    if (any(proc_log$script == "04.fuzzy_create_full_db_training")) {
        already_processed <- TRUE
        prev_ts <- tail(proc_log$timestamp[proc_log$script == "04.fuzzy_create_full_db_training"], 1)
        cli_alert_warning(" Database already processed by script 04 on {prev_ts}")
    }
}
if (already_processed && !any(arguments == "--force")) {
    dbDisconnect(con, shutdown=TRUE)
    cli_alert_info(" Skipping. Use --force to re-run script 04 on this database.")
    quit(save="no", status=0)
}

# Convert the other cp_train variables from database instead!  
# - Subset only the actives though - but using the metadata from this script as the ground truth!
tbls_to_load <- c("metadata", "mesh", "db", "features_1pass_3", "sample_type", "phenotypic_active")
cli_progress_bar("Loading database tables", total=length(tbls_to_load), clear=FALSE)
cp_compiled <- list()
for (tbl in tbls_to_load) {
    if (tbl %in% avail_tbls) cp_compiled[[tbl]] <- as.matrix(dbReadTable(con, tbl))
    cli_progress_update()
}
dbDisconnect(con, shutdown=TRUE)
cp_train <- subsetRows(cp_compiled, rows=which(cp_compiled$sample_type=="training" & cp_compiled$phenotypic_active[,2]=="active"))
cp_train <- global_rownames(cp_train)
cli_alert_success(" Data loaded!")

# Go through and get the pubmed CIDs, and match them to the drugbank ones
cli_alert_info(" Initial CID matching")
dbmatch <- c()
cli_progress_bar("Matching PubChem CIDs", total=nrow(cp_train$mesh), clear=FALSE)
for(i in 1:nrow(cp_train$mesh)) {
    cli_progress_update()
    dbt <- which(db_cid[,2]==as.numeric(cp_train$mesh[i,"pubchem_cid"]))
    if(length(dbt)==0) {
        dbt <- which(db_cid[,2]==cp_train$mesh[i,16])
        if(length(dbt)==0) {
            # If the CID and the SMILES match CID are not applicable, then it's not in the current DrugBank rip
            dbmatch[i] <- NA
        } else {
            dbmatch[i] <- dbt[1]
        }      
    } else {
        dbmatch[i] <- dbt[1]
    }
}
cli_progress_done()

# Target exists in db_all$cett$targets$general_information as "name" column
cli_alert_info(" Matching with DrugBank targets")
tar_info <- unlist(db_all$cett$targets$general_information[,6])
tar_match <- list()

# Updated version
cli_progress_bar("Matching DrugBank targets", total=nrow(cp_train$mesh), clear=FALSE)
for(i in 1:nrow(cp_train$mesh)) {
    cli_progress_update()
    if(!is.na(dbmatch[i])) {
        tar_t <- which(tar_info==unlist(db_cid[dbmatch[i],3]))
        if(length(tar_t)!=0) {
            tar_match[[i]] <- db_all$cett$targets$general_information[tar_t,1:2]
        } else {
            tar_match[[i]] <- NA
        }
    } else {
        # We may have manually annotated this though, so make sure manual annotations are ok
        if(!is.na(cp_train$db[i,"name"])) {
            # This means there is a drugbank annotation added manually, so use this!
            t0 <- which(grepl(cp_train$db[i,"name"],unlist(db_all$cett$targets$general_information[,2])))[1]
            tar_match[[i]] <- db_all$cett$targets$general_information[t0,1:2]
        } else {
            # Otherwise it's an NA all around
            tar_match[[i]] <- NA
        }
    }
    # Update functional class if --useremaps is indicated
}
cli_progress_done()

# More accurate ATC codes?
cli_alert_info(" Repairing ATC codes")
atc_info <- unlist(db_all$drugs$atc_codes[,10])
atc_match <- list()
cli_progress_bar("Matching ATC codes", total=nrow(cp_train$mesh), clear=FALSE)
for(i in 1:nrow(cp_train$mesh)) {
    cli_progress_update()
    if(!is.na(dbmatch[i])) {
        atc_t <- which(atc_info==unlist(db_cid[dbmatch[i],3]))
        if(length(atc_t)!=0) {
            atc_match[[i]] <- db_all$drugs$atc_codes[atc_t,1:9]
        } else {
            atc_match[[i]] <- NA
        }
    } else {
        atc_match[[i]] <- NA
    }
}
cli_progress_done()

# Create unique training sets for target and atc levels 1-4
full_db_training <- list()
tag <<- generateTag(hqp="Anonymous")
full_db_training$tag <- tag
cli_alert_info(" Creating training sets")
cli_progress_bar("Building training sets", total=nrow(cp_train$mesh)*5, clear=FALSE)
for(i in 1:5) {
    counter <- 1
    if(i<5) { # If it's *not* target
        full_db_training[[paste0("level_",i)]] <- i
        for(ii in 1:length(atc_match)) {
            if(any(is.null(atc_match[[ii]]))) atc_match[[ii]] <- NA 
            if(is.na(as.matrix(atc_match[[ii]])[1,1])) {
            # ii is an index number, don't change it, it's important
                wline <- ii 
            } else {
                wline <- cbind(ii,unique(atc_match[[ii]][,seq(2,8,by=2)[i]]))
            }
            # Fix if it's empty
            if(length(wline)>1) {
                wmed <- matrix(data=cp_train$features_1pass_3[rep(ii,nrow(wline)),],nrow=nrow(wline),ncol=ncol(cp_train$features_1pass_3))
                rownames(wmed) <- rownames(cp_train$metadata)[rep(ii,nrow(wline))]
                colnames(wmed) <- colnames(cp_train$features_1pass_3)
                if(counter==1) {
                    counter <- counter + 1
                    full_db_training[[paste0("level_",i)]]$db <- wline
                    full_db_training[[paste0("level_",i)]]$median <- wmed
                } else {
                    full_db_training[[paste0("level_",i)]]$db <- rbind(full_db_training[[paste0("level_",i)]]$db,wline)
                    full_db_training[[paste0("level_",i)]]$median <- rbind(full_db_training[[paste0("level_",i)]]$median,wmed)
                }
            } 
            cli_progress_update()
        }
    } else { # Otherwise it's target
        tartemp <- c()
        full_db_training[["target"]] <- i
        for(ii in 1:length(tar_match)) {
            if(any(is.null(tar_match[[ii]]) | length(tar_match[[ii]])==0)) tar_match[[ii]] <- NA 
            if(is.na(as.matrix(tar_match[[ii]])[1,1])) {
                wline <- ii
            } else {
                wline <- cbind(ii,unique(tar_match[[ii]][,2]))
                if(length(wline[,2])>8) {
                    wline <- wline[1:8,]
                    # print(wline) # Just for QC purposes
                }
                tartemp[ii] <- paste(wline[,2],collapse=";")
            }
            
            if(length(wline)>1) {
                wmed <- matrix(data=cp_train$features_1pass_3[rep(ii,nrow(wline)),],nrow=nrow(wline),ncol=ncol(cp_train$features_1pass_3))
                rownames(wmed) <- rownames(cp_train$metadata)[rep(ii,nrow(wline))]
                colnames(wmed) <- colnames(cp_train$features_1pass_3)
                if(counter==1) {
                    counter <- counter + 1
                    full_db_training[["target"]]$db <- wline
                    full_db_training[["target"]]$median <- wmed
                } else {
                    full_db_training[["target"]]$db <- rbind(full_db_training[["target"]]$db,wline)
                    full_db_training[["target"]]$median <- rbind(full_db_training[["target"]]$median,wmed)
                }
            } 
            cli_progress_update()
        }
        full_db_training[["target"]]$target_mesh_match <- tar_match
        full_db_training[["target"]]$collapsed <- tartemp
    }
    gc(verbose=FALSE)
}
cli_progress_done()

# Save the list to speed things up
options(echo=TRUE)
save(full_db_training,file="data/output/full_db_training_scale_dba_remap.Rdata",version=2)
options(echo=FALSE)

# Log this run so re-runs are skipped unless --force is passed
con_log <- dbConnect(duckdb(), dbdir=db_path, read_only=FALSE)
dbWriteTable(con_log, "_processing_log",
             data.frame(
                 script    = "04.fuzzy_create_full_db_training",
                 timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                 row.names = NULL
             ),
             append=TRUE, overwrite=FALSE)
dbDisconnect(con_log, shutdown=TRUE)
cli_alert_success(" Processing log updated")

