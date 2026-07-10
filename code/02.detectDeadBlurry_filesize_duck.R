# Initial option declarations
options(stringsAsFactors=FALSE) # Otherwise we need to force them as strings repeatedly
options(echo=FALSE) # Rscript needs options(echo=TRUE) to make an output file
options(cli.progress_show_after=0) # For progress bars
options(cli.condition="always") # For progress bars
options(warn=(-1)) # Ignore NA warnings
arguments <- commandArgs(trailingOnly=FALSE) # Will capture all arguments, so can search these later

library(cli)     # For progress bars with an ETA
library(DBI)     # Database interface
library(duckdb)  # Embedded analytical database
library(rootSolve)
library(stringi)
library(stringr)
library(openxlsx)

# Function to update lists by subsetting the same number of members from each matrix/vector
subsetRows <- function(listIn,subsetThese) {
    listOut <- listIn
    for(iii in 1:length(listIn)) {
        if(names(listIn)[iii]=="tag") {
            listOut[[iii]] <- listIn[[iii]]
        } else {
            if(!is.null(nrow(listIn[[iii]]))) {
                listOut[[iii]] <- listIn[[iii]][subsetThese,]
            } else {
                listOut[[iii]] <- listIn[[iii]][subsetThese]
            }
        }
    }
    return(listOut)
}

findBlank <- function(codeIn,dead_cutoff=0.1) {
    dummy <- c(NA,NA); names(dummy) <- c("Blank.Status","Blank.Cell.Density")
    if(any(as.numeric(substr(codeIn,nchar(codeIn)-1,nchar(codeIn))) %in% c(1:4,21:24))) {
        return(dummy)
    } 
    if(any(substr(codeIn,nchar(codeIn)-2,nchar(codeIn)-2) %in% c("O","P"))) {
        return(dummy)
    } 
    if(any(substr(codeIn,4,4)=="F")) {
        return(dummy)
    } 
    # Get the blank and check if it has low cell density
    w_wells <- both_data[which(fplates==substr(codeIn,1,17)),]
    codecol <- substr(codeIn,nchar(codeIn)-1,nchar(codeIn))
    platecol <- substr(w_wells[,1],nchar(w_wells[,1])-1,nchar(w_wells[,1]))
    fromcol <- w_wells[which(platecol==codecol),]
    blankwell <- dummy
    if(nrow(fromcol)<5) {
        return(dummy)
    }
    if(substr(codeIn,nchar(codeIn)-2,nchar(codeIn)-2) %in% LETTERS[seq(1,16,by=2)]) {
        blankwell <- fromcol[which(substr(fromcol[,1],nchar(fromcol[,1])-2,nchar(fromcol[,1])-2)=="O"),3]
        if(length(blankwell)==0) return(dummy)
        blankOut <- ifelse(blankwell<=dead_cutoff,"dead","alive")
        codeOut <- c(blankOut,blankwell); names(codeOut) <- c("Blank.Status","Blank.Cell.Density")
    } else {
        blankwell <- fromcol[which(substr(fromcol[,1],nchar(fromcol[,1])-2,nchar(fromcol[,1])-2)=="P"),3]
        if(length(blankwell)==0) return(dummy)
        blankOut <- ifelse(blankwell<=dead_cutoff,"dead","alive")
        codeOut <- c(blankOut,blankwell); names(codeOut) <- c("Blank.Status","Blank.Cell.Density")
    }
    return(codeOut)
}

# Then, match the bad images to their actual data so we can remove them
# "25.01.20_KCDF_250nL.015.R3_Int_R" is the format of the bad variable plates
# Ignore the date, take the KCDX format, the plate, and the rep.  Make a metadata-friendly plate
convertName <- function(nameIn) {
    tn <- nameIn
    tempsplit <- strsplit(tn,"/")[[1]]
    tempsplit <- tempsplit[which(nchar(tempsplit)==max(nchar(tempsplit)))]
    tempsplit <- strsplit(tempsplit,"_")[[1]]

    # In case it's a KCDX plate which has odd filename conventions
    if(substr(tn,nchar(tn)-1,nchar(tn)-1)=="_") {
        # This will be library, concentration, plate, and rep
        nameOut <- paste0(tempsplit[2],"-",tempsplit[3])
    } else {
        # This will be library, concentration, plate, and rep
        nameOut <- paste0(tempsplit[1],"-",tempsplit[3],"-",tempsplit[2],"-",tempsplit[4])
    }

    return(nameOut)
}

# Density normalization using maxima
minmaxNorm <- function(listItemIn=filesizes_raw[[119]],type="dapi") {
    # Can be "merged","dapi","gfp","txred","cy5"
    listItemIn <- listItemIn[[type]]

    dens <- density(listItemIn)
    minmax <- uniroot.all(approxfun(dens$x[-1],diff(dens$y)),interval=range(dens$x[-1]))
    if(length(minmax)==1) {
        dens <- density(listItemIn,adjust=0.75)
        minmax <- uniroot.all(approxfun(dens$x[-1],diff(dens$y)),interval=range(dens$x[-1]))
    }
    if(length(minmax)==1) {
        minmax <- c((min(listItemIn)+minmax)/2.125,minmax)
        # Convert to actual x vals
        maxima <- unlist(sapply(minmax,function(xx) { 
            closest_t <- (abs(xx-dens$x))+xx
            which(closest_t==min(closest_t))
        }))
    } else {
        # Convert to actual x vals
        minmax_x <- unlist(sapply(minmax,function(xx) { 
            closest_t <- (abs(xx-dens$x))+xx
            which(closest_t==min(closest_t))
        }))
        # which are the minima and maxima?
        minima <- c(); maxima <- c()
        for(ii in 1:length(minmax_x)) {
            v0 <- dens$y[minmax_x[ii]]
            v1 <- dens$y[minmax_x[ii]+1]
            if(v1>v0) {
                minima <- c(minima,minmax_x[ii])
            } else {
                maxima <- c(maxima,minmax_x[ii])
            }
        }
    }

    # Of the maxima, what are the phenotypic maxima?
    if(length(maxima)>2) {
        # If there are more than 2 maxima, just take the top 2
        for(q in 1:(length(maxima)-2)) {
           maxima <- maxima[-which(dens$y[maxima]==min(dens$y[maxima]))[1]]
        }
        listItemOut <- (listItemIn-dens$x[maxima[1]])/(dens$x[maxima[2]]-dens$x[maxima[1]])
    } else if(length(maxima)==1) {
        # If there is just the 1 main (normal cells) maximum, just set that to 1
        listItemOut <- (listItemIn)/(dens$x[maxima])
    } else {
        listItemOut <- (listItemIn-dens$x[maxima[1]])/(dens$x[maxima[2]]-dens$x[maxima[1]])
    }
    # Then normalize the plate to these peaks
    # listItemOut <- (listItemIn-dens$x[maxima[1]])/(dens$x[maxima[2]]-dens$x[maxima[1]])
    return(listItemOut)
}

# Roll across fractions to find active adjacent fractions, to prioritize ones with active neighbours
rollingCheck <- function(itemNumber) {
    indicesIn <- bad$indices_means[[itemNumber]]
    filesizeIn <- bad$filesizes_means[[itemNumber]]
    prioritize <- c(); prioritize[1] <- "low"
    for(ii in 2:length(indicesIn)) {
        if(indicesIn[ii]==(indicesIn[ii-1]+2)) {
            # Examine both fractions and see which is deader
            if(filesizeIn[ii-1]>filesizeIn[ii]) {
                prioritize[ii] <- "high"
                prioritize[ii-1] <- "med"
            } else {
                prioritize[ii] <- "med" 
                prioritize[ii-1] <- "high"
            }
        } else {
            prioritize[ii] <- "low"            
        }
    }
    return(prioritize)
}

# Header line
cli_h1("Flagging dead/empty wells")

# let's set --wd as the working directory, defaulting at the directory that it's run in
if(any(arguments=="--wd")) {
    arg_ind <- which(arguments=="--wd")
    trailing_arg <- arguments[arg_ind[1]+1]
    if(length(arg_ind)>=1) {
        invisible(tryCatch(setwd(trailing_arg),error=function(e) e, finally=function(x) setwd(trailing_arg)))
    }
    cli_alert_success(" Set working directory to {trailing_arg}")
}

# Ensure required directories exist (create them silently if not present)
required_dirs <- c(
    "code",         # Output: filesizes scan
    "data/output",  # Output: flat-file exports (deadwells CSV, dead_wells TSV)
    "data/db",      # Input/Output: DuckDB files from script 01 and this script
    "data/input",   # Input: filesizes_all_fluorphores.Rds
    "legends",      # Input: KCB_Fractions_Compiled.xlsx
    "functions",    # Source: functions.R, justGetFileSizes.R
    "drugbank"      # Input: drugbank rip
)
for (d in required_dirs) {
    if (!dir.exists(d)) {
        dir.create(d, recursive=TRUE, showWarnings=FALSE)
        cli_alert_info(" Created missing directory: {d}")
    }
}

# Run the 'justGetFileSizes.R' first to generate filesizes
source("functions/functions.R")
cli_alert_success(" Loaded dependencies")

# Read a cp_compiled-style list back from a DuckDB file written by script 01.
# Tables that exist in the db are read and assembled; missing tables are silently skipped.
read_from_duckdb <- function(db_path) {
    cli_alert_info(" Opening DuckDB at [{db_path}]")
    con <- dbConnect(duckdb(), dbdir=db_path, read_only=TRUE)
    on.exit(dbDisconnect(con, shutdown=TRUE), add=TRUE)
    available <- dbListTables(con)
    matrix_tables <- c("metadata", "mesh", "db", "counts", "sd", "url", "method",
                       "features_raw", "features_1pass_1", "features_1pass_2", "features_1pass_3")
    result <- list()
    for (tbl in matrix_tables) {
        if (!tbl %in% available) next
        df <- dbReadTable(con, tbl)
        result[[tbl]] <- as.matrix(df)
    }
    if ("sample_type" %in% available) {
        st <- dbReadTable(con, "sample_type")
        result[["type"]] <- st$sample_type
    }
    if ("_tag" %in% available) {
        tag_df <- dbReadTable(con, "_tag")
        result[["tag"]] <- setNames(tag_df$value, tag_df$key)
    }
    cli_alert_success(" Loaded [{length(result)}] components from DuckDB")
    result
}

# Write a cp_compiled-style list into a DuckDB file (appending to existing tables).
# dataset_label is used as a suffix on the table names so noDBA and justDBA
# can coexist in the same database file.
write_compiled_to_duckdb <- function(compiled_list, db_path, dataset_label, overwrite=FALSE) {
    cli_alert_info(" Writing [{dataset_label}] to DuckDB at [{db_path}]")
    con <- dbConnect(duckdb(), dbdir=db_path, read_only=FALSE)
    on.exit(dbDisconnect(con, shutdown=TRUE), add=TRUE)
    cp_index <- compiled_list$metadata[, "CP_Index"]
    matrix_tables <- c("metadata", "mesh", "db", "counts", "sd", "url", "method",
                       "features_raw", "features_1pass_1", "features_1pass_2", "features_1pass_3")
    for (tbl in matrix_tables) {
        if (!tbl %in% names(compiled_list)) next
        df <- as.data.frame(compiled_list[[tbl]], stringsAsFactors=FALSE)
        if (!"CP_Index" %in% colnames(df)) df <- cbind(CP_Index=cp_index, df, stringsAsFactors=FALSE)
        tbl_name <- paste0(tbl, "__", dataset_label)
        dbWriteTable(con, tbl_name, df, overwrite=overwrite, append=!overwrite)
    }
    if ("type" %in% names(compiled_list)) {
        type_df <- data.frame(CP_Index=cp_index, sample_type=compiled_list$type, stringsAsFactors=FALSE)
        dbWriteTable(con, paste0("sample_type__", dataset_label), type_df, overwrite=overwrite, append=!overwrite)
    }
    if (!is.null(compiled_list$tag)) {
        tag_val <- if (is.list(compiled_list$tag)) unlist(compiled_list$tag) else compiled_list$tag
        tag_df  <- data.frame(key=names(tag_val), value=as.character(tag_val), stringsAsFactors=FALSE)
        dbWriteTable(con, paste0("_tag__", dataset_label), tag_df, overwrite=TRUE)
    }
    cli_alert_success(" [{dataset_label}] written to DuckDB")
}

# Check if script 02 has already been applied to this database.
# Pass --force on the command line to bypass and re-run anyway.
already_processed_02 <- FALSE
db_path_02 <- "data/db/cellpainting.duckdb"
if (file.exists(db_path_02)) {
    con_check_02 <- dbConnect(duckdb(), dbdir=db_path_02, read_only=TRUE)
    avail_tbls_02 <- dbListTables(con_check_02)
    if ("_processing_log" %in% avail_tbls_02) {
        proc_log_02 <- dbReadTable(con_check_02, "_processing_log")
        if (any(proc_log_02$script == "02.detectDeadBlurry_filesize")) {
            already_processed_02 <- TRUE
            prev_ts_02 <- tail(proc_log_02$timestamp[proc_log_02$script == "02.detectDeadBlurry_filesize"], 1)
            cli_alert_warning(" Database already processed by script 02 on {prev_ts_02}")
        }
    }
    dbDisconnect(con_check_02, shutdown=TRUE)
}

if (already_processed_02 && !any(arguments == "--force")) {
    cli_alert_info(" Skipping. Use --force to re-run script 02 on this database.")
    quit(save="no", status=0)
}

useMethod <- "filesize"
if(any(arguments=="--method")) {
    arg_ind <- which(arguments=="--method")
    trailing_arg <- arguments[arg_ind[1]+1]
    if(trailing_arg!="filesize") {
        useMethod <- "transformer"
    }
    cli_alert_success(" Set detection method to [{trailing_arg}]")
} else {
    cli_alert_success(" Set detection method to [filesize]")    
}

if(useMethod=="filesize") {
    if(file.exists("data/input/filesizes_all_fluorphores.Rds")) {
        cli_alert_success(" Existing filesize scan present at [data/input/filesizes_all_fluorphores.Rds]")
        if(any(arguments=="--force")) {
            cli_alert_success(" Forcing re-scan of file sizes")
            source("functions/justGetFileSizes.R")
        }
    } else {
        cli_alert_success(" No existing filesize scan present at [data/input/filesizes_all_fluorphores.Rds], scanning paths")
        source("functions/justGetFileSizes.R")
    }

    # Once screened, time to check it out
    filesizes_raw <- readRDS("data/input/filesizes_all_fluorphores.Rds")

    # Be sure to take the newest versions of the files!!!
    cli_alert_info(" Purging duplicate (older/incomplete) runs")
    repairedNames <- sapply(names(filesizes_raw),function(x) {
        t0 <- strsplit(x,"/")[[1]]
        if(!any(grepl("ontam",t0))) {
            t1 <- strsplit(t0[length(t0)],"_")[[1]]
            library0 <- t1[1]
            if(nchar(library0)>4) {
                library0 <- t1[which(grepl("KC",t1))]
            }
            plate0 <- sprintf("%03d",as.numeric(t1[which(!is.na(as.numeric(t1)))]))
            if(length(plate0)==0) {
                t00 <- strsplit(t1[3],"-")[[1]]
                plate0 <- t00[2]
                conc0 <- t00[1]
                rep0 <- t00[3]
            } else {
                conc0 <- t1[which(!is.na(as.numeric(substr(t1,1,1))) & is.na(as.numeric(substr(t1,nchar(t1),nchar(t1)))) & nchar(t1)>2 & nchar(t1)<8 & substr(t1,nchar(t1)-1,nchar(t1)-1)!="s")]
                if(length(conc0)==0) {
                    conc0 <- t1[which(grepl("serial",t1))]
                }
                rep0 <- t1[which(nchar(t1)==2 & substr(t1,1,1)=="R")]
            }
            paste(library0,plate0,conc0,rep0,sep="-")
        }
    }); repairedNames <- unlist(repairedNames)

    dups <- repairedNames[which(repairedNames %in% repairedNames[which(duplicated(repairedNames))])]
    unq <- unique(dups); removeThese <- c()
    cli_progress_bar("Purging duplicate runs", total=length(unq), clear=FALSE)
    for(i in 1:length(unq)) {
        cli_progress_update()
        newest <- NA
        working0 <- dups[which(dups==unq[i])]
        s0 <- sapply(names(working0),function(x) strsplit(x,"_"))
        s1 <- sapply(s0,function(x) which(grepl("202.-",x)))
        s1repair <- which(sapply(s1,function(x) length(as.numeric(x)))==0)
        if(length(s1repair)>0) {
            s0 <- s0[-s1repair]
            s1 <- s1[-s1repair]
        }
        if(sum(sapply(s1,length))==1) {
            newest <- names(sapply(s1,length))[which(sapply(s1,length)==1)]
            removeThese <- c(removeThese,newest)
        } else {
            dates0 <- c()
            for(k in 1:length(s1)) {
                dates0[k] <- s0[[k]][as.numeric(s1[k])]
            }
            s2 <- sapply(dates0,function(x) strsplit(x,"-")[[1]])
            s3 <- apply(s2,1,max)
            # Check which is newer now
            if(length(which(s2[1,]==s3[1]))>1) {
                if(length(which(s2[2,]==s3[2]))>1) {
                    newest <- s2[,which(s2[3,]==s3[3])]
                } else {
                    newest <- s2[,which(s2[2,]==s3[2])]
                }
            } else {
                newest <- s2[,which(s2[1,]==s3[1])]
            }
            newest <- paste(newest,collapse="-")
            removeThese <- c(removeThese,names(working0)[which(is.na(as.numeric(sapply(s0,function(x) which(x==newest)))))])
        }
    }
    cli_progress_done()
    remove_from_masters <- which(names(repairedNames) %in% removeThese)
    f1 <- length(filesizes_raw)
    f0 <- length(remove_from_masters); 
    filesizes_raw <- filesizes_raw[-remove_from_masters]
    cli_alert_success(" Purged [{f0}] duplicate entries across [{f1}] experiments")

    filesizes_merged <- lapply(filesizes_raw,function(x) minmaxNorm(x,type="merged"))
    filesizes_dapi <- lapply(filesizes_raw,function(x) minmaxNorm(x,type="dapi"))
    filesizes_gfp <- lapply(filesizes_raw,function(x) minmaxNorm(x,type="gfp"))
    filesizes_txred <- lapply(filesizes_raw,function(x) minmaxNorm(x,type="txred"))
    filesizes_cy5 <- lapply(filesizes_raw,function(x) minmaxNorm(x,type="cy5"))

    wells384 <- paste0(sort(rep(LETTERS[1:16],24)),rep(sprintf("%02d",1:24),16))
    wells3456 <- sort(rep(paste0(sort(rep(LETTERS[1:16],24)),rep(sprintf("%02d",1:24),16)),9))

    # Get the means
    filesizes_means_merged <- lapply(filesizes_merged,function(x) {
        sapply(wells384,function(y) {
            mean(x[which(wells3456==y)])
        })
    })
    filesizes_means_dapi <- lapply(filesizes_dapi,function(x) {
        sapply(wells384,function(y) {
            mean(x[which(wells3456==y)])
        })
    })
    filesizes_means_gfp <- lapply(filesizes_gfp,function(x) {
        sapply(wells384,function(y) {
            mean(x[which(wells3456==y)])
        })
    })
    filesizes_means_txred <- lapply(filesizes_txred,function(x) {
        sapply(wells384,function(y) {
            mean(x[which(wells3456==y)])
        })
    })
    filesizes_means_cy5 <- lapply(filesizes_cy5,function(x) {
        sapply(wells384,function(y) {
            mean(x[which(wells3456==y)])
        })
    })

    filesizes <- list()
    filesizes$merged <- filesizes_means_merged
    filesizes$dapi <- filesizes_means_dapi
    filesizes$gfp <- filesizes_means_gfp
    filesizes$txred <- filesizes_means_txred
    filesizes$cy5 <- filesizes_means_cy5

    # Check for permutations of file intensities to see if we can best get a coefficient to provide cytotoxicity data
    # Try first dapi/txred or something like that, see if when proportions are off we can flag these.  penalizedLDA?
    cutoff <- list()
    cutoff$merged <- 0.2
    filesizes_means <- list()
    well_dba_indices <- list(); well_dba_values <- list(); well_dba_names <- list()
    cli_progress_bar("Computing well cutoffs", total=length(filesizes_means_dapi), clear=FALSE)
    for(v in 1:length(filesizes_means_dapi)) {
        cli_progress_update()
        # v <- 393
        # pass1 is the merged cutoff
        filesizes_means[[v]] <- filesizes_means_merged[[v]]
        pass1 <- which(filesizes_means_merged[[v]]<cutoff$merged)
        # pass2 is a dapi check
        cutoff$dapi[v] <- iqm(filesizes_means_dapi[[v]])-(sd(filesizes_means_dapi[[v]])*2.1)
        pass2 <- which(filesizes_means_dapi[[v]]<cutoff$dapi[v])
        # pass3 is a txred check
        cutoff$txred[v] <- iqm(filesizes_means_txred[[v]])-(sd(filesizes_means_txred[[v]])*2.1)
        pass3 <- which(filesizes_means_txred[[v]]<cutoff$txred[v])
        # pass4 is mitotracker
        cutoff$cy5[v] <- iqm(filesizes_means_cy5[[v]])+(sd(filesizes_means_dapi[[v]])*3)
        pass4 <- which(filesizes_means_cy5[[v]]>cutoff$cy5[v] & filesizes_means_merged[[v]]<cutoff$merged)
        well_dba_indices[[v]] <- sort(unique(c(pass4,intersect(pass1,pass2),intersect(pass1,pass3))))
        well_dba_values[[v]] <- filesizes_means_merged[[v]][well_dba_indices[[v]]]
        well_dba_names[[v]] <- wells384[well_dba_indices[[v]]]
    }
    cli_progress_done()
    names(filesizes_means) <- names(filesizes_means_dapi)

    # Convert the bad list name to a common format
    workingList <- filesizes_means
    newNames_workingList <- sapply(names(workingList),convertName)
    names(workingList) <- newNames_workingList

    bad <- list()
    bad$cutoff <- cutoff
    bad$filesizes_raw <- filesizes_raw
    bad$filesizes_norm <- filesizes
    bad$filesizes_means <- filesizes_means
    bad$filesizes_means_codes <- newNames_workingList
    bad$well_dba_indices <- well_dba_indices
    bad$well_dba_values <- well_dba_values
    bad$well_dba_names <- well_dba_names

    if(!any(arguments=="--annotateonly")) {
        if(any(arguments=="--hqp")) {
            arg_ind <- which(arguments=="--hqp")
            arg_trailing <- arguments[arg_ind[1]+1]
            if(length(arg_ind)>=1) {
                tag <<- invisible(tryCatch(generateTag(hqp=arg_trailing),error=function(e) e, finally=function(x) generateTag(hqp=arg_trailing)))
            }
            cli_alert_success(" HQP set as {arg_trailing}")
        }
        # If no hqp exists in the arguments
        if(!any(ls()=="tag")) {
            cli_alert_warning(" No HQP set in arguements so using [Anonymous] to generate info tag")
            tag <<- generateTag(hqp="Anonymous")
        }
        options(echo=TRUE)
        saveRDS(bad,"data/output/dba_summary.Rds",version=2)
        options(echo=FALSE)
        cli_alert_success(" Writing file [data/output/dba_summary.Rds]")

        # Export dead-well summary to CSV for reporting
        allwells <- c(); allvalues <- c(); justallwells <- c()
        for(i in 1:length(bad$filesizes_means)) {
            wplate <- bad$well_dba_names[[i]]   # named: well positions flagged as dead
            if(length(wplate>0)) {
                allwells <- c(allwells,paste(bad$filesizes_means_codes[i],wplate,sep="-"))
                allvalues <- c(allvalues,bad$well_dba_values[[i]])  # named: scaled cell-density values
                justallwells <- c(justallwells,wplate)
            }
        }

        alllibs <- sapply(allwells,function(x) strsplit(x,"-")[[1]][1])
        allplates <- sapply(allwells,function(x) strsplit(x,"-")[[1]][3])
        allreps <- sapply(allwells,function(x) strsplit(x,"-")[[1]][4])
        export_this <- cbind(allwells,alllibs,allplates,allreps,justallwells,allvalues)
        colnames(export_this) <- c("code","library","plate","rep","well","scaled.cell.density")
        options(echo=TRUE)
        write.table(export_this,paste0("data/output/deadwells_all.csv"),sep=",",row.names=FALSE,col.names=TRUE,quote=FALSE)
        options(echo=FALSE)
        fname <- paste0("data/output/deadwells_all.csv")
        cli_alert_success(" Writing file [{fname}]")

        cli_alert_success(" Loading dependencies")
        source("functions/functions.R")
        options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx16192m"))
        cli_alert_info(" Clearing RAM")
        gc()

        # Only the metadata table is needed here to build sample codes and row count.
        # Keep the connection open so we can write the quality annotation back without
        # a second open/close cycle.  The feature matrices are never loaded.
        cli_alert_info(" Connecting to [data/db/cellpainting.duckdb] (metadata only)")
        db_con    <- dbConnect(duckdb(), dbdir="data/db/cellpainting.duckdb", read_only=FALSE)
        cc_meta   <- dbReadTable(db_con, "metadata")  # fast: one small table
        bad <- readRDS("data/output/dba_summary.Rds")

        cc_comp_codes <- apply(cc_meta, 1, function(x) paste(x[7], x[2], paste0("R", x[11]), x[5], sep="-"))

        # Create a second set that converts the old names to the new ones
        cc_comp_codes2 <- cc_comp_codes
        cc_comp_codes2[which(substr(cc_comp_codes2,1,4)=="ANTI")] <- gsub("ANTI","ATCA",cc_comp_codes2[which(substr(cc_comp_codes2,1,4)=="ANTI")])

        # Then swap conc and plate
        cc_comp_codes3 <- sapply(cc_comp_codes2,function(x) {
            st <- strsplit(x,"-")[[1]]
            paste(st[1],st[3],st[2],st[4],st[5],sep="-")
        }); names(cc_comp_codes3) <- NULL

        # Bind the bad wells to their codes
        badcodes <- c(); badnames <- c(); badcombo <- c()
        cli_progress_bar("Building bad-well codes", total=length(bad$well_dba_names), clear=FALSE)
        for(i in 1:length(bad$well_dba_names)) { 
            cli_progress_update()
            tname <- bad$filesizes_means_codes[i]
            badnames_t <- strsplit(names(bad[[2]])[i],"/")[[1]][length(strsplit(names(bad[[2]])[i],"/")[[1]])]
            if(length(bad$well_dba_names[[i]])>0) {
                tcode <- paste0(tname,"-",bad$well_dba_names[[i]])
                badnames <- c(badnames,badnames_t)
                badcodes <- c(badcodes,tcode)
                badcombo <- c(badcombo,paste0(badnames_t,"_",bad$well_dba_names[[i]]))
            }
        }
        cli_progress_done()

        # Then finally match to cc_comp_codes.  These will be the indices to remove
        allbad <- c()
        cli_progress_bar(" Processing ",total=length(badcodes),clear=FALSE)
        for(i in 1:length(badcodes)) {
            cli_progress_update()
            temp <- which(cc_comp_codes3==badcodes[i])
            if(length(temp)>0) {
                allbad <- c(allbad,temp)
            }
        }

        # Build quality vector and write it as a new table in the open connection
        imageQuality <- rep("good", nrow(cc_meta))
        imageQuality[allbad] <- "bad"

        if(any(arguments=="--hqp")) {
            arg_ind <- which(arguments=="--hqp")
            arg_trailing <- arguments[arg_ind[1]+1]
            if(length(arg_ind)>=1) {
                tag <<- invisible(tryCatch(generateTag(hqp=arg_trailing),error=function(e) e, finally=function(x) generateTag(hqp=arg_trailing)))
            }
            cli_alert_success(" HQP set as {arg_trailing}")
        }
        # If no hqp exists in the arguments
        if(!any(ls()=="tag")) {
            cli_alert_warning(" No HQP set in arguments so using [Anonymous] to generate info tag")
            tag <<- generateTag(hqp="Anonymous")
        }

        cli_alert_success(" Adding [deadwells_filesize] table to [data/db/cellpainting.duckdb]")
        options(echo=TRUE)
        dw_df <- data.frame(CP_Index=cc_meta$CP_Index,
                            image_quality=imageQuality,
                            stringsAsFactors=FALSE)
        dbWriteTable(db_con, "deadwells_filesize", dw_df, overwrite=TRUE)

        # Log this run so re-runs are skipped unless --force is passed
        dbWriteTable(db_con, "_processing_log",
                     data.frame(
                         script    = "02.detectDeadBlurry_filesize",
                         timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                         row.names = NULL
                     ),
                     append=TRUE, overwrite=FALSE)
        dbDisconnect(db_con, shutdown=TRUE)
        options(echo=FALSE)
        cli_alert_success(" Processing log updated")
    } else {
        cli_alert_success(" Only annotating data")        
    }
} else {
    print("Currently no other method is available, set --method filesize to continue")
}
