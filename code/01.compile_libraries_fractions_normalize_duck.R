# Initial option declarations
options(stringsAsFactors=FALSE) # Otherwise we need to force them as strings repeatedly
options(echo=FALSE) # Rscript needs options(echo=TRUE) to make an output file
options(cli.progress_show_after=0) # For progress bars
options(cli.condition="always") # For progress bars
options(warn=(-1)) # Ignore NA warnings
arguments <- commandArgs(trailingOnly=FALSE) # Will capture all arguments, so can search these later

library(cli)          # For progress bars with an ETA
library(DBI)          # Database interface
library(duckdb)       # Embedded analytical database
library(matrixStats)  # Fast colMedians, colSds etc. (avoids apply overhead)

# Header line
cli_h1("Compiling and normalizing input data")

# Get the peak of the density distribution
denPeak <- function(dataIn,adjval=0.7) {
    denT <- density(dataIn,adjust=adjval,na.rm=TRUE)
    denP <- denT$x[which(denT$y==max(denT$y,na.rm=TRUE))]
    return(denP)
}

# For arguments:
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
    "data/input/training",  # Input: training .RData files
    "data/output",          # Output: legacy .Rds files (fallback)
    "data/db",              # Output: DuckDB database file
    "logs/images",          # Output: QC PNG plots
    "functions",            # Source: functions.R
    "drugbank"             # Input: drugbank rip
)
for (d in required_dirs) {
    if (!dir.exists(d)) {
        dir.create(d, recursive=TRUE, showWarnings=FALSE)
        cli_alert_info(" Created missing directory: {d}")
    }
}

# DuckDB helper: write a named list of matrices/data.frames into a single DuckDB
# file, one table per element. A CP_Index column is prepended to every matrix
# so rows can be joined back together later.
write_to_duckdb <- function(compiled_list, db_path, overwrite=FALSE) {
    cli_h2("Writing to DuckDB")
    cli_alert_info(" Opening database at [{db_path}]")
    con <- dbConnect(duckdb(), dbdir=db_path, read_only=FALSE)
    on.exit(dbDisconnect(con, shutdown=TRUE), add=TRUE)

    # Tables that are 2-D (matrix or data.frame) and get a CP_Index key column
    matrix_tables <- c("metadata", "mesh", "db", "counts", "sd", "url", "method",
                       "features_raw", "features_1pass_1", "features_1pass_2", "features_1pass_3")

    # The row-key to attach to every table so tables can be joined
    cp_index <- compiled_list$metadata[, "CP_Index"]

    manifest_rows <- list()

    cli_progress_bar("Writing tables", total=length(matrix_tables), clear=FALSE)
    for (tbl in matrix_tables) {
        if (!tbl %in% names(compiled_list)) {
            cli_alert_warning(" Skipping [{tbl}] — not found in compiled list")
            cli_progress_update()
            next
        }
        df <- as.data.frame(compiled_list[[tbl]], stringsAsFactors=FALSE)
        # Prepend the CP_Index key so every table is joinable
        if (!"CP_Index" %in% colnames(df)) {
            df <- cbind(CP_Index=cp_index, df, stringsAsFactors=FALSE)
        }
        dbWriteTable(con, tbl, df,
                     overwrite=overwrite, append=!overwrite)
        manifest_rows[[tbl]] <- data.frame(
            table_name   = tbl,
            nrow         = nrow(df),
            ncol         = ncol(df),
            written_at   = as.character(Sys.time()),
            stringsAsFactors = FALSE
        )
        cli_progress_update()
    }

    # Store the type vector as a small lookup table
    if ("type" %in% names(compiled_list)) {
        type_df <- data.frame(CP_Index=cp_index,
                              sample_type=compiled_list$type,
                              stringsAsFactors=FALSE)
        dbWriteTable(con, "sample_type", type_df, overwrite=overwrite, append=!overwrite)
        manifest_rows[["sample_type"]] <- data.frame(
            table_name="sample_type", nrow=nrow(type_df), ncol=2,
            written_at=as.character(Sys.time()), stringsAsFactors=FALSE)
    }

    # Store the tag as a key-value metadata table
    if ("tag" %in% names(compiled_list) && !is.null(compiled_list$tag)) {
        tag_val <- compiled_list$tag
        if (is.list(tag_val)) tag_val <- unlist(tag_val)
        tag_df <- data.frame(
            key   = names(tag_val),
            value = as.character(tag_val),
            stringsAsFactors = FALSE
        )
        dbWriteTable(con, "_tag", tag_df, overwrite=TRUE)
    }

    # Write the manifest last
    manifest_df <- do.call(rbind, manifest_rows)
    dbWriteTable(con, "_manifest", manifest_df, overwrite=TRUE)

    all_tables <- dbListTables(con)
    cli_alert_success(" DuckDB contains [{length(all_tables)}] tables: {paste(all_tables, collapse=', ')}")
    invisible(db_path)
}

# Confirm the normalization method
norm_method <- "density"
if(any(arguments=="--norm")) {
    arg_ind <- which(arguments=="--norm")
    trailing_arg <- arguments[arg_ind[1]+1]
    if(length(arg_ind)>=1) {
        cli_alert_success(" Set normalization method {trailing_arg}")
        norm_method <- trailing_arg
    } else {
        cli_alert_warning(" Normalization change indicated but no normalization method provided.  Defaulting to [density]")
    }
}

# Functions and library loading
source("functions/functions.R")
cli_alert_success(" Loaded dependencies")

# Check if script 01 has already been applied to this database.
# Pass --force on the command line to bypass and re-run anyway.
already_processed_01 <- FALSE
db_path_check <- "data/db/cellpainting.duckdb"
if (file.exists(db_path_check)) {
    con_check <- dbConnect(duckdb(), dbdir=db_path_check, read_only=TRUE)
    avail_tbls_check <- dbListTables(con_check)
    if ("_processing_log" %in% avail_tbls_check) {
        proc_log_check <- dbReadTable(con_check, "_processing_log")
        if (any(proc_log_check$script == "01.compile_libraries_fractions_normalize")) {
            already_processed_01 <- TRUE
            prev_ts_01 <- tail(proc_log_check$timestamp[proc_log_check$script == "01.compile_libraries_fractions_normalize"], 1)
            cli_alert_info(" Database already processed by script 01 on {prev_ts_01}")
        }
    }
    dbDisconnect(con_check, shutdown=TRUE)
}
if (already_processed_01 && !any(arguments == "--force")) {
    cli_alert_info(" Skipping. Use --force to re-run script 01 on this database.")
    quit(save="no", status=0)
}

if(any(arguments=="--emerg")) {
    cli_alert_info(" Emergency mode activated, re-emerging database from Rds")
    cp_compiled <- readRDS("data/output/cp_compiled_all_norm.Rds")
    cp_compiled$metadata[,1] <- vapply(seq_len(nrow(cp_compiled$metadata)),
                                       function(i) cli::hash_sha1(paste(cp_compiled$metadata[i,], collapse="-")),
                                       character(1))
    colnames(cp_compiled$metadata)[1] <- "CP_Index"
    db_path <- "data/db/cellpainting.duckdb"
    options(echo=TRUE)
    write_to_duckdb(cp_compiled, db_path=db_path, overwrite=TRUE)
    options(echo=FALSE)
    cli_alert_success(" DuckDB written to [{db_path}]")

    # Log this run so re-runs are skipped unless --force is passed
    con_log <- dbConnect(duckdb(), dbdir=db_path, read_only=FALSE)
    dbWriteTable(con_log, "_processing_log",
                    data.frame(
                        script    = "01.compile_libraries_fractions_normalize",
                        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                        row.names = NULL
                    ),
                    append=TRUE, overwrite=FALSE)
    dbDisconnect(con_log, shutdown=TRUE)
    cli_alert_success(" Processing log updated")
    quit(save="no", status=0)
}

# Mendeley download URLs — one per library (in order: ATCA, BIO2, CTNS, PHKN, PSYC, TMOL, TRKB)
mendeley_urls <- c(
    ATCA = "https://data.mendeley.com/api/datasets/3kw7jpbc5h/draft/files/3c3a3941-0cfa-40f2-9696-311085bf60f9",
    BIO2 = "https://data.mendeley.com/api/datasets/3kw7jpbc5h/draft/files/8d60a435-bbda-4bfb-93cb-d5ac600f1d6a",
    CTNS = "https://data.mendeley.com/api/datasets/3kw7jpbc5h/draft/files/8c695fcc-950f-4ac2-9707-d83f7c6edf22",
    PHKN = "https://data.mendeley.com/api/datasets/3kw7jpbc5h/draft/files/88e6fff7-f494-4622-87aa-5ad0d7994045",
    PSYC = "https://data.mendeley.com/api/datasets/3kw7jpbc5h/draft/files/d53fe3ef-0988-44e7-bd27-ae6df5e040aa",
    TMOL = "https://data.mendeley.com/api/datasets/3kw7jpbc5h/draft/files/976fa53b-37e4-474a-9986-28f9f7ffbcf1",
    TRKB = "https://data.mendeley.com/api/datasets/3kw7jpbc5h/draft/files/3e333beb-bac8-4049-8fe2-77fdd9d00142"
)

# Check if any .RData files exist in the training directory
fileList <- list.files("data/input/training", pattern="\\.RData$|\\.Rdata$",
                       recursive=FALSE, full.names=TRUE)

if (length(fileList) == 0) {
    cli_alert_warning(" No training .RData files found in [data/input/training] — downloading from Mendeley")
    cli_progress_bar("Downloading libraries", total=length(mendeley_urls), clear=FALSE)
    for (lib in names(mendeley_urls)) {
        dest <- file.path("data/input/training", paste0(lib, ".RData"))
        tryCatch(
            download.file(mendeley_urls[lib], destfile=dest, mode="wb", quiet=TRUE),
            error = function(e) cli_alert_danger(" Failed to download {lib}: {conditionMessage(e)}")
        )
        cli_progress_update()
    }
    # Re-scan after download
    fileList <- list.files("data/input/training", pattern="\\.RData$|\\.Rdata$",
                           recursive=FALSE, full.names=TRUE)
    if (length(fileList) == 0) {
        cli_alert_danger(" Download failed — no training files available. Check your internet connection and try again.")
        quit(save="no", status=1)
    }
    cli_alert_success(" Downloaded {length(fileList)} library files to [data/input/training]")
}

# Load all training files
cli_progress_bar("Loading training data", total=length(fileList), clear=FALSE)
for (i in 1:length(fileList)) {
    load(fileList[i])
    cli_progress_update()
}
cli_alert_success(" Training set loaded")

# Repair bioactives2
if(any(colnames(cp_bioa2$metadata) %in% c("Image_CPD_MMOL_CONC","Image_CPD_SMILES"))) {
    cli_alert_info(" Performing some cleanup")
    cp_bioa2$metadata[,"Image_Metadata_CPD_SMILES"] <- cp_bioa2$metadata[,"Image_CPD_SMILES"]
    cp_bioa2$metadata <- cp_bioa2$metadata[,-which(colnames(cp_bioa2$metadata) %in% c("Image_CPD_MMOL_CONC","Image_CPD_SMILES"))]
    if(any(arguments=="--nogpt")) {
        removedannots <- round((length(which(is.na(cp_bioa2$drugbank_old[,12])))/nrow(cp_bioa2$drugbank_old))*100,0)
        cp_bioa2$drugbank <- cp_bioa2$drugbank_old
        cli_alert_info(" Removed the GPT-annotated drugbank target annotations ({removedannots}% of Bioactives2 library)")
    } else {
        cli_alert_info(" Using GPT-boosted drugbank target annotations (annotated additional {removedannots} of Bioactives2 library)")
    }
}

# Compile the training compounds first
cli_alert_info(" Binding the training and test sets")
cp_complete <- list()
cp_complete$metadata <- rbind(cp_bioa2$metadata,cp_pharmakon$metadata,cp_anticancer$metadata,cp_psycho$metadata,cp_caithness$metadata,cp_trkb$metadata,cp_tmol$metadata)
cp_complete$counts <- rbind(cp_bioa2$counts,cp_pharmakon$counts,cp_anticancer$counts,cp_psycho$counts,cp_caithness$counts,cp_trkb$counts,cp_tmol$counts)
cp_complete$median <- rbind(cp_bioa2$median,cp_pharmakon$median,cp_anticancer$median,cp_psycho$median,cp_caithness$median,cp_trkb$median,cp_tmol$median)
cp_complete$texture <- rbind(cp_bioa2$texture,cp_pharmakon$texture,cp_anticancer$texture,cp_psycho$texture,cp_caithness$texture,cp_trkb$texture,cp_tmol$texture)
cp_complete$sd <- rbind(cp_bioa2$sd,cp_pharmakon$sd,cp_anticancer$sd,cp_psycho$sd,cp_caithness$sd,cp_trkb$sd,cp_tmol$sd)
cp_complete$url <- rbind(cp_bioa2$url,cp_pharmakon$url,cp_anticancer$url,cp_psycho$url,cp_caithness$url,cp_trkb$url,cp_tmol$url)
cp_complete$method <- rbind(cp_bioa2$method,cp_pharmakon$method,cp_anticancer$method,cp_psycho$method,cp_caithness$method,cp_trkb$method,cp_tmol$method)
cp_complete$mesh <- rbind(cp_bioa2$mesh,cp_pharmakon$mesh,cp_anticancer$mesh,cp_psycho$mesh,cp_caithness$mesh,cp_trkb$mesh,cp_tmol$mesh)
cp_complete$drugbank <- rbind(cp_bioa2$drugbank,cp_pharmakon$drugbank,cp_anticancer$drugbank,cp_psycho$drugbank,cp_caithness$drugbank,cp_trkb$drugbank,cp_tmol$drugbank)
cp_complete$type <- rep("training",nrow(cp_complete$metadata))

# Update the ArbPlate IDs
newStarts <- which(cp_complete$metadata[,12]==1 & cp_complete$metadata[,5]=="A01")
ArbPlate <- cp_complete$metadata[,12]
for(p in 2:length(newStarts)) {
    oldmax <- cp_complete$metadata[newStarts[p]-1,12]
    ArbPlate[newStarts[p]:nrow(cp_complete$metadata)] <- ArbPlate[newStarts[p]:nrow(cp_complete$metadata)]+oldmax
}
cp_complete$metadata[,12] <- ArbPlate


# Duplicate a working list here and remove the individual variables that were in the before
cp_working <- cp_complete
# matrix(as.numeric()) is faster than apply(..., as.numeric) because apply
# coerces to a list internally before simplifying back to a matrix
cp_working$median  <- matrix(as.numeric(cp_complete$median[,-1]),
                                nrow=nrow(cp_complete$median),
                                dimnames=list(NULL, colnames(cp_complete$median)[-1]))
cp_working$texture <- matrix(as.numeric(cp_complete$texture[,-1]),
                                nrow=nrow(cp_complete$texture),
                                dimnames=list(NULL, colnames(cp_complete$texture)[-1]))
cp_working$features <- cbind(cp_working$median, cp_working$texture)
rownames(cp_working$features) <- cp_complete$metadata[,"CP_Index"]
rm(cp_complete)
gc(verbose=FALSE)

# Vectorised plate-ID repair: use sub() + sprintf() instead of row-wise sapply/strsplit
# This is ~10-50x faster on large metadata matrices
fixPlateCode <- function(vec) {
    # Normalise separator, extract numeric part, zero-pad to 3 digits
    normed  <- gsub("_", "-", vec)
    prefix  <- sub("-.*", "", normed)
    numpart <- regmatches(normed, regexpr("[0-9]+", normed))
    ifelse(is.na(vec) | nchar(vec) <= 1, NA,
            paste(prefix, sprintf("%03d", as.integer(numpart)), sep="-"))
}
c3 <- fixPlateCode(cp_working$metadata[,3])
c7 <- fixPlateCode(cp_working$metadata[,7])

# ifelse() on vectors is faster than apply(cbind(), 1, ...)
fixedCols <- ifelse(is.na(c3), c7, c3)
cp_working$metadata[,3] <- paste(fixedCols, cp_working$metadata[,5], sep="-")
cp_working$metadata[,7] <- fixedCols

# NA imputation: matrixStats::colMedians is a compiled C loop, far faster than
# apply()-ing median() column-by-column from R
col_medians <- colMedians(cp_working$features, na.rm=TRUE)
na_counts   <- colSums(is.na(cp_working$features))
cp_features_1pass_1 <- cp_working$features
# Only impute columns with a manageable number of NAs (same threshold as before)
impute_cols <- which(na_counts > 0 & na_counts < 5000)
for (j in impute_cols) {
    cp_features_1pass_1[is.na(cp_features_1pass_1[, j]), j] <- col_medians[j]
}
# Columns with >= 5000 NAs get wiped to all-NA (preserving original logic)
wipe_cols <- which(na_counts >= 5000)
if (length(wipe_cols) > 0) cp_features_1pass_1[, wipe_cols] <- NA

# QC plot for raw data
cli_alert_info(" Generating plot of raw data in [logs/images]")
options(echo=TRUE)
png(filename="logs/images/feature_raw.png",width=10,height=7,units="in",bg="transparent",res=600)
    plot(cp_working$features[,1110],xlab="Sample Index",ylab=paste0(colnames(cp_working$features)[1110]," values"),pch=19,col=rgb(0,0,0,0.2),cex=0.5)
dev.off()
options(echo=FALSE)

unq_pids <- unique(cp_working$metadata[,12])
rownames(cp_features_1pass_1) <- seq(1,nrow(cp_features_1pass_1)); rownames(cp_working$features) <- rownames(cp_features_1pass_1)
cli_alert_info(" Normalizing data")
cli_progress_bar("First pass normalization        ",total=length(unq_pids),clear=FALSE)
for(i in 1:length(unq_pids)) {
    s_working <- which(cp_working$metadata[,12]==unq_pids[i])
    # Pass the matrix slice directly - avoids rebuilding a named matrix each iteration
    cp_features_1pass_1[s_working,] <- rescale(cp_features_1pass_1[s_working,,drop=FALSE],
                                                cp_features2=cp_working$features,
                                                adjval=0.7, meth=norm_method)
    cli_progress_update()
}
# QC plot for 1-pass normalized data
cli_alert_info(" Generating plot of 1-pass normalized data in [logs/images]")
options(echo=TRUE)
png(filename="logs/images/feature_norm_1pass.png",width=10,height=7,units="in",bg="transparent",res=600)
    plot(cp_features_1pass_1[,1110],xlab="Sample Index",ylab=paste0(colnames(cp_features_1pass_1)[1110]," values"),pch=19,col=rgb(0,0,0,0.2),cex=0.5)
dev.off()
options(echo=FALSE)

cp_features_1pass_2 <- cp_features_1pass_1
cli_progress_bar("Second pass normalization       ",total=length(unq_pids),clear=FALSE)
for(i in 1:length(unq_pids)) {
    s_working <- which(cp_working$metadata[,12]==unq_pids[i])
    # colMedians on the slice is faster than apply(..., median) - same compiled C path
    plate_medians <- colMedians(cp_features_1pass_2[s_working,,drop=FALSE], na.rm=TRUE)
    cp_features_1pass_2[s_working,] <- sweep(cp_features_1pass_2[s_working,,drop=FALSE], 2, plate_medians)
    cli_progress_update()
}
# QC plot for 2-pass normalized data
cli_alert_info(" Generating plot of 2-pass normalized data in [logs/images]")
options(echo=TRUE)
png(filename="logs/images/feature_norm_2pass.png",width=10,height=7,units="in",bg="transparent",res=600)
    plot(cp_features_1pass_2[,1110],xlab="Sample Index",ylab=paste0(colnames(cp_features_1pass_2)[1110]," values"),pch=19,col=rgb(0,0,0,0.2),cex=0.5)
dev.off()
options(echo=FALSE)

cp_features_1pass_3 <- cp_features_1pass_2
cli_progress_bar("Third (final) pass normalization",total=length(unq_pids),clear=FALSE)
for(i in 1:length(unq_pids)) {
    s_working <- which(cp_working$metadata[,12]==unq_pids[i])
    # Compute denPeak once per column for this plate, then subtract as a vector sweep
    # avoids calling apply() which has per-column function-call overhead
    peaks <- apply(cp_features_1pass_3[s_working,,drop=FALSE], 2,
                    function(x) denPeak(x, adjval=0.7)[1])
    cp_features_1pass_3[s_working,] <- sweep(cp_features_1pass_3[s_working,,drop=FALSE], 2, peaks)
    cli_progress_update()
}
# QC plot for 3-pass normalized data
cli_alert_info(" Generating plot of 3-pass normalized data in [logs/images]")
options(echo=TRUE)
png(filename="logs/images/feature_norm_3pass.png",width=10,height=7,units="in",bg="transparent",res=600)
    plot(cp_features_1pass_3[,1110],xlab="Sample Index",ylab=paste0(colnames(cp_features_1pass_3)[1110]," values"),pch=19,col=rgb(0,0,0,0.2),cex=0.5)
dev.off()
options(echo=FALSE)

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

cp_compiled <- list()
cp_compiled$metadata <- cp_working$metadata
cp_compiled$mesh <- cp_working$mesh
cp_compiled$db <- cp_working$drugbank
cp_compiled$counts <- cp_working$counts
cp_compiled$sd <- cp_working$sd
cp_compiled$url <- cp_working$url
cp_compiled$method <- cp_working$method
cp_compiled$features_raw <- cp_working$features
cp_compiled$features_1pass_1 <- cp_features_1pass_1
cp_compiled$features_1pass_2 <- cp_features_1pass_2
cp_compiled$features_1pass_3 <- cp_features_1pass_3
cp_compiled$type <- cp_working$type
cp_compiled$tag <- tag

# vapply is faster than apply for a scalar-returning function on rows:
# it pre-allocates the result vector and avoids the overhead of apply's list path
cp_compiled$metadata[,1] <- vapply(seq_len(nrow(cp_compiled$metadata)),
                                    function(i) cli::hash_sha1(paste(cp_compiled$metadata[i,], collapse="-")),
                                    character(1))
colnames(cp_compiled$metadata)[1] <- "CP_Index"
cp_compiled <- global_rownames(cp_compiled)

# Write to DuckDB (primary output)
db_path <- "data/db/cellpainting.duckdb"
overwrite_db <- any(arguments=="--force")

if (file.exists(db_path) && !overwrite_db) {
    cli_alert_warning(" DuckDB file already exists at [{db_path}]. Use --force to overwrite.")
} else {
    if (file.exists(db_path) && overwrite_db) {
        cli_alert_info(" --force flag set: removing existing DuckDB file")
        file.remove(db_path)
    }
    options(echo=TRUE)
    write_to_duckdb(cp_compiled, db_path=db_path, overwrite=TRUE)
    options(echo=FALSE)
    cli_alert_success(" DuckDB written to [{db_path}]")

    # Log this run so re-runs are skipped unless --force is passed
    con_log <- dbConnect(duckdb(), dbdir=db_path, read_only=FALSE)
    dbWriteTable(con_log, "_processing_log",
                    data.frame(
                        script    = "01.compile_libraries_fractions_normalize",
                        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                        row.names = NULL
                    ),
                    append=TRUE, overwrite=FALSE)
    dbDisconnect(con_log, shutdown=TRUE)
    cli_alert_success(" Processing log updated")
}

# Fallback: also write the legacy .Rds so nothing breaks
options(echo=TRUE)
saveRDS(cp_compiled, "data/output/cp_compiled_all_norm.Rds", version=2)
options(echo=FALSE)
cli_alert_success(" Fallback [cp_compiled_all_norm.Rds] also written to data/output")
cat("\n")





