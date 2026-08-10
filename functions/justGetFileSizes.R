# =============================================================================
# functions/justGetFileSizes.R
#
# Purpose:  Scan raw Cell Painting image directories and record TIFF file sizes
#           for each fluorescence channel (merged, DAPI, GFP, TxRed, Cy5) on a
#           per-plate basis.  File size is used as a proxy for cell density to
#           flag dead or empty wells.
#
#           This script is sourced by 02.detectDeadBlurry_filesize_duck.R and
#           is NOT intended to be run standalone.
#
# Inputs:
#   basedir1, basedir2  root paths to the microscopy storage drives
#                       (defaults: Z:/ and Y:/ on Windows; /media/galactica/
#                       and /media/battlestar/ on Linux).  UPDATE THESE to match
#                       your local storage mount points before running.
#
# Outputs:
#   data/input/filesizes_all_fluorphores.Rds   named list of per-plate file sizes
#
# Notes:
#   - On first run the scan can take a long time for large image archives.
#     Subsequent runs skip plates already recorded in the existing .Rds file.
#   - Pass --force in script 02 to force a full re-scan.
# =============================================================================

options(stringsAsFactors=FALSE)
options(echo=FALSE)
options(cli.progress_show_after=0)
options(cli.condition="always")
options(warn=(-1))

library(cli) # For progress bars with an ETA

# ── Storage mount points ──────────────────────────────────────────────────────
# Update these paths to match your local microscopy storage drives.
if(tolower(names(which(ps::ps_os_type())))[length(tolower(names(which(ps::ps_os_type()))))] == "windows") {
    basedir1 <- "Z:/"; basedir2 <- "Y:/"
} else {
    basedir1 <- "/media/galactica/"; basedir2 <- "/media/battlestar/"
}

# ── Discover image directories ────────────────────────────────────────────────
cli_alert_info(" Discovering image directories on {basedir1} and {basedir2}")

allDirs  <- list.dirs(paste0(basedir1, "cell_painting/raw_images"), recursive=FALSE)
allNames1 <- allDirs
if(length(allNames1) > 0 && substr(allNames1[1], 1, 2) == "..") {
    allNames1 <- substr(allNames1, 3, nchar(allNames1))
}

allDirs2  <- list.dirs(paste0(basedir2, "cell_painting/raw_images"), recursive=FALSE)
allNames2 <- allDirs2
if(length(allNames2) > 0 && substr(allNames2[1], 1, 2) == "..") {
    allNames2 <- substr(allNames2, 3, nchar(allNames2))
}

allDirs3  <- paste0(basedir2, "cell_painting/to_be_processed")
allNames3 <- allDirs3
if(length(allNames3) > 0 && substr(allNames3[1], 1, 2) == "..") {
    allNames3 <- substr(allNames3, 3, nchar(allNames3))
}

# If a plate appears in both drives, prefer the copy in the second drive
names1  <- substr(allNames1, nchar(allNames1)-3, nchar(allNames1))
names2  <- substr(allNames2, nchar(allNames2)-3, nchar(allNames2))
rmNames <- which(names1 %in% names2)
if(length(rmNames) > 0) {
    allDirs <- c(allNames1[-rmNames], allNames2, allNames3)
}
cli_alert_success(" Found {length(allDirs)} top-level image director{?y/ies} to scan")

# ── Load existing file-size cache ─────────────────────────────────────────────
rds_path <- "data/input/filesizes_all_fluorphores.Rds"
filesizes <- list()
if(file.exists(rds_path)) {
    cli_alert_info(" Loading existing file-size cache from [{rds_path}]")
    filesizes_raw <- readRDS(rds_path)
    existo <- sapply(names(filesizes_raw), function(x) strsplit(x, "/")[[1]][length(strsplit(x, "/")[[1]])])
    names(existo) <- NULL
    cli_alert_success(" Cache loaded — {length(existo)} plate{?s} already recorded; these will be skipped")
} else {
    cli_alert_info(" No existing cache found — starting a fresh scan")
    filesizes_raw <- list()
    existo <- NULL
}

# ── Scan directories ──────────────────────────────────────────────────────────
for(k in 1:length(allDirs)) {
    cli_alert_info(" [{k}/{length(allDirs)}] Scanning directory: {allDirs[k]}")
    imgDirs <- list.dirs(allDirs[k], recursive=FALSE)

    if(length(imgDirs) == 0) {
        cli_alert_warning(" Path empty — nothing to scan")
        next
    }

    allNames <- imgDirs
    if(substr(allNames[1], 1, 2) == "..") {
        allNames <- substr(allNames, 3, nchar(allNames))
    }
    allNames_existo <- sapply(allNames, function(x) strsplit(x, "/")[[1]][length(strsplit(x, "/")[[1]])])
    names(allNames_existo) <- NULL

    # Skip plates already in the cache
    if(length(which(allNames_existo %in% existo)) > 0) {
        theseAreDone <- which(allNames_existo %in% existo)
        allNames <- allNames[-theseAreDone]
        imgDirs  <- imgDirs[-theseAreDone]
    }

    if(length(imgDirs) == 0) {
        cli_alert_success(" All plates in this directory are already cached")
        next
    }

    cli_progress_bar(
        paste0("Scanning plates in ", basename(allDirs[k])),
        total = length(imgDirs),
        clear = FALSE
    )
    for(i in 1:length(imgDirs)) {
        plate_name <- basename(imgDirs[i])

        # list.files() calls can be slow on large network drives — status is
        # shown inside the progress bar format string via the name token.
        cli_progress_update(set=i-1, status=paste0("listing files: ", plate_name))

        mergeFiles <- list.files(imgDirs[i], recursive=TRUE, full.names=TRUE, pattern="Plate_M")
        cli_progress_update(set=i-1, status=paste0("DAPI: ",   plate_name))
        dapiFiles  <- list.files(imgDirs[i], recursive=TRUE, full.names=TRUE, pattern=glob2rx("*Plate_R*d0*"))
        cli_progress_update(set=i-1, status=paste0("GFP: ",    plate_name))
        gfpFiles   <- list.files(imgDirs[i], recursive=TRUE, full.names=TRUE, pattern=glob2rx("*Plate_R*d1*"))
        cli_progress_update(set=i-1, status=paste0("TxRed: ",  plate_name))
        txredFiles <- list.files(imgDirs[i], recursive=TRUE, full.names=TRUE, pattern=glob2rx("*Plate_R*d2*"))
        cli_progress_update(set=i-1, status=paste0("Cy5: ",    plate_name))
        cy5Files   <- list.files(imgDirs[i], recursive=TRUE, full.names=TRUE, pattern=glob2rx("*Plate_R*d3*"))

        if(length(mergeFiles) > 10) {
            filesizes[[allNames[i]]]$merged <- file.size(mergeFiles)
            filesizes[[allNames[i]]]$dapi   <- file.size(dapiFiles)
            filesizes[[allNames[i]]]$gfp    <- file.size(gfpFiles)
            filesizes[[allNames[i]]]$txred  <- file.size(txredFiles)
            filesizes[[allNames[i]]]$cy5    <- file.size(cy5Files)
        }
        cli_progress_update()
    }
    cli_alert_success(" Directory complete — {length(imgDirs)} new plate{?s} scanned")
}

# ── Save results ──────────────────────────────────────────────────────────────
options(echo=TRUE)
saveRDS(append(filesizes_raw, filesizes), rds_path, version=2)
options(echo=FALSE)
cli_alert_success(" File-size cache saved to [{rds_path}] ({length(filesizes)} new plate{?s} added)")
