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

# Initial option declarations
options(stringsAsFactors=FALSE) # Otherwise we need to force them as strings repeatedly
options(echo=FALSE) # Rscript needs options(echo=TRUE) to make an output file
options(cli.progress_show_after=0) # For progress bars
options(cli.condition="always") # For progress bars
options(warn=(-1)) # Ignore NA warnings
arguments <- commandArgs(trailingOnly=FALSE) # Will capture all arguments, so can search these later

library(cli) # For progress bars with an ETA

if(tolower(names(which(ps::ps_os_type())))[length(tolower(names(which(ps::ps_os_type()))))]=="windows") {
    basedir1 <- "Z:/"; basedir2 <- "Y:/"
} else {
    basedir1 <- "/media/galactica/"; basedir2 <- "/media/battlestar/"
}

# Run this on the server to screen the training library for DBA plates
allDirs <- list.dirs(paste0(basedir1,"cell_painting/raw_images"),recursive=FALSE)
allNames1 <- allDirs
if(substr(allNames1[1],1,2)=="..") allNames1 <- substr(allNames1,3,nchar(allNames1))

allDirs2 <- list.dirs(paste0(basedir2,"cell_painting/raw_images"),recursive=FALSE)
allNames2 <- allDirs2
if(substr(allNames2[1],1,2)=="..") allNames2 <- substr(allNames2,3,nchar(allNames2))

allDirs3 <- paste0(basedir2,"cell_painting/to_be_processed")
allNames3 <- allDirs3
if(substr(allNames3[1],1,2)=="..") allNames3 <- substr(allNames3,3,nchar(allNames3))

# Check if the directories are duplicated, and if they are, choose the ones in the unknowns directory
names1 <- substr(allNames1,nchar(allNames1)-3,nchar(allNames1))
names2 <- substr(allNames2,nchar(allNames2)-3,nchar(allNames2))
rmNames <- which(names1 %in% names2)
if(length(rmNames)>0) {
    allDirs <- c(allNames1[-rmNames],allNames2,allNames3)
}

filesizes <- list()
if(file.exists("code/filesizes_all_fluorphores.Rds")) {
    filesizes_raw <- readRDS("code/filesizes_all_fluorphores.Rds")
    existo <- sapply(names(filesizes_raw),function(x) strsplit(x,"/")[[1]][length(strsplit(x,"/")[[1]])])
    names(existo) <- NULL
} else {
    existo <- NULL
}

for(k in 1:length(allDirs)) {
    print(paste0("Scanning ",allDirs[k]))
    imgDirs <- list.dirs(allDirs[k],recursive=FALSE)
    if(length(imgDirs)>0) {
        allNames <- imgDirs
        if(substr(allNames[1],1,2)=="..") {
            allNames <- substr(allNames,3,nchar(allNames))
        }
        allNames_existo <- sapply(allNames,function(x) strsplit(x,"/")[[1]][length(strsplit(x,"/")[[1]])])
        names(allNames_existo) <- NULL

        if(length(which(allNames_existo %in% existo))>0) {
            theseAreDone <- which(allNames_existo %in% existo)
            allNames <- allNames[-theseAreDone]
            imgDirs <- imgDirs[-theseAreDone]
        }

        if(length(imgDirs)>0) {
            cli_progress_bar(" Processing ",total=length(imgDirs),clear=FALSE)
            for(i in 1:length(imgDirs)) {
                mergeFiles <- list.files(imgDirs[i],recursive=TRUE,full.names=TRUE,pattern="Plate_M")
                dapiFiles <- list.files(imgDirs[i],recursive=TRUE,full.names=TRUE,pattern=glob2rx("*Plate_R*d0*"))
                gfpFiles <- list.files(imgDirs[i],recursive=TRUE,full.names=TRUE,pattern=glob2rx("*Plate_R*d1*"))
                txredFiles <- list.files(imgDirs[i],recursive=TRUE,full.names=TRUE,pattern=glob2rx("*Plate_R*d2*"))
                cy5Files <- list.files(imgDirs[i],recursive=TRUE,full.names=TRUE,pattern=glob2rx("*Plate_R*d3*"))
                if(length(mergeFiles)>10) {
                    filesizes[[allNames[i]]]$merged <- file.size(mergeFiles)
                    filesizes[[allNames[i]]]$dapi <- file.size(dapiFiles)
                    filesizes[[allNames[i]]]$gfp <- file.size(gfpFiles)
                    filesizes[[allNames[i]]]$txred <- file.size(txredFiles)
                    filesizes[[allNames[i]]]$cy5 <- file.size(cy5Files)
                }
                cli_progress_update()
            }
            
            print(paste0("---> Up to date as of ",Sys.time()))
        } else {
            print(paste0("---> Up to date as of ",Sys.time()))
        }
    } else {
        print(paste0("---> Path empty as of ",Sys.time()))
    }
}

# filesizes1 <- readRDS("D:/Dropbox/Data/McMaster/Cell Painting/bad_images/filesizes_all_fluorphores.Rds")
# filesizes2 <- readRDS("D:/Dropbox/Data/McMaster/Cell Painting/bad_images/filesizes_all_fluorphores_22.Rds")
# # filesizes2 <- filesizes
# filesizes3 <- append(filesizes1,filesizes2)
# View(names(filesizes3))

saveRDS(append(filesizes_raw,filesizes),"code/filesizes_all_fluorphores.Rds",version=2)


