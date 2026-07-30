options(stringsAsFactors=FALSE)

fancy_barplot <- function(listIn) {
    vals <- c(); vals_sd <- c()
    for(ii in 1:length(listIn)) {
        vals[ii] <- mean(listIn[[ii]])
        vals_sd[ii] <- sd(listIn[[ii]])
    }
    names(vals) <- names(allData)
    t1 <- barplot(vals,ylim=c(0,max(vals)*1.2),ylab="2HAX Intensity/cell (a.u.)",cex.names=0.7)
    t1
    arrows(x0=t1,y0=vals,y1=vals+vals_sd,angle=90,length=0.075)
    arrows(x0=t1,y0=vals,y1=vals-vals_sd,angle=90,length=0.075)
}

# Analyze the output now
allFiles <- list.files("analyzed/",full.names=TRUE)
data_median <- read.csv(allFiles[grep("Image",allFiles)],sep=",",header=TRUE)
metadata <- data_median[,which(grepl("Metadata",colnames(data_median)))]

# Test plot
H2AX_intensity <- data_median[,"Intensity_TotalIntensity_H2AX_Nuclei"]
nuclei_counts <- data_median[,"Count_Nuclei"]
# conc_subset <- which(data_median[,"Image_Metadata_Code"] %in% c("KCB25-11B","DMSO"))
# plot((H2AX_intensity/nuclei_counts),pch=19,col=as.numeric(as.factor(metadata[,2]))) # Everything just as a test
# plot((H2AX_intensity/nuclei_counts)[conc_subset],pch=19,col=as.numeric(as.factor(metadata[,2]))[conc_subset]) # Just KCB25-11B and DMSO

# Get all the data into a nice list
nameo <- apply(metadata[,1:2],1,paste,collapse="_")
nameo <- gsub(" ","",nameo)
unq <- unique(metadata[,1])
unq <- unq[c(6,1,4,3,2,5)]
allData <- list()
for(i in 1:length(unq)) {
    these <- which(metadata[,1]==unq[i])
    these_data <- (H2AX_intensity/nuclei_counts)[these]
    unq2 <- unique(metadata[these,2])
    for(k in 1:length(unq2)) {
        allData[[unique(nameo[these][which(metadata[these,2]==unq2[k])])]] <- these_data[which(metadata[these,2]==unq2[k])]
    }
}
# Reorder by drug and concentration
allData <- allData[c(1,4,3,2,7,6,5,9,8,14,13,12,11,10,19,18,17,16,15)]

# Make a fancy barplot with error bars, perfect for svg output and opening in Inkscape/Illustrator
fancy_barplot(lapply(allData,function(x) log(x,10)))

# Make a table for GraphPad also
gp_table <- matrix(data=NA,nrow=48,ncol=length(allData),dimnames=list(c(),names(allData)))
meta_table <- matrix(data=NA,nrow=2,ncol=length(allData),dimnames=list(c(),names(allData)))
for(i in 1:length(allData)) {
    working <- allData[[i]]
    if(length(working)>48) {
        working <- sample(working,48)
    } 
    gp_table[,i] <- working
    meta_table[,i] <- strsplit(names(allData)[i],"_")[[1]]
}
gp_final <- rbind(meta_table,gp_table)
write.table(gp_final,"graphpad_table.csv",sep=",",row.names=FALSE,col.names=TRUE,quote=FALSE)

