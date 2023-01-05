library(dplyr)
library(xcms)
library(CAMERA)
library(here)
here()

register(SerialParam(), default=TRUE)

#ReadData
combined <- "COMBINED_132"
samples <- list.files(paste0("/Volumes/GoogleDrive-108763975232397414829/My Drive/Ecuador_XCMS_Projects/Duroia/data/",combined,"/"),pattern = ".mzML",full.names = TRUE) 
blank <-  list.files("/Volumes/GoogleDrive-108763975232397414829/My Drive/Ecuador_XCMS_Projects/Duroia/data/Blank/",pattern = ".mzML",full.names = TRUE)
RTI <-  list.files("/Volumes/GoogleDrive-108763975232397414829/My Drive/Ecuador_XCMS_Projects/Duroia/data/Standard/",pattern = ".mzML",full.names = TRUE)


MSE_level_1 <- list.files("/Volumes/GoogleDrive-108763975232397414829/My Drive/Ecuador_XCMS_Projects/Duroia/data/MSE/",pattern = "*MSE.mzML",full.names = TRUE)


MSE_deconvoluted<- list.files("/Volumes/GoogleDrive-108763975232397414829/My Drive/Ecuador_XCMS_Projects/Duroia/data/MSE/",pattern = "*_Spectrum.mzML",full.names = TRUE)

#datafiles <- c(samples,blank,RTI[1],MSE_level_1,MSE_deconvoluted)
#s_groups <-  c(rep("Samples",length(samples)),rep("Blank",length(blank)),rep("Standard",length(RTI[1])),rep("MSE_L1",length(MSE_level_1)),rep("MSE_L2",length(MSE_deconvoluted)))

datafiles <- c(samples,blank,RTI[1],MSE_level_1)
s_groups <-  c(rep("Samples",length(samples)),rep("Blank",length(blank)),rep("Standard",length(RTI[1])),rep("MSE_L1",length(MSE_level_1)))

pd <- data.frame(file = basename(datafiles), injection_idx = 1:length(datafiles), sample = gsub(".mzML","",basename(datafiles)), group =s_groups)
rawData <- readMSData(datafiles, pdata = new("NAnnotatedDataFrame", pd), centroided. = TRUE,mode= "onDisk", verbose=T)

cwp <- CentWaveParam(peakwidth = c(2, 10), ppm=15, integrate = 2, noise=500, snthresh=0)
processedData <- findChromPeaks(rawData, param = cwp)

pdp <- PeakDensityParam(sampleGroups = pData(processedData)$group, bw=4, binSize=0.5, minSamples=1, minFraction = 0.001)

processedData <- groupChromPeaks(processedData, param = pdp)

medWidth <- median(chromPeaks(processedData)[, "rtmax"] -
                     chromPeaks(processedData)[, "rtmin"])

#pgp <- PeakGroupsParam(minFraction = 0.45, span = 0.6,subset=which(pData(processedData)$group == "Standard"),subsetAdjust = "previous")

#processedData <- adjustRtime(processedData, param = pgp) 

processedData <- fillChromPeaks(processedData,
                                param = FillChromPeaksParam(fixedRt = medWidth))



xset <- as(filterMsLevel(processedData, msLevel = 1L), "xcmsSet")
sampclass(xset) <-pData(processedData)$group

xsa <- xsAnnotate(xset, polarity = "negative")
xsaF <- groupFWHM(xsa, sigma = 6, perfwhm = 1)
xsaC <- groupCorr(xsaF, cor_eic_th = 0.6, pval = 0.05, graphMethod = "hcs",
                  calcCiS = TRUE, calcCaS = FALSE, calcIso = FALSE)
xsaFI <- findIsotopes(xsaC, maxcharge = 2, maxiso = 3, minfrac = 0.5, ppm = 10, intval = "maxo")

xsaFA <- findAdducts(xsaFI, polarity = "negative", 
                     max_peaks = 100, multiplier = 3, ppm = 10)
#edgelist <- getEdgelist(xsaFA)

xset5 <- getPeaklist(xsaFA)

#Remove Peaks Found in the Blank
colnames(xset5)[13:25] <- pData(processedData)$sample
xset6 <- xset5[which(xset5$BN50_57 < 1000),]

#Shorten Retention Time Window
xset7 <- xset6[which(xset6$rt > 25 & xset6$rt < 660),]

#remove anything with max value less than 5000

abund_df <- xset7[,c(13:25)]
abund_df[is.na(abund_df)] <- 0

xset7$max_value <- as.numeric(apply(abund_df, 1, max))

xset8 <- xset7[which(as.numeric(apply(abund_df, 1, max)) > 5000),c(1,4,13:25,28,29)]

max_pc <- aggregate(xset8$max_value, by = list(xset8$pcgroup), max)


xset9 <- as.data.frame(xset8 %>%
                         group_by(pcgroup) %>%
                         filter(max_value == max(max_value, na.rm=TRUE)))

xset9 <- xset9[order(xset9$max_value,decreasing = T),]








#### Match Features between MS1 and MS2

MSE_deconvoluted<- list.files("/Volumes/GoogleDrive-108763975232397414829/My Drive/Ecuador_XCMS_Projects/Duroia/data/MSE/",pattern = ".mgf",full.names = TRUE)


library(MergeION)

library(stringr)
str_split_fixed(before$type, "_and_", 2)

mse_DF<- data.frame()
for(file in MSE_deconvoluted){
mse_mgf2 <- readMGF2(file)
mse_DF_single <- data.frame(scans= mse_mgf2$metadata$TITLE, precursorMZ_Intensity = mse_mgf2$metadata$PEPMASS, retentionTime= as.numeric(mse_mgf2$metadata$RTINSECONDS))
mse_DF_single$precursorMZ <- as.numeric(str_split_fixed(mse_DF_single$precursorMZ_Intensity," ",2)[,1])
mse_DF_single$precursorMZ_TIC <- as.numeric(str_split_fixed(mse_DF_single$precursorMZ_Intensity," ",2)[,2])
mse_DF_single$scans <- as.numeric(gsub("Spectrum ","",mse_DF_single$scans))
mse_DF_single$file <- gsub(".mgf","",basename(file))
mse_DF <- rbind(mse_DF,mse_DF_single)}

table(mse_DF$file)

head(mse_DF)


#features_spec <- featureDefinitions(processedData)
features_spec <- xset5

FT231= features_spec[which.min(abs(features_spec$mz - 231.07)), ]

x=203

x = mse_DF$scans[which.min(abs(as.numeric(mse_DF$precursorMZ) - 231.07)) ]

match_features <- function(x){
  ms_feature <- features_spec[row.names(features_spec)  == x,]
  mz <- ms_feature$mz
  mzmin <- mz - 0.5
  mzmax <- mz + 0.5
  rt <- ms_feature$rt
  rtmin <- rt-60
  rtmax <- rt+60
  pot_feature <- mse_DF[mse_DF$precursorMZ >= mzmin & mse_DF$precursorMZ <= mzmax & mse_DF$retentionTime >= rtmin & mse_DF$retentionTime <= rtmax,]
  
  if(nrow(pot_feature) == 0){return(data.frame(FeatID = x,msms_scan = "no_match",file="no_match"))}
  if(nrow(pot_feature) == 1){return(data.frame(FeatID = x ,msms_scan = pot_feature$scans, file=pot_feature$file))}
  if(nrow(pot_feature) > 1) {
    closest_feature <- pot_feature[which.min(abs(pot_feature$precursorMZ - mz)+abs(pot_feature$retentionTime - rt)),]
    return(data.frame(FeatID = x ,msms_scan = closest_feature$scans,file=closest_feature$file))}
}

#match_features <- function(x){
#  msms_scan <- mse_DF[mse_DF$scans == x,]
#  mz <- as.numeric(msms_scan$precursorMZ)
#  mzmin <- mz - 0.02
#  mzmax <- mz + 0.02
#  rt <- as.numeric(msms_scan$retentionTime)*60
#  rtmin <- rt-60
#  rtmax <- rt+60
#  pot_feature <- features_spec[features_spec$mzmed >= mzmin & features_spec$mzmed <= mzmax & features_spec$rtmed >= rtmin & features_spec$rtmed <= rtmax,]
  
#  if(nrow(pot_feature) == 0){return(data.frame(scan = x,feature_number = "no_match"))}
#  if(nrow(pot_feature) == 1){return(data.frame(scan = x,feature_number = as.character(row.names(pot_feature))))}
#  if(nrow(pot_feature) > 1) {
#    closest_feature <- pot_feature[which.min(abs(pot_feature$mzmed - mz)+abs(pot_feature$rtmed*60 - rt)),]
#    return(data.frame(scan = x ,feature_number = as.character(row.names(closest_feature))))}
#}


match_features(x)

ms_ms_match_results <- lapply(row.names(features_spec),FUN = match_features)
ms_ms_match_results_df <- do.call("rbind", ms_ms_match_results)

table(ms_ms_match_results_df$msms_scan =="no_match")
table(ms_ms_match_results_df$file)

#merge onto xset5

xset_Matched <-merge(as.data.frame(xset5),ms_ms_match_results_df,by.x=0,by.y="FeatID")
names(xset_Matched)


abund_df <- xset_Matched[,c(14:26)]
abund_df[is.na(abund_df)] <- 0

xset_Matched$max_value <- as.numeric(apply(abund_df, 1, max))
xset_Matched_simple <- xset_Matched[order(c(as.numeric(xset_Matched$pcgroup),xset_Matched$max_value)),c("mz","rt","pcgroup","max_value","isotopes","adduct","msms_scan","file")]

table(unique(xset_Matched_simple$pcgroup) %in% unique(xset_Matched_simple$pcgroup[xset_Matched_simple$msms_scan !="no_match"]))

#### Need to take a look at matches and set real filter
#### Need to only look at major TIC amounts.... and see what portion of TIC has ms/ms matches...



