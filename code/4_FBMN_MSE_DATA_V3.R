library(dplyr)
library(xcms)
library(CAMERA)
library(here)
library(MergeION)
library(stringr)

here()

source("https://raw.githubusercontent.com/jorainer/xcms-gnps-tools/master/customFunctions.R")

#Serial vs multicore. It appears that paraellizing often causes errors. It's slower to run in serial but seems to work better. 
register(SerialParam(), default=TRUE)
#alternatively on a PC we can:
#register(bpstart(SnowParam(4)))


#set working directory to the google drive space.
#The following will determine which folder based on if it's a PC or a Mac. edit the line setWD() for the path where google data are stored. 

if (.Platform$OS.type == "unix") {
  register(SerialParam(), default=TRUE)
  setwd("/Volumes/GoogleDrive-108763975232397414829/My Drive/Ecuador_XCMS_Projects/")
  #register(bpstart(MulticoreParam(4)))
} else {
  register(bpstart(SnowParam(4)))
  setwd("E:/My Drive/Ecuador_XCMS_Projects/")
}
register(SerialParam())


### The UPLC Master file has all of the information required for loading data and the order that files were run in, which gets used later on for standards 
UPLC_Results_Master <- read.csv("./UPLC_RESULTS_AMAZON_PROJECT_JUNE2022_V6.csv")

#First we  subset to a specifc project....
project <- "DUROIA"

UPLC_Results <- UPLC_Results_Master[UPLC_Results_Master[,project] == "YES",]

#get all of the information from the uplc_results table for building the metatadata tables associated with the project and store in pd dataframe

### Get the file paths for each of the datatypes.
#MS level 1 all
mzml_all <- UPLC_Results$Converted_Check



#define sample groups (samples, blanks and standards)
s_groups <-  UPLC_Results$sample_type
id <- UPLC_Results$X.id.
sample <- gsub(".mzML","",basename(mzml_all))


pd <- data.frame(injection_idx = id, sample=sample, file = basename(mzml_all), path= mzml_all, group =s_groups)

#Read all of the data as an XCMS project.
#this slow but can't be spead up with par processing. Best on fast hard drives.
rawData <- readMSData(pd$path, pdata = new("NAnnotatedDataFrame", pd), centroided. = TRUE,mode= "onDisk", verbose=T)


#Peak chromatographic peaks
#slowest step, this in theory can be sped up with par processing but might break
cwp <- CentWaveParam(peakwidth = c(2, 10), ppm=15, integrate = 2, noise=500, snthresh=0)

processedData <- findChromPeaks(rawData, param = cwp)

#group peaks
pdp <- PeakDensityParam(sampleGroups = pData(processedData)$group, bw=4, binSize=0.5, minSamples=1, minFraction = 0.001)

processedData <- groupChromPeaks(processedData, param = pdp)


#RT correction
medWidth <- median(chromPeaks(processedData)[, "rtmax"] -
                     chromPeaks(processedData)[, "rtmin"])

pgp <- PeakGroupsParam(minFraction = 0.45, span = 0.6,subset=which(pData(processedData)$group == "Standard"),subsetAdjust = "previous")

processedData <- adjustRtime(processedData, param = pgp) 

#regroup peaks
processedData <- groupChromPeaks(processedData, param = pdp)

#now that we have a set of features - go back and fill the blanks
processedData <- fillChromPeaks(processedData,
                                param = FillChromPeaksParam(fixedRt = medWidth))



#convert xcms object into a feature matrix so that it can be compared with MS2 data

xset <- as(filterMsLevel(processedData, msLevel = 1L), "xcmsSet")
sampclass(xset) <-pData(processedData)$group


#USE CAMERA for feature grouping into compounds 
xsa <- xsAnnotate(xset, polarity = "negative")
xsaF <- groupFWHM(xsa, sigma = 6, perfwhm = 1)
xsaC <- groupCorr(xsaF, cor_eic_th = 0.6, pval = 0.05, graphMethod = "hcs",
                  calcCiS = TRUE, calcCaS = FALSE, calcIso = FALSE)
xsaFI <- findIsotopes(xsaC, maxcharge = 2, maxiso = 3, minfrac = 0.5, ppm = 10, intval = "maxo")

xsaFA <- findAdducts(xsaFI, polarity = "negative", 
                     max_peaks = 100, multiplier = 3, ppm = 10)
edgelist <- getEdgelist(xsaFA)

xset5 <- getPeaklist(xsaFA)
xset5$feature_id <- row.names(processedData@msFeatureData$featureDefinitions)
xset5$row_number <- row.names(xset5)


#Remove Peaks Found in the Blank
nsamples <- length(pData(processedData)$sample)
colnames(xset5)[seq(13,12+nsamples,1)] <- pData(processedData)$sample

#convert NA abundances to 0

abund_df <- xset5[,seq(13,12+nsamples,1)]
abund_df[is.na(abund_df)] <- 0

xset5$Blank_max_value <- as.numeric(apply(abund_df[,which(pData(processedData)$group == "Blank")], 1, max))
xset5$max_value <- as.numeric(apply(abund_df, 1, max))

xset6 <- xset5[which(xset5$Blank_max_value < 1000),]


#Shorten Retention Time Window
xset7 <- xset6[which(xset6$rt > 25 & xset6$rt < 660),]

#remove anything with max value less than 5000
xset8 <- xset7[which(xset7$max_value > 5000),]

max_pc <- aggregate(xset8$max_value, by = list(xset8$pcgroup), max)


#xset9 <- as.data.frame(xset8 %>%
                         group_by(pcgroup) %>%
                         filter(max_value == max(max_value, na.rm=TRUE)))

xset8 <- xset8[order(xset8$max_value,decreasing = T),]



#### Match Features between MS1 and MS2

#MS level 2 
MSE_deconvoluted <-  UPLC_Results$MSE_Check[UPLC_Results$MSE_Check != ""]
MSE_deconvoluted <- paste0("./Duroia/data/MSE/msdial/MS_Dial_100_075/",list.files(path = "./Duroia/data/MSE/msdial/MS_Dial_100_075/",pattern = ".mgf"))
MSE_deconvoluted <- paste0("./Duroia/data/MSE/msdial/MS_Dial_100_05/",list.files(path = "./Duroia/data/MSE/msdial/MS_Dial_100_05//",pattern = ".mgf"))



mse_DF<- data.frame()
for(file in MSE_deconvoluted){
mse_mgf2 <- readMGF2(file)
mse_DF_single <- data.frame(scans= as.numeric(mse_mgf2$metadata$SCANS), precursorMZ= as.numeric(mse_mgf2$metadata$PEPMASS), retentionTime= as.numeric(mse_mgf2$metadata$RTINMINUTES)*60)
mse_DF_single$file <- gsub(".mgf","",basename(file))
mse_DF <- rbind(mse_DF,mse_DF_single)}

table(mse_DF$file)

head(mse_DF)


#features_spec <- featureDefinitions(processedData)
features_spec <- xset8


features_spec[features_spec$row_number %in% c(273:276),]
FT231= features_spec[which.min(abs(features_spec$mz - 231.07)), ]

x=features_spec$row_number[which.min(abs(features_spec$mz - 231.03)) ]
x=features_spec$row_number[which.min(abs(features_spec$mz - 825.3)) ]
x=features_spec$row_number[which.min(abs(features_spec$mz - 609.14)) ]
x=features_spec$row_number[which.min(abs(features_spec$mz - 301.)) ]
x=features_spec$row_number[which.min(abs(features_spec$mz - 203.08)) ]

pot_features_feature_spec <- features_spec$row_number[which(features_spec$mz > 609.1 & features_spec$mz < 609.2 )]
ms_feature <- features_spec[features_spec$row_number %in% pot_features_feature_spec,]
 x=1284
 x=1285

 "mz","rt","feature_id"
 
match_features <- function(x){
  ms_feature <- features_spec[features_spec$row_number  == x,]
  mz <- ms_feature$mz
  mzmin <- mz - 0.5
  mzmax <- mz + 0.5
  rt <- ms_feature$rt
  rtmin <- rt-45
  rtmax <- rt+45
  pot_feature <- mse_DF[mse_DF$precursorMZ >= mzmin & mse_DF$precursorMZ <= mzmax & mse_DF$retentionTime >= rtmin & mse_DF$retentionTime <= rtmax,]
  
  if(nrow(pot_feature) == 0){return(data.frame(feature_id = x, mz = mz,rt=rt,msms_scan = "no_match",file="no_match"))}
  if(nrow(pot_feature) == 1){return(data.frame(feature_id = x, mz = mz,rt=rt,msms_scan = pot_feature$scans, file=pot_feature$file))}
  if(nrow(pot_feature) > 1) {
    if(length(unique(pot_feature$file))==1){
    closest_feature <- pot_feature[which.min(abs(pot_feature$precursorMZ - mz)+abs(pot_feature$retentionTime - rt)),]
    return(data.frame(feature_id = x,mz = mz,rt=rt,msms_scan = closest_feature$scans,file=closest_feature$file))}
    if(length(unique(pot_feature$file))>1){
      closest_feature <- data.frame()
      for(mse_file in unique(pot_feature$file)){
        pot_feature2 <- pot_feature[pot_feature$file==mse_file,]
        closest_feature <- rbind(closest_feature,pot_feature2[which.min(abs(pot_feature2$precursorMZ - mz)+abs(pot_feature2$retentionTime - rt)),])
      }
      return(data.frame(feature_id = x,mz = mz,rt=rt,msms_scan = closest_feature$scans,file=closest_feature$file))}
    }
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

#### now I itterate over the above and not over the feature table because I have multiple matches...
table(ms_ms_match_results_df$msms_scan =="no_match")
table(ms_ms_match_results_df$file)


#merge onto xset8
#xset_Matched <-merge(as.data.frame(xset8),ms_ms_match_results_df,by.x=0,by.y="FeatID")
#row.names(xset_Matched) <- xset_Matched$Row.names
#names(xset_Matched)



#### Need to take a look at matches and set real filter
#### Need to only look at major TIC amounts.... and see what portion of TIC has ms/ms matches...


#write mgf file 

####### Write .mgf containing MSMS spectrum for each compound #######

PEPMASS=609.1465




MSE_deconvoluted
mse_spec_DF <- data.frame()
for(file in MSE_deconvoluted){
  mse_mgf2 <- readMGF2(file)
  spec_list <- mse_mgf2$sp
  for(spec in 1:length(spec_list)){
    mse_spec_DF <<- rbind(mse_spec_DF,data.frame(file=gsub(".mgf","",basename(file)),scan=spec,mz=unlist(spec_list[[spec]])[,1],inten=unlist(spec_list[[spec]])[,2]))
    }}
  
mse_spec_DF
  
### add that it will write a spectra for every unique aquisition from a file...

text_to_write <- c()
feature_list_to_write_mgf <- ms_ms_match_results_df[ms_ms_match_results_df$file != "no_match",][order(as.numeric(ms_ms_match_results_df[ms_ms_match_results_df$file != "no_match",]$feature_id)),]
#feature_order <- as.numeric(xset_Matched$Row.names)[order(as.numeric(xset_Matched$Row.names))]
for(k in 1:nrow(feature_list_to_write_mgf)) {
  if (k== 1) {text_to_write = c(text_to_write,paste("COM=Experimentexported on", Sys.time(),"\n",sep = "\t"))}
  #### fix 
  feature_data <- feature_list_to_write_mgf[k,c("mz","rt","feature_id","msms_scan","file")]
  #parentabund <- xset_Matched[which(xset_Matched$row_number==k),xset_Matched$file[xset_Matched$row_number==k]]
  #if(is.null(parentabund)) { parentabund = 1}
  
    #ind_spec <- mse_spec_DF[which(feature_data$file == mse_spec_DF$file & as.numeric(feature_data$msms_scan) == mse_spec_DF$scan),]

####have to add iterate over each file! otherwise it's going to pull the scan from only one file!!!
    ind_spec <- library_manager(mse_mgf2, query = c(paste0("SCANS=",feature_data$msms_scan)),logical="AND")$SELECTED$sp[[1]]
    ind_spec_peaks <- c()
    for(n in 1:nrow(ind_spec)){ind_spec_peaks <- c(ind_spec_peaks,paste(ind_spec[n,1],ind_spec[n,2],sep="\t"))}
    text_to_write = c(text_to_write, c(c("BEGIN IONS",
                                         paste("SCANS=",k,sep = ""),
                                         paste("TITLE=msLevel 2; retentionTime",round(feature_data$rt,5),"; scanNum", k,"; precMz", round(feature_data$mz,5),"; precCharge 1-",sep = " "),
                                         paste("RTINSECONDS=",round(feature_data$rt,5), sep=""),
                                         paste("RTINMINUTES=",round(feature_data$rt,5)/60, sep=""),
                                         paste("PEPMASS=",round(feature_data$mz,5), sep=""),
                                         "ION=[M-H]-",
                                         "CHARGE=1-",
                                         paste("FEATURE_ID=",feature_data$feature_id, sep=""),
                                         paste("PEAK_ID=",feature_data$pcgroup, sep="")),
                                         ind_spec_peaks,
                                         "END IONS \n"))  }}


write(text_to_write, file=here("results","./mgf_test_V4.mgf"))


### export abundance table


#row.names(xset9) <- xset9$feature_id 

xset_Matched <- xset_Matched[order(as.numeric(xset_Matched$Row.names)),]
names(xset_Matched)[1:8] <- c("Row.names","mzmed","mzmin","mzmax","rtmed","rtmin","rtmax","npeaks")	
gnps_feature_abund_table <- xset_Matched[c(1:12,seq(14,13+nsamples,1))]
write.table(gnps_feature_abund_table, file = here("results","feauter_abundance_table_all_features.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t", na = "NA")


##try what they did in the XCMS code...


## get data
featuresDef <- featureDefinitions(processedData)
featuresIntensities <- featureValues(processedData, value = "into")

## generate data table
dataTable <- merge(featuresDef, featuresIntensities, by=0, all=TRUE)
dataTable <- dataTable[, !(names(dataTable) %in% c("peakidx"))]

###XMCS has not been working

#export edge list

edgelist <- getEdgelist(xsaFA)

edgelist_sub <- edgelist[edgelist$Annotation != "", ]
write.csv(edgelist_sub, file = here("results","camera_iin_edgelist.csv"), row.names = FALSE,quote = FALSE, na = "")
