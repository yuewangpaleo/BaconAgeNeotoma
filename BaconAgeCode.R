###################################################################
########  The code is to 
########    1) choose qualified cores to build Bacon age-depth model;
########    2) build Bacon age-depth model using 16 combination values of accumulation rate and subdivision thickness;
########    3) choose the best Bacon age.
########  by Yue Wang, Nov 2018
#####################################################################

library(Bchron)
library(dplyr)
library(neotoma)
setwd("D:/ClimateRefugia/Age")


##################========full core==================############################################

####  Step 1: choose cores based on age controls & samples
####    criteria: 
####      1. age controls > 2;
####      2. maximum interval < 3000 years;
####      3. pollen samples > 3;

## get all the pollen records that are older than 100 BP in Neotoma 
record_all <- get_dataset(loc = c(-172.25,10.25,-48.25,79.75),ageyoung = 100,datasettype = 'pollen')  ### set young age boundary as 100 BP to exclude the modern surface samples
pollen_all <- get_download(record_all)

datasetid <- vector()
for (i in 1:length(pollen_all)) {
  
  agecontrol <- get_chroncontrol(pollen_all[[i]]$chronologies[[1]]$chronology.id[1])
  number <- nrow(agecontrol$chron.control)
  number_sample <- nrow(pollen_all[[i]]$counts)
  maxinterval <- max(diff(agecontrol$chron.control$age))
  
  if (is.na(maxinterval)==FALSE)
    if (number > 2 & maxinterval < 3000)
      if (number_sample > 3)
        datasetid <- c(datasetid,pollen_all[[i]]$dataset$dataset.meta$dataset.id)
  
}
pollen <- get_download(datasetid)
saveRDS(pollen,"pollen.RDS")

rm(list = ls(all=TRUE))


####  Step 2: use bulk-baconizing, preparing parameter files and age files
####    - bulk-baconizing folder is named by combinations of accumulation rate value and thickness value,
####      as "bulk-baconizing_accumulation_thickness". For example, if the accumulation rate is 5 cm/yr and 
####      and the thickness is 5 cm, then the folder's name is "bulk-baconizing_5_5".
####    - run bulk-baconizing.rmd from chunk 1 to chunk 8, before Bacon age model running.
####    - then run the lines below to check missing age file cores.

####    Some cores cannot write age files using bulk-baconizing.rmd because they don't have any geological
####    chronology age controls, and most of their chronology age controls are not trustable but based on 
####    estimation or interpolation. Only one trustable age controls exists in those cores. We need to 
####    delete those cores.

####    Some cores cannot write age files using bulk-baconizing.rmd because in Neotoma, these cores'
####    geological chronology file does not have any age controls, but in fact the geological chronology 
####    age controls are stored in chronology file not geological chronology file. Bulk-baconizing.rmd 
####    cannot write age files for those cores. We need to write age files for those cores using the 
####    lines below.

datasetid_all <- read.csv("./bulk-baconizing_5_5/data/params/bacon_params_v4.csv")
datasetid_noagefile <- subset(datasetid_all,is.na(suitable)==TRUE)

## delete the cores that only one trustable age control exists
datasetid_noagefile <- datasetid_noagefile[-grep("only one age",datasetid_noagefile$notes),]

## get the pollen records
pollen_noagefile <- get_download(datasetid_noagefile$datasetid)

## get age file example from successful cores
agecontrolsample <- read.csv("./bulk-baconizing_5_5/Cores/ADELINE/ADELINE.csv")

## produce age files
parameter <- expand.grid(accprior=c(5,10,20,50),thickness=c(5,10,15,20))
for (i in 1:length(pollen_noagefile)) {
  agecontrol <- get_chroncontrol(pollen_noagefile[[i]])$chron.control
  agecontrol <- agecontrol[!is.na(agecontrol$depth),]
  agecontrol <- agecontrol[order(agecontrol$depth),]
  
  agecontrolfile <- as.data.frame(matrix(NA,nrow = nrow(agecontrol),ncol = ncol(agecontrolsample)))
  colnames(agecontrolfile) <- colnames(agecontrolsample)  
  agecontrolfile$age <- agecontrol$age
  agecontrolfile$error <- (agecontrol$age.old-agecontrol$age.young)/2
  agecontrolfile$error[is.na(agecontrolfile$error)] <- 2
  agecontrolfile$depth <- agecontrol$depth
  agecontrolfile$cc[agecontrol$control.type=="Radiocarbon"] <- 1
  agecontrolfile$cc[agecontrol$control.type!="Radiocarbon"] <- 0
  agecontrolfile$labid <- seq(1,nrow(agecontrolfile))
  
  depth <- pollen_fail_full[[i]]$sample.meta$depth
  
  handle <- pollen_fail_full[[i]]$dataset$dataset.meta$collection.handle
  
  for (j in 1:16) {
    dir.create(file.path(paste0("./bulk-baconizing_",parameter$accprior[j],"_",parameter$thickness[j],"/Cores"),handle))
    write.csv(agecontrolfile,paste0("./bulk-baconizing_",parameter$accprior[j],"_",parameter$thickness[j],"/Cores/",handle,"/",handle,".csv"),row.names = FALSE)
    write.table(depth,paste0("./bulk-baconizing_",parameter$accprior[j],"_",parameter$thickness[j],"/Cores/",handle,"/",handle,"_depths.txt"),row.names = FALSE,col.names = FALSE)
  }
}

rm(list = ls(all=TRUE))


####  Step 3: run bulk-baconizing.rmd for all the 16 folders.


####  Step 4: get the best Bacon age 
datasetid_all <- read.csv("./bulk-baconizing_5_5/data/params/bacon_params_v4.csv")

##  There are cores that have varve years. We need to discard them.
datasetid <- subset(datasetid_all,success==1)

parameter <- expand.grid(accprior=c(5,10,20,50),thickness=c(5,10,15,20))
notes_selection <- vector(length = nrow(datasetid))
bestacc <- vector(length = nrow(datasetid))
bestthick <- vector(length = nrow(datasetid))
siteID <- vector(length = nrow(datasetid))
lon <- vector(length = nrow(datasetid))
lat <- vector(length = nrow(datasetid))
sitename <- vector(length = nrow(datasetid))
description <- vector(length = nrow(datasetid))

for (k in 1:nrow(datasetid)) {
  
  ## get age controls and get calendar age controls
  corename <- as.character(datasetid$handle[k])
  ## some handles are dates, different from folder names. 
  ## change those handles respectively 
  if (corename == "5/2/2013")
    corename <- "05-02-13"
  if (corename == "5/20/2014")
    corename <- "05-20-14"
  if (corename == "5/21/2002")
    corename <- "05-21-2"
  if (corename == "5/21/2004")
    corename <- "05-21-4"
  if (corename == "5/21/2005")
    corename <- "05-21-5"
  agecontrol <- read.csv(paste0("./bulk-baconizing_5_5/Cores/",corename,"/",corename,".csv"))
  ## drop radiocarbon dates that are not qualified
  agecontrol_fail <- subset(agecontrol,cc==1 & age<71 & age>46401)
  if (nrow(agecontrol_fail)!=0)
    agecontrol <- agecontrol[-as.numeric(rownames(agecontrol_fail)),]
  rownames(agecontrol) <- seq(1,nrow(agecontrol),by=1)
  ## calibrate radiocarbon age
  agecontrol_RC <- subset(agecontrol,cc==1 & age >=71 & age <=46401)
  if (nrow(agecontrol_RC) > 0) {
    agecontrol_RC_cal <- BchronCalibrate(agecontrol_RC$age, agecontrol_RC$error, calCurves = rep('intcal13',nrow(agecontrol_RC)))
    agesamples <- sampleAges(agecontrol_RC_cal)
    agecontrol[row.names(agecontrol_RC),2]  <- matrix(apply(agesamples,2,quantile,probs=0.5))
  }
  
  ## get Bacon results
  goodage <- vector()
  orderage <- vector()
  bestage <- 0
  distance <- 1000000000000000
  for (i in 1:16) {
    agefilename <- list.files(path = paste0("./bulk-baconizing_",parameter$accprior[i],"_",parameter$thickness[i],"/Cores/",corename),pattern = "/*_ages.txt")
    if (length(agefilename)>0) {
      BaconAge_i <- read.delim(paste0("./bulk-baconizing_",parameter$accprior[i],"_",parameter$thickness[i],"/Cores/",corename,"/",agefilename))[,c(1,5)]
      
      ## sometimes the top depths are extrapolated too young. We need to get rid of these outputs.
      if (min(BaconAge_i$mean) >= -70) {
        goodage <- c(goodage,i)
        
        ## criteria 1: ages are with depth order
        if (min(diff(BaconAge_i$mean))>=0){
          orderage <- c(orderage,i)
          
          ## criteria 2: residuals are smallest
          lo <- loess(BaconAge_i$mean~BaconAge_i$depth,span = 0.5)
          residual <- (agecontrol$age-predict(lo,agecontrol$depth))^2
          distance_i <- sum(residual,na.rm = TRUE)
          if (distance_i < distance) {
            distance <- distance_i
            bestage <- i
          }
        }
      }
    }
  }
  
  if (length(goodage)==0) 
    notes_selection[k] <- paste0("no good age file.")
  
  if (length(goodage)!=0 & length(orderage)==0) 
    notes_selection[k] <- paste0("no order age file.")
  
  if (length(goodage)!=0 & length(orderage)!=0) {
    ## accumulation rate is increased to 100 cm/yr if high accumulation rate is needed
    if (datasetid$acc.mean.old[k]==100) {
      bestacc[k] <- 100
      bestthick[k] <- parameter$thickness[bestage]
    }
    if (datasetid$acc.mean.old[k]!=100) {
      bestacc[k] <- parameter$accprior[bestage]
      bestthick[k] <- parameter$thickness[bestage]
    }
    dir.create(file.path("./bulk-baconizing_result/Cores_full",corename))
    filelist <- list.files(paste0("./bulk-baconizing_",parameter$accprior[bestage],"_",parameter$thickness[bestage],"/Cores/",corename,"/"))
    file.copy(from = paste0("./bulk-baconizing_",parameter$accprior[bestage],"_",parameter$thickness[bestage],"/Cores/",corename,"/",filelist),to = paste0("./bulk-baconizing_result/Cores_full/",corename,"/"))
  }
   
  pollen <- get_download(datasetid$datasetid[k])
  siteID[k] <- pollen[[1]]$dataset$site.data$site.id
  sitename[k] <- pollen[[1]]$dataset$site.data$site.name
  lon[k] <- pollen[[1]]$dataset$site.data$long
  lat[k] <- pollen[[1]]$dataset$site.data$lat
  description[k] <- pollen[[1]]$dataset$site.data$description
}

### analyze the cores that are not qualified
nogoodage <- which(notes_selection=="no good age file.")
for (k in nogoodage) {
  ## get age controls and get calendar age controls
  corename <- as.character(datasetid$handle[k])
  ## some handles are dates, different from folder names. 
  ## change those handles respectively 
  if (corename == "5/2/2013")
    corename <- "05-02-13"
  if (corename == "5/20/2014")
    corename <- "05-20-14"
  if (corename == "5/21/2002")
    corename <- "05-21-2"
  if (corename == "5/21/2004")
    corename <- "05-21-4"
  if (corename == "5/21/2005")
    corename <- "05-21-5"
  agecontrol <- read.csv(paste0("./bulk-baconizing_5_5/Cores/",corename,"/",corename,".csv"))
  ## drop radiocarbon dates that are not qualified
  agecontrol_fail <- subset(agecontrol,cc==1 & age<71 & age>46401)
  if (nrow(agecontrol_fail)!=0)
    agecontrol <- agecontrol[-as.numeric(rownames(agecontrol_fail)),]
  rownames(agecontrol) <- seq(1,nrow(agecontrol),by=1)
  ## calibrate radiocarbon age
  agecontrol_RC <- subset(agecontrol,cc==1 & age >=71 & age <=46401)
  if (nrow(agecontrol_RC) > 0) {
    agecontrol_RC_cal <- BchronCalibrate(agecontrol_RC$age, agecontrol_RC$error, calCurves = rep('intcal13',nrow(agecontrol_RC)))
    agesamples <- sampleAges(agecontrol_RC_cal)
    agecontrol[row.names(agecontrol_RC),2]  <- matrix(apply(agesamples,2,quantile,probs=0.5))
  }
  
  ## get Bacon results
  orderage <- vector()
  bestage <- 0
  distance <- 1000000000000000
  for (i in 1:16) {
    agefilename <- list.files(path = paste0("./bulk-baconizing_",parameter$accprior[i],"_",parameter$thickness[i],"/Cores/",corename),pattern = "/*_ages.txt")
    if (length(agefilename)>0) {
      BaconAge_i <- read.delim(paste0("./bulk-baconizing_",parameter$accprior[i],"_",parameter$thickness[i],"/Cores/",corename,"/",agefilename))[,c(1,5)]
      
      ## criteria 1: ages are with depth order
      if (min(diff(BaconAge_i$mean))>=0){
        orderage <- c(orderage,i)
        
        ## criteria 2: residuals are smallest
        lo <- loess(BaconAge_i$mean~BaconAge_i$depth,span = 0.5)
        residual <- (agecontrol$age-predict(lo,agecontrol$depth))^2
        distance_i <- sum(residual,na.rm = TRUE)
        if (distance_i < distance) {
          distance <- distance_i
          bestage <- i
        }
      }
    }
  }
  
  if (length(orderage)==0) 
    notes_selection[k] <- paste0("no order age file.")
  
  if (length(orderage)!=0) {
    ## accumulation rate is increased to 100 cm/yr if high accumulation rate is needed
    if (datasetid$acc.mean.old[k]==100) {
      bestacc[k] <- 100
      bestthick[k] <- parameter$thickness[bestage]
    }
    if (datasetid$acc.mean.old[k]!=100) {
      bestacc[k] <- parameter$accprior[bestage]
      bestthick[k] <- parameter$thickness[bestage]
    }
    dir.create(file.path("./bulk-baconizing_result/Cores_full",corename))
    filelist <- list.files(paste0("./bulk-baconizing_",parameter$accprior[bestage],"_",parameter$thickness[bestage],"/Cores/",corename,"/"))
    file.copy(from = paste0("./bulk-baconizing_",parameter$accprior[bestage],"_",parameter$thickness[bestage],"/Cores/",corename,"/",filelist),to = paste0("./bulk-baconizing_result/Cores_full/",corename,"/"))
  }
  
  pollen <- get_download(datasetid$datasetid[k])
  siteID[k] <- pollen[[1]]$dataset$site.data$site.id
  sitename[k] <- pollen[[1]]$dataset$site.data$site.name
  lon[k] <- pollen[[1]]$dataset$site.data$long
  lat[k] <- pollen[[1]]$dataset$site.data$lat
  description[k] <- pollen[[1]]$dataset$site.data$description
}

noorderage <- which(notes_selection=="no order age file.")
for (k in noorderage) {
  ## get age controls and get calendar age controls
  corename <- as.character(datasetid$handle[k])
  ## some handles are dates, different from folder names. 
  ## change those handles respectively 
  if (corename == "5/2/2013")
    corename <- "05-02-13"
  if (corename == "5/20/2014")
    corename <- "05-20-14"
  if (corename == "5/21/2002")
    corename <- "05-21-2"
  if (corename == "5/21/2004")
    corename <- "05-21-4"
  if (corename == "5/21/2005")
    corename <- "05-21-5"
  agecontrol <- read.csv(paste0("./bulk-baconizing_5_5/Cores/",corename,"/",corename,".csv"))
  ## drop radiocarbon dates that are not qualified
  agecontrol_fail <- subset(agecontrol,cc==1 & age<71 & age>46401)
  if (nrow(agecontrol_fail)!=0)
    agecontrol <- agecontrol[-as.numeric(rownames(agecontrol_fail)),]
  rownames(agecontrol) <- seq(1,nrow(agecontrol),by=1)
  ## calibrate radiocarbon age
  agecontrol_RC <- subset(agecontrol,cc==1 & age >=71 & age <=46401)
  if (nrow(agecontrol_RC) > 0) {
    agecontrol_RC_cal <- BchronCalibrate(agecontrol_RC$age, agecontrol_RC$error, calCurves = rep('intcal13',nrow(agecontrol_RC)))
    agesamples <- sampleAges(agecontrol_RC_cal)
    agecontrol[row.names(agecontrol_RC),2]  <- matrix(apply(agesamples,2,quantile,probs=0.5))
  }
  
  ## get Bacon results
  goodage <- 0
  bestage <- 0
  distance <- 1000000000000000
  for (i in 1:16) {
    agefilename <- list.files(path = paste0("./bulk-baconizing_",parameter$accprior[i],"_",parameter$thickness[i],"/Cores/",corename),pattern = "/*_ages.txt")
    if (length(agefilename)>0) {
      BaconAge_i <- read.delim(paste0("./bulk-baconizing_",parameter$accprior[i],"_",parameter$thickness[i],"/Cores/",corename,"/",agefilename))[,c(1,5)]
      
      ## sometimes the top depths are extrapolated too young. We need to get rid of these outputs.
      if (min(BaconAge_i$mean) >= -70) {
        goodage <- c(goodage,i)
        
        ## criteria 2: residuals are smallest
        lo <- loess(BaconAge_i$mean~BaconAge_i$depth,span = 0.5)
        residual <- (agecontrol$age-predict(lo,agecontrol$depth))^2
        distance_i <- sum(residual,na.rm = TRUE)
        if (distance_i < distance) {
          distance <- distance_i
          bestage <- i
        }
      }
    }
  }
  
  if (length(goodage)==0) 
    notes_selection[k] <- paste0("no good age file.")
  
  if (length(goodage)!=0) {
    ## accumulation rate is increased to 100 cm/yr if high accumulation rate is needed
    if (datasetid$acc.mean.old[k]==100) {
      bestacc[k] <- 100
      bestthick[k] <- parameter$thickness[bestage]
    }
    if (datasetid$acc.mean.old[k]!=100) {
      bestacc[k] <- parameter$accprior[bestage]
      bestthick[k] <- parameter$thickness[bestage]
    }
    dir.create(file.path("./bulk-baconizing_result/Cores_full",corename))
    filelist <- list.files(paste0("./bulk-baconizing_",parameter$accprior[bestage],"_",parameter$thickness[bestage],"/Cores/",corename,"/"))
    file.copy(from = paste0("./bulk-baconizing_",parameter$accprior[bestage],"_",parameter$thickness[bestage],"/Cores/",corename,"/",filelist),to = paste0("./bulk-baconizing_result/Cores_full/",corename,"/"))
  }
  
  pollen <- get_download(datasetid$datasetid[k])
  siteID[k] <- pollen[[1]]$dataset$site.data$site.id
  sitename[k] <- pollen[[1]]$dataset$site.data$site.name
  lon[k] <- pollen[[1]]$dataset$site.data$long
  lat[k] <- pollen[[1]]$dataset$site.data$lat
  description[k] <- pollen[[1]]$dataset$site.data$description
}

datasetid <- cbind(datasetid,notes_selection,bestacc,bestthick,siteID,sitename,lon,lat,description)
siteinfo_full <- datasetid[,c(20:23,2,1,18,6,19,7,8,9,3,5,17,24)]
write.csv(siteinfo_full,"./bulk-baconizing_result/SiteInfo_fullcore.csv",row.names = FALSE)

rm(list = ls(all=TRUE))

##############============End: full core==================############################################


##################===========part core==================############################################

####  Step 1: choose and build Bacon ages for cores that only some sections meet the criteria
####    criteria: 
####      age controls that 
####                  conitinuous;
####                  maximum interval < 3000 years;
####                  in the section pollen samples > 3;
####     > 3

## get all the pollen records that are older than 100 BP in Neotoma 
record_all <- get_dataset(loc = c(-172.25,10.25,-48.25,79.75),ageyoung = 100,datasettype = 'pollen')  ### set young age boundary as 100 BP to exclude the modern surface samples
pollen_all <- get_download(record_all)

datasetid <- vector()
for (i in 1:length(pollen_all)) {
  
  agecontrol <- get_chroncontrol(pollen_all[[i]]$chronologies[[1]]$chronology.id[1])
  
  ages <- agecontrol$chron.control$age[which(agecontrol$chron.control$age > 1000)]
  agediff <- diff(ages)
  locations <- which(agediff < 3000)
  result <- rle(diff(locations))
  
  number <- nrow(agecontrol$chron.control)
  number_sample <- nrow(pollen_all[[i]]$counts)
  maxinterval <- max(diff(agecontrol$chron.control$age[which(is.na(agecontrol$chron.control$age)==FALSE)]))
  
  if (maxinterval > 3000)
    if (any(result$lengths > 2 & result$values==1)==TRUE)
      if (number_sample > 3) 
        datasetid <- c(datasetid,pollen_all[[i]]$dataset$dataset.meta$dataset.id)
  
}

datasetid_full <- read.csv("./bulk-baconizing_result/SiteInfo_fullcore.csv")$datasetID
datasetid_part <- datasetid[!datasetid %in% datasetid_full]
pollen <- get_download(datasetid_part)
saveRDS(pollen,"pollen_partcore.RDS")

rm(list = ls(all=TRUE))


####  Step 2: produce top depth and bottom depth for each core
pollen <- readRDS("pollen_partcore.RDS")
siteID <- vector(length = length(pollen))
sitename <- vector(length = length(pollen))
lon <- vector(length = length(pollen))
lat <- vector(length = length(pollen))
description <- vector(length = length(pollen))
datasetID <- vector(length = length(pollen))
handle <- vector(length = length(pollen))
topdepth <- vector(length = length(pollen))
botdepth <- vector(length = length(pollen))

for (i in 1:length(pollen)) {
  
  agefile <- get_chroncontrol(pollen[[i]])$chron.control
  agefile$diff <- c(diff(agefile$age),10000)
  agegood <- subset(agefile, age > 100 & diff < 3000)
  row_name <- as.numeric(row.names(agegood))
  result_section <- split(row_name,cummax(c(1,diff(row_name))))
  r1 <- as.vector(rapply(result_section,length,how = "unlist"))
  m <- which(r1==max(r1))
  section <- c(result_section[[m]],result_section[[m]][length(result_section[[m]])]+1)
  agefile_new <- agefile[section,]

  depth <- pollen[[i]]$sample.meta$depth
  depth_new <- depth[which(depth>=agefile_new$depth[1] & depth<=agefile_new$depth[nrow(agefile_new)])]

  topdepth[i] <- depth_new[i]
  botdepth[i] <- depth_new[length(depth_new)]
  siteID[i] <- pollen[[i]]$dataset$site.data$site.id
  sitename[i] <- pollen[[i]]$dataset$site.data$site.name
  lon[i] <- pollen[[i]]$dataset$site.data$long
  lat[i] <- pollen[[i]]$dataset$site.data$lat
  description[i] <- pollen[[i]]$dataset$site.data$description
  datasetID[i] <- pollen[[i]]$dataset$dataset.meta$dataset.id
  handle[i] <- pollen[[i]]$dataset$dataset.meta$collection.handle
}
siteinfo_partcore <- as.data.frame(cbind(siteID,sitename,lon,lat,datasetID,handle,topdepth,botdepth,description))
write.csv(siteinfo_partcore,"./bulk-baconizing_result/SiteInfo_partcore.csv",row.names = FALSE)

rm(list = ls(all=TRUE))


####  Step 3: use bulk-baconizing, preparing parameter files and age files
####    - similarly to full core part, bulk-baconizing folder is named as 
####      "bulk-baconizing_partcore_accumulation_thickness". 
####    - run bulk-baconizing.rmd from chunk 1 to chunk 8, before Bacon age model running.
####    - then run the lines below to check missing age file cores.

datasetid_all <- read.csv("./bulk-baconizing_partcore_5_5/data/params/bacon_params_v4.csv")
datasetid_noagefile <- subset(datasetid_all,is.na(suitable)==TRUE)
siteinfo_partcore <- read.csv("./bulk-baconizing_result/SiteInfo_partcore.csv")
siteinfo_partcore_noagefile <- subset(siteinfo_partcore,datasetID %in% datasetid_noagefile$datasetid)

## get the pollen records
pollen_noagefile <- get_download(datasetid_noagefile$datasetid)

## get age file example from successful cores
agecontrolsample <- read.csv("./bulk-baconizing_5_5/Cores/ADELINE/ADELINE.csv")

## produce age files
parameter <- expand.grid(accprior=c(5,10,20,50),thickness=c(5,10,15,20))
for (i in 1:length(pollen_noagefile)) {
  agecontrol <- get_chroncontrol(pollen_noagefile[[i]])$chron.control
  agecontrol <- agecontrol[!is.na(agecontrol$depth),]
  agecontrol <- agecontrol[order(agecontrol$depth),]
  
  agecontrolfile <- as.data.frame(matrix(NA,nrow = nrow(agecontrol),ncol = ncol(agecontrolsample)))
  colnames(agecontrolfile) <- colnames(agecontrolsample)  
  agecontrolfile$age <- agecontrol$age
  agecontrolfile$error <- (agecontrol$age.old-agecontrol$age.young)/2
  agecontrolfile$error[is.na(agecontrolfile$error)] <- 2
  agecontrolfile$depth <- agecontrol$depth
  agecontrolfile$cc[agecontrol$control.type=="Radiocarbon"] <- 1
  agecontrolfile$cc[agecontrol$control.type!="Radiocarbon"] <- 0
  agecontrolfile$labid <- seq(1,nrow(agecontrolfile))
  
  depth <- pollen_noagefile[[i]]$sample.meta$depth
  depth_new <- depth[which(depth>=siteinfo_partcore_noagefile$topdepth[i] & depth<=siteinfo_partcore_noagefile$botdepth[i])]
  
  handle <- pollen_noagefile[[i]]$dataset$dataset.meta$collection.handle
  
  for (j in 1:16) {
    dir.create(file.path(paste0("./bulk-baconizing_partcore_",parameter$accprior[j],"_",parameter$thickness[j],"/Cores"),handle))
    write.csv(agecontrolfile,paste0("./bulk-baconizing_partcore_",parameter$accprior[j],"_",parameter$thickness[j],"/Cores/",handle,"/",handle,".csv"),row.names = FALSE)
    write.table(depth_new,paste0("./bulk-baconizing_partcore_",parameter$accprior[j],"_",parameter$thickness[j],"/Cores/",handle,"/",handle,"_depths.txt"),row.names = FALSE,col.names = FALSE)
    params <- read.csv(paste0("./bulk-baconizing_partcore_",parameter$accprior[j],"_",parameter$thickness[j],"/data/params/bacon_params_v4.csv"))
    params$suitable[is.na(params$suitable)==TRUE] <- 1
    write.csv(params,paste0("./bulk-baconizing_partcore_",parameter$accprior[j],"_",parameter$thickness[j],"/data/params/bacon_params_v4.csv"),row.names = FALSE)
  }
}

rm(list = ls(all=TRUE))


####  Step 4: run bulk-baconizing.rmd for all the 16 folders.


####  Step 5: get the best Bacon age priors
parameter <- expand.grid(accprior=c(5,10,20,50),thickness=c(5,10,15,20))
datasetid <- read.csv("./bulk-baconizing_partcore_5_5/data/params/bacon_params_v4.csv")
siteinfo_partcore_old <- read.csv("./bulk-baconizing_partcore_5_5/data/SiteInfo_partcore.csv")
siteinfo_datasetidrow <- match(datasetid$datasetid,siteinfo_partcore_old$datasetID)
siteinfo_partcore <- siteinfo_partcore_old[siteinfo_datasetidrow,]

notes_selection <- vector(length = nrow(datasetid))
bestacc <- vector(length = nrow(datasetid))
bestthick <- vector(length = nrow(datasetid))

for (k in 1:nrow(datasetid)) {
  
  ## get age controls and get calendar age controls
  corename <- as.character(datasetid$handle[k])
  agecontrol <- read.csv(paste0("./bulk-baconizing_partcore_5_5/Cores/",corename,"/",corename,".csv"))
  ## drop radiocarbon dates that are not qualified
  agecontrol_fail <- subset(agecontrol,cc==1 & age<71 & age>46401)
  if (nrow(agecontrol_fail)!=0)
    agecontrol <- agecontrol[-as.numeric(rownames(agecontrol_fail)),]
  ## get age controls that are qualified
  agecontrol_position <- which(agecontrol$depth>=siteinfo_partcore$topdepth[k] & agecontrol$depth<=siteinfo_partcore$botdepth[k])
  agecontrol_top <- max(agecontrol_position[1]-1,1)
  agecontrol_bot <- min(agecontrol_position[length(agecontrol_position)]+1,nrow(agecontrol))
  agecontrol <- agecontrol[agecontrol_top:agecontrol_bot,]
  agecontrol <- subset(agecontrol,depth>=siteinfo_partcore$topdepth[k] & depth<=siteinfo_partcore$botdepth[k])
  rownames(agecontrol) <- seq(1,nrow(agecontrol),by=1)
  agecontrol_RC <- subset(agecontrol,cc==1 & age >=71 & age <=46401)
  if (nrow(agecontrol_RC) > 0) {
    agecontrol_RC_cal <- BchronCalibrate(agecontrol_RC$age, agecontrol_RC$error, calCurves = rep('intcal13',nrow(agecontrol_RC)))
    agesamples <- sampleAges(agecontrol_RC_cal)
    agecontrol[row.names(agecontrol_RC),2]  <- matrix(apply(agesamples,2,quantile,probs=0.5))
  }
  
  ## get Bacon results
  goodage <- vector()
  orderage <- vector()
  bestage <- 0
  distance <- 1000000000000000
  for (i in 1:16) {
    agefilename <- list.files(path = paste0("./bulk-baconizing_partcore_",parameter$accprior[i],"_",parameter$thickness[i],"/Cores/",corename),pattern = "/*_ages.txt")
    if (length(agefilename)>0){
      BaconAge_i <- read.delim(paste0("./bulk-baconizing_partcore_",parameter$accprior[i],"_",parameter$thickness[i],"/Cores/",corename,"/",agefilename))[,c(1,5)]
      
      ## sometimes the top depths are extrapolated too young. We need to get rid of these outputs.
      if (min(BaconAge_i$mean) >= -70) {
        goodage <- c(goodage,i)
        
        ## criteria 1: ages are with depth order
        if (min(diff(BaconAge_i$mean))>=0){
          orderage <- c(orderage,i)
          
          ## criteria 2: residuals are smallest
          lo <- loess(BaconAge_i$mean~BaconAge_i$depth,span = 0.5)
          residual <- (agecontrol$age-predict(lo,agecontrol$depth))^2
          distance_i <- sum(residual,na.rm = TRUE)
          if (distance_i < distance) {
            distance <- distance_i
            bestage <- i
          }
        }
      }
    }
  }
  
  if (length(goodage)==0) 
    notes_selection[k] <- paste0("no good age file.")
  
  if (length(goodage)!=0 & length(orderage)==0) 
    notes_selection[k] <- paste0("no order age file.")
  
  if (length(goodage)!=0 & length(orderage)!=0) {
    ## accumulation rate is increased to 100 cm/yr if high accumulation rate is needed
    if (datasetid$acc.mean.old[k]==100) {
      bestacc[k] <- 100
      bestthick[k] <- parameter$thickness[bestage]
    }
    else {
      bestacc[k] <- parameter$accprior[bestage]
      bestthick[k] <- parameter$thickness[bestage]
    }
    
    dir.create(file.path("./bulk-baconizing_result/Cores_part",corename))
    filelist <- list.files(paste0("./bulk-baconizing_partcore_",parameter$accprior[bestage],"_",parameter$thickness[bestage],"/Cores/",corename,"/"))
    file.copy(from = paste0("./bulk-baconizing_partcore_",parameter$accprior[bestage],"_",parameter$thickness[bestage],"/Cores/",corename,"/",filelist),to = paste0("./bulk-baconizing_result/Cores_part/",corename,"/"))
  }
}

siteinfo_partcore <- cbind(siteinfo_partcore,datasetid[,c(1,3,5:9)],bestacc,bestthick)
siteinfo_partcore <- siteinfo_partcore[,c(1:5,11,6:7,18,14,19,15:17,12:13,10)]
write.csv(siteinfo_partcore,"./bulk-baconizing_result/SiteInfo_partcore.csv",row.names = FALSE)

rm(list = ls(all=TRUE))

##############============End: part core==================############################################
