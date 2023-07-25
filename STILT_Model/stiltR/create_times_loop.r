#!/usr/bin/env Rscript

create_times <- function(stationfile, index)
{
   print(paste('Creating times for', stationfile$station[index]))
   
   #.libPaths("/home/mmolen/R/x86_64-pc-linux-gnu-library/3.4/")
   #require(base) # for julian
   #library('fields')
   #library('foreign')
   #library('methods')
   #library('methods')
   
   #script to create object "Times" w/ starting times for HF
   #output is a matrix with 4 columns: 
   #-fjul (fractional julian day since 1/1/1960)
   #-lat (deg N)
   #-lon (deg E)
   #-altitude (meters above ground)
   #
   #  $Id: create_times.r,v 1.2 2007-06-27 11:54:03 skoerner Exp $
   #---------------------------------------------------------------------------------------------------
   
   ##############
   #path to store receptor information (should be the same as in setStiltparam.r)
#   path    <- "/home/dkivits/STILT/STILT_Model/Output/" 
   path<-"/scratch-shared/dkivits/STILT/footprints/"
   outname <- paste('.RDataTimes_',stationfile$code[index],'.hf',sep='') #name for object with receptor information
   ##############
   # julian(mm,dd,yyyy)
   
   y1 <- 2021
   m1 <-    1
   d1 <-    5
   h1 <-    0
   y2 <- 2021
   m2 <-    1
   d2 <-    6
   h2 <-    23
   dt <- 1/24

   fjul <- seq(julian(m1,d1,y1)+h1/24,julian(m2,d2,y2)+h2/24,dt)

   hours <- rep(seq(1,24,1), times = length(fjul)/24)
   l <- list(jul = fjul, hours = hours)
   df <- data.frame(l)
      
   # # Filter mountain type sites
   # if (stationfile$type[index] == "2") {
   #    print("This is a mountain site")

   #    if (stationfile$utc2lst_offset[index] >= 1) {
   #       # Select only the hours between 23 and 04 local time for mountain sites
   #       dfsubset <- subset(df, hours >= abs(22 + stationfile$utc2lst_offset[index] - 24) & hours <= (5 + stationfile$utc2lst_offset[index]))
   #       fjul <- dfsubset$jul

   #    } else {
   #       # Select only the hours between 23 and 04 local time for mountain sites
   #       dfsubset <- subset(df, hours >= (22 + stationfile$utc2lst_offset[index]) & hours <= (5 + stationfile$utc2lst_offset[index]))
   #       fjul <- dfsubset$jul
   #    }
   # } else {
   #    print("This is a lowland site")

   #    # Select only the hours between 11 and 16 local time for lowland sites
   #    dfsubset <- subset(df, hours >= (10 + stationfile$utc2lst_offset[index]) & hours <= (17 + stationfile$utc2lst_offset[index]))
   #    fjul <- dfsubset$jul
   # }

   fjul <- round(fjul,6)

   lat = stationfile$lat[index]
   lon = stationfile$lon[index]
   agl = stationfile$alt[index]

   assignr(outname,cbind(fjul,lat,lon,agl),path=path,printTF=T)
}  