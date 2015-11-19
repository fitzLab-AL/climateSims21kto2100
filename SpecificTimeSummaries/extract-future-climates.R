
# FUNCTIONS ---------------------------------------------------------------

ncdfToFutDecade = function(ncdf, outDirectory, model, interval, years){
    ## Testing Lines: for development
    # ncdf = netCDF
    # outDirectory = netCDFlistHist[[1]]
    # model = netCDFlistHist[1]
    # interval = 12
    # years=c(5,15,25,35,45,55,65,75,85,95)
    
    ##needed libraries
    require(raster)
    require(ncdf4)
    
    ##runs though each variable
    for(vn in 1:ncdf$nvars){
        ##extract easily "grabbed" data to be used later in the function
        var = ncdf$var[[vn]]
        varName = gsub(" ", "_", var$name)
        varAdd = var$addOffset
        varScale = var$scaleFact
        months = 1:var$varsize[3]
        quarters = list(c(1,2,3), c(2,3,4), c(3,4,5), c(4,5,6), c(5,6,7), c(6,7,8), c(7,8,9), c(8,9,10), c(9,10,11), c(10,11,12), c(11,12,1), c(12,1,2))
        days = var$dim[[3]]$vals
        longVals = var$dim[[1]]$vals
        latVals = var$dim[[2]]$vals
        
        januarySeq <- seq(from=1, to=length(months), by=interval)
        januarySeq <- januarySeq[years]
        
        ##extracts data, 200 years at a time, 100 years around each interval
        for(i in 1:length(januarySeq)){
            ##sets year for extraction
            startDay = days[januarySeq[i]]
            endDay = days[januarySeq[i]+interval-1]
            
            ##gets BP year values
            BP = which(days >= startDay & days <= endDay)
            
            if(length(BP)>0){
                ##get year data as matrix, by month
                exData = lapply(BP, function(x){ncvar_get(ncdf, var, start=c(1,1,x), count=c(250,140,1))})
                
                ##create raster, and corrects spatial reference location
                ##projection from www.spatialreference.org
                rastList = lapply(exData, function(x){raster(apply(x,1,rev), xmn=round(min(longVals)), xmx=round(max(longVals)), ymn=round(min(latVals)), ymx=round(max(latVals)), crs='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')})
                
                quarterRastList = lapply(quarters, function(qu){rastList[qu]})
                
                ##calculates the desired values of the rasters
                ##by quarter
                if(varName!="tmax" & varName!="tmin"){
                    quartVal = lapply(quarterRastList, function(quar){sum(stack(quar))})
                    
                }else{
                    quartVal = lapply(quarterRastList, function(quar){mean(stack(quar))})
                }
                
                ##max and min values of each cell from the quarters
                quartMin = calc(x=stack(quartVal), fun=min)
                quartMax = calc(x=stack(quartVal), fun=max)
                
                ##by month
                monthMin = calc(x=stack(rastList), fun=min)
                monthMax = calc(x=stack(rastList), fun=max)
                
                ##by year
                if(varName!="tmax" & varName!="tmin"){
                    rastCalc = sum(stack(rastList))
                }else{
                    rastCalc = mean(stack(rastList))
                }
                
                ##Calculates seasonal Variations of data
                if(varName!="tmax" & varName!="tmin" & varName!="gdd0" & varName!="gdd5"){
                    rastVar = calc(x=stack(rastList), fun=sd)/mean(stack(rastList))
                }else{
                    rastVar = calc(x=stack(rastList), fun=sd)
                }
                
                ##get time period start and end BP
                BPend = days[BP[length(BP)]]
                BPstart = days[BP[1]]
                
                ##sets up output directory
                modOutDir = paste(outDirectory, as.character(BPstart), "_", as.character(BPend), "/", sep="")
                
                ##creates output directory if it does not already exist
                dir.create(modOutDir, showWarnings=TRUE, recursive=TRUE)
                
                ##writes out pertinent rasters
                if(varName!="tmax"){
                    writeRaster(quartMin, filename=paste(modOutDir, "qt_lwr_", toupper(varName), ".tif", sep=""), overwrite=T)
                }
                
                if(varName!="tmin"){
                    writeRaster(quartMax, filename=paste(modOutDir, "qt_hgr_", toupper(varName), ".tif", sep=""), overwrite=T)    
                }
                
                if(varName!="tmax"){
                    writeRaster(monthMin, filename=paste(modOutDir, "mo_lwr_", toupper(varName), ".tif", sep=""), overwrite=T)
                }
                
                if(varName!="tmin"){
                    writeRaster(monthMax, filename=paste(modOutDir, "mo_hgr_", toupper(varName), ".tif", sep=""), overwrite=T)
                }
                
                if(varName!="tmin" & varName!="tmax"){
                    writeRaster(rastCalc, filename=paste(modOutDir, "an_sum_", toupper(varName), ".tif", sep=""), overwrite=T)
                }else{
                    writeRaster(rastCalc, filename=paste(modOutDir, "an_avg_", toupper(varName), ".tif", sep=""), overwrite=T)
                }

                if(varName!="tmax" & varName!="tmin" & varName!="gdd0" & varName!="gdd5"){
                    writeRaster(rastVar, filename=paste(modOutDir, "an_cv_", toupper(varName), ".tif", sep=""), overwrite=T)
                }else{
                    writeRaster(rastVar, filename=paste(modOutDir, "an_sd_", toupper(varName), ".tif", sep=""), overwrite=T)
                }
            } 
        }  
    }
}

createFutureETR = function(ncdf, outDirectory, model, interval, years){
    ## Testing Lines: for development
    # ncdf = netCDF
    # outDirectory = outETDir45[[1]]
    # model = etScenList45[1]
    # interval=12 
    # years=c(5,15,25,35,45,55,65,75,85,95)
    
    ##needed libraries
    require(raster)
    require(ncdf4)
    
    ##extract easily "grabbed" data to be used later in the function
    petVar = ncdf$var[[1]]
    aetVar = ncdf$var[[2]]  
    varName = "etr"
    months = 1:petVar$varsize[3]  
    quarters = list(c(1,2,3), c(2,3,4), c(3,4,5), c(4,5,6), c(5,6,7), c(6,7,8), c(7,8,9), c(8,9,10), c(9,10,11), c(10,11,12), c(11,12,1), c(12,1,2))
    days = petVar$dim[[3]]$vals
    longVals = petVar$dim[[1]]$vals
    latVals = petVar$dim[[2]]$vals
    
    januarySeq <- seq(from=1, to=length(months), by=interval)
    januarySeq <- januarySeq[years]
    
    ##extracts data, 200 years at a time, 100 years around each interval
    for(i in 1:length(januarySeq)){
        ##sets year for extraction
        startDay = days[januarySeq[i]]
        endDay = days[januarySeq[i]+interval-1]
        
        ##gets BP year values
        BP = which(days >= startDay & days <= endDay)
        
        if(length(BP)>0){
            ##get year data as matrix, by month
            exPetData = lapply(BP, function(x){ncvar_get(ncdf, petVar, start=c(1,1,x), count=c(250,140,1))})
            exAetData = lapply(BP, function(x){ncvar_get(ncdf, aetVar, start=c(1,1,x), count=c(250,140,1))})

            ##create raster, and corrects spatial reference location
            ##projection from www.spatialreference.org
            rastPetList = lapply(exPetData, function(x){raster(apply(x,1,rev), xmn=round(min(longVals)), xmx=round(max(longVals)), ymn=round(min(latVals)), ymx=round(max(latVals)), crs='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')})
            rastAetList = lapply(exAetData, function(x){raster(apply(x,1,rev), xmn=round(min(longVals)), xmx=round(max(longVals)), ymn=round(min(latVals)), ymx=round(max(latVals)), crs='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')})
            
            ##calculate ETR
            rastEtrList = mapply(function(aet, pet){aet/pet}, aet=rastAetList, pet=rastPetList)
            
            ##calculates the desired values of the rasters
            ##by quarter
            quarterRastList = lapply(quarters, function(qu){rastEtrList[qu]})
            quartVal = lapply(quarterRastList, function(quar){mean(stack(quar))})
            
            ##max and min values of each cell from the quarters
            quartMin = calc(x=stack(quartVal), fun=min)
            quartMax = calc(x=stack(quartVal), fun=max)
            
            ##by month
            monthMin = calc(x=stack(rastEtrList), fun=min)
            monthMax = calc(x=stack(rastEtrList), fun=max)
            
            ##by year
            rastSum = mean(stack(rastEtrList), na.rm=T)
            
            ##Calculates seasonal Variations of data
            rastVar = calc(x=stack(rastEtrList), fun=sd)/mean(stack(rastEtrList), na.rm=T)
            
            ##get time period start and end BP
            BPend = days[BP[length(BP)]]
            BPstart = days[BP[1]]
            
            ##sets up output directory
            modOutDir = paste(outDirectory, as.character(BPstart), "_", as.character(BPend), "/", sep="")
            
            ##creates output directory if it does not already exist
            dir.create(modOutDir, showWarnings=FALSE)
            
            ##writes out pertinent rasters
            writeRaster(quartMin, filename=paste(modOutDir, "qt_lwr_", toupper(varName), ".tif", sep=""), overwrite=T)
            writeRaster(quartMax, filename=paste(modOutDir, "qt_hgr_", toupper(varName), ".tif", sep=""), overwrite=T)
            writeRaster(monthMin, filename=paste(modOutDir, "mo_lwr_", toupper(varName), ".tif", sep=""), overwrite=T)
            writeRaster(monthMax, filename=paste(modOutDir, "mo_hgr_", toupper(varName), ".tif", sep=""), overwrite=T)
            writeRaster(rastSum, filename=paste(modOutDir, "an_avg_", toupper(varName), ".tif", sep=""), overwrite=T)
            writeRaster(rastVar, filename=paste(modOutDir, "an_cv_", toupper(varName), ".tif", sep=""), overwrite=T)
        } 
    }  
}

createFutureWDEI = function(etncdf, precipncdf, outDirectory, model, interval, years){
    ## Testing Lines: for development
    # etncdf = netCDFlist45[1]
    # precipncdf = netCDFlist45[3]
    # outDirectory = outDirList45[1]
    # model = scenerioList45[1]
    # interval = 12
    # years=c(5,15,25,35,45,55,65,75,85,95)
    
    ##needed libraries
    require(raster)
    require(ncdf4)
    
    ##opens etNCDF container, and 
    ##extract easily "grabbed" data to be used later in the function
    etOpen = nc_open(etncdf)
    petVar = etOpen$var[[1]]
    petVarName = petVar$name
    nc_close(etOpen)
    
    ##opens the percip container
    precipOpen = nc_open(precipncdf)
    precipVar = precipOpen$var[[1]]
    precipVarName = precipVar$name
    nc_close(precipOpen)
    
    varName = "wdi"
    months = 1:petVar$varsize[3]  
    quarters = list(c(1,2,3), c(2,3,4), c(3,4,5), c(4,5,6), c(5,6,7), c(6,7,8), c(7,8,9), c(8,9,10), c(9,10,11), c(10,11,12), c(11,12,1), c(12,1,2))
    days = petVar$dim[[3]]$vals
    longVals = petVar$dim[[1]]$vals
    latVals = petVar$dim[[2]]$vals
    
    januarySeq <- seq(from=1, to=length(months), by=interval)
    januarySeq <- januarySeq[years]
    
    ##extracts data
    for(i in 1:length(januarySeq)){
        ##sets year for extraction
        startDay = days[januarySeq[i]]
        endDay = days[januarySeq[i]+interval-1]
        
        ##gets BP year values
        BP = which(days >= startDay & days <= endDay)
        
        if(length(BP)>0){
            ##get year data as matrix, by month
            etOpen = nc_open(etncdf)
            exPetData = lapply(BP, function(x){ncvar_get(etOpen, petVar, start=c(1,1,x), count=c(250,140,1))})
            nc_close(etOpen)
            
            precipOpen = nc_open(precipncdf)
            exPrecipData = lapply(BP, function(x){ncvar_get(precipOpen, precipVar, start=c(1,1,x), count=c(250,140,1))})
            nc_close(precipOpen)
            
            ##create raster, and corrects spatial reference location
            ##projection from www.spatialreference.org
            rastPetList = lapply(exPetData, function(x){raster(apply(x,1,rev), xmn=round(min(longVals)), xmx=round(max(longVals)), ymn=round(min(latVals)), ymx=round(max(latVals)), crs='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')})
            rastPrecipList = lapply(exPrecipData, function(x){raster(apply(x,1,rev), xmn=round(min(longVals)), xmx=round(max(longVals)), ymn=round(min(latVals)), ymx=round(max(latVals)), crs='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')})
            
            ##calculate WDEI
            rastWdeiList = mapply(function(precip, pet){precip-pet}, precip=rastPrecipList, pet=rastPetList)
            
            ##calculates the desired values of the rasters
            ##by quarter
            quarterRastList = lapply(quarters, function(qu){rastWdeiList[qu]})
            quartVal = lapply(quarterRastList, function(quar){sum(stack(quar))})
            
            ##max and min values of each cell from the quarters
            quartMin = calc(x=stack(quartVal), fun=min)
            quartMax = calc(x=stack(quartVal), fun=max)
            
            ##by month
            monthMin = calc(x=stack(rastWdeiList), fun=min)
            monthMax = calc(x=stack(rastWdeiList), fun=max)
            
            ##by year
            rastSum = sum(stack(rastWdeiList))
            
            ##Calculates seasonal Variations of data
            rastVar = calc(x=stack(rastWdeiList), fun=sd)/mean(stack(rastWdeiList))
            
            ##get time period start and end BP
            BPend = days[BP[length(BP)]]
            BPstart = days[BP[1]]
            
            ##sets up output directory
            modOutDir = paste(outDirectory, as.character(BPstart), "_", as.character(BPend), "/", sep="")
            
            ##creates output directory if it does not already exist
            dir.create(modOutDir, showWarnings=TRUE)
            
            ##writes out pertinent rasters
            writeRaster(quartMin, filename=paste(modOutDir, "qt_lwr_", toupper(varName), ".tif", sep=""), overwrite=T)
            writeRaster(quartMax, filename=paste(modOutDir, "qt_hgr_", toupper(varName), ".tif", sep=""), overwrite=T)
            writeRaster(monthMin, filename=paste(modOutDir, "mo_lwr_", toupper(varName), ".tif", sep=""), overwrite=T)
            writeRaster(monthMax, filename=paste(modOutDir, "mo_hgr_", toupper(varName), ".tif", sep=""), overwrite=T)
            writeRaster(rastSum, filename=paste(modOutDir, "an_sum_", toupper(varName), ".tif", sep=""), overwrite=T)
            writeRaster(rastVar, filename=paste(modOutDir, "an_cv_", toupper(varName), ".tif", sep=""), overwrite=T)
        }  
    }    
}


# LOAD LIBRARIES ----------------------------------------------------------

library(ncdf4)
library(raster)


# SETUP DIRECTORIES AND FILES ---------------------------------------------

climMods <- c("ACCESS1-3", "CanESM2", "CESM1-CAM5", "CNRM-CM5", "CSIRO-Mk3-6-0", "GFDL-CM3", "GISS-E2-R", "HadGEM2-ES", "inmcm4", "IPSL-CM5A-MR", "MIROC5", "MRI-CGCM3")

## Define directories
# inDir45 = "M:/ClimateData/future/cmip5_UofWiscData/rcp45/"
# inDir85 = "M:/ClimateData/future/cmip5_UofWiscData/rcp85/"
# outDir45 = "M:/paleoCLMs/data/derived/extractedClimate/futureClims/rcp45/"
# outDir85 = "M:/paleoCLMs/data/derived/extractedClimate/futureClims/rcp85/"
inDir45 = "E:/00-PaleoCLM/Data/Climate/Future/DavidData/cmip5/rcp45"
outDir45 = "E:/00-PaleoCLM/Data/Climate/Future/newExtraction/rcp45"
inDir85 = "E:/00-PaleoCLM/Data/Climate/Future/DavidData/cmip5/rcp85"
outDir85 = "E:/00-PaleoCLM/Data/Climate/Future/newExtraction/rcp85"

## Lists the desired .nc files to extract data from, some of the primary climate variables were not used
netCDFlist45 = list.files(inDir45, pattern=".nc", recursive=T, full.names=T)
netCDFlist45 = netCDFlist45[grep("ET.nc|gdd.nc|prcp.nc|temp.nc", netCDFlist45)]
netCDFlist85 = list.files(inDir85, pattern=".nc", recursive=T, full.names=T)
netCDFlist85 = netCDFlist85[grep("ET.nc|gdd.nc|prcp.nc|temp.nc", netCDFlist85)]

## Sets up output locations, based on using subdirectories of both inDir and outDir, all named the same
## Copied to have the same length as netCDF files
scenerioList45 <- sapply(strsplit(netCDFlist45, "/"), function(x, y){x[[which(x %in% y)]]}, climMods)
scenerioList85 <- sapply(strsplit(netCDFlist85, "/"), function(x, y){x[[which(x %in% y)]]}, climMods)
outDirList45 = lapply(scenerioList45, function(scen){paste(outDir45, "/", scen, "/", sep="")})
outDirList85 = lapply(scenerioList85, function(scen){paste(outDir85, "/", scen, "/", sep="")})


# RUN FUNCTIONS -----------------------------------------------------------

## Extraction of future climates RCP4.5
for(ind in 1:length(netCDFlist45)){
    netCDF <- nc_open(netCDFlist45[ind])
    ncdfToFutDecade(netCDF, outDirList45[ind], scenerioList45[ind], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))
    nc_close(netCDF)
}

##extraction of future climates RCP8.5
for(ind in 1:length(netCDFlist85)){
    netCDF <- nc_open(netCDFlist85[ind])
    ncdfToFutDecade(netCDF, outDirList85[ind], scenerioList85[ind], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))
    nc_close(netCDF)
}

##Runs createFutureETR for future climates RCP4.5
etCDFlist45 = netCDFlist45[grep("ET.nc", netCDFlist45)]
outETDir45 = outDirList45[grep("ET.nc", netCDFlist45)]
etScenList45 = scenerioList45[grep("ET.nc", netCDFlist45)]
for(ind in 1:length(etCDFlist45)){
    netCDF = nc_open(etCDFlist45[ind])
    createFutureETR(netCDF, outETDir45[ind], etScenList45[ind], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))
    nc_close(netCDF)
}

##Runs createFutureETR for future climates RCP8.5
etCDFlist85 = netCDFlist85[grep("ET.nc", netCDFlist85)]
outETDir85 = outDirList85[grep("ET.nc", netCDFlist85)]
etScenList85 = scenerioList85[grep("ET.nc", netCDFlist85)]
for(ind in 1:length(etCDFlist85)){
    netCDF = nc_open(etCDFlist85[ind])
    createFutureETR(netCDF, outETDir85[ind], etScenList85[ind], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))
    nc_close(netCDF)
}

## Runs createWDEI for future climates RCP4.5
##ACCESS1-3_45
createFutureWDEI(netCDFlist45[1], netCDFlist45[3], outDirList45[1], scenerioList45[1], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))
##CanESM2_45
createFutureWDEI(netCDFlist45[5], netCDFlist45[7], outDirList45[5], scenerioList45[5], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))
##CESM1-CAM5_45
createFutureWDEI(netCDFlist45[9], netCDFlist45[11], outDirList45[9], scenerioList45[9], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))
##CNRM-CM5_45
createFutureWDEI(netCDFlist45[13], netCDFlist45[15], outDirList45[13], scenerioList45[13], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))
##CSIRO-MK3-6-0_45
createFutureWDEI(netCDFlist45[17], netCDFlist45[19], outDirList45[17], scenerioList45[17], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))
##GFDL-cm3_45
createFutureWDEI(netCDFlist45[21], netCDFlist45[23], outDirList45[21], scenerioList45[21], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))
##GISS-E2-R_45
createFutureWDEI(netCDFlist45[25], netCDFlist45[27], outDirList45[25], scenerioList45[25], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))
##HadGEM2-ES_45
createFutureWDEI(netCDFlist45[29], netCDFlist45[31], outDirList45[29], scenerioList45[29], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))
##inmcm4_45
createFutureWDEI(netCDFlist45[33], netCDFlist45[35], outDirList45[33], scenerioList45[33], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))
##IPSL-CM5A-MR_45
createFutureWDEI(netCDFlist45[37], netCDFlist45[39], outDirList45[37], scenerioList45[37], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))
##MIROC5_45
createFutureWDEI(netCDFlist45[41], netCDFlist45[43], outDirList45[41], scenerioList45[41], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))
##MRI-CGCM3_45
createFutureWDEI(netCDFlist45[45], netCDFlist45[47], outDirList45[45], scenerioList45[45], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))

## Runs createWDEI for future climates RCP8.5
##ACCESS1-3_85
createFutureWDEI(netCDFlist85[1], netCDFlist85[3], outDirList85[1], scenerioList85[1], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))
##CanESM2_85
createFutureWDEI(netCDFlist85[5], netCDFlist85[7], outDirList85[5], scenerioList85[5], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))
##CESM1-CAM5_85
createFutureWDEI(netCDFlist85[9], netCDFlist85[11], outDirList85[9], scenerioList85[9], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))
##CNRM-CM5_85
createFutureWDEI(netCDFlist85[13], netCDFlist85[15], outDirList85[13], scenerioList85[13], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))
##CSIRO-MK3-6-0_85
createFutureWDEI(netCDFlist85[17], netCDFlist85[19], outDirList85[17], scenerioList85[17], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))
##GFDL-cm3_85
createFutureWDEI(netCDFlist85[21], netCDFlist85[23], outDirList85[21], scenerioList85[21], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))
##GISS-E2-R_85
createFutureWDEI(netCDFlist85[25], netCDFlist85[27], outDirList85[25], scenerioList85[25], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))
##HadGEM2-ES_85
createFutureWDEI(netCDFlist85[29], netCDFlist85[31], outDirList85[29], scenerioList85[29], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))
##inmcm4_85
createFutureWDEI(netCDFlist85[33], netCDFlist85[35], outDirList85[33], scenerioList85[33], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))
##IPSL-CM5A-MR_85
createFutureWDEI(netCDFlist85[37], netCDFlist85[39], outDirList85[37], scenerioList85[37], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))
##MIROC5_85
createFutureWDEI(netCDFlist85[41], netCDFlist85[43], outDirList85[41], scenerioList85[41], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))
##MRI-CGCM3_85
createFutureWDEI(netCDFlist85[45], netCDFlist85[47], outDirList85[45], scenerioList85[45], interval=12, years=c(5,15,25,35,45,55,65,75,85,95))


# CHANGE FOLDER NAMES -----------------------------------------------------

climMods <- c("ACCESS1-3", "CanESM2", "CESM1-CAM5", "CNRM-CM5", "CSIRO-Mk3-6-0", "GFDL-CM3", "GISS-E2-R", "HadGEM2-ES", "inmcm4", "IPSL-CM5A-MR", "MIROC5", "MRI-CGCM3")

period <- seq(2010, 2100, by=10)

yStart <- c(1461, 5113, 8766, 12418, 16071, 19723, 23376, 27028, 30681, 34333)
yEnd <- c(1795, 5448, 9100, 12753, 16405, 20058, 23710, 27363, 31015, 34667)

# RCP4.5
dir <- paste(outDir45, "/", climMods, sep="")

oldNames <- unlist(lapply(dir, function(x, i, j){paste(x, "/", i, "_", j, sep="")}, yStart, yEnd))
newNames <- unlist(lapply(dir, function(x, i){paste(x, "/", i, sep="")}, period))

mapply(file.rename, oldNames, newNames)

# RCP8.5
dir <- paste(outDir85, "/", climMods, sep="")

oldNames <- unlist(lapply(dir, function(x, i, j){paste(x, "/", i, "_", j, sep="")}, yStart, yEnd))
newNames <- unlist(lapply(dir, function(x, i){paste(x, "/", i, sep="")}, period))

mapply(file.rename, oldNames, newNames)
