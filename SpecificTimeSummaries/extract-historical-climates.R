
# FUNCTIONS ---------------------------------------------------------------

ncdfToHistDecade = function(ncdf, outDirectory, interval){
##    Testing Lines: for development
#     ncdf = netCDF
#     outDirectory = outDirListHist[[2]]
#     interval = 12
    
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
        
        ##creates empty lists, to be filled below
        quartMinList = list()
        quartMaxList = list()
        monthMinList = list()
        monthMaxList = list()
        yrList = list()
        varYrList = list()
        
        if(length(januarySeq)>0){
            ##extracts data by month, and sumerizes the data by year, returns one year worth of data
            for(i in 1:length(januarySeq)){
                ##sets year for extraction
                startDay = days[januarySeq[i]]
                endDay = days[januarySeq[i]+interval-1]
                
                ##gets BP year values
                BP = which(days >= startDay & days <= endDay)
                
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
                quartMinList[i] = quartMin
                quartMaxList[i] = quartMax
                
                ##by month
                monthMin = calc(x=stack(rastList), fun=min)
                monthMax = calc(x=stack(rastList), fun=max)
                monthMinList[i] = monthMin
                monthMaxList[i] = monthMax
                
                ##by year
                if(varName!="tmax" & varName!="tmin"){
                    rastCalc = sum(stack(rastList))
                    yrList[i] = rastCalc
                }else{
                    rastCalc = mean(stack(rastList))
                    yrList[i] = rastCalc
                }
                
                ##Calculates seasonal Variations of data
                if(varName!="tmax" & varName!="tmin" & varName!="gdd0" & varName!="gdd5"){
                    rastVar = calc(x=stack(rastList), fun=sd)/mean(stack(rastList))
                    varYrList[i] = rastVar
                }else{
                    rastVar = calc(x=stack(rastList), fun=sd)
                    varYrList[i] = rastVar
                }
            }
            
            ##averages the quarters together
            quartMinAve = mean(stack(quartMinList))
            quartMaxAve = mean(stack(quartMaxList))
            
            ##averages the months together
            monthMinAve = mean(stack(monthMinList))
            monthMaxAve = mean(stack(monthMaxList))
            
            ##averages the years together
            yrAve = mean(stack(yrList))
            varYrAve = mean(stack(varYrList))
            
            ##creates output directory if it does not already exist
            dir.create(outDirectory, showWarnings=TRUE, recursive=TRUE)

            ##writes out pertinent rasters
            if(varName!="tmax"){
                writeRaster(quartMinAve, filename=paste(outDirectory, "qt_lwr_", toupper(varName), ".tif", sep=""), overwrite=T)
            }
            
            if(varName!="tmin"){
                writeRaster(quartMaxAve, filename=paste(outDirectory, "qt_hgr_", toupper(varName), ".tif", sep=""), overwrite=T)
            }
            
            if(varName!="tmax"){
                writeRaster(monthMinAve, filename=paste(outDirectory, "mo_lwr_", toupper(varName), ".tif", sep=""), overwrite=T)
            }

            if(varName!="tmin"){
                writeRaster(monthMaxAve, filename=paste(outDirectory, "mo_hgr_", toupper(varName), ".tif", sep=""), overwrite=T)
            }
            
            if(varName!="tmin" & varName!="tmax"){
                writeRaster(yrAve, filename=paste(outDirectory, "an_sum_", toupper(varName), ".tif", sep=""), overwrite=T)
            }else{
                writeRaster(yrAve, filename=paste(outDirectory, "an_avg_", toupper(varName), ".tif", sep=""), overwrite=T)
            }
            
            if(varName!="tmax" & varName!="tmin" & varName!="gdd0" & varName!="gdd5"){
                writeRaster(varYrAve, filename=paste(outDirectory, "an_cv_", toupper(varName), ".tif", sep=""), overwrite=T)
            }else{
                writeRaster(varYrAve, filename=paste(outDirectory, "an_sd_", toupper(varName), ".tif", sep=""), overwrite=T)
            }
        }  
    }
}

createHistETR = function(ncdf, outDirectory, interval){
##    Testing Lines: for development
#     ncdf = netCDF
#     outDirectory = outETDirHist[[1]]
#     model = etCDFlistHist[1]
#     interval=12
    
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
    
    ##creates empty lists, to be filled below
    quartMinList = list()
    quartMaxList = list()
    monthMinList = list()
    monthMaxList = list()
    yrList = list()
    varYrList = list()
    
    if(length(januarySeq)>0){
        ##extracts data by month, and sumerizes the data by year, returns one year worth of data
        for(i in 1:length(januarySeq)){
            ##sets year for extraction
            startDay = days[januarySeq[i]]
            endDay = days[januarySeq[i]+interval-1]
            
            ##gets BP year values
            BP = which(days >= startDay & days <= endDay)
            
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
            quartMinList[i] = quartMin
            quartMaxList[i] = quartMax
            
            ##by month
            monthMin = calc(x=stack(rastEtrList), fun=min)
            monthMax = calc(x=stack(rastEtrList), fun=max)
            monthMinList[i] = monthMin
            monthMaxList[i] = monthMax
            
            ##by year
            rastSum = mean(stack(rastEtrList), na.rm=T)
            yrList[i] = rastSum
            
            ##Calculates seasonal Variations of data
            rastVar = calc(x=stack(rastEtrList), fun=sd)/mean(stack(rastEtrList), na.rm=T)
            varYrList[i] = rastVar
        }
        
        ##averages the quarters together
        quartMinAve = mean(stack(quartMinList))
        quartMaxAve = mean(stack(quartMaxList))
        
        ##averages the months together
        monthMinAve = mean(stack(monthMinList))
        monthMaxAve = mean(stack(monthMaxList))
        
        ##averages the years together
        centAve = mean(stack(yrList))
        varCentAve = mean(stack(varYrList))
        
        ##creates output directory if it does not already exist
        dir.create(outDirectory, showWarnings=FALSE)
        
        ##writes out pertinent rasters
        writeRaster(quartMinAve, filename=paste(outDirectory, "qt_lwr_", toupper(varName), ".tif", sep=""), overwrite=T)
        writeRaster(quartMaxAve, filename=paste(outDirectory, "qt_hgr_", toupper(varName), ".tif", sep=""), overwrite=T)
        writeRaster(monthMinAve, filename=paste(outDirectory, "mo_lwr_", toupper(varName), ".tif", sep=""), overwrite=T)
        writeRaster(monthMaxAve, filename=paste(outDirectory, "mo_hgr_", toupper(varName), ".tif", sep=""), overwrite=T)
        writeRaster(centAve, filename=paste(outDirectory, "an_avg_", toupper(varName), ".tif", sep=""), overwrite=T)
        writeRaster(varCentAve, filename=paste(outDirectory, "an_cv_", toupper(varName), ".tif", sep=""), overwrite=T)
    }  
}

createHistWDEI = function(etncdf, precipncdf, outDirectory, model, interval){
##    Testing Lines: for development
#     etncdf = netCDFlistHist[1]
#     precipncdf = netCDFlistHist[3]
#     outDirectory = outDirListHist[[1]]
#     model = scenerioListHist[1]
#     interval = 12
    
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
    
    ##creates empty lists, to be filled below
    quartMinList = list()
    quartMaxList = list()
    monthMinList = list()
    monthMaxList = list()
    yrList = list()
    varYrList = list()
    
    if(length(januarySeq)>0){
        
        ##extracts data
        for(i in 1:length(januarySeq)){
            ##sets year for extraction
            startDay = days[januarySeq[i]]
            endDay = days[januarySeq[i]+interval-1]
            
            ##gets BP year values
            BP = which(days >= startDay & days <= endDay)
            
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
            quartMinList[i] = quartMin
            quartMaxList[i] = quartMax
            
            ##by month
            monthMin = calc(x=stack(rastWdeiList), fun=min)
            monthMax = calc(x=stack(rastWdeiList), fun=max)
            monthMinList[i] = monthMin
            monthMaxList[i] = monthMax
            
            ##by year
            rastSum = sum(stack(rastWdeiList))
            yrList[i] = rastSum
            
            ##Calculates seasonal Variations of data
            rastVar = calc(x=stack(rastWdeiList), fun=sd)/mean(stack(rastWdeiList))
            varYrList[i] = rastVar
            
        }
        ##averages the quarters together
        quartMinAve = mean(stack(quartMinList))
        quartMaxAve = mean(stack(quartMaxList))
        
        ##averages the months together
        monthMinAve = mean(stack(monthMinList))
        monthMaxAve = mean(stack(monthMaxList))
        
        ##averages the years together
        centAve = mean(stack(yrList))
        varCentAve = mean(stack(varYrList))
        
        ##get time period start and end BP
        BPend = days[BP[length(BP)]]
        BPstart = days[BP[1]]
        
        ##creates output directory if it does not already exist
        dir.create(outDirectory, showWarnings=FALSE)
        
        ##writes out pertinent rasters
        writeRaster(quartMinAve, filename=paste(outDirectory, "qt_lwr_", toupper(varName), ".tif", sep=""), overwrite=T)
        writeRaster(quartMaxAve, filename=paste(outDirectory, "qt_hgr_", toupper(varName), ".tif", sep=""), overwrite=T)
        writeRaster(monthMinAve, filename=paste(outDirectory, "mo_lwr_", toupper(varName), ".tif", sep=""), overwrite=T)
        writeRaster(monthMaxAve, filename=paste(outDirectory, "mo_hgr_", toupper(varName), ".tif", sep=""), overwrite=T)
        writeRaster(centAve, filename=paste(outDirectory, "an_sum_", toupper(varName), ".tif", sep=""), overwrite=T)
        writeRaster(varCentAve, filename=paste(outDirectory, "an_cv_", toupper(varName), ".tif", sep=""), overwrite=T)
    }
}


# LOAD LIBRARIES ----------------------------------------------------------

library(ncdf4)
library(raster)


# SETUP DIRECTORIES AND FILES ---------------------------------------------
climMods <- c("ACCESS1-3", "CanESM2", "CESM1-CAM5", "CNRM-CM5", "CSIRO-Mk3-6-0", "GFDL-CM3", "GISS-E2-R", "HadGEM2-ES", "inmcm4", "IPSL-CM5A-MR", "MIROC5", "MRI-CGCM3")

#inDirHist = "M:/ClimateData/future/cmip5_UofWiscData/historical/"
#outDirHist = "M:/paleoCLMs/data/derived/extractedClimate/futureClims/historical/"
inDirHist = "E:/00-PaleoCLM/Data/Climate/Future/DavidData/cmip5/historical"
outDirHist = "E:/00-PaleoCLM/Data/Climate/Future/newExtraction/historical"

netCDFlistHist <- list.files(inDirHist, pattern=".nc", recursive=T, full.names=T)
netCDFlistHist <- netCDFlistHist[grep("ET.nc|gdd.nc|prcp.nc|temp.nc", netCDFlistHist)]

scenerioListHist <- sapply(strsplit(netCDFlistHist, "/"), function(x, y){x[[which(x %in% y)]]}, climMods)
outDirListHist <- lapply(scenerioListHist, function(scen){paste(outDirHist, "/", scen, "/", sep="")})


# RUN FUNCTIONS -----------------------------------------------------------

## Extraction of historical
for(ind in 1:length(netCDFlistHist)){
    netCDF <- nc_open(netCDFlistHist[ind])
    ncdfToHistDecade(netCDF, outDirListHist[[ind]],  interval=12)
    nc_close(netCDF)
}

## ETR of historical
etCDFlistHist = netCDFlistHist[grep("ET.nc", netCDFlistHist)]
outETDirHist = outDirListHist[grep("ET.nc", netCDFlistHist)]
for(ind in 1:length(etCDFlistHist)){
    netCDF <- nc_open(etCDFlistHist[ind])
    createHistETR(netCDF, outETDirHist[[ind]], interval=12)
    nc_close(netCDF)
}

## WDI of historical
##ACCESS1-3
createHistWDEI(netCDFlistHist[1], netCDFlistHist[3], outDirListHist[[1]], scenerioListHist[1], interval=12)
##CanESM2
createHistWDEI(netCDFlistHist[5], netCDFlistHist[7], outDirListHist[[5]], scenerioListHist[5], interval=12)
##CESM1-CAM5
createHistWDEI(netCDFlistHist[9], netCDFlistHist[11], outDirListHist[[9]], scenerioListHist[9], interval=12)
##CNRM-CM5
createHistWDEI(netCDFlistHist[13], netCDFlistHist[15], outDirListHist[[13]], scenerioListHist[13], interval=12)
##CSIRO-MK3-6-0
createHistWDEI(netCDFlistHist[17], netCDFlistHist[19], outDirListHist[[17]], scenerioListHist[17], interval=12)
##GFDL-cm3
createHistWDEI(netCDFlistHist[21], netCDFlistHist[23], outDirListHist[[21]], scenerioListHist[21], interval=12)
##GISS-E2-R
createHistWDEI(netCDFlistHist[25], netCDFlistHist[27], outDirListHist[[25]], scenerioListHist[25], interval=12)
##HadGEM2-ES
createHistWDEI(netCDFlistHist[29], netCDFlistHist[31], outDirListHist[[29]], scenerioListHist[29], interval=12)
##inmcm4
createHistWDEI(netCDFlistHist[33], netCDFlistHist[35], outDirListHist[[33]], scenerioListHist[33], interval=12)
##IPSL-CM5A-MR
createHistWDEI(netCDFlistHist[37], netCDFlistHist[39], outDirListHist[[37]], scenerioListHist[37], interval=12)
##MIROC5
createHistWDEI(netCDFlistHist[41], netCDFlistHist[43], outDirListHist[[41]], scenerioListHist[41], interval=12)
##MRI-CGCM3
createHistWDEI(netCDFlistHist[45], netCDFlistHist[47], outDirListHist[[45]], scenerioListHist[45], interval=12)

