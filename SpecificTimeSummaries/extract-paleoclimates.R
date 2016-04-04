
# FUNCTIONS ---------------------------------------------------------------

ncdfToSumerizationByQM = function(ncdf, outDirectory, model, interval, brack){
    ## Function to convert all of the variables in the given netCDF4 file into 
    ## a series of rasters, summerized month and quarter data by 200 years at 500
    ## year intervals
    ##
    ## Input Variables:
    ## ncdf = the ncdf file
    ## outDirectory = where to write the rasters 
    ## model = character string of which climate model being used
    ## interval = negative integer, the number of decades to set up the year sequence, ex: -50 = 500 years
    ## brack = possitive integer, the number of decade to bracket the each interval year by for the summarization
    ##
    ## Output Variables:
    ## None, rasters written to disk
    ##
    ################################
    ## Testing Lines: for development
    # ncdf = netCDF
    # outDirectory = outDirList[[1]]
    # model = scenerioList[1]
    # interval = -50
    # brack = 10
    ################################

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
        years = var$dim[[4]]$vals
        longVals = var$dim[[1]]$vals
        latVals = var$dim[[2]]$vals
        
        ##sets up 500 year intervals 
        if(model=="ECBilt"){
            years = years + 100
            yearSeq = seq(from=0, to=-length(years), by=interval)
        }else{
            yearSeq = seq(from=0, to=-length(years), by=interval)
        }
        
        ##extracts data, 200 years at a time, 100 years around each interval
        for(i in 1:length(yearSeq)){
            ##sets year for extraction
            endYearClose = yearSeq[i]+brack
            endYearFar = yearSeq[i]-brack
            
            ##creates empty lists, to be filled below
            quartMinList = list()
            quartMaxList = list()
            monthMinList = list()
            monthMaxList = list()
            yrList = list()
            varYrList = list()
            ##gets BP year values
            BP = which(years >= endYearFar & years <= endYearClose)
            
            if(length(BP)>0){
                ##extracts data by month, and sumerizes the data by year, returns one year worth of data
                for(yr in 1:length(BP)){
                    r=BP[yr]
                    
                    ##get year data as matrix, by month
                    exData = lapply(months, function(x){ncvar_get(ncdf, var, start=c(1,1,x,r), count=c(250,140,1,1))})
                    
                    ##create raster, and corrects spatial reference location
                    ##projection from www.spatialreference.org
                    rastList = lapply(exData, function(x){raster(apply(x,1,rev), xmn=round(min(longVals)), xmx=round(max(longVals)), ymn=round(min(latVals)), ymx=round(max(latVals)), crs='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')})
                    
                    quarterRastList = lapply(quarters, function(qu){rastList[qu]})
                    ##calculates the desired values of the rasters by quarter
                    if(varName!="tmax" & varName!="tmin"){
                        quartVal = lapply(quarterRastList, function(quar){sum(stack(quar))})                       
                    }else{
                        quartVal = lapply(quarterRastList, function(quar){mean(stack(quar))})
                    }

                    ##max and min values of each cell from the quarters
                    quartMin = calc(x=stack(quartVal), fun=min)
                    quartMax = calc(x=stack(quartVal), fun=max)
                    quartMinList[yr] = quartMin
                    quartMaxList[yr] = quartMax
                    
                    ##by month
                    monthMin = calc(x=stack(rastList), fun=min)
                    monthMax = calc(x=stack(rastList), fun=max)
                    monthMinList[yr] = monthMin
                    monthMaxList[yr] = monthMax
                    
                    ##by year
                    if(varName!="tmax" & varName!="tmin"){
                        rastSum = sum(stack(rastList))
                        yrList[yr] = rastSum
                    }else{
                        rastMean = mean(stack(rastList))
                        yrList[yr] = rastMean
                    }
                    
                    ##Calculates seasonal Variations of data
                    if(varName!="tmax" & varName!="tmin" & varName!="gdd0" & varName!="gdd5"){
                        rastVar = calc(x=stack(rastList), fun=sd)/mean(stack(rastList))
                        varYrList[yr] = rastVar
                    }else{
                        rastVar = calc(x=stack(rastList), fun=sd)
                        varYrList[yr] = rastVar
                    }
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
                BPstart = years[BP[length(BP)]]
                BPend = years[BP[1]]
                
                ##sets up output directory
                modOutDir = paste(outDirectory, as.character(BPstart), "_", as.character(BPend), "/", sep="")
                
                ##creates output directory if it does not already exist
                dir.create(modOutDir, recursive=TRUE, showWarnings=FALSE)
                
                ##writes out pertinent rasters
                if(varName!="tmax"){
                    writeRaster(quartMinAve, filename=paste(modOutDir, "qt_lwr_", toupper(varName), ".tif", sep=""), overwrite=T)
                }
                
                if(varName!="tmin"){
                    writeRaster(quartMaxAve, filename=paste(modOutDir, "qt_hgr_", toupper(varName), ".tif", sep=""), overwrite=T)
                }
                
                if(varName!="tmax"){
                    writeRaster(monthMinAve, filename=paste(modOutDir, "mo_lwr_", toupper(varName), ".tif", sep=""), overwrite=T)
                }
                
                if(varName!="tmin"){
                    writeRaster(monthMaxAve, filename=paste(modOutDir, "mo_hgr_", toupper(varName), ".tif", sep=""), overwrite=T)
                }
                
                if(varName!="tmin" & varName!="tmax"){
                    writeRaster(centAve, filename=paste(modOutDir, "an_sum_", toupper(varName), ".tif", sep=""), overwrite=T)
                }else{
                    writeRaster(centAve, filename=paste(modOutDir, "an_avg_", toupper(varName), ".tif", sep=""), overwrite=T)
                    
                }
                
                if(varName!="tmax" & varName!="tmin" & varName!="gdd0" & varName!="gdd5"){
                    writeRaster(varCentAve, filename=paste(modOutDir, "an_cv_", toupper(varName), ".tif", sep=""), overwrite=T)
                }else{
                    writeRaster(varCentAve, filename=paste(modOutDir, "an_sd_", toupper(varName), ".tif", sep=""), overwrite=T)
                }
                
            } 
        }  
    }
}

createETR = function(ncdf, outDirectory, model, interval, brack){
    ## A function to calculate the etr from the available aet and pet data
    ## located in a netCDF file, uses much of the same code as ncdfToSumerizationByQM
    ## 
    ## Input Variables:
    ## ncdf = the ncdf file
    ## outDirectory = where to write the rasters 
    ## model = character string of which climate model being used
    ## interval = negative integer, the number of decades to set up the year sequence, ex: -50 = 500 years
    ## brack = possitive integer, the number of decade to bracket the each interval year by for the summarization
    ##
    ## Output Variables:
    ## None, rasters written to disk
    ##
    ################################
    ## Testing Lines: for development
    # ncdf = netCDF
    # outDirectory = outETDir[[2]]
    # model = etScenList[2]
    # interval = -50
    # brack = 10
    ################################
    
    ##needed libraries
    require(raster)
    require(ncdf4)
    
    ##extract easily "grabbed" data to be used later in the function
    petVar = ncdf$var[[1]]
    aetVar = ncdf$var[[2]]  
    varName = "etr"
    months = 1:petVar$varsize[3]  
    quarters = list(c(1,2,3), c(2,3,4), c(3,4,5), c(4,5,6), c(5,6,7), c(6,7,8), c(7,8,9), c(8,9,10), c(9,10,11), c(10,11,12), c(11,12,1), c(12,1,2))
    years = petVar$dim[[4]]$vals
    longVals = petVar$dim[[1]]$vals
    latVals = petVar$dim[[2]]$vals
    
    ##sets up 500 year intervals 
    if(model=="ECBilt"){
        years = years + 100
        yearSeq = seq(from=0, to=-length(years), by=interval)
    }else{
        yearSeq = seq(from=0, to=-length(years), by=interval)
    }
    
    ##extracts data, 200 years at a time, 100 years around each interval
    for(i in 1:length(yearSeq)){
        ##sets year for extraction
        endYearClose = yearSeq[i]+brack
        endYearFar = yearSeq[i]-brack
        
        ##creates empty lists, to be filled below
        quartMinList = list()
        quartMaxList = list()
        monthMinList = list()
        monthMaxList = list()
        yrList = list()
        varYrList = list()

        ##gets BP year values
        BP = which(years >= endYearFar & years <= endYearClose)

        if(length(BP)>0){
            ##extracts data by month, and sumerizes the data by year, returns one year worth of data
            for(yr in 1:length(BP)){
                r=BP[yr]

                ##get year data as matrix, by month
                exPetData = lapply(months, function(x){ncvar_get(ncdf, petVar, start=c(1,1,x,r), count=c(250,140,1,1))})
                exAetData = lapply(months, function(x){ncvar_get(ncdf, aetVar, start=c(1,1,x,r), count=c(250,140,1,1))})
                
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
                quartMinList[yr] = quartMin
                quartMaxList[yr] = quartMax
                
                ##by month
                monthMin = calc(x=stack(rastEtrList), fun=min)
                monthMax = calc(x=stack(rastEtrList), fun=max)
                monthMinList[yr] = monthMin
                monthMaxList[yr] = monthMax
                
                ##by year
                rastSum = mean(stack(rastEtrList), na.rm=T)
                yrList[yr] = rastSum
                
                ##Calculates seasonal Variations of data
                rastVar = calc(x=stack(rastEtrList), fun=sd)/mean(stack(rastEtrList), na.rm=T)
                varYrList[yr] = rastVar
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
            BPstart = years[BP[length(BP)]]
            BPend = years[BP[1]]

            ##sets up output directory
            modOutDir = paste(outDirectory, as.character(BPstart), "_", as.character(BPend), "/", sep="")

            ##creates output directory if it does not already exist
            dir.create(modOutDir, showWarnings=FALSE)

            ##writes out pertinent rasters
            writeRaster(quartMinAve, filename=paste(modOutDir, "qt_lwr_", toupper(varName), ".tif", sep=""), overwrite=T)
            writeRaster(quartMaxAve, filename=paste(modOutDir, "qt_hgr_", toupper(varName), ".tif", sep=""), overwrite=T)
            writeRaster(monthMinAve, filename=paste(modOutDir, "mo_lwr_", toupper(varName), ".tif", sep=""), overwrite=T)
            writeRaster(monthMaxAve, filename=paste(modOutDir, "mo_hgr_", toupper(varName), ".tif", sep=""), overwrite=T)
            writeRaster(centAve, filename=paste(modOutDir, "an_avg_", toupper(varName), ".tif", sep=""), overwrite=T)
            writeRaster(varCentAve, filename=paste(modOutDir, "an_cv_", toupper(varName), ".tif", sep=""), overwrite=T)
        } 
    }  
}

createWDEI = function(etncdf, precipncdf, outDirectory, model, interval, brack){
    ## A function to calculate the wdei from the available precip and pet data
    ## located in a netCDF file, uses much of the same code as ncdfToSumerizationByQM
    ##
    ## Input Variables:
    ## ncdf = the ncdf file
    ## outDirectory = where to write the rasters 
    ## model = character string of which climate model being used
    ## interval = negative integer, the number of decades to set up the year sequence, ex: -50 = 500 years
    ## brack = possitive integer, the number of decade to bracket the each interval year by for the summarization
    ##
    ## Output Variables:
    ## None, rasters written to disk
    ################################
    ## Testing Lines: for development
    # etncdf = netCDFlist[1]
    # precipncdf = netCDFlist[3]
    # outDirectory = outDirList[[1]]
    # tempDirectory = tempDir
    # model = scenerioList[1]
    # interval = -50
    # brack = 10
    ################################

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
    years = petVar$dim[[4]]$vals
    longVals = petVar$dim[[1]]$vals
    latVals = petVar$dim[[2]]$vals
    
    ##sets up 500 year intervals 
    if(model=="ECBilt"){
        years = years + 100
        yearSeq = seq(from=0, to=-length(years), by=interval)
    }else{
        yearSeq = seq(from=0, to=-length(years), by=interval)
    }
    
    ##extracts data, 200 years at a time, 100 years around each interval
    for(i in 1:length(yearSeq)){
        ##sets year for extraction
        endYearClose = yearSeq[i]+brack
        endYearFar = yearSeq[i]-brack
        
        ##creates empty lists, to be filled below
        quartMinList = list()
        quartMaxList = list()
        monthMinList = list()
        monthMaxList = list()
        yrList = list()
        varYrList = list()
        
        ##gets BP year values
        BP = which(years >= endYearFar & years <= endYearClose)
        
        if(length(BP)>0){
            ##extracts data by month, and sumerizes the data by year, returns one year worth of data
            for(yr in 1:length(BP)){
                r=BP[yr]

                ##get year data as matrix, by month
                etOpen = nc_open(etncdf)
                exPetData = lapply(months, function(x){ncvar_get(etOpen, petVar, start=c(1,1,x,r), count=c(250,140,1,1))})
                nc_close(etOpen)
                precipOpen = nc_open(precipncdf)
                exPrecipData = lapply(months, function(x){ncvar_get(precipOpen, precipVar, start=c(1,1,x,r), count=c(250,140,1,1))})
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
                quartMinList[yr] = quartMin
                quartMaxList[yr] = quartMax
                
                ##by month
                monthMin = calc(x=stack(rastWdeiList), fun=min)
                monthMax = calc(x=stack(rastWdeiList), fun=max)
                monthMinList[yr] = monthMin
                monthMaxList[yr] = monthMax
                
                ##by year
                rastSum = sum(stack(rastWdeiList))
                yrList[yr] = rastSum
                
                ##Calculates seasonal Variations of data
                rastVar = calc(x=stack(rastWdeiList), fun=sd)/mean(stack(rastWdeiList))
                varYrList[yr] = rastVar
                
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
                BPstart = years[BP[length(BP)]]
                BPend = years[BP[1]]
                
                ##sets up output directory
                modOutDir = paste(outDirectory, as.character(BPstart), "_", as.character(BPend), "/", sep="")
                
                ##creates output directory if it does not already exist
                dir.create(modOutDir, showWarnings=FALSE)
                
                ##writes out pertinent rasters
                writeRaster(quartMinAve, filename=paste(modOutDir, "qt_lwr_", toupper(varName), ".tif", sep=""), overwrite=T)
                writeRaster(quartMaxAve, filename=paste(modOutDir, "qt_hgr_", toupper(varName), ".tif", sep=""), overwrite=T)
                writeRaster(monthMinAve, filename=paste(modOutDir, "mo_lwr_", toupper(varName), ".tif", sep=""), overwrite=T)
                writeRaster(monthMaxAve, filename=paste(modOutDir, "mo_hgr_", toupper(varName), ".tif", sep=""), overwrite=T)
                writeRaster(centAve, filename=paste(modOutDir, "an_sum_", toupper(varName), ".tif", sep=""), overwrite=T)
                writeRaster(varCentAve, filename=paste(modOutDir, "an_cv_", toupper(varName), ".tif", sep=""), overwrite=T)
            } 
        }  
    }    
}      


# LOAD LIBRARIES ----------------------------------------------------------

library(ncdf4)
library(raster)


# SETUP DIRECTORIES AND FILES ---------------------------------------------

climMods <- c("CCSM", "ECBilt")

## Define directories
inDir = "/home/diego/Datos/Proyectos/00-PaleoCLM/Data/Climate/Climate-NetCDF"
outDir = "/home/diego/Datos/Proyectos/00-PaleoCLM/08-NatureData-ClimateData/R-project/Data-v2"

## Lists the desired .nc files to extract data from, some of the primary climate variables were not used
netCDFlist = list.files(inDir, pattern=".nc", recursive=T, full.names=T)
netCDFlist = netCDFlist[grep("ET.nc|gdd.nc|prcp.nc|temp.nc", netCDFlist)]

## Sets up output locations, based on using subdirectories of both inDir and outDir, all named the same
## Copied to have the same length as netCDF files
scenerioList = sapply(strsplit(netCDFlist, "/"), function(x, y){x[[which(x %in% y)]]}, climMods)
outDirList = lapply(scenerioList, function(scen){paste(outDir, "/", scen, "/", sep="")})


# RUN FUNCTIONS -----------------------------------------------------------

## Runs ncdfToSumerizationByQM
## extraction and summerization
for(ind in 1:length(netCDFlist)){
    netCDF = nc_open(netCDFlist[ind])
    ncdfToSumerizationByQM(netCDF, outDirList[ind], scenerioList[ind], interval=-50, brack=10)
    nc_close(netCDF)
}

##Runs createETR
etCDFlist = netCDFlist[c(1,5)]
outETDir = outDirList[c(1,5)]
etScenList = scenerioList[c(1,5)]
for(ind in 1:length(etCDFlist)){
    netCDF = nc_open(etCDFlist[ind])
    createETR(netCDF, outETDir[ind], etScenList[ind], interval=-50, brack=10)
    nc_close(netCDF)
}

##Runs createWDEI
createWDEI(netCDFlist[1], netCDFlist[3], outDirList[1], scenerioList[1], interval=-50, brack=10)
createWDEI(netCDFlist[5], netCDFlist[8], outDirList[5], scenerioList[5], interval=-50, brack=10)


# CHANGE FOLDER NAMES -----------------------------------------------------

# CCSM
dir <- paste(outDir, "/CCSM", sep="")
period <- seq(0, 22000, by=500)

yStart <- seq(10, -2190, by=-50)
yStart[1] <- 3

yEnd <- seq(-10, -2210, -50)
yEnd[45] <- -2200

oldNames <- paste(dir, "/", yStart, "_", yEnd, sep="")
newNames <- paste(dir, "/", period, "BP", sep="")

mapply(file.rename, oldNames, newNames)


# ECBilt
dir <- paste(outDir, "/ECBilt", sep="")

period <- seq(0, 21000, by=500)

yStart <- seq(10, -2090, by=-50)
yStart[1] <- -1

yEnd <- seq(-10, -2110, -50)
yEnd[43] <- -2100

oldNames <- paste(dir, "/", yStart, "_", yEnd, sep="")
newNames <- paste(dir, "/", period, "BP", sep="")

mapply(file.rename, oldNames, newNames)



