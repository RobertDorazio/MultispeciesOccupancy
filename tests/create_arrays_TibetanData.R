### Include packages and source code

path_wd = getwd()
path_data = sub(pattern="tests", replacement="data/TibetanPlateau", path_wd)



### OBJECTIVE:  Read data from files and create data arrays to fit a multi-species occupancy model


### Read list of data arrays from rds file

dlist = readRDS(paste(path_data, 'realdata.rds', sep='/'))

dlist$Y = dlist$Y[,, -c(4,8)]  #  omit detections of Brown bear and Pallas cat (too few detections for analysis)

## dlist$Y = dlist$Y[,, c(5, 10, 12)]  # retain only Grey wolf, stone marten, White lipped deer

## dlist$Y = dlist$Y[,, c(5, 10)]  # retain only Grey wolf, stone marten



### Save each array in the list as a flat csv file, which will be read and reshaped by Julia

write(as.vector(dlist$Y), ncolumns=1, file=paste("Y_TibetanData.csv", sep="/") )
write(as.vector(dlist$X), ncolumns=1, file=paste("X_TibetanData.csv", sep="/") )
write(as.vector(dlist$V), ncolumns=1, file=paste("V_TibetanData.csv", sep="/") )
