setwd("J:/Bioinformatics-R-Projects/Heatmaps/Pinebio") # or where ever you copy the files
rawdata <- read.delim2('137.txt', header = T)
nuniqdf <- rawdata[!duplicated(rawdata[ , c("GeneSymbol")]),]
raws <- nuniqdf[,1]

data <- data.matrix(nuniqdf[1:20000,2:10]) # remove the 20000 
data_scale <- t(scale(t(data)))
data_scale <- na.omit(data_scale)
heatmap(data_scale)
