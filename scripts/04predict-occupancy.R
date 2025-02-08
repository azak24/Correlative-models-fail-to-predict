################################################################################
# Using top model to predict occupancy across study area
# Outputs .csv files of by-cell occupancy probability
# Includes both 2023 and 2004 code
# Amanda Zak
# August 2023
################################################################################

### 2023 Model 

# import grid (forested grid cells only)
grid <- readRDS("output/grid_forested.RData")

# add covariate for the percent of forest(landscape) in the early-successional phase
grid$PercESpF4500 <- grid$PercES4500/grid$PercF4500

# create new matrix to hold predicted occupancy values
pred <- grid[,c(1,2,3)]
pred$Occ <- NA

# import standardization values
stnd <- read.csv("output/standardization2023.csv",row.names = NULL)

# model (these values come from top model output in Program MARK)
a <- -2.1594519 # intercept
sfb <- 0.8779362 # snowfall beta
fb <- 0.7056529 # forest in buffer beta
esb <- 1.6499307 # es of forest in buffer beta

# put variables on standardized scale
sfStnd <- (grid$sf2-stnd[1,1])/stnd[2,1]
fbStnd <- (grid$PercF4500-stnd[1,7])/stnd[2,7]
esbStnd <- (grid$PercESpF4500-stnd[1,10])/stnd[2,10]

# calculate occupancy probability
pred$Occ <- a + sfb*sfStnd + fb*fbStnd + esb*esbStnd

# back-transform
pred$OccReal <- exp(pred$Occ)/(1+exp(pred$Occ))

# average occupancy is close to what was predicted in MARK - 10%
mean(pred$OccReal)

# export as csv file to create raster in GIS
write.csv(pred,"output/PredictedOccupancy2023.csv", row.names = FALSE)

###############################################################################

# 2004 predicted occupancy (forested grid cells only)
grid04 <- readRDS("output/grid04_forested.RData")

# add covariate for the percent of forest(landscape) in the early-successional phase
grid04$PercESpF4500 <- grid04$PercES4500/grid04$PercF4500

# create new matrix to hold predicted occupancy values
pred04 <- grid04[,c(1,2,3)]
pred04$Occ <- NA

# import standardization values
stnd04 <- read.csv("output/standardization2004.csv")

# model (from Program MARK output)
a <- -0.0985304 # intercept
sfb <- 0.4419514  # snowfall beta
fb <- -0.0156494 # forest in buffer beta
esb <- 0.7554851 # es of forest in buffer beta

# put variables on standardized scale
sfStnd <- (grid04$sf1-stnd04[1,1])/stnd04[2,1]
fbStnd <- (grid04$PercF4500-stnd04[1,5])/stnd04[2,5]
esbStnd <- (grid04$PercESpF4500-stnd04[1,8])/stnd04[2,8]

# calculate occupancy probability
pred04$Occ <- a + sfb*sfStnd + fb*fbStnd + esb*esbStnd

# back-transform
pred04$OccReal <- exp(pred04$Occ)/(1+exp(pred04$Occ))

# average occupancy is close to what was predicted in MARK - 45%
mean(pred04$OccReal)

# export as csv file to create raster in GIS
write.csv(pred04,"output/PredictedOccupancy2004.csv", row.names = FALSE)

