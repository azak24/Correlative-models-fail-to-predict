################################################################################
# Predicting 2023 occupancy using 2004 model under two prediction scenarios
  # 1. keeping land cover constant at 2004 level but using 2023 climate
  # 2. predicting with both 2023 climate and land cover values
# Outputs .csv files of by-cell occupancy probability
# Amanda Zak
# August 2023
################################################################################

# load data
grid <- read.csv("data/StudyAreaGrid.csv")
grid04 <- read.csv("data/StudyAreaGrid2004.csv")
stnd <- read.csv("output/standardization2023.csv",row.names = NULL)
stnd04 <- read.csv("output/standardization2004.csv")
for23_keep <- readRDS("output/forested23.RData") # will use this to subset for grid cells that are forested in 2023

# add covariate for the percent of forest(landscape) in the early-successional phase
grid$PercESpF4500 <- grid$PercES4500/grid$PercF4500
grid04$PercESpF4500 <- grid04$PercES4500/grid04$PercF4500

################################################################################
# Prediction 1:
# Predicting occupancy in 2023 using 2004 model, 2004 land data, but 2023 climate data
  # Assumes no land cover change over time- only climate change

# create data frame to hold covariates values and predictions
pred1 <- grid04[,c("OBJECTID","PercF4500","PercESpF4500")] # 2004 land cover values
sf2 <- grid[,c("sf2")] # 2023 climate values
pred1 <- cbind(pred1,sf2)
pred1$Occ <- NA

# keep only the grid cells that are in the 2023 grid
pred1 <- pred1[which(pred1$OBJECTID %in% for23_keep$OBJECTID),]

# load model parameters
a <- -0.0985304 # intercept
sfb <- 0.4419514  # snowfall beta
fb <- -0.0156494 # forest in buffer beta
esb <- 0.7554851 # es of forest in buffer beta

# put variables on standardized scale
sfStnd <- (pred1$sf2-stnd04[1,1])/stnd04[2,1] # 2023 snowfall data
fbStnd <- (pred1$PercF4500-stnd04[1,5])/stnd04[2,5]
esbStnd <- (pred1$PercESpF4500-stnd04[1,8])/stnd04[2,8]

# calculate occupancy probability
pred1$Occ <- a + sfb*sfStnd + fb*fbStnd + esb*esbStnd

# back-transform
pred1$OccReal <- exp(pred1$Occ)/(1+exp(pred1$Occ))
mean(pred1$OccReal)

# save
write.csv(pred1,"output/PredOcc_2004model2023sf.csv", row.names = FALSE)


################################################################################
# Prediction 2
# Using 2004 model and 2023 covariate values for climate and land cover
  # Allows for land cover change over time

# create data frame to hold covariates values and predictions
pred2 <- grid[,c("OBJECTID","PercF4500","PercESpF4500","sf2")]
pred2$Occ <- NA

# keep only the grid cells that are in the 2023 grid to account for forest loss
pred2 <- pred2[which(pred2$OBJECTID %in% for23_keep$OBJECTID),]

# load model parameters
a <- -0.0985304 # intercept
sfb <- 0.4419514  # snowfall beta
fb <- -0.0156494 # forest in buffer beta
esb <- 0.7554851 # es of forest in buffer beta

# put 2023 variables on 2004 standardized scale
sfStnd <- (pred2$sf2-stnd04[1,1])/stnd04[2,1]
fbStnd <- (pred2$PercF4500-stnd04[1,5])/stnd04[2,5]
esbStnd <- (pred2$PercESpF4500-stnd04[1,8])/stnd04[2,8]

# calculate occupancy probability
pred2$Occ <- a + sfb*sfStnd + fb*fbStnd + esb*esbStnd

# back-transform
pred2$OccReal <- exp(pred2$Occ)/(1+exp(pred2$Occ))

# average occupancy
mean(pred2$OccReal)

# save
write.csv(pred2,"output/PredOcc_2004model2023covs.csv", row.names = FALSE)
