################################################################################
# This script saves subsets of the 2023 and 2004 grid cells with >90% forest cover
  # so that estimates of study area occupancy probabilities do not include
  # non-forested cover

# In this script:
  # Identify forested (>90% forest) and nonforested (<90% forest) grid cells
  # Save these subsets

# Amanda Zak
# February 2024
################################################################################

# Import .csv files of area (sq m) of forest cover in each grid cell in each year
for23 <- read.csv("data/Forest2023.csv")
for04 <- read.csv("data/Forest2004.csv")

# calculate % forested cover of each grid cell
for23$percFor <- (for23$FORES_1)/129600
for04$percFor <- (for04$FORES_1)/129600

# extract >=90% forested cells
for23_keep <- subset(for23,for23$percFor >= 0.90)
for04_keep <- subset(for04,for04$percFor >= 0.90)

# save these subsets to refer to later
saveRDS(for23_keep,"output/forested23.RData")
saveRDS(for04_keep,"output/forested04.RData")

# load grid cell .csv files -these contain the covariate values for each grid cell
grid <- read.csv("data/StudyAreaGrid.csv")
grid04 <- read.csv("data/StudyAreaGrid2004.csv")

# subset the grid cells to keep only the forested ones
grid_keep <- subset(grid,grid$OBJECTID %in% for23_keep$OBJECTID_1)
grid04_keep <- subset(grid04,grid04$OBJECTID %in% for04_keep$OBJECTID_1)
saveRDS(grid_keep,"output/grid_forested.RData")
saveRDS(grid04_keep,"output/grid04_forested.RData")
