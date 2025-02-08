###############################################################################
# New capt hist script that uses the most recent 2004 and 2023 land covariates
# Amanda Zak
# September 2023
###############################################################################
library(dplyr)
library(tidyr)

# Read in DNA results for each sample
dna <- read.csv("data/DNA.csv")

# Remove the data flagged "remove" (censored for various reasons)
dna <- subset(dna, dna$Note != "Remove")

# Read in list of completed surveys for each visit to each transect
visits <- read.csv("data/Surveys.csv")

# Add detection covariates to each survey
snowcover <- read.csv("data/SnowCover.csv")
snowcover$ID <- paste0(snowcover$transectid,"-",snowcover$survey_visit_no)
snowcover$scYN <- ifelse(snowcover$snow_condition == "0 - 25%",0,1)
snowcover$scPerc1 <- ifelse(snowcover$snow_condition == "25 - 75%",1,0)
snowcover$scPerc2 <- ifelse(snowcover$snow_condition == "75 - 100%",1,0)
visits <- left_join(visits,snowcover,by="ID",keep=FALSE)
# note: some surveys were incomplete (logging on transect, encountered body of water, etc.)
  # so 'visits' has fewer records than 'snowcover'

# Create snowshoe hare capt hist
sh <- matrix(data = NA, nrow = length(unique(visits$Transect)), ncol = 10)
colnames(sh) <- c("Transect","V1","V2","V3","sc1","sc2","sc3","scPerc1","scPerc2","scPerc3")
sh[,"Transect"] <- c(unique(visits$Transect))
sh <- as.data.frame(sh)
for (i in 1:nrow(sh)) { # for each transect
  for (j in 1:3) { # for each visit
    if (j %in% visits[which(visits$Transect == sh$Transect[i]),"Visit"]) { # if the transect was visited
      # if sh detected, matrix gets a 1
      # if not detected, matrix gets a 0
      sh[i,j+1] <- ifelse("Snowshoe Hare" %in% dna[which(dna$Transect == sh$Transect[i] & dna$Visit == j),"Species"],
                          "1","0")
      # fill in snow cover covs
      sh[i,j+4] <- visits[which(visits$Transect == sh$Transect[i] & visits$Visit == j),9]
      sh[i,j+7] <- paste0(visits[which(visits$Transect == sh$Transect[i] & visits$Visit == j),10]," ",
                          visits[which(visits$Transect == sh$Transect[i] & visits$Visit == j),11])
    } else {
      sh[i,j+1] <- "." # if not visited, gets placeholder
      sh[i,j+4] <- "9"
      sh[i,j+7] <- "9 9"
    }
  }
}

# Import site coordinates and covariates
coords <- read.csv("data/AllTransectCoords.csv") # coordinates for the 325 sample units (not all were surveyed)
grid <- read.csv("data/StudyAreaGrid.csv") # all grid cells in study area and their covariate values

# Add a covariate to grid: Calculate percent of forest(landscape) that is early-successional(landscape)
grid$PercESpF4500 <- grid$PercES4500/grid$PercF4500

# Join coords with transects, then join covariates
sh <- inner_join(sh,coords,by=join_by("Transect"=="TransectID"))
sh <- sh[,c("Transect","V1","V2","V3","sc1","sc2","sc3","scPerc1","scPerc2","scPerc3","CoordX","CoordY")]
colnames(sh)[11:12] <- c("Longitude","Latitude")
ch_sh <- inner_join(sh,grid,by=join_by("Latitude","Longitude"))

# Reduce columns
ch_sh <- ch_sh[,c(1:12,17:20,23:28)]

# Standardize covariates
for (i in 13:22) {
  var <- colnames(ch_sh)[i]
  eval(parse(text=paste0("ch_sh$",var,"Stnd <- (ch_sh$",var," - mean(ch_sh$",
                         var,"))/sd(ch_sh$",var,")")))
}

# Save mean & sd for back-transforming later
stnd <- matrix(data = NA, ncol = 10, nrow = 2)
colnames(stnd) <- colnames(ch_sh)[13:22]
for (i in 1:ncol(stnd)) {
  eval(parse(text=paste0("stnd[1,",i,"] <- mean(ch_sh$",
                         colnames(ch_sh)[i+12],")")))
  eval(parse(text=paste0("stnd[2,",i,"] <- sd(ch_sh$",colnames(ch_sh)[i+12],")")))
}

write.csv(stnd,"output/standardization2023.csv",row.names=FALSE)

# save capture history object for use in model fit script
saveRDS(ch_sh,"output/2023capthist.RData")

# Write final text file for use in Program MARK
fn <- "output/2023capthist.inp"
ch_sh$tag <- paste("/*", ch_sh$Transect, "*/")
ch_sh$capthist <- paste0(ch_sh$tag," ",ch_sh$V1,ch_sh$V2,ch_sh$V3," ","1"," ",
                         ch_sh$sc1," ",ch_sh$sc2," ",ch_sh$sc3," ",ch_sh$scPerc1," ",
                         ch_sh$scPerc2," ",ch_sh$scPerc3," ",ch_sh$PercESStnd," ",
                         ch_sh$PercentConStnd," ",ch_sh$PercF4500Stnd," ",
                         ch_sh$PercES4500Stnd," ",
                         ch_sh$PercESpF4500Stnd," ",ch_sh$PercHD4500Stnd,
                         " ",ch_sh$sf2Stnd," ",ch_sh$sd2Stnd," ",
                         ch_sh$dd2Stnd," ",ch_sh$bf2Stnd,";")

write(paste("/* 2023 Snowshoe Hare Capture History Northern PA",
            "\n  Transect ID",
            "\n  3 capture occasions ('.' for occasions not surveyed)",
            "\n  Covariates (standardized): ",
            "\n   Snow cover present (0 = <=25% ground cover, 1 = >25% ground cover) for each capture occasion ('9' for occasions not surveyed)",
            "\n   Percent snow cover present (00 = <25% ground cover, 10 = 25 - 75% ground cover, 01 = >75% ground cover) for each capture occasion ('9' for occasions not surveyed)",
            "\n   Proportion early-successional cover within sample unit",
            "\n   Proportion conifer cover within sample unit",
            "\n   Proportion forested cover within a square-shaped buffer around the sample unit with a side length of 9 km",
            "\n   Proportion early-successional cover within a square-shaped buffer around the sample unit with a side length of 9 km",
            "\n   Proportion of forested cover that is early-successional within a square-shaped buffer around the sample unit with a side length of 9 km",
            "\n   Proportion human-developed cover within a square-shaped buffer around the sample unit with a side length of 9 km",
            "\n   Average snowfall during the period 2013-2023, in cm",
            "\n   Average number of days per year with snow cover during the period 2013-2023",
            "\n   Average degree days (sum of maximum daily temperatures below 0 Celsius) during the period 2013-2023, in Celsius",
            "\n   Average number of days per year with a maximum temperature below 0 Celsius during the period 2013-2023",
            "*/"),
      file=fn)
write(ch_sh$capthist, fn, append=T)
