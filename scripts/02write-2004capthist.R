################################################################################
# 2004 data capture history script using newest land covariates
# Amanda Zak
# September 2023
################################################################################
library(stringr)
library(dplyr)

# Import 2004 transect coordinates
coords <- read.csv("data/2004Transects.csv", skipNul = TRUE)

# Import 2004 capture history
inp <- read.csv("data/2004CaptHist.csv", skipNul = TRUE)

# Need to match up coords$SiteID and inp$Site so that they can be joined
# First remove duplicate rows
coords <- na.omit(distinct(coords[,c("SiteID","SiteID1","Latitude1","Longitude1")]))
# Re-format SiteID and Site so that they match
coords$SiteID2 <- paste0(coords$SiteID,ifelse(coords$SiteID1 == "t","t","r"))
inp$SiteID2 <- str_sub_all(inp$Site,3)
inp$SiteID2 <- str_sub_all(inp$SiteID2,end=-4)
inp$SiteID2 <- as.character(inp$SiteID2)
# Now join
coordsJoin <- inner_join(coords,inp,by="SiteID2",relationship="many-to-one")
# However, there are 2 transects with multiple different coordinates - 81r and 147r
  # Likely because the surveyor started at different locations on different visits
# To solve, we will use the midpoint for each duplicated survey
  # Then drop the duplicates
coordsJoin[83,3] <- (coordsJoin[83,3] + coordsJoin[84,3])/2
coordsJoin[83,4] <- (coordsJoin[83,4] + coordsJoin[84,4])/2
coordsJoin <- coordsJoin[-84,]
coordsJoin[138,3] <- (coordsJoin[138,3] + coordsJoin[139,3])/2
coordsJoin[138,4] <- (coordsJoin[138,4] + coordsJoin[139,4])/2
coordsJoin <- coordsJoin[-139,]

# Import 2004 transect covariates
trans <- read.csv("data/2004TransectCovariates.csv")

# Add a covariate: Calculate percent of forest(landscape) that is early-successional(landscape)
trans$PercESpF4500 <- trans$PercES4500/trans$PercF4500

# Join the site covariates to the capture history and pare down to essential columns
captHist <- inner_join(coordsJoin,trans,by="SiteID2")
captHist <- captHist[,c(6:9,33:40,29:30)] # keep site, obs, env covariates, and lat/long

# Standardize variables and save values for back-transforming
for (i in 5:12) {
  var <- colnames(captHist)[i]
  eval(parse(text=paste0("captHist$",var,"Stnd <- (captHist$",var," - mean(captHist$",
                         var,"))/sd(captHist$",var,")")))
}

stnd04 <- matrix(data = NA, ncol = 8, nrow = 2)
colnames(stnd04) <- colnames(captHist)[5:12]
for (i in 1:ncol(stnd04)) {
  eval(parse(text=paste0("stnd04[1,",i,"] <- mean(captHist$",
                         colnames(captHist)[i+4],")")))
  eval(parse(text=paste0("stnd04[2,",i,"] <- sd(captHist$",colnames(captHist)[i+4],")")))
}

write.csv(stnd04,"output/standardization2004.csv",row.names = FALSE)

# save capture history object for use in model fit script
saveRDS(captHist,"output/2004capthist.RData")

# Write text file for use in Program MARK
fn <- "output/2004capthist.inp"
captHist$capthist <- paste0(captHist$Site," ",captHist$Obs1,captHist$Obs2,captHist$Obs3,
                            " ","1"," ",captHist$PercESStnd," ",captHist$PercES4500Stnd," ",
                            captHist$PercESpF4500Stnd," ",captHist$PercF4500Stnd," ",
                            captHist$sf1Stnd,";")

write(paste("/* 2004 Snowshoe Hare Capture History Northern PA",
            "\n  Transect ID",
            "\n  3 capture occasions ('.' for occasions not surveyed)",
            "\n  Covariates (standardized): ",
            "\n   Proportion of early-successional cover within the survey unit",
            "\n   Proportion early-successional cover within a square-shaped buffer around the sample unit with a side length of 9 km",
            "\n   Proportion of forested cover that is early-successional within a square-shaped buffer around the sample unit with a side length of 9 km",
            "\n   Proportion forested cover within a square-shaped buffer around the sample unit with a side length of 9 km",
            "\n   Average snowfall during the period 1994-2004, in cm",
            "*/"),
      file=fn)
write(captHist$capthist, fn, append=T)
