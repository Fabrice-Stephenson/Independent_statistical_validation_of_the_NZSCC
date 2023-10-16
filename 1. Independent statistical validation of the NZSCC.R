##===================  INDEPENDENT STATISTICAL VALIDATION OF THE NZSCC                     ===================##

##------------------------------  Author: Fabrice Stephenson
##------------------------------  Aim: Using independent evaluation data, test the classification strength  
##------------------------------       and underlying statistical model of the NZSCC (considering the heterogeneity 
##------------------------------       in environmental coverage and statistical uncertainty)
##------------------------------  Start date : 04/01/2023
##------------------------------  End date : 08/05/2023

##=============================================================================================================##
# list of species and env preds
rm(list = ls())
gc()

####=======================        OVERVIEW                                     ==============================####

# Here we provide an example of the workflow and R code used in Stephenson, Tablada, Rowden, Bulmer, Bowden & Geange
# "Independent statistical validation of the New Zealand Seafloor Community Classification" 

# We provide independent evaluation data for demersal fish, to exemplify the validation of the NZSCC which 
# includes assessing the classification strength and underlying statistical model of the NZSCC (considering the 
# heterogeneity in environmental coverage and statistical uncertainty).


####=======================        LOAD PACKAGES AND DATA FILES                ==============================####
# load packages needed
require(vegan); require(ecodist); require(parallel); library(stringr)
library(raster); require(devEMF); require(reshape); require(ggplot2)

home <- "~" # change to file directory
setwd(home)

# load custom function to run pairwise permanova tests to test for between groups diffs
source("Pairwise_adonis.R")

# Albers equal area projection - European Petroleum Survey Group (EPSG):9191
MPI.proj <- "+proj=aea +lat_0=-40 +lon_0=175 +lat_1=-30 +lat_2=-50 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

#### NZ SCC FILES
# Load the NZSCC Classification (from Stephenson et al., 2022)
NZSCC <- raster("NZSCC/CMB.GF75_EEZ_TS-1km_Final.tif")
# plot(NZSCC)

# Environmental Coverage (from Stephenson et al., 2022)
Env.Cov <- raster("NZSCC/Env_Cov.tif")
# plot(Env.Cov)

# Uncertainty in compositional turnover (from Stephenson et al., 2022)
SD <- raster("NZSCC/SD.ras_CMB.tif")
# plot(SD)

# Bathymetry
bathy <- raster("NZSCC/Bathy.tif")
# plot(bathy)

#### ENVIRONMENTAL DISTANCE
# Turnover (environmental transformation) object for multi-level classification from Stephenson et al., 2022
load("NZSCC/CMB.R.source")

#### BIOLOGICAL DATA
# DF
# DF <- read.csv("DF/Catch_F2.csv", sep = ";", stringsAsFactors=F)
# points(DF$X, DF$Y) #Check location of projected points
#
# # Tag locations with spatial information
# DF$NZSCC <- raster::extract(NZSCC, DF[,c("X","Y")])
# DF <- DF[DF$NZSCC>0,] # remove areas outside of study area (646 rows for DF)
# DF$Env.Cov <- extract(Env.Cov, DF[,c("X","Y")])
# DF$SD <- extract(SD, DF[,c("X","Y")])

# Split Env Cov into categories based on Stephenson et al., 2021 - Presence-only habitat suitability models for vulnerable marine
# ecosystem indicator taxa in the South Pacific have reached their predictive limit
# subjectively categorised into three classes representing areas of: low (< 0.1), moderate (0.1-0.5), and high 
# (0.51-1.0) environmental coverage
# DF$Env.Cov.C <- 0
# DF <- na.omit(DF)
# DF[DF$Env.Cov <= 0.1,]$Env.Cov.C <- "Low"
# DF[DF$Env.Cov > 0.1,]$Env.Cov.C <- "Mod"
# DF[DF$Env.Cov > 0.5,]$Env.Cov.C <- "High"


#### FINAL FILES
# SAVE
# save(DF, file = "R_files/DF.source")

# LOAD 
load("DF.source")

# DF is a dataframe with data on demersal fish with 33,948 records (rows) and 15 columns:
# ID: code representing unique trawl tows; weight:  biomass of catch (unstandardised, kg); tow.area: area (km2 towed);
# Std.wgt: standardised biomass (kg.km-2); common_name: common name attributed to taxa; scientific_name: latin name
# LAT: latitude (degrees) of trawl; LON: longitude of trawl (degrees); X: projected X coordinate of trawl;
# Y: projected Y coordinate of trawl; NZSCC: NZSCC group number for each trawl location; Env.Cov: Environmetnal coverage
# calculated for the NZSCC for each trawl location ; SD: Standard deviation of the mean compositional turnover 
# calculated for the NZSCC for each trawl location ; Env.Cov.C: categorical classification of the Environmetnal coverage
# calculated for the NZSCC for each trawl location 

# 1. OVERVIEW OF DATA   -------------------------------------------------------------------------------------####
# extract metadata and only retain trawl locations (rather than species records)   
meta <- DF[,c("ID","X","Y", "NZSCC","SD","Env.Cov", "Env.Cov.C")]
meta <- unique(meta) # 4099 unique locations

# Transfrom species data from long to wide format based on trawl locations (i.e., produce a speices matrix)
DF.matrix <- cast(DF, ID ~ scientific_name, value = "weight", sum) # 269 species at 4099 unique locations
DF.matrix[is.na(DF.matrix)] <- 0 # remove NAs
ID <- DF.matrix[,1:2]
DF.matrix[DF.matrix > 0] <- 1 # transform abundance into presence-absence
DF.matrix$ID <- ID$ID; rm(ID)

# merge with metadata
DF.matrix <- merge(meta, DF.matrix, by = "ID")

# FINAL DATASET FOR FURTHER ANALYSIS
DF <- DF.matrix

# Investigate location of samples in relation to environmental coverage and uncertainty
table(DF$NZSCC) # number of samples within each group
unique(DF$NZSCC); # groups represented by the data
length(unique(DF$NZSCC)) # 49 out of 75 groups represented

# Distribution of samples acros bathymetric gradient
Bathy <- extract(bathy, DF[,c("X","Y")])*-1
Bathy[Bathy <0] <- 0
hist(Bathy, freq=T, breaks=seq(0,1400,10), xlab = "Water depth (m)")

# stats on min; mean and max environmental coverage of new samples
min(na.omit(DF$Env.Cov));mean(na.omit(DF$Env.Cov)); max(na.omit(DF$Env.Cov))

# number of samples in each category of environmental coverage
table(DF$Env.Cov.C) 

# stats on min; mean and max SD of compositional turnover of new samples
min(na.omit(DF$SD));mean(na.omit(DF$SD)); max(na.omit(DF$SD))


# 2. PRESENCE/ABSENCE --------------------------------------------------------------------------------------------####
### ANOSIM of each group - pairwise comparison between each group ###
# Only select samples for NZ SCC groups where there are more than X samples within a group 
grps.sum <- as.data.frame(table(DF.matrix$NZSCC)) # count number of samples in each group
grps.sum <- grps.sum[grps.sum[,2] >= 5, 1] # only keep groups with more than 5 samples

# remove samples for groups where there are fewer than X samples
DF.cut <- DF[DF[,"NZSCC"] %in% grps.sum,] # remove all samples where groups did not have 5 samples min

# create species only dataframe for dissimilarity calculation
DF.spe <- DF.cut[,c(8:(ncol(DF.cut)))]
DF.spe <- DF.spe[,colSums(DF.spe[,1:length(DF.spe)]) > 0]  
DF.spe <- DF.spe[rowSums(DF.spe[1:nrow(DF.spe),]) > 0,] # remove samples where no longer have all species

# extract NZ SCC groups for trimmed down dataframe
DF.class <- DF.cut[,c("NZSCC","ID")]
DF.class <- DF.class[rownames(DF.class) %in% rownames(DF.spe), ]

# distance matrix for species 
DF.dist <- vegdist(DF.spe, distance="bray", na.rm = T, upper = T)

# setup dataframe for results of ANOSIM
Summ.table <- data.frame(Num.Classes = character(), Prop.total.class=double(),Prop.inter.class.diff=double(),
                         R2 = double(), p.val = double(), stringsAsFactors=FALSE)

# ANOSIM test consistent with Stephenson et al., 2022
perm <- anosim(DF.dist, DF.class$NZSCC, permutations = 200, parallel = 5)

Summ.table[1,1] <- "NZ SCC (75 groups)"
Summ.table[,2] <- length(grps.sum)/75
  
Perm.pair <- pairwise.adonis(DF.dist, DF.class$NZSCC,perm = 99)
gc()
Summ.table[,3] <- nrow(Perm.pair[Perm.pair$p.value <0.05,]) / nrow(Perm.pair)
  
Summ.table[,4] <- perm$statistic
Summ.table[,5] <- perm$signif

#### COMPARE RESULTS TO A RANDOM CLASSIFICATION BOOTSTRAPPED 25 TIMES
Summ.table.rd <- data.frame(Num.Classes = character(), Prop.total.class=double(),Prop.inter.class.diff=double(),
                         R2 = double(), p.val = double(), stringsAsFactors=FALSE)
set.seed(1987)

for(j in 1:25){
  # generate ramdom group attribution based on number of groups within DF dataframe
  DF.class <- sample(1:49, nrow(DF.spe), replace=TRUE) # random group assignment
  perm <- anosim(DF.dist, DF.class, permutations = 100, parallel = 5)
  Summ.table.rd[j,1] <- "Randomly assigned"
  Summ.table.rd[j,2] <- length(grps.sum)/75
  Summ.table.rd[j,4] <- perm$statistic
  Summ.table.rd[j,5] <- perm$signif
  print(paste("Finished iteration: ",j, " of 25", sep = ""))
}

Summ.table.rd$Prop.inter.class.diff <- 0

setwd(home)
save(Summ.table.rd, file = "R_files/Summ.table.rd.DF.source")
# load("R_files/Summ.table.rd.source")

Summ.table.F <- rbind(Summ.table, c("Randomly assigned", colMeans(Summ.table.rd[,2:5])))
Summ.table.F
Summ.table.F[,c(2:5)] <- round(Summ.table.F[,c(2:5)],3) 

save(Summ.table.F, file = "R_files/Summ.table.F.DF.source")
load("R_files/Summ.table.F.DF.source")

# SIGNIFICANCE OF PAIRWISE GROUP DIFFERENCES
Perm.pair[c('Grp1','vs', 'Grp2')] <- str_split_fixed(Perm.pair$pairs, ' ', 3)
cor.pair <- Perm.pair[,c("Grp1", "Grp2", "p.value")]
n.grp <- 75
mat <- matrix(nrow = n.grp, ncol = n.grp)
rownames(mat) <- c(1:n.grp)
colnames(mat) <- c(1:n.grp)

for (i in 1:n.grp){
  for (j in 1:n.grp){
    v <- cor.pair[cor.pair$Grp1 == i & cor.pair$Grp2 == j,]
    
    if (is.na(v[1,"p.value"])){
      next} 
    if (v[1,"p.value"] <= 0.05){
      mat[i,j] <- "***"
    } else {
      mat[i,j] <- "ns"
    }
  }
}

write.csv(mat, file = "R_files/DF_sig_PA.csv")

### COMPARISON BETWEEN BIOLOGICAL AND ENVIRONMENTAL DISTANCE ###
# Comparison between biological and environmental distance for more than X samples (start at 5; 50 ; 200) - different colours for environmental coverage (to see if relationship)
load("R_files/DF.source")
# aggregate data into wide format (species matrix)
meta <- DF[,c("ID","X","Y", "NZSCC","SD","Env.Cov", "Env.Cov.C")]
meta <- unique(meta) # 4099 unique locations
DF.matrix <- cast(DF, ID ~ scientific_name, value = "weight", sum) # 269 species at 4099 unique locations
DF.matrix[is.na(DF.matrix)] <- 0
ID <- DF.matrix[,1:2]
DF.matrix[DF.matrix > 0] <- 1
DF.matrix$ID <- ID$ID; rm(ID)

DF.matrix <- merge(meta, DF.matrix, by = "ID")
DF <- DF.matrix

grps.sum <- as.data.frame(table(DF.matrix$NZSCC)) # count number of samples in each group
grps.sum <- grps.sum[grps.sum[,2] >= 5, 1] # only keep groups with more than 5 samples
DF.cut <- DF[DF[,"NZSCC"] %in% grps.sum,] # remove all samples where groups did not have 5 samples min

DF.spe <- DF.cut[,c(8:(ncol(DF.cut)))]
DF.spe <- DF.spe[,colSums(DF.spe[,1:length(DF.spe)]) > 0]  
DF.spe <- DF.spe[rowSums(DF.spe[1:nrow(DF.spe),]) > 0,] # remove samples where no longer have all species

DF.Env <- DF.cut[,c("X","Y")]
DF.Env <- DF.Env[rownames(DF.Env) %in% rownames(DF.spe), ]
DF.Env <- extract(CMB.R, DF.Env)

#### DISTANCES
DF.dist <- vegdist(DF.spe, distance="bray", na.rm = T, upper = T)
# Using extended bray curtis distance as per Stephenson et al., 2018 (to account for samples
# with no species in common - i.e., 100$ dissimilarity it estimates beyond this based on other samples) 
DF.dist.Ext <- stepacross(DF.dist, path = "extended", toolong = 0.8)
DF.Env.dist <- vegdist(DF.Env, distance="euclidean", na.rm = T, upper = T)

# vegan::mantel(DF.dist, DF.Env.dist, permutations = 500, parallel = 4)
mtel <- vegan::mantel(DF.dist.Ext, DF.Env.dist, na.rm = T,permutations = 500, parallel = 4)

# extracting environmental coverage and standard deviation for subset of samples
DF.cut <- DF.cut[rownames(DF.cut) %in% rownames(DF.spe), ]
env.cov <- DF.cut$Env.Cov
SD <- DF.cut$SD

# reformat distance matrices to vectors for plotting
# aa <- as.vector(DF.dist)
sp <- as.vector(DF.dist.Ext) # species
ev <- as.vector(DF.Env.dist) # turnover
ec <- as.vector(rep(env.cov, length(env.cov)/2)) # related env coverage
sd <- as.vector(rep(SD, length(SD)/2)) # related standard deviation of turnover

#new data frame with vectorized distance matrices
mat <- data.frame(sp, ev, ec, sd)
# subsample randomly for plotting
mat <- mat[sample(1:50000, replace = F),]
mat <- na.omit(mat)

# extended bray similarity vs environmental distance with colors representing env cov
mm <- ggplot(mat, aes(y = sp, x = ev)) + 
  geom_point(size = 2, aes(color = ec))+ 
  labs(x = "Environmental distance", y = "Extended Bray-Curtis Dissimilarity", colour = "Environmental coverage")+
  # geom_smooth(method = "loess", colour = "black", alpha = 0.2) +
  annotate("text", x=0.47, y=0.6, label= paste("Mantel statistic r = ", round(mtel$statistic, 2), " (sig: ***)", sep = "")) +
  theme(axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
        axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
        axis.title= element_text(face = "bold", size = 14, colour = "black"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.position = c(.75, .07),
        legend.direction="horizontal",
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 11, face = "bold"))+
  scale_colour_gradientn(colours = c("red3","orange2","dodgerblue3"), limits = c(0.0, 1.0), values = c(0.0,0.1,0.5,1.0), breaks = c(0.0,0.1,0.5,1.0), labels = c(0.0,0.1,0.5,1.0))
mm
ggsave("DF_PA_EnvCov.jpeg",dpi = 600)

# extended bray similarity vs environmental distance with colors representing SD
mm2 <- ggplot(mat, aes(y = sp, x = ev)) + 
  geom_point(size = 2, aes(color = sd))+ 
  labs(x = "Environmental distance", y = "Extended Bray-Curtis Dissimilarity", fill = "SD") + 
  # geom_smooth(method = "loess", colour = "black", alpha = 0.2) +
  annotate("text", x=0.47, y=0.6, label= paste("Mantel statistic r = ", round(mtel$statistic, 2), " (sig: ***)", sep = "")) +
  theme(axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
        axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
        axis.title= element_text(face = "bold", size = 14, colour = "black"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.position = c(.75, .07),
        legend.direction="horizontal",
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 11, face = "bold")) +
  scale_colour_gradientn(colours = c("red3","orange2","dodgerblue3"),limits = c(0.002, 0.004), breaks = c(0.002,0.003,0.004), labels = c(0.002,0.003,0.004))
mm2
ggsave("DF_PA_SD.jpeg",dpi = 600)


# 3. ABUNDANCE  -----------------------------------------------------------------------------------####
load("DF.source")

meta <- DF[,c("ID","X","Y", "NZSCC","SD","Env.Cov", "Env.Cov.C")]
meta <- unique(meta) # 4099 unique locations
DF.matrix <- cast(DF, ID ~ scientific_name, value = "Std.wgt", sum) # 269 species at 4099 unique locations
DF.matrix[is.na(DF.matrix)] <- 0
ID <- DF.matrix[,1:2]

DF.matrix <- merge(meta, DF.matrix, by = "ID")
DF <- DF.matrix

### ANOSIM of each group - pairwise comparison between each group ###
# Only select samples for NZ SCC groups where there are more than X samples within a group 
grps.sum <- as.data.frame(table(DF.matrix$NZSCC)) # count number of samples in each group
grps.sum <- grps.sum[grps.sum[,2] > 4,1] # only keep groups with more than 5 samples

# remove samples for groups where there are fewer than X samples
DF.cut <- DF[DF[,"NZSCC"] %in% grps.sum,] # remove all samples where groups did not have 5 samples min

# create species only dataframe for dissimilarity calculation
DF.spe <- DF.cut[,c(8:(ncol(DF.cut)))]
DF.spe <- DF.spe[,colSums(DF.spe[,1:length(DF.spe)]) > 0]  
DF.spe <- DF.spe[rowSums(DF.spe[1:nrow(DF.spe),]) > 0,] # remove samples where no longer have all species

# extract NZ SCC groups for trimmed down dataframe
DF.class <- DF.cut[,c("NZSCC","ID")]
DF.class <- DF.class[rownames(DF.class) %in% rownames(DF.spe), ]

# distance matrix for species 
DF.dist <- vegdist(DF.spe, distance="bray", na.rm = T, upper = T)

# setup dataframe for results of ANOSIM
Summ.table <- data.frame(Num.Classes = character(), Prop.total.class=double(),Prop.inter.class.diff=double(),
                         R2 = double(), p.val = double(), stringsAsFactors=FALSE)

# ANOSIM test consistent with Stephenson et al., 2022
perm <- anosim(DF.dist, DF.class$NZSCC, permutations = 200, parallel = 5)
# Permanova test to compare
# perm2 <- adonis2(DF.dist ~ DF.class$NZSCC, permutations = 200, parallel = 5)

Summ.table[1,1] <- "NZ SCC (75 groups)"
Summ.table[,2] <- length(grps.sum)/75

Perm.pair <- pairwise.adonis(DF.dist, DF.class$NZSCC,perm = 99)
gc()
Summ.table[,3] <- nrow(Perm.pair[Perm.pair$p.value <0.05,]) / nrow(Perm.pair)

Summ.table[,4] <- perm$statistic
Summ.table[,5] <- perm$signif

#### COMPARE RESULTS TO A RANDOM CLASSIFICATION BOOTSTRAPPED 25 TIMES
Summ.table.rd <- data.frame(Num.Classes = character(), Prop.total.class=double(),Prop.inter.class.diff=double(),
                            R2 = double(), p.val = double(), stringsAsFactors=FALSE)
set.seed(1987)

for(j in 1:25){
  # generate ramdom group attribution based on number of groups within DF dataframe
  DF.class <- sample(1:49, nrow(DF.spe), replace=TRUE) # random group assignment
  perm <- anosim(DF.dist, DF.class, permutations = 100, parallel = 5)
  Summ.table.rd[j,1] <- "Randomly assigned"
  Summ.table.rd[j,2] <- length(grps.sum)/75
  Summ.table.rd[j,4] <- perm$statistic
  Summ.table.rd[j,5] <- perm$signif
}

setwd(home)
save(Summ.table.rd, file = "R_files/Summ.table.rd.DF.Abu.source")
load("Summ.table.rd.source")

Summ.table.F <- rbind(Summ.table, c("Randomly assigned", colMeans(Summ.table.rd[,2:5])))
Summ.table.F
Summ.table.F[,c(2:5)] <- round(Summ.table.F[,c(2:5)],3) 

save(Summ.table.F, file = "Summ.table.F.DF.Abu.source")

# SIGNIFICANCE OF PAIRWISE GROUP DIFFERENCES
Perm.pair[c('Grp1','vs', 'Grp2')] <- str_split_fixed(Perm.pair$pairs, ' ', 3)
cor.pair <- Perm.pair[,c("Grp1", "Grp2", "p.value")]
n.grp <- 75
mat <- matrix(nrow = n.grp, ncol = n.grp)
rownames(mat) <- c(1:n.grp)
colnames(mat) <- c(1:n.grp)

for (i in 1:n.grp){
  for (j in 1:n.grp){
    v <- cor.pair[cor.pair$Grp1 == i & cor.pair$Grp2 == j,]
    
    if (is.na(v[1,"p.value"])){
      next} 
    if (v[1,"p.value"] <= 0.05){
      mat[i,j] <- "***"
    } else {
      mat[i,j] <- "ns"
    }
  }
}

write.csv(mat, file = "DF_sig_ABU.csv")

### COMPARISON BETWEEN BIOLOGICAL AND ENVIRONMENTAL DISTANCE ###
load("DF.source")

meta <- DF[,c("ID","X","Y", "NZSCC","SD","Env.Cov", "Env.Cov.C")]
meta <- unique(meta) # 4099 unique locations
DF.matrix <- cast(DF, ID ~ scientific_name, value = "Std.wgt", sum) # 269 species at 4099 unique locations
DF.matrix[is.na(DF.matrix)] <- 0
ID <- DF.matrix[,1:2]

DF.matrix <- merge(meta, DF.matrix, by = "ID")
# Comparison between biological and environmental distance for more than X samples (start at 5; 50 ; 200) - different colours for environmental coverage (to see if relationship)
DF <- DF.matrix

grps.sum <- as.data.frame(table(DF.matrix$NZSCC)) # count number of samples in each group
grps.sum <- grps.sum[grps.sum[,2] > 4,1] # only keep groups with more than 5 samples
DF.cut <- DF[DF[,"NZSCC"] %in% grps.sum,] # remove all samples where groups did not have 5 samples min

DF.spe <- DF.cut[,c(8:(ncol(DF.cut)))]
DF.spe <- DF.spe[,colSums(DF.spe[,1:length(DF.spe)]) > 0]  
DF.spe <- DF.spe[rowSums(DF.spe[1:nrow(DF.spe),]) > 0,] # remove samples where no longer have all species

DF.Env <- DF.cut[,c("X","Y")]
DF.Env <- DF.Env[rownames(DF.Env) %in% rownames(DF.spe), ]
DF.Env <- extract(CMB.R, DF.Env)

#### DISTANCES
DF.dist <- vegdist(DF.spe, distance="bray", na.rm = T, upper = T)
# Using extended bray curtis distance as per Stephenson et al., 2018 (to account for samples
# with no species in common - i.e., 100$ dissimilarity it estimates beyond this based on other samples) 
DF.dist.Ext <- stepacross(DF.dist, path = "extended", toolong = 0.8)
DF.Env.dist <- vegdist(DF.Env, distance="euclidean", na.rm = T, upper = T)

# vegan::mantel(DF.dist, DF.Env.dist, permutations = 500, parallel = 4)
mtel <- vegan::mantel(DF.dist.Ext, DF.Env.dist, na.rm = T,permutations = 500, parallel = 4)
mtel

# extracting environmental coverage and standard deviation for subset of samples
DF.cut <- DF.cut[rownames(DF.cut) %in% rownames(DF.spe), ]
env.cov <- DF.cut$Env.Cov
SD <- DF.cut$SD

# reformat distance matrices to vectors for plotting
# sp <- as.vector(DF.dist)
sp <- as.vector(DF.dist.Ext) # species
ev <- as.vector(DF.Env.dist) # turnover
ec <- as.vector(rep(env.cov, length(env.cov)/2)) # related env coverage
ec <- ec[1:7463316]
sd <- as.vector(rep(SD, length(SD)/2)) # related standard deviation of turnover
sd <- sd[1:7463316]

#new data frame with vectorized distance matrices
mat <- data.frame(sp, ev, ec, sd)
# subsample randomly for plotting
mat <- mat[sample(1:10000, replace = F),]
mat <- na.omit(mat)

# extended bray similarity vs environmental distance with colors representing env cov
mm <- ggplot(mat, aes(y = sp, x = ev)) + 
  geom_point(size = 2, aes(color = ec))+ 
  labs(x = "Environmental distance", y = "Extended Bray-Curtis Dissimilarity", colour = "Environmental coverage")+
  # geom_smooth(method = "loess", colour = "black", alpha = 0.2) +
  annotate("text", x=0.47, y=1, label= paste("Mantel statistic r = ", round(mtel$statistic, 2), " (sig: ***)", sep = "")) +
  theme(axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
        axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
        axis.title= element_text(face = "bold", size = 14, colour = "black"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.position = c(.75, .075),
        legend.direction="horizontal",
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 11, face = "bold"))+
  scale_colour_gradientn(colours = c("red3","orange2","dodgerblue3"), limits = c(0.0, 1.0), values = c(0.0,0.1,0.5,1.0), breaks = c(0.0,0.1,0.5,1.0), labels = c(0.0,0.1,0.5,1.0))
mm
ggsave("DF_ABU_EnvCov.jpeg",dpi = 600)

# extended bray similarity vs environmental distance with colors representing SD
mm <- ggplot(mat, aes(y = sp, x = ev)) + 
  geom_point(size = 2, aes(color = sd))+ 
  labs(x = "Environmental distance", y = "Extended Bray-Curtis Dissimilarity", fill = "SD") + 
  # geom_smooth(method = "loess", colour = "black", alpha = 0.2) +
  annotate("text", x=0.47, y=1, label= paste("Mantel statistic r = ", round(mtel$statistic, 2), " (sig: ***)", sep = "")) +
  theme(axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
        axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
        axis.title= element_text(face = "bold", size = 14, colour = "black"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.position = c(.85, .075),
        legend.direction="horizontal",
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 11, face = "bold")) +
  scale_colour_gradientn(colours = c("red3","orange2","dodgerblue3"),limits = c(0.002, 0.004), breaks = c(0.002,0.003,0.004), labels = c(0.002,0.003,0.004))
mm
ggsave("DF_ABU_SD.jpeg",dpi = 600)