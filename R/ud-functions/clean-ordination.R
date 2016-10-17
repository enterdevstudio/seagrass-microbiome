## define location metadata in data frame
nmds.mat$site <- rep(NA)
nmds.mat$above.below <- rep(NA)
nmds.mat$host.env <- rep(NA)
nmds.mat$is.leaf <- rep(FALSE)
nmds.mat$is.roots <- rep(FALSE)
nmds.mat$is.sediment <- rep(FALSE)
nmds.mat$is.water <- rep(FALSE)
nmds.mat$ocean <- rep(NA)
nmds.mat$lat <- rep(NA)
nmds.mat$lon <- rep(NA)
nmds.mat$sample.type <- rep(NA)
nmds.mat$sample.loc <- rep(NA)
nmds.mat$plot.loc <- rep(NA)
nmds.mat$subsite.code <- rep(NA)

## abiotic environment
nmds.mat$temp.C <- rep(NA)
nmds.mat$salinity <- rep(NA)
# nmds.mat$day.length.h <- rep(NA)

## host community genetics
nmds.mat$genotype.richness <- rep(NA)
nmds.mat$allele.richness <- rep(NA)
# nmds.mat$inbreeding <- rep(NA)

## leaf characters
nmds.mat$mean.sheath.width <- rep(NA)
# nmds.mat$se.sheath.width <- rep(NA)
nmds.mat$mean.sheath.length <- rep(NA)
# nmds.mat$se.sheath.length <- rep(NA)
nmds.mat$longest.leaf.cm <- rep(NA)
# nmds.mat$se.longest.leaf.cm <- rep(NA)
nmds.mat$leaf.C <- rep(NA)
# nmds.mat$se.leaf.C <- rep(NA)
nmds.mat$leaf.N <- rep(NA)
# nmds.mat$se.leaf.N <- rep(NA)

## host community structure
nmds.mat$zmarina.above.biomass <- rep(NA)
# nmds.mat$se.zmarina.above.biomass <- rep(NA)
nmds.mat$zmarina.below.biomass <- rep(NA)
# nmds.mat$se.zmarina.below.biomass <- rep(NA)
nmds.mat$mean.zmarina.shoots.m2 <- rep(NA)
# nmds.mat$se.zmarina.shoots.m2 <- rep(NA)
# nmds.mat$mean.periphyton <- rep(NA)
# nmds.mat$se.periphyton <- rep(NA)
nmds.mat$mean.macroalgae <- rep(NA)
# nmds.mat$se.macroalgae <- rep(NA)

## animal community structure
nmds.mat$mean.mesograzer.b <- rep(NA)
# nmds.mat$se.mesograzer.b <- rep(NA)
# nmds.mat$mean.grazer.richness <- rep(NA)
# nmds.mat$se.grazer.richness <- rep(NA)
# nmds.mat$mean.std.epibiota <- rep(NA)
# nmds.mat$se.std.epibiota <- rep(NA)
# nmds.mat$std.mollusc.b <- rep(NA)
# nmds.mat$se.std.mollusc.b <- rep(NA)
nmds.mat$std.crustacean.b <- rep(NA)
# nmds.mat$se.std.crustacean.b <- rep(NA)

## for loop
for(f in 1:length(rownames(nmds.mat))){
  tryCatch({
  curr.site <- as.character(rownames(nmds.mat)[f])
  curr.dat <- subset(meta, SampleID == curr.site)
  curr.bio <- subset(biotic, Site == as.character(curr.dat$sub.code))
  
  ## pull in data
  ## define location metadata in data frame
  nmds.mat$site[f] <- as.character(curr.bio$Site.Code[1])
  nmds.mat$ocean[f] <- as.character(curr.bio$Ocean[1])
  nmds.mat$lat[f] <- as.numeric(as.character(curr.bio$Latitude[1]))
  nmds.mat$lon[f] <- as.numeric(as.character(curr.bio$Longitude[1]))
  nmds.mat$sample.type[f] <- as.character(curr.dat$SampleType[1])
  if(nmds.mat$sample.type[f] == 'Leaf'){
    nmds.mat$above.below[f] <- 'Above'
    nmds.mat$host.env[f] <- 'Host'
    nmds.mat$is.leaf[f] <- TRUE
  }
  if(nmds.mat$sample.type[f] == 'Roots'){
    nmds.mat$above.below[f] <- 'Below'
    nmds.mat$host.env[f] <- 'Host'
    nmds.mat$is.roots[f] <- TRUE
  }
  if(nmds.mat$sample.type[f] == 'Sediment'){
    nmds.mat$above.below[f] <- 'Below'
    nmds.mat$host.env[f] <- 'Environment'
    nmds.mat$is.sediment[f] <- TRUE
  }
  if(nmds.mat$sample.type[f] == 'Water'){
    nmds.mat$above.below[f] <- 'Above'
    nmds.mat$host.env[f] <- 'Environment'
    nmds.mat$is.water[f] <- TRUE
  }
  nmds.mat$sample.loc[f] <- as.character(curr.bio$Depth.Categorical[1])
  nmds.mat$subsite.code[f] <- as.character(curr.dat$SubsiteNumber[1])
  nmds.mat$plot.loc[f] <- as.character(curr.dat$PlotLocation[1])
  
  ## abiotic environment
  nmds.mat$temp.C[f] <- as.numeric(as.character(curr.bio$Temperature.C[1]))
  nmds.mat$salinity[f] <- as.numeric(as.character(curr.bio$Salinity.ppt[1]))
#   nmds.mat$day.length.h[f] <- as.numeric(as.character(curr.bio$Day.length.hours[1]))
  
  ## host community genetics
  nmds.mat$genotype.richness[f] <- as.numeric(as.character(curr.bio$GenotypicRichness[1]))
  nmds.mat$allele.richness[f] <- as.numeric(as.character(curr.bio$AllelicRichness[1]))
#   nmds.mat$inbreeding[f] <- as.numeric(as.character(curr.bio$Inbreeding[1]))
  
  ## leaf characters
  nmds.mat$mean.sheath.width[f] <- as.numeric(as.character(curr.bio$Sheath.Width.cm.[1]))
#   nmds.mat$se.sheath.width[f] <- as.numeric(as.character(curr.bio$SE.Sheath.Width.cm.[1]))
  nmds.mat$mean.sheath.length[f] <- as.numeric(as.character(curr.bio$Shealth.Length.cm.[1]))
#   nmds.mat$se.sheath.length[f] <- as.numeric(as.character(curr.bio$SE.Shealth.Length.cm.[1]))
  nmds.mat$longest.leaf.cm[f] <- as.numeric(as.character(curr.bio$Longest.Leaf.Length.cm.[1]))
#   nmds.mat$se.longest.leaf.cm[f] <- as.numeric(as.character(curr.bio$SE.Longest.Leaf.Length.cm.[1]))
  nmds.mat$leaf.C[f] <- as.numeric(as.character(curr.bio$Mean.Leaf.PercC[1]))
#   nmds.mat$se.leaf.C[f] <- as.numeric(as.character(curr.bio$SE.Mean.Leaf.PercC[1]))
  nmds.mat$leaf.N[f] <- as.numeric(as.character(curr.bio$Mean.Leaf.PercN[1]))
#   nmds.mat$se.leaf.N[f] <- as.numeric(as.character(curr.bio$SE.Mean.Leaf.PercN[1]))
  
  ## host community structure
  nmds.mat$zmarina.above.biomass[f] <- as.numeric(as.character(curr.bio$Above.Zmarina.g[1]))
#   nmds.mat$se.zmarina.above.biomass[f] <- as.numeric(as.character(curr.bio$SE.Mean.Above.Zmarina.g[1]))
  nmds.mat$zmarina.below.biomass[f] <- as.numeric(as.character(curr.bio$Below.Zmarina.g[1]))
#   nmds.mat$se.zmarina.below.biomass[f] <- as.numeric(as.character(curr.bio$SE.Below.Zmarina.g[1]))
  nmds.mat$mean.zmarina.shoots.m2[f] <- as.numeric(as.character(curr.bio$Mean.Shoots.Zmarina.per.m2[1]))
#   nmds.mat$se.zmarina.shoots.m2[f] <- as.numeric(as.character(curr.bio$SE.Shoots.Zmarina.per.m2[1]))
#   nmds.mat$mean.periphyton[f] <- as.numeric(as.character(curr.bio$Mean.Site.Std.Periphyton[1]))
#   nmds.mat$se.periphyton[f] <- as.numeric(as.character(curr.bio$SE.Site.Std.Periphyton[1]))
  nmds.mat$mean.macroalgae[f] <- as.numeric(as.character(curr.bio$Mean.Macroalgae.g.m2[1]))
#   nmds.mat$se.macroalgae[f] <- as.numeric(as.character(curr.bio$SE.Mean.Macroalgae.g.m2[1]))
  
  ## animal community structure
  nmds.mat$mean.mesograzer.b[f] <- as.numeric(as.character(curr.bio$Std.Total.Biomass.Mesograzers1[1]))
#   nmds.mat$se.mesograzer.b[f] <- as.numeric(as.character(curr.bio$SE.Std.Total.Biomass.Mesograzers[1]))
#   nmds.mat$mean.grazer.richness[f] <- as.numeric(as.character(curr.bio$Mean.Grazer.Richness[1]))
#   nmds.mat$se.grazer.richness[f] <- as.numeric(as.character(curr.bio$SE.Grazer.Richness[1]))
#   nmds.mat$mean.std.epibiota[f] <- as.numeric(as.character(curr.bio$Mean.Site.Std.TotalEpibiota[1]))
#   nmds.mat$se.std.epibiota[f] <- as.numeric(as.character(curr.bio$SE.Site.Std.TotalEpibiota[1]))
#   nmds.mat$std.mollusc.b[f] <- as.numeric(as.character(curr.bio$Std.Total.Biomass.Molluscs1[1]))
#   nmds.mat$se.std.mollusc.b[f] <- as.numeric(as.character(curr.bio$SE.Std.Total.Biomass.Molluscs[1]))
  nmds.mat$std.crustacean.b[f] <- as.numeric(as.character(curr.bio$Std.Total.Biomass.Crustacean1[1]))
#   nmds.mat$se.std.crustacean.b[f] <- as.numeric(as.character(curr.bio$SE.Std.Total.Biomass.Crustacean[1]))
  
  }, error = function(e){})
  
  ## drop us a line
  if(f == 1){bar.vec <- c(na.omit(seq(1:length(rownames(nmds.mat)))[1:length(rownames(nmds.mat)) * round(length(rownames(nmds.mat)) / 10)]))
             cat('|')}
  if(f %in% bar.vec == TRUE){cat('=====|')}
}

nmds.mat$subsite.tag <- paste(nmds.mat$site, nmds.mat$subsite.code, sep = '.')
nmds.mat$lump.tag <- paste(nmds.mat$subsite.tag, nmds.mat$sample.type, sep = '.')
nmds.mat$plant.tag <- paste(nmds.mat$subsite.tag, nmds.mat$plot.loc, sep = '.')
nmds.mat$C.N.ratio <- nmds.mat$leaf.C / nmds.mat$leaf.N
nmds.mat <- subset(nmds.mat, sample.type != 'Epiphyte')

# ## mark re-ran samples
# reruns <- read.csv('/Users/Ashkaan/Dropbox/SMP/data/rerun_500_numseqs_fix.csv', header = TRUE)
# nmds.mat$rerun <- rownames(nmds.mat) %in% reruns$sample

names(nmds.mat)[c(1, 2)] <- c('X1', 'X2')

## na.omit if you wanna be super clean
nmds.mat.na <- na.omit(nmds.mat)
# nmds.mat.na <- nmds.mat[!is.na(nmds.mat$sample.type), ]

## add scores from first pca axis describing each broad data type
# pca.abiotic <- princomp(nmds.mat.na[, c(5,10:12)]) ## loads mostly on latitude and salinity
# nmds.mat.na$abiotic.char <- pca.abiotic$scores[, 1]
# 
# pca.leaf.char <- princomp(nmds.mat.na[, c(16:25)]) ## loads mostly on longest leaf (cm)
# nmds.mat.na$leaf.char <- pca.leaf.char$scores[, 1]
# 
# pca.host.comm <- princomp(nmds.mat.na[, c(13,14,26:34)]) ## loads on mean shoots m2 and below ground biomass
# nmds.mat.na$host.comm <- pca.host.comm$scores[, 1]
# 
# pca.animal <- princomp(nmds.mat.na[, c(36:45)]) ## loads on mesograzer b and crustacean b
# nmds.mat.na$animal.comm <- pca.animal$scores[, 1]


## end of biotic data cleanup

# pca.all <- princomp(nmds.mat.na[, c(12:14, 18:47)])
# loadings(pca.all)

# longest.leaf.cm, zmarina.above.biomass, zmarina.below.biomass, mean.zmarina.shoots.m2, 
# mean.mesograzer.b, mean.macroalgae, std.crustacean.b


