## define location metadata in data frame
nmds.mat.uw$site <- rep(NA)
nmds.mat.uw$above.below <- rep(NA)
nmds.mat.uw$host.env <- rep(NA)
nmds.mat.uw$is.leaf <- rep(FALSE)
nmds.mat.uw$is.roots <- rep(FALSE)
nmds.mat.uw$is.sediment <- rep(FALSE)
nmds.mat.uw$is.water <- rep(FALSE)
nmds.mat.uw$ocean <- rep(NA)
nmds.mat.uw$lat <- rep(NA)
nmds.mat.uw$lon <- rep(NA)
nmds.mat.uw$sample.type <- rep(NA)
nmds.mat.uw$sample.loc <- rep(NA)
nmds.mat.uw$plot.loc <- rep(NA)
nmds.mat.uw$subsite.code <- rep(NA)

## abiotic environment
nmds.mat.uw$temp.C <- rep(NA)
nmds.mat.uw$salinity <- rep(NA)
# nmds.mat.uw$day.length.h <- rep(NA)

## host community genetics
nmds.mat.uw$genotype.richness <- rep(NA)
nmds.mat.uw$allele.richness <- rep(NA)
# nmds.mat.uw$inbreeding <- rep(NA)

## leaf characters
nmds.mat.uw$mean.sheath.width <- rep(NA)
# nmds.mat.uw$se.sheath.width <- rep(NA)
nmds.mat.uw$mean.sheath.length <- rep(NA)
# nmds.mat.uw$se.sheath.length <- rep(NA)
nmds.mat.uw$longest.leaf.cm <- rep(NA)
# nmds.mat.uw$se.longest.leaf.cm <- rep(NA)
nmds.mat.uw$leaf.C <- rep(NA)
# nmds.mat.uw$se.leaf.C <- rep(NA)
nmds.mat.uw$leaf.N <- rep(NA)
# nmds.mat.uw$se.leaf.N <- rep(NA)

## host community structure
nmds.mat.uw$zmarina.above.biomass <- rep(NA)
# nmds.mat.uw$se.zmarina.above.biomass <- rep(NA)
nmds.mat.uw$zmarina.below.biomass <- rep(NA)
# nmds.mat.uw$se.zmarina.below.biomass <- rep(NA)
nmds.mat.uw$mean.zmarina.shoots.m2 <- rep(NA)
# nmds.mat.uw$se.zmarina.shoots.m2 <- rep(NA)
# nmds.mat.uw$mean.periphyton <- rep(NA)
# nmds.mat.uw$se.periphyton <- rep(NA)
nmds.mat.uw$mean.macroalgae <- rep(NA)
# nmds.mat.uw$se.macroalgae <- rep(NA)

## animal community structure
nmds.mat.uw$mean.mesograzer.b <- rep(NA)
# nmds.mat.uw$se.mesograzer.b <- rep(NA)
# nmds.mat.uw$mean.grazer.richness <- rep(NA)
# nmds.mat.uw$se.grazer.richness <- rep(NA)
# nmds.mat.uw$mean.std.epibiota <- rep(NA)
# nmds.mat.uw$se.std.epibiota <- rep(NA)
# nmds.mat.uw$std.mollusc.b <- rep(NA)
# nmds.mat.uw$se.std.mollusc.b <- rep(NA)
nmds.mat.uw$std.crustacean.b <- rep(NA)
# nmds.mat.uw$se.std.crustacean.b <- rep(NA)

## for loop
for(f in 1:length(rownames(nmds.mat.uw))){
  tryCatch({
  curr.site <- as.character(rownames(nmds.mat.uw)[f])
  curr.dat <- subset(meta, SampleID == curr.site)
  curr.bio <- subset(biotic, Site == as.character(curr.dat$sub.code))
  
  ## pull in data
  ## define location metadata in data frame
  nmds.mat.uw$site[f] <- as.character(curr.bio$Site.Code[1])
  nmds.mat.uw$ocean[f] <- as.character(curr.bio$Ocean[1])
  nmds.mat.uw$lat[f] <- as.numeric(as.character(curr.bio$Latitude[1]))
  nmds.mat.uw$lon[f] <- as.numeric(as.character(curr.bio$Longitude[1]))
  nmds.mat.uw$sample.type[f] <- as.character(curr.dat$SampleType[1])
  if(nmds.mat.uw$sample.type[f] == 'Leaf'){
    nmds.mat.uw$above.below[f] <- 'Above'
    nmds.mat.uw$host.env[f] <- 'Host'
    nmds.mat.uw$is.leaf[f] <- TRUE
  }
  if(nmds.mat.uw$sample.type[f] == 'Roots'){
    nmds.mat.uw$above.below[f] <- 'Below'
    nmds.mat.uw$host.env[f] <- 'Host'
    nmds.mat.uw$is.roots[f] <- TRUE
  }
  if(nmds.mat.uw$sample.type[f] == 'Sediment'){
    nmds.mat.uw$above.below[f] <- 'Below'
    nmds.mat.uw$host.env[f] <- 'Environment'
    nmds.mat.uw$is.sediment[f] <- TRUE
  }
  if(nmds.mat.uw$sample.type[f] == 'Water'){
    nmds.mat.uw$above.below[f] <- 'Above'
    nmds.mat.uw$host.env[f] <- 'Environment'
    nmds.mat.uw$is.water[f] <- TRUE
  }
  nmds.mat.uw$sample.loc[f] <- as.character(curr.bio$Depth.Categorical[1])
  nmds.mat.uw$subsite.code[f] <- as.character(curr.dat$SubsiteNumber[1])
  nmds.mat.uw$plot.loc[f] <- as.character(curr.dat$PlotLocation[1])
  
  ## abiotic environment
  nmds.mat.uw$temp.C[f] <- as.numeric(as.character(curr.bio$Temperature.C[1]))
  nmds.mat.uw$salinity[f] <- as.numeric(as.character(curr.bio$Salinity.ppt[1]))
#   nmds.mat.uw$day.length.h[f] <- as.numeric(as.character(curr.bio$Day.length.hours[1]))
  
  ## host community genetics
  nmds.mat.uw$genotype.richness[f] <- as.numeric(as.character(curr.bio$GenotypicRichness[1]))
  nmds.mat.uw$allele.richness[f] <- as.numeric(as.character(curr.bio$AllelicRichness[1]))
#   nmds.mat.uw$inbreeding[f] <- as.numeric(as.character(curr.bio$Inbreeding[1]))
  
  ## leaf characters
  nmds.mat.uw$mean.sheath.width[f] <- as.numeric(as.character(curr.bio$Sheath.Width.cm.[1]))
#   nmds.mat.uw$se.sheath.width[f] <- as.numeric(as.character(curr.bio$SE.Sheath.Width.cm.[1]))
  nmds.mat.uw$mean.sheath.length[f] <- as.numeric(as.character(curr.bio$Shealth.Length.cm.[1]))
#   nmds.mat.uw$se.sheath.length[f] <- as.numeric(as.character(curr.bio$SE.Shealth.Length.cm.[1]))
  nmds.mat.uw$longest.leaf.cm[f] <- as.numeric(as.character(curr.bio$Longest.Leaf.Length.cm.[1]))
#   nmds.mat.uw$se.longest.leaf.cm[f] <- as.numeric(as.character(curr.bio$SE.Longest.Leaf.Length.cm.[1]))
  nmds.mat.uw$leaf.C[f] <- as.numeric(as.character(curr.bio$Mean.Leaf.PercC[1]))
#   nmds.mat.uw$se.leaf.C[f] <- as.numeric(as.character(curr.bio$SE.Mean.Leaf.PercC[1]))
  nmds.mat.uw$leaf.N[f] <- as.numeric(as.character(curr.bio$Mean.Leaf.PercN[1]))
#   nmds.mat.uw$se.leaf.N[f] <- as.numeric(as.character(curr.bio$SE.Mean.Leaf.PercN[1]))
  
  ## host community structure
  nmds.mat.uw$zmarina.above.biomass[f] <- as.numeric(as.character(curr.bio$Above.Zmarina.g[1]))
#   nmds.mat.uw$se.zmarina.above.biomass[f] <- as.numeric(as.character(curr.bio$SE.Mean.Above.Zmarina.g[1]))
  nmds.mat.uw$zmarina.below.biomass[f] <- as.numeric(as.character(curr.bio$Below.Zmarina.g[1]))
#   nmds.mat.uw$se.zmarina.below.biomass[f] <- as.numeric(as.character(curr.bio$SE.Below.Zmarina.g[1]))
  nmds.mat.uw$mean.zmarina.shoots.m2[f] <- as.numeric(as.character(curr.bio$Mean.Shoots.Zmarina.per.m2[1]))
#   nmds.mat.uw$se.zmarina.shoots.m2[f] <- as.numeric(as.character(curr.bio$SE.Shoots.Zmarina.per.m2[1]))
#   nmds.mat.uw$mean.periphyton[f] <- as.numeric(as.character(curr.bio$Mean.Site.Std.Periphyton[1]))
#   nmds.mat.uw$se.periphyton[f] <- as.numeric(as.character(curr.bio$SE.Site.Std.Periphyton[1]))
  nmds.mat.uw$mean.macroalgae[f] <- as.numeric(as.character(curr.bio$Mean.Macroalgae.g.m2[1]))
#   nmds.mat.uw$se.macroalgae[f] <- as.numeric(as.character(curr.bio$SE.Mean.Macroalgae.g.m2[1]))
  
  ## animal community structure
  nmds.mat.uw$mean.mesograzer.b[f] <- as.numeric(as.character(curr.bio$Std.Total.Biomass.Mesograzers1[1]))
#   nmds.mat.uw$se.mesograzer.b[f] <- as.numeric(as.character(curr.bio$SE.Std.Total.Biomass.Mesograzers[1]))
#   nmds.mat.uw$mean.grazer.richness[f] <- as.numeric(as.character(curr.bio$Mean.Grazer.Richness[1]))
#   nmds.mat.uw$se.grazer.richness[f] <- as.numeric(as.character(curr.bio$SE.Grazer.Richness[1]))
#   nmds.mat.uw$mean.std.epibiota[f] <- as.numeric(as.character(curr.bio$Mean.Site.Std.TotalEpibiota[1]))
#   nmds.mat.uw$se.std.epibiota[f] <- as.numeric(as.character(curr.bio$SE.Site.Std.TotalEpibiota[1]))
#   nmds.mat.uw$std.mollusc.b[f] <- as.numeric(as.character(curr.bio$Std.Total.Biomass.Molluscs1[1]))
#   nmds.mat.uw$se.std.mollusc.b[f] <- as.numeric(as.character(curr.bio$SE.Std.Total.Biomass.Molluscs[1]))
  nmds.mat.uw$std.crustacean.b[f] <- as.numeric(as.character(curr.bio$Std.Total.Biomass.Crustacean1[1]))
#   nmds.mat.uw$se.std.crustacean.b[f] <- as.numeric(as.character(curr.bio$SE.Std.Total.Biomass.Crustacean[1]))
  
  }, error = function(e){})
  
  ## drop us a line
  if(f == 1){bar.vec <- c(na.omit(seq(1:length(rownames(nmds.mat.uw)))[1:length(rownames(nmds.mat.uw)) * round(length(rownames(nmds.mat.uw)) / 10)]))
             cat('|')}
  if(f %in% bar.vec == TRUE){cat('=====|')}
}

nmds.mat.uw$subsite.tag <- paste(nmds.mat.uw$site, nmds.mat.uw$subsite.code, sep = '.')
nmds.mat.uw$lump.tag <- paste(nmds.mat.uw$subsite.tag, nmds.mat.uw$sample.type, sep = '.')
nmds.mat.uw$plant.tag <- paste(nmds.mat.uw$subsite.tag, nmds.mat.uw$plot.loc, sep = '.')
nmds.mat.uw$C.N.ratio <- nmds.mat.uw$leaf.C / nmds.mat.uw$leaf.N
nmds.mat.uw <- subset(nmds.mat.uw, sample.type != 'Epiphyte')

# ## mark re-ran samples
# reruns <- read.csv('/Users/Ashkaan/Dropbox/SMP/data/rerun_500_numseqs_fix.csv', header = TRUE)
# nmds.mat.uw$rerun <- rownames(nmds.mat.uw) %in% reruns$sample

names(nmds.mat.uw)[c(1, 2)] <- c('X1', 'X2')

## na.omit if you wanna be super clean
nmds.mat.uw.na <- na.omit(nmds.mat.uw)
# nmds.mat.uw.na <- nmds.mat.uw[!is.na(nmds.mat.uw$sample.type), ]

## add scores from first pca axis describing each broad data type
# pca.abiotic <- princomp(nmds.mat.uw.na[, c(5,10:12)]) ## loads mostly on latitude and salinity
# nmds.mat.uw.na$abiotic.char <- pca.abiotic$scores[, 1]
# 
# pca.leaf.char <- princomp(nmds.mat.uw.na[, c(16:25)]) ## loads mostly on longest leaf (cm)
# nmds.mat.uw.na$leaf.char <- pca.leaf.char$scores[, 1]
# 
# pca.host.comm <- princomp(nmds.mat.uw.na[, c(13,14,26:34)]) ## loads on mean shoots m2 and below ground biomass
# nmds.mat.uw.na$host.comm <- pca.host.comm$scores[, 1]
# 
# pca.animal <- princomp(nmds.mat.uw.na[, c(36:45)]) ## loads on mesograzer b and crustacean b
# nmds.mat.uw.na$animal.comm <- pca.animal$scores[, 1]


## end of biotic data cleanup

# pca.all <- princomp(nmds.mat.uw.na[, c(12:14, 18:47)])
# loadings(pca.all)

# longest.leaf.cm, zmarina.above.biomass, zmarina.below.biomass, mean.zmarina.shoots.m2, 
# mean.mesograzer.b, mean.macroalgae, std.crustacean.b


