##--------------------------------------------------
## this script reproduces analyses from a study
## on the microbiomes of a globally distributed 
## seagrass species , Z. marina, collected by the 
## Zostera Ecological Network (ZEN).
##
## - Ashkaan K. Fahimipour
##--------------------------------------------------

##--------------------------------------------------
## libraries and user-defined functions
##--------------------------------------------------
setwd('/Users/Ashkaan/Dropbox/ZEN-microbiome-project/')
source('./R/ud-functions/microbiome-functions.R')
source('./R/ud-functions/summary-SE.R')
set.seed(777)

##--------------------------------------------------
## import data
##--------------------------------------------------
## (1) non-rarefied biom table
biom <- read_biom('./data/pick-otus/otu_table_filt_JSON.biom')

## (2) rarefied biom table to depth of 1000
# biom <- read_biom('./data/pick-otus/otu_table_filt_rare1000_JSON.biom')

## coerce to matrix
dat <- as.data.frame(as(biom_data(biom), 'matrix'), header = TRUE)

## drop NOPNA (control) samples
dat <- dat[, -which(lapply(strsplit(colnames(dat), split = 'NOPNA'), length) == 2)]

## (2) import metadata
tax <- as.matrix(observation_metadata(biom), sep = ',')
meta <- as.data.frame(read.csv('./data/meta.csv'))
biotic <- read.csv('./data/meta-biotic.csv', header = TRUE)

## clean up metadata
meta$SubsiteNumber[which(meta$SubsiteNumber == 1)] <- 'A'
meta$SubsiteNumber[which(meta$SubsiteNumber == 2)] <- 'B'
meta$sub.code <- paste(meta$ZenSite, meta$SubsiteNumber, sep = '.')
meta$plant.code <- paste(meta$ZenSite, meta$SubsiteNumber, meta$PlotLocation, sep = '.')

## Exclude epiphyte samples
dont.use <- meta$SampleID[which(meta$SampleType == 'Epiphyte')]
dat <- dat[, -which(names(dat) %in% dont.use)]

## remove sites with less than k reads
k.filter <- 1000
if(length(which(colSums(dat) < k.filter)) > 0){dat <- dat[, -c(which(colSums(dat) <= k.filter))]}

## remove OTUs that occur less than k times
dat <- dat[-which(rowSums(dat != 0) < 1), ]

## trim taxonomy list to OTUs that are present in dat
tax.trim <- as.matrix(tax[which(rownames(tax) %in% rownames(dat)), ])
tax.2 <- tax.to.long(tax.trim)

# ##-----------------------------
# ## normalize biom table by TMM
# ##-----------------------------
# counts <- DGEList(dat)
# norm.dat <- calcNormFactors(counts, method = 'TMM')
# norm.dat.2 <- cpm(norm.dat, normalized.lib.sizes = TRUE)
# dat <- norm.dat.2
# ##-----------------------------

## dividing elements by sample sums
dat.rel <- t(t(dat) / rowSums(t(dat)))

## sqrt transform normalized biom table for downstream analyses
dat.rel.trans <- sqrt(dat.rel)

##--------------------------------------------------
## SECTION 1: ORDINATION
##--------------------------------------------------

## (1) vegan distances or...
dist.metric.veg <- vegdist(t(dat.rel.trans), method = 'canberra', binary = FALSE)

## (2) UW unifrac distance computed in phyloseq (need to load ps code snippet below, first)
# dist.metric.veg.uw <- UniFrac(ps, weighted = FALSE)

## (3) W unifrac distance computed in phyloseq
# dist.metric.veg.w <- UniFrac(ps, weighted = TRUE)

## (a) nmds ordination
nmds <- metaMDS(dist.metric.veg, k = 2, trymax = 25, pca = TRUE)
nmds.mat <- as.data.frame(nmds$points)
nmds$stress

## (b) pcoa with vegan
# pcoa <- cmdscale(dist.metric.veg, k = 2, eig = FALSE)
# nmds.mat <- data.frame(pcoa)

## coerce dist into matrix object for your viewing pleasure
dist.metric <- as.matrix(dist.metric.veg)

## custom script to clean it all up and rename variables to ones I like
source('./R/ud-functions/clean-ordination.R')

##-----------------------------
##-----------------------------
# ## bring in metadata
# metadata <- nmds.mat
# metadata$coast <- biotic$Coast[match(metadata$subsite.tag, biotic$Site)]
# 
# ## otu table
# otus <- dat
# 
# ## drop rare otus
# otus <- otus[-which(rowSums(otus != 0) < 5), ]
# 
# ## pick leaves only by LEAVING OUT roots
# otus <- otus[, match(rownames(subset(metadata, sample.type != 'Roots')), colnames(otus))]
# 
# ## convert to integer counts
# otus.2 <- apply(otus, 2, function(x) as.integer(x))
# rownames(otus.2) <- rownames(otus)
# 
# ## extract only those samples in common between the two tables
# common.sample.ids <- intersect(rownames(metadata), colnames(otus.2))
# otus.2 <- otus.2[, common.sample.ids]
# otus.2 <- t(otus.2)
# metadata <- metadata[common.sample.ids, ]
# 
# ## add source/sink data to data frame
# metadata$SourceSink <- rep('sink')
# metadata$SourceSink[which(with(metadata, sample.type == 'Water' | sample.type == 'Sediment'))] <- rep('source')
# 
# ## extract the source environments and source/sink indices
# unique(metadata$coast)
# train.ix <- which(metadata$SourceSink == 'source' & metadata$coast == 'East.Atlantic')
# test.ix <- which(metadata$SourceSink == 'sink'& metadata$coast == 'East.Atlantic')
# envs <- metadata$sample.type
# metadata$Description <- metadata$site
# desc <- metadata$Description
# 
# ## load SourceTracker code
# source('./R/ud-functions/source-tracker.r')
# 
# ## to skip tuning, run this
# alpha1 <- 2e-3
# alpha2 <- 1e-2
# 
# ## train SourceTracker object on training data
# st <- sourcetracker(otus.2[train.ix,], envs[train.ix], rarefaction_depth = min(rowSums(otus.2)))
# 
# ## estimate source proportions in test data
# results <- predict(st, otus.2[test.ix, ], alpha1 = alpha1, alpha2 = alpha2, rarefaction_depth = min(rowSums(otus.2)))
# 
# ## write to disk
# write.csv(results$proportions, './data/sourcetracker/roots-EastAtlantic.csv')
# write.csv(results$proportions_sd, './data/sourcetracker/roots-EastAtlantic-sd.csv')
# 
# ## load
# # source.tracker <- read.csv('./data/source-tracker-lumped.csv', header = TRUE, row.names = 1)
# source.tracker <- as.data.frame(results$proportions)
# source.tracker$sample.type <- nmds.mat$sample.type[match(rownames(source.tracker), rownames(nmds.mat))]
# 
# # source.tracker.sd <- read.csv('./data/source-tracker-lumped-sd.csv', header = TRUE, row.names = 1)
# source.tracker.sd <- as.data.frame(results$proportions_sd)
# names(source.tracker.sd) <- c('Sediment.sd', 'Water.sd', 'Unknown.sd')
# source.tracker <- cbind(source.tracker, source.tracker.sd)

##-----------------------------
## load pre-ran instances
##-----------------------------
source.tracker <- read.csv('./data/sourcetracker/sourcetracker.csv', header = TRUE, row.names = 1)

## cube root function
cubrt_trans <- function() trans_new('cubrt', function(x) 10*x^(1/3), function(x) x^(1/3))

main.l <- ggplot(subset(source.tracker, type == 'Leaf'), aes(x = Sediment, y = Water)) +
  stat_density2d(aes(fill = ..density..), geom = 'tile', contour = FALSE, n = 200) +
  # geom_errorbar(aes(ymin = Water - Water.sd, ymax = Water + Water.sd), size = 0.25) +
  # geom_errorbarh(aes(xmin = Sediment - Sediment.sd, xmax = Sediment + Sediment.sd), size = 0.25) +
  scale_fill_gradient2(low = '#f7f7f7', mid = '#f7fcf5', high = '#00441b', na.value = '#f7f7f7') +
  geom_point(shape = 21, size = 2, alpha = 0.8, fill = '#74c476', stroke = 0.1) +
  xlim(-0.1, 1) +
  ylim(-0.1, 1) +
  # scale_x_sqrt(limits = c(0, 1.1)) +
  # scale_y_sqrt(limits = c(0, 1.05)) +
  xlab(' ') +
  ylab('Proportion of Water Taxa') +
  expand_limits(x = c(-0.1, 1), y = c(-0.1, 1)) +
  theme_classic() +
  theme(axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5)) +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
        legend.position = 'none', axis.text.x = element_blank())

main.r <- ggplot(subset(source.tracker, type == 'Roots'), aes(x = Sediment, y = Water)) +
  stat_density2d(aes(fill = ..density..), geom = 'tile', contour = FALSE, n = 200) +
  # geom_errorbar(aes(ymin = Water - Water.sd, ymax = Water + Water.sd), size = 0.25) +
  # geom_errorbarh(aes(xmin = Sediment - Sediment.sd, xmax = Sediment + Sediment.sd), size = 0.25) +
  scale_fill_gradient2(low = '#f7f7f7', mid = '#fcfbfd', high = '#3f007d', na.value = '#f7f7f7') +
  geom_point(shape = 21, size = 2, alpha = 0.8, fill = '#807dba', stroke = 0.1) +
  xlim(-0.1, 1) +
  ylim(-0.1, 1) +
  # scale_x_sqrt(limits = c(0, 1)) +
  # scale_y_sqrt(limits = c(0, 1)) +
  xlab('Proportion of Sediment Taxa') +
  ylab('Proportion of Water Taxa') +
  expand_limits(x = c(-0.1, 1), y = c(-0.1, 1)) +
  theme_classic() +
  theme(axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5)) +
  theme(strip.background = element_blank(), strip.text.x = element_blank(), legend.position = 'none')

## plot together on a grid
dev.off()
gridExtra::grid.arrange(main.l, main.r, nrow = 2)

## END

