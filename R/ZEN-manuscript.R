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

##-----------------------------
## normalize biom table by TMM
##-----------------------------
counts <- DGEList(dat)
norm.dat <- calcNormFactors(counts, method = 'TMM')
norm.dat.2 <- cpm(norm.dat, normalized.lib.sizes = TRUE)
dat <- norm.dat.2
##-----------------------------

## dividing elements by sample sums
dat.rel <- t(t(dat) / rowSums(t(dat)))

## sqrt transform normalized biom table for downstream analyses
dat.rel.trans <- sqrt(dat.rel)

##----------------------------
## make phyloseq object with tree
##----------------------------
## note: reading the tree spits out a warning due to Newick formatting,
tree <- read_tree_greengenes('./data/pick-otus/pruned_fast.tre') 
tree$edge.length <- compute.brlen(tree, 1)$edge.length

## prep biom table to coerce to phyloseq object
ps.dat <- dat
ps.otu <- otu_table(ps.dat, taxa_are_rows = TRUE)

## prep metadata
ps.meta <- sample_data(nmds.mat[which(rownames(nmds.mat) %in% colnames(ps.dat)), ])

## prep taxonomy
tax.3 <- as.matrix(tax.2)
rownames(tax.3) <- tax.3[, 'otu']
tax.3 <- tax.3[which(rownames(tax.3) %in% rownames(ps.dat)), ]
ps.tax <- tax_table(tax.3)
colnames(ps.tax) <- colnames(tax.3)
taxa_names(ps.tax) <- rownames(tax.3)

## coerce into phyloseq object
ps <- phyloseq(ps.otu, ps.meta, ps.tax, tree)

##--------------------------------------------------
## SECTION 1: ORDINATION
##--------------------------------------------------

## (1) vegan distances or...
dist.metric.veg <- vegdist(t(dat.rel.trans), method = 'canberra', binary = FALSE)

## (2) UW unifrac distance computed in phyloseq (need to load ps code snippet below, first)
dist.metric.veg.uw <- UniFrac(ps, weighted = FALSE)

# ## (a) nmds ordination for canberra
# nmds <- metaMDS(dist.metric.veg, k = 2, trymax = 30, pca = TRUE)
# nmds.mat <- as.data.frame(nmds$points)
# nmds$stress
# 
# ## (a) nmds ordination for UniFrac
# nmds.uw <- metaMDS(dist.metric.veg.uw, k = 2, trymax = 30, pca = TRUE)
# nmds.mat.uw <- as.data.frame(nmds.uw$points)
# nmds.uw$stress

## (b) pcoa with vegan
pcoa <- cmdscale(dist.metric.veg, k = 2, eig = FALSE)
nmds.mat <- data.frame(pcoa)

## (b) pcoa with vegan UniFrac
pcoa.uw <- cmdscale(dist.metric.veg.uw, k = 2, eig = FALSE)
nmds.mat.uw <- data.frame(pcoa.uw)

## coerce dist into matrix object for your viewing pleasure
dist.metric <- as.matrix(dist.metric.veg)
dist.metric.uw <- as.matrix(dist.metric.veg.uw)

## custom script to clean it all up and rename variables to ones I like
source('./R/ud-functions/clean-ordination.R')
source('./R/ud-functions/clean-ordination-uw.R')

## choose color palette
colorpal <- c('#78c679', '#045a8d', '#9970ab', '#cb181d')

## reorder factors to be more sensible
nmds.mat$sample.type <- factor(nmds.mat$sample.type, levels = c('Leaf', 'Water', 'Roots', 'Sediment'))
nmds.mat.uw$sample.type <- factor(nmds.mat.uw$sample.type, levels = c('Leaf', 'Water', 'Roots', 'Sediment'))

## ggplot function
ord.1 <- ggplot(data = nmds.mat[!is.na(nmds.mat$sample.type), ], 
                aes(x = X1, y = X2, 
                    fill = as.factor(sample.type),
                    shape = as.factor(host.env),
                    colour = as.factor(sample.type))) +
  stat_ellipse(data = nmds.mat[!is.na(nmds.mat$sample.type), ],
               aes(x = X1, y = X2, colour = as.factor(sample.type)),
               alpha = 0.15, geom = 'polygon', size = 0.0, linetype = 1,
               type = 't', level = 0.95) +
  # geom_path(aes(group = plant.tag), alpha = 0.8, size = 0.5, linetype = 1, colour = '#bdbdbd') +
  # geom_path(aes(group = plant.tag)) +
  geom_point(alpha = 0.75, colour = '#252525', size = 2, 
             shape = 21,
             stroke = 0.2) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        panel.border = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(color = 'black')) +
  xlab(' ') +
  ylab('Axis 2') +
  scale_shape_manual(values = c(21, 22), guide = guide_legend(title = NULL)) +
  scale_colour_manual(values = colorpal, guide = guide_legend(title = NULL)) +
  scale_fill_manual(values = colorpal, guide = guide_legend(title = NULL)) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.5),
        legend.position = 'none',
        axis.line.y = element_line(color = "black", size = 0.5))

## ggplot function
ord.2 <- ggplot(data = nmds.mat.uw[!is.na(nmds.mat.uw$sample.type), ], 
                aes(x = X1, y = X2, 
                    fill = as.factor(sample.type),
                    shape = as.factor(host.env),
                    colour = as.factor(sample.type))) +
  stat_ellipse(data = nmds.mat.uw[!is.na(nmds.mat.uw$sample.type), ],
               aes(x = X1, y = X2, colour = as.factor(sample.type)),
               alpha = 0.15, geom = 'polygon', size = 0.0, linetype = 1,
               type = 't', level = 0.95) +
  # geom_path(aes(group = plant.tag), alpha = 0.8, size = 0.5, linetype = 1, colour = '#bdbdbd') +
  # geom_path(aes(group = plant.tag)) +
  geom_point(alpha = 0.75, colour = '#252525', size = 2, 
             shape = 21,
             stroke = 0.2) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        panel.border = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(color = 'black')) +
  xlab('Axis 1') +
  ylab('Axis 2') +
  scale_shape_manual(values = c(21, 22), guide = guide_legend(title = NULL)) +
  scale_colour_manual(values = colorpal, guide = guide_legend(title = NULL)) +
  scale_fill_manual(values = colorpal, guide = guide_legend(title = NULL)) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.5),
        legend.position = 'none',
        axis.line.y = element_line(color = "black", size = 0.5))

##---------------------
## centroid analysis canberra
##---------------------
## subsite scale
nmds.mat$cluster <- with(nmds.mat, paste(subsite.tag, sample.type, sep = '.'))
cent.subsite <- aggregate(nmds.mat[ ,c('X1', 'X2')], list(group = nmds.mat$cluster), mean)
cent.subsite <- cbind(cent.subsite, nmds.mat[match(cent.subsite$group, nmds.mat$cluster), -c(1, 2)])
rownames(cent.subsite) <- cent.subsite$group

## subsite level centroids
cent.dist <- 1 - as.matrix(vegdist(cent.subsite[, 2:3], method = 'euclidean'))

## pare for above/belowground
cent.dist <- cent.dist[rownames(cent.dist) %in% subset(cent.subsite, sample.type == 'Leaf')$group, ]
cent.dist <- cent.dist[, colnames(cent.dist) %in% subset(cent.subsite, sample.type == 'Water')$group]

## melt it
cent.melt <- melt(cent.dist)
cent.melt$tag.1 <- cent.subsite$subsite.tag[match(cent.melt$Var1, cent.subsite$group)]
cent.melt$tag.2 <- cent.subsite$subsite.tag[match(cent.melt$Var2, cent.subsite$group)]

## add type factor
cent.melt$same <- rep(FALSE)
cent.melt$same[which(cent.melt$tag.1 == cent.melt$tag.2)] <- rep(TRUE)

## summary stats
sum.cent <- summarySE(cent.melt, measurevar = 'value', groupvars = 'same')
sum.cent$type <- rep('Leaf')

## stats
currentStat <- diffMeans(cent.melt$value, cent.melt$same)
cat("Initial Test Statistic: ", currentStat)
p.val <- permutationTest(cent.melt$value, cent.melt$same, testStat = diffMeans)
cat("p-value: ", p.val)

##--------
## subsite scale
nmds.mat$cluster <- with(nmds.mat, paste(subsite.tag, sample.type, sep = '.'))
cent.subsite <- aggregate(nmds.mat[ ,c('X1', 'X2')], list(group = nmds.mat$cluster), mean)
cent.subsite <- cbind(cent.subsite, nmds.mat[match(cent.subsite$group, nmds.mat$cluster), -c(1, 2)])
rownames(cent.subsite) <- cent.subsite$group

## subsite level centroids
cent.dist <- 1 - as.matrix(dist(cent.subsite[, 2:3]))

## pare for above/belowground
cent.dist <- cent.dist[rownames(cent.dist) %in% subset(cent.subsite, sample.type == 'Roots')$group, ]
cent.dist <- cent.dist[, colnames(cent.dist) %in% subset(cent.subsite, sample.type == 'Sediment')$group]

## melt it
cent.melt <- melt(cent.dist)
cent.melt$tag.1 <- cent.subsite$subsite.tag[match(cent.melt$Var1, cent.subsite$group)]
cent.melt$tag.2 <- cent.subsite$subsite.tag[match(cent.melt$Var2, cent.subsite$group)]

## add type factor
cent.melt$same <- rep(FALSE)
cent.melt$same[which(cent.melt$tag.1 == cent.melt$tag.2)] <- rep(TRUE)

## summary stats
sum.cent.2 <- summarySE(cent.melt, measurevar = 'value', groupvars = 'same')
sum.cent.2$type <- rep('Roots')

## stats
currentStat <- diffMeans(cent.melt$value, cent.melt$same)
cat("Initial Test Statistic: ", currentStat)
p.val <- permutationTest(cent.melt$value, cent.melt$same, testStat = diffMeans)
cat("p-value: ", p.val)

## rbind
sum.cent <- rbind(sum.cent, sum.cent.2)

##---------------------
## centroid analysis UniFrac
##---------------------
## subsite scale
nmds.mat.uw$cluster <- with(nmds.mat.uw, paste(subsite.tag, sample.type, sep = '.'))
cent.subsite.uw <- aggregate(nmds.mat.uw[ ,c('X1', 'X2')], list(group = nmds.mat.uw$cluster), mean)
cent.subsite.uw <- cbind(cent.subsite.uw, nmds.mat.uw[match(cent.subsite.uw$group, nmds.mat.uw$cluster), -c(1, 2)])
rownames(cent.subsite.uw) <- cent.subsite.uw$group

## subsite level centroids
cent.dist.uw <- 1 - as.matrix(vegdist(cent.subsite.uw[, 2:3], method = 'euclidean'))

## pare for above/belowground
cent.dist.uw <- cent.dist.uw[rownames(cent.dist.uw) %in% subset(cent.subsite.uw, sample.type == 'Leaf')$group, ]
cent.dist.uw <- cent.dist.uw[, colnames(cent.dist.uw) %in% subset(cent.subsite.uw, sample.type == 'Water')$group]

## melt it
cent.melt.uw <- melt(cent.dist.uw)
cent.melt.uw$tag.1 <- cent.subsite.uw$subsite.tag[match(cent.melt.uw$Var1, cent.subsite.uw$group)]
cent.melt.uw$tag.2 <- cent.subsite.uw$subsite.tag[match(cent.melt.uw$Var2, cent.subsite.uw$group)]

## add type factor
cent.melt.uw$same <- rep(FALSE)
cent.melt.uw$same[which(cent.melt.uw$tag.1 == cent.melt.uw$tag.2)] <- rep(TRUE)

## summary stats
sum.cent.uw <- summarySE(cent.melt.uw, measurevar = 'value', groupvars = 'same')
sum.cent.uw$type <- rep('Leaf')

## stats
currentStat <- diffMeans(cent.melt.uw$value, cent.melt.uw$same)
cat("Initial Test Statistic: ", currentStat)
p.val <- permutationTest(cent.melt.uw$value, cent.melt.uw$same, testStat = diffMeans)
cat("p-value: ", p.val)

##--------
## subsite scale
nmds.mat.uw$cluster <- with(nmds.mat.uw, paste(subsite.tag, sample.type, sep = '.'))
cent.subsite.uw <- aggregate(nmds.mat.uw[ ,c('X1', 'X2')], list(group = nmds.mat.uw$cluster), mean)
cent.subsite.uw <- cbind(cent.subsite.uw, nmds.mat.uw[match(cent.subsite.uw$group, nmds.mat.uw$cluster), -c(1, 2)])
rownames(cent.subsite.uw) <- cent.subsite.uw$group

## subsite level centroids
cent.dist.uw <- 1 - as.matrix(dist(cent.subsite.uw[, 2:3]))

## pare for above/belowground
cent.dist.uw <- cent.dist.uw[rownames(cent.dist.uw) %in% subset(cent.subsite.uw, sample.type == 'Roots')$group, ]
cent.dist.uw <- cent.dist.uw[, colnames(cent.dist.uw) %in% subset(cent.subsite.uw, sample.type == 'Sediment')$group]

## melt it
cent.melt.uw <- melt(cent.dist.uw)
cent.melt.uw$tag.1 <- cent.subsite.uw$subsite.tag[match(cent.melt.uw$Var1, cent.subsite.uw$group)]
cent.melt.uw$tag.2 <- cent.subsite.uw$subsite.tag[match(cent.melt.uw$Var2, cent.subsite.uw$group)]

## add 'same' factor to denote same site as binary variable
cent.melt.uw$same <- rep(FALSE)
cent.melt.uw$same[which(cent.melt.uw$tag.1 == cent.melt.uw$tag.2)] <- rep(TRUE)

## summary stats
sum.cent.uw.2 <- summarySE(cent.melt.uw, measurevar = 'value', groupvars = 'same')
sum.cent.uw.2$type <- rep('Roots')

## stats
currentStat <- diffMeans(cent.melt.uw$value, cent.melt.uw$same)
cat("Initial Test Statistic: ", currentStat)
p.val <- permutationTest(cent.melt.uw$value, cent.melt.uw$same, testStat = diffMeans)
cat("p-value: ", p.val)

## rbind
sum.cent.uw <- rbind(sum.cent.uw, sum.cent.uw.2)

## ggplot it
same.l <- ggplot(data = sum.cent, aes(x = (same), y = value, group = type, fill = type)) +
  geom_errorbar(aes(ymin = value - 1*se, ymax = value + 1*se), size = 0.25, width = 0.02) +
  geom_path(size = 0.5, linetype = 2) +
  geom_point(alpha = 0.95, colour = '#252525', size = 2.5,
             shape = 21,
             stroke = 0.2) +
  # stat_smooth(method = 'lm', size = 0.75, linetype = 2, formula = y ~ x, colour = '#252525') +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        panel.border = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(color = 'black')) +
  xlab(' ') +
  ylab('Similarity') +
  # scale_colour_manual(values = colorpal, guide = guide_legend(title = NULL)) +
  scale_fill_manual(values = colorpal[c(1, 3)], guide = guide_legend(title = NULL)) +
  theme(axis.line.x = element_line(color = "black", size = 0.5),
        axis.text.x = element_blank(),
        legend.position = 'none',
        axis.line.y = element_line(color = "black", size = 0.5))

same.r <- ggplot(data = sum.cent.uw, aes(x = as.factor(same), y = value, group = type, fill = type)) +
  geom_errorbar(aes(ymin = value - 1*se, ymax = value + 1*se), size = 0.25, width = 0.02) +
  geom_path(size = 0.5, linetype = 2) +
  geom_point(alpha = 0.95, colour = '#252525', size = 2.5,
             shape = 21,
             stroke = 0.2) +
  # stat_smooth(method = 'lm', size = 0.75, linetype = 2, formula = y ~ x, colour = '#252525') +
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        panel.border = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(color = 'black')) +
  xlab(' ') +
  ylab('Similarity') +
  # scale_colour_manual(values = colorpal, guide = guide_legend(title = NULL)) +
  scale_fill_manual(values = colorpal[c(1, 3)], guide = guide_legend(title = NULL)) +
  theme(axis.line.x = element_line(color = "black", size = 0.5),
        axis.text.x = element_blank(),
        legend.position = 'none',
        axis.line.y = element_line(color = "black", size = 0.5))

## plot together on a grid
dev.off()
gridExtra::grid.arrange(ord.1, same.l, ord.2, same.r, nrow = 2)





