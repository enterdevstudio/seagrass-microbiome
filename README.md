# <img src="figures/plant_vec.jpg" width=50, height=61)> The ZEN Microbiome Project

> May 10, 2016

### 1. Introduction 
The first step in understanding the symbiotic relationship between seagrass microbes and their host is to characterize the baseline microbiota and the differences that are associated with host genotype and environment. In the Moore Foundation proposal, we proposed to address the following set of questions:

- Is there a ‘core microbiome’ or set of commonly occurring microorganisms that appear in all assemblages associated with seagrasses?
- How much is microbial community composition influenced by genetic and ecological variation in the host, and is this different for different host tissues?
- Are there significant co-occurrence or co-exclusion relationships between pairs of microbial taxa within the seagrass microbiome?

### 2. Analyses
These analyses depend heavily on Python scripts bundled with [*macQIIME* v1.9](http://www.wernerlab.org/software/macqiime) and the *biom*, *ggplot2*, and *vegan* packages in the statistical programming environment `R`.

#### Importing, sorting and cleaning the data

Let's import all of our data. These data are in the *BIOM* format, generated from [these](http://edhar.genomecenter.ucdavis.edu/~gjospin/Seagrass/) libraries using *QIIME* 1.9. OTUs were picked against the most recent *GreenGenes* database at a cutoff similarity of 97%. I rarefied samples to a depth of 1000, and normalized OTU read counts by taxon copy number using [*PICRUSt*](https://picrust.github.io/picrust/). Since, R has a tough time with HDF5 formatted files, I converted the BIOM table to JSON using [*biom convert*](http://biom-format.org/documentation/biom_conversion.html) prior to these analyses. Bash scripts that generated data are a modified version of James Meadows' scripts found [here](https://github.com/jfmeadow/ReproducibleDemo/blob/master/QIIME/pickTheseOTUs.sh); these are updated to work with *QIIME* 1.9 and tweaked to do some extra bits like filtering out chloroplasts/mitochondria and tree pruning. I'll get these posted soon.

I'll load a bunch of libraries and functions I've written in R over the last few months that automate various useful microbiome things and data reshaping. I'll get these posted soon, too.

```r
setwd('/Users/Ashkaan/Dropbox/')
source('./General Functions/microbiome_functions.R')
set.seed(34543)
```

Now I'll import the data. I also want to clean it up a little. This includes renaming some of the metadata columns and excluding samples we don't want to include in our analyses (e.g., NOPNA samples, epiphytes).

```r
## import data
biom <- read_biom('./ZEN-microbiome/data/otu_table_ZEN2_filt_rare1000_JSON_copynum.biom')
dat <- as.data.frame(as(biom_data(biom), "matrix"), header = TRUE)

## drop NOPNA samples
dat <- dat[, -which(lapply(strsplit(colnames(dat), split = 'NOPNA'), length) == 2)]

## import metadata
tax <- as.matrix(observation_metadata(biom), sep = ',')
meta <- as.data.frame(read.csv('./SMP/data/ZEN_1/97_percent/map_for_R.csv'))
biotic <- read.csv('./ZEN-microbiome/data/ZEN_biotic_data_2016_04_26.csv', header = TRUE)

## clean up meta file
meta$SubsiteNumber[which(meta$SubsiteNumber == 1)] <- 'A'
meta$SubsiteNumber[which(meta$SubsiteNumber == 2)] <- 'B'
meta$sub.code <- paste(meta$ZenSite, meta$SubsiteNumber, sep = '.')

## trim taxonomy list to OTUs that are present
tax.trim <- as.matrix(tax[which(rownames(tax) %in% rownames(dat)), ])
tax.2 <- tax.to.long(tax.trim)

## Exclude epiphyte samples
dont.use <- meta$SampleID[which(meta$SampleType == 'Epiphyte')]
dat <- dat[, -which(names(dat) %in% dont.use)]

## remove OTUs that are never sampled after rarefaction
dat <- dat[-which(rowSums(dat != 0) == 0), ]

## remove sites with no reads
if(length(which(colSums(dat != 0) < 1)) > 0){
  dat <- dat[, -c(which(colSums(dat != 0) < 1))] ## removes sites with less than k reads
}
```

Finally, I'll normalize and `cube-root` transform my BIOM table for analysis downstream. I found that this transformation yields the most reliable ordination.

```r
## add relative abundances to new biom table
dat.rel <- t(t(dat)/rowSums(t(dat)))

## transform data
dat.rel.trans <- dat.rel^(1/3)
```

#### 2a. Community Ordination
I want to analyze patterns in microbial community composition among our different sample types (Leaf, Root, Sediment and Water). To do this, I'll calculate Canberra *distances* between community pairs using `cube-root` transformed relative abundances. I'll then visualize these in 2-dimensions using nonmetric multidimensional scaling (NMDS; Fig. 1).

```r
dist.metric.veg <- vegdist(t(dat.rel.trans), method = 'canberra')

## nmds ordination
nmds <- metaMDS(dist.metric.veg, k = 2, maxit = 99, pc = TRUE)
tsne.mat <- as.data.frame(nmds$points)

## reconvert dist matrices into matrix objects
dist.metric <- as.matrix(dist.metric.veg)

## custom script to clean it up and rename variables to ones I like
source('./SMP/code/clean_ZEN_biotic.R')

## ggplots
colorpal <- c('#4daf4a', '#e41a1c', '#984ea3', '#377eb8')
ggplot(data = tsne.mat[!is.na(tsne.mat$sample.type), ], 
       aes(x = X1, y = X2, fill = as.factor(sample.type),
#            size = comm.coop,
           colour = as.factor(sample.type))) +
  stat_ellipse(data = tsne.mat[!is.na(tsne.mat$sample.type), ], 
               aes(x = X1, y = X2, colour = as.factor(sample.type)), 
               alpha = 0.25, geom = 'polygon',
               size = 0, linetype = 3, type = 't', level = 0.95) +
#   stat_ellipse(data = tsne.mat[!is.na(tsne.mat$sample.type), ], 
#                aes(x = X1, y = X2, group = as.factor(k.mean)), 
#                alpha = 0, geom = 'polygon', colour = 'gray40',
#                size = 0.4, linetype = 2, type = 't', level = 0.95) +
  geom_point(alpha = 0.75, 
#              colour = '#f7f7f7', 
             size = 3,          
             shape = 21) +
  theme_bw() + 
  theme_classic() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        panel.border = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  xlab('Axis 1') +
  ylab('Axis 2') +
#   scale_radius(range = c(2, 6.5)) +
  scale_colour_manual(values = colorpal, guide = guide_legend(title = NULL)) +
#   scale_shape_manual(values = c(21, 23, 22, 24), guide = guide_legend(title = NULL)) +
  scale_fill_manual(values = colorpal, guide = guide_legend(title = NULL)) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5))
```


#### Figure 1
<img src="figures/fig-1a.jpg" width=500, height=343)>

The points in Fig. 1 represent microbial communities, colored by sample type. Points that are closer together in NMDS space have more similar bacterial communities, or lower beta-diversity; ellipses represent 95% confidence intervals. The `cube-root` transformation yielded an ordination stress of `0.196`, indicating that an acceptable solution was found.

There is a clear pattern of differentiation among microbial communities found on different parts of the plant and environment. The composition of microbiomes from seagrass roots, soil and ocean water are significantly different (stats), whereas the composition of leaf microbiomes overlap with all three. 

Roughly, above- and belowground microbiomes are separated along the *x*-axis, while the *y*-axis separates plant-associated communities from environmental ones. Moreover, root and sediment microbiomes have higher beta-diversity than do leaves and water. Assuming that water and sediment are two potential sources of colonists for leaves and roots, this would be worth exploring further.

The results in Fig. 1b look qualitatively similar to the results of the Canberra ordination, which reveals about the drivers of community composition in the seagrass microbiome. Namely, it suggests that the OTUs that really differ in presence and relative abundance between our four sample types also exhibit a low level of shared phylogenetic branch length on average. If the differences between the four groups were driven primarily by differences in the relative abundance of close relatives, then we would not expect to see such comparable (to Canberra) beta-diversity patterns in NMDS space when the unweighted unifrac metric was applied.

#### Summary of section 2a
I analyzed beta-diversity in community composition using two different metrics: `Canberra`, which computes community dissimilarities using differences in the presence and relative abundance of OTUs and `unweighted Unifrac`, which weighs community dissimilarity based on the phylogenetic distance between OTUs that are present in the system. Because they yield similar beta-diversity patterns, we can conclude that the major differences between microbiomes sampled from different habitats (i.e., leaf, root, water or sediment) were due to differences in the presence and abundance of phylogenetically distinct OTUs. Moreover, these differences were large between soil and root communities but not between water and leaf communities. It is therefore conceivable that different assembly mechanisms operate on the *Z. marina* leaf and root microbiomes. I would like to evaluate this possibility. Do seagrass hosts really *'select'* advantageous bacteria, is their relative abundance directly linked to their availability in the environment and is this different for leaves and roots?

#### 2b. Defining a *'core'* microbiome
Here, I fit generalized linear models with binomial error distributions and logit link functions to the data in order to determine which OTUs were significantly associated with different sample types. I'm going to fit one model for each sample type as a function of the scaled relative abundance of each OTU using a one-vs-all scheme. I'll then adjust p-values for multiple comparisons using the Benjamini-Hochberg procedure, and select the model with a q-value < 0.01 and the largest positive regression coefficient as the strongest association. In other words, I'm modeling the probability that the community comes from each of the 4 sample types as a function of the relative abundance of each OTU. The analysis reveals about 1100 OTUs with non-random habitat associations. I propose that the OTUs significantly associated with Leaf and Root communities comprise their *'core'* microbiomes.

```r
pres.dat <- t(dat)

## drop OTUs appearing in < 10 sites for being underpowered
if(length(which(colSums(pres.dat != 0) < 10)) > 0){
  pres.dat <- pres.dat[, -c(which(colSums(pres.dat != 0) < 10))]
}

## save for distance-decay and plant-dis
pre.pre.pres.dat <- as.matrix(pres.dat)
pre.pre.pres.dat <- t(t(pre.pre.pres.dat)/rowSums(t(pre.pre.pres.dat)))

# ## turn into pres/abs
# pres.dat[pres.dat > 0] <- 1

## save for distance-decay and plant-dis
pre.pres.dat <- as.matrix(pres.dat)

## add metadata ### NEED TO CHANGE TO tsne.mat AND RERUN EVERYTHING INCLUDING MMINTE!!!
pres.dat <- cbind(pres.dat, tsne.mat[match(rownames(pres.dat), rownames(tsne.mat)), ])
# pres.frame <- data.frame('otu' = names(pres.dat[, -c((ncol(pres.dat) - ncol(tsne.mat.na)):ncol(pres.dat))]), 'coef' = rep(NA), 'p.adjust' = rep(NA), 'effect' = rep(NA))

## list
glm.list <- list()

## one vs all iterations
for(u in 1:(ncol(pres.dat) - ncol(tsne.mat))){
  tryCatch({
    
    pres.dat$glm.otu <- pres.dat[, u]
    
    ## one-vs-all scheme
    leaf.glm <- glm(as.factor(is.leaf) ~ scale(glm.otu), data = pres.dat, family = 'binomial')
    roots.glm <- glm(as.factor(is.roots) ~ scale(glm.otu), data = pres.dat, family = 'binomial')
    sediment.glm <- glm(as.factor(is.sediment) ~ scale(glm.otu), data = pres.dat, family = 'binomial')
    water.glm <- glm(as.factor(is.water) ~ scale(glm.otu), data = pres.dat, family = 'binomial')
    
    ## pick one with highest positive coefficient with sig p-value
    glm.frame <- data.frame(rbind(summary(leaf.glm)$coefficients[-1, ],
                       summary(roots.glm)$coefficients[-1, ],
                       summary(sediment.glm)$coefficients[-1, ],
                       summary(water.glm)$coefficients[-1, ]),
                       'otu' = as.character(rep(colnames(pres.dat)[u])))
    glm.frame$class <- c('Leaf', 'Roots', 'Sediment', 'Water')

    glm.list[[u]] <- glm.frame

#     ## organize into data frame
#     pres.frame$p.adjust[u] <- glm.frame$p.adjust
#     pres.frame$coef[u] <- glm.frame$Estimate
#     pres.frame$effect[u] <- glm.frame$class
    
  }, error = function(e){})
  if(u == 1){bar.vec <- c(na.omit(seq(1:(ncol(pres.dat)))[1:(ncol(pres.dat)) * round((ncol(pres.dat)) / 10)]))
             cat('|')}
  if(u %in% bar.vec == TRUE){cat('=====|')}
}

## summarize results
# pres.frame <- na.omit(pres.frame)
# pres.frame <- subset(pres.frame, p.adjust <= 0.01)

pres.frame <- do.call('rbind', glm.list)
names(pres.frame)[c(4, 6)] <- c('p', 'effect')
pres.frame$p.adjust <- p.adjust(pres.frame$p, 'fdr')
pres.frame <- subset(pres.frame, p.adjust <= 0.01)

## add taxonomy
pres.frame$phylum <- tax.2$phylum[match(pres.frame$otu, tax.2$otu)]
pres.frame$class <- tax.2$class[match(pres.frame$otu, tax.2$otu)]
pres.frame$order <- tax.2$order[match(pres.frame$otu, tax.2$otu)]
pres.frame$family <- tax.2$family[match(pres.frame$otu, tax.2$otu)]
pres.frame$genus <- tax.2$genus[match(pres.frame$otu, tax.2$otu)]
pres.frame$sp <- tax.2$sp[match(pres.frame$otu, tax.2$otu)]

## abs frame
abs.frame <- subset(pres.frame, Estimate < 0)
pres.frame <- subset(pres.frame, Estimate > 0)

## pull out argmaxes
pres.frame <- do.call('rbind', as.list(by(pres.frame, pres.frame$otu, 
                             function(x){
                               subset(x, otu == otu)[which.max(
                                 subset(x, otu == otu)$Estimate), ]
                               }))
)

abs.frame <- do.call('rbind', as.list(by(abs.frame, abs.frame$otu, 
                                          function(x){
                                            subset(x, otu == otu)[which.min(
                                              subset(x, otu == otu)$Estimate), ]
                                          }))
)
```

#### Fig. 2
<img src="figures/fig-2.jpg" width=300, height=250)>






