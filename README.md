# The ZEN microbiome
#### Ashkaan K Fahimipour
February 29, 2016

## Importing, sorting and cleaning the data

First, let's load a bunch of libraries and functions I've written over the last few months that do various microbiome things. These analyses are performed exclusively in R.

```
setwd('/Users/Ashkaan/Dropbox/')
source('./General Functions/microbiome_functions.R')
```
Let's import all of our data. These data are in the BIOM format, generated from merged libraries using macQIIME 1.9. OTUs were picked against the most recent *greengenes* database at a cutoff similarity of 97%. I rarefied samples to a depth of 100, and normalized OTU read counts by taxon copy number using the PICRUST software. Rarefying to this depth isn't ideal, but yields a decent sample of the data.  I find that the results presented below do not depend on the choice of rarefaction depth up to 500 or whether read counts were corrected for copy number. Since, R has a tough time with HDF5 formatted files, I converted the BIOM table into JSON using *biom convert* in macQIIME prior to these analyses. Contact me for bash scripts that take care of BIOM conversions.

```
## import data
biom <- read_biom('./SMP/data/ZEN_1/97_percent/rarefied_100/otu_table_json_even100_norm_by_copy.biom')
dat <- as.data.frame(as(biom_data(biom), "matrix"), header = TRUE)

## remove OTUs with zero abundance
dat <- dat[-which(rowSums(dat.rel) == 0), ]

## import metadata
tax <- as.matrix(observation_metadata(biom), sep = ',')
meta <- read.csv('./SMP/data/ZEN_1/97_percent/map_for_R.csv')
meta <- as.data.frame(meta)
biotic <- read.csv('./SMP/data/ZEN_biotic_data_2016_01_20.csv', header = TRUE)
```

I want to clean up the OTU table and metadata a little. Most of this is filtering out data that aren't useful to us, and renaming some of the metadata columns.

```
## clean up meta file
meta$SubsiteNumber[which(meta$SubsiteNumber == 1)] <- 'A'
meta$SubsiteNumber[which(meta$SubsiteNumber == 2)] <- 'B'
meta$sub.code <- paste(meta$ZenSite, meta$SubsiteNumber, sep = '.')

## trim taxonomy list to OTUs that are present
tax.trim <- as.matrix(tax[which(rownames(tax) %in% rownames(dat)), ])
tax.2 <- tax.to.long(tax.trim)

## exclude epiphyte samples
dont.use <- meta$SampleID[which(meta$SampleType == 'Epiphyte')]
dat <- dat[, -which(names(dat) %in% dont.use)]

## filter out chloroplasts
chloros <- tax.2$otu[which(tax.2$class == 'Chloroplast')]
dat <- dat[-which(rownames(dat) %in% chloros), ]
```

Finally, the last thing I want to do before analysis is to normalize my BIOM table, and to transform the relative abundances for more well-behaved residuals downstream. I find log(1 + x) and root transforms are good for these sorts of data.

```
## add relative abundances to new biom table
## computing row-wise is so much faster!
dat.rel <- t(t(dat)/rowSums(t(dat)))

## transform data to improve residuals
dat.rel.trans <- log(1 + dat.rel)
```

## 1. Community Ordination
I want to get an initial sense of bacterial community similarity across our different sample types (Leaf, Root, Sediment and Water). To do this, I'm first going to compute Bray-Curtis *distances* of log(1 + x) transformed relative abundances between community pairs. I'll then project these into 2-dimensions using t-distributed stochastic neighbor embedding, or t-SNE. t-SNE is sometimes thought of as a non-linear analogue of nMDS.

```
## distance matrix of transformed relative abundances
dist.metric <- vegdist(t(dat.rel.trans), method = 'bray')
t.sne <- tsne(dist.metric, k = 2, perplexity = 70, initial_dims = NA, max_iter = 1000, whiten = 0)
dist.metric <- as.matrix(dist.metric)

## make distance matrix
tsne.mat = data.frame(t.sne)
rownames(tsne.mat) <- rownames(dist.metric)

## custom script to clean up meta data and rename variables to ones that I like
source('./SMP/code/clean_ZEN_biotic.R')

```
Plotting our results tells us a lot about the structure of these data.

```
colorpal <- c('#006837', '#de2d26', '#8856a7', '#3182bd')
ggplot(data = tsne.mat.na, aes(x = X1, y = X2, fill = as.factor(sample.type))) + 
  stat_ellipse(alpha = 0.2, geom = 'polygon', size = .1, linetype = 1, type = 't', level = 0.95) +
  geom_point(shape = 21, alpha = 0.65, size = 5) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        panel.border = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  xlab('t-SNE 1') +
  ylab('t-SNE 2') +
  scale_radius(range = c(1.5, 7.5)) +
  scale_colour_manual(values = colorpal, guide = guide_legend(title = NULL)) +
  scale_fill_manual(values = colorpal, guide = guide_legend(title = NULL)) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
  	axis.ticks = element_blank())
```

#### Figure 1

![Fig. 1](figures/tSNE.pdf "tSNE")

In the Fig. 1, points represent bacterial communities, colored by sample type. Points that are closer together in t-SNE space have more similar bacterial communities. The ellipses represent 95% confidence intervals. The first thing I notice is that there is pretty clear differentiation among bacterial communities found on different parts of the plant and environment. Roughly, the *x*-axis appears to differentiate above- and belowground communities (Fig. 2), while the *y*-axis appears to seperate plant-associated communities from environmental ones (Fig. 3). Moreover, it is clear that roots and sediment seperate more than leaves and water. It would be great to know more about why this is.


#### Figure 2
![Fig. 2](figures/above_tSNE.pdf "tSNE")

#### Figure 3
![Fig. 3](figures/host_tSNE.pdf "tSNE")

### 1a. Permutational MANOVA
Something we might want to know right away is whether any of our (a)biotic variables correlate with axes scores from our community ordination. There are a couple of ways to do this. The first is with a PERMANOVA test. Although I think there are some limitations associated with this sort of analyses, it could be worth doing and will at least give us some intuition about the data. I'll rely on the *adonis* function in the *vegan* package for PERMANOVAs. First we need to merge our metadata with our OTU table. To pare down the number of possible covariates in our analyses, I performed a quick PCA on the metadata, and selected variables that loaded on the first 3 components.

```
pca.all <- princomp(tsne.mat.na[, c(10:14, 16:45)])
loadings(pca.all)
```

The first 3 components load onto 7 variables. These are:

1. Longest Leaf (cm)
2. *Z. marina* above ground biomass
3. *Z. marina* below ground biomass
4. Mean *Z. marina* shoots density
5. Mean mesograzer biomass
6. Mean macroalgae biomass
7. Standardized crustacean biomass

We'll be including these in our PERMANOVA along with latitude, longitude and site ID...

```
## first add sample types to transposed biom table
adonis.dat <- as.data.frame(t(dat.rel.trans))
adonis.dat$sample.type <- rep(0)
adonis.dat$site <- rep(0)
adonis.dat$lat <- rep(0)
adonis.dat$lon <- rep(0)
adonis.dat$longest.leaf.cm <- rep(0)
adonis.dat$zmarina.above.biomass <- rep(0)
adonis.dat$zmarina.below.biomass <- rep(0)
adonis.dat$mean.zmarina.shoots.m2 <- rep(0)
adonis.dat$mean.mesograzer.b <- rep(0)
adonis.dat$std.crustacean.b <- rep(0)
adonis.dat$mean.macroalgae <- rep(0)
```

Now, let's loop through sites to make sure everyone gets matched with their proper metadata.

```

## for loop
for(q in 1:length(rownames(adonis.dat))){
  curr.site <- as.character(rownames(adonis.dat)[q])
  curr.dat <- subset(tsne.mat.na, rownames(tsne.mat.na) == curr.site)
  adonis.dat$sample.type[q] <- as.character(curr.dat$sample.type[1])
  adonis.dat$site[q] <- as.character(curr.dat$site[1])
  adonis.dat$lat[q] <- as.numeric(as.character(curr.dat$lat[1]))
  adonis.dat$lon[q] <- as.numeric(as.character(curr.dat$lon[1]))
  adonis.dat$longest.leaf.cm[q] <- as.numeric(curr.dat$longest.leaf.cm[1])
  adonis.dat$zmarina.above.biomass[q] <- as.numeric(curr.dat$zmarina.above.biomass[1])
  adonis.dat$zmarina.below.biomass[q] <- as.numeric(curr.dat$zmarina.below.biomass[1])
  adonis.dat$mean.zmarina.shoots.m2[q] <- as.numeric(curr.dat$mean.zmarina.shoots.m2[1])
  adonis.dat$mean.mesograzer.b[q] <- as.numeric(curr.dat$mean.mesograzer.b[1])
  adonis.dat$std.crustacean.b[q] <- as.numeric(curr.dat$std.crustacean.b[1])
  adonis.dat$mean.macroalgae[q] <- as.numeric(curr.dat$mean.macroalgae[1])
  if(q == 1){bar.vec <- c(na.omit(seq(1:length(rownames(adonis.dat)))[1:length(rownames(adonis.dat)) * round(length(rownames(adonis.dat)) / 10)]))
             cat('|')}
  if(q %in% bar.vec == TRUE){cat('=====|')}
}
adonis.dat <- na.omit(adonis.dat)

```

Now that our data are in the right format, let's actually perform the PERMANOVA. I'm going to limit interactions in our model to *sample type* by *variable* interactions, since that is of primary interest right now; i.e., we want to know how community composition depends on different covariates, and how these relationships depend on *sample type*. We're going to constrain permutations to within sites by passing the *site* variable to the *strata* argument.


```
ad.mod <- adonis(adonis.dat[, -c(as.numeric(ncol(adonis.dat) - 11):as.numeric(ncol(adonis.dat)))] ~ 
                   adonis.dat$sample.type * adonis.dat$lat +
                   adonis.dat$sample.type * adonis.dat$lon +
                   adonis.dat$sample.type * adonis.dat$site +
                   adonis.dat$sample.type * adonis.dat$longest.leaf.cm +
                   adonis.dat$sample.type * adonis.dat$zmarina.above.biomass +
                   adonis.dat$sample.type * adonis.dat$zmarina.below.biomass +
                   adonis.dat$sample.type * adonis.dat$mean.zmarina.shoots.m2 +
                   adonis.dat$sample.type * adonis.dat$mean.mesograzer.b +
                   adonis.dat$sample.type * adonis.dat$std.crustacean.b +
                   adonis.dat$sample.type * adonis.dat$mean.macroalgae, 
                 method = 'bray', strata = adonis.dat$site, nperm = 999)
ad.mod

```

![Fig. 4](figures/permanova_table.png =x400 "PERMANOVA")

### Interpretation
A lot of covariates are significantly correlated with community composition; essentially all of them except *mean macroalgae* are correlated with community composition in all bacterial samples. Most of the variation is explained by *site* and the *sample type x site* interaction, indicating strong among-site variation in microbiome composition. While highly significant, the other covariates explain less than 24% of the residual variation all together. We can also see that microbiomes of different *sample types* are correlated to these data in different ways (i.e., *sample type by x* interactions ). But, the PERMANOVA framework doesn't permit any sort of post-hoc analyses to tell us *how* they differ. For this reason, I'll take an alternative approach similar to Kembel et al. (2011) PNAS and look for correlations between these data and axes scores resulting from t-SNE ordination of microbial community compositions.
 
### 1b. Linear models vs. axes scores
Another way to explore the relationships between covariates and community compositions are to regress our covariates directly onto our axes scores from the t-SNE. Since there is a reasonable expectation that communities within the same site are non-independent, I will fit a Maximum Likelihood random intercept model using the *nlme* package.

```
library(nlme)

## 1st axis
mod.x1 <- lme(X1 ~ 
                as.factor(sample.type) * lat +
                as.factor(sample.type) * longest.leaf.cm +
                as.factor(sample.type) * zmarina.above.biomass +
                as.factor(sample.type) * zmarina.below.biomass +
                as.factor(sample.type) * mean.zmarina.shoots.m2 +
                as.factor(sample.type) * mean.mesograzer.b +
                as.factor(sample.type) * std.crustacean.b,
              random = ~1|site, data = tsne.mat.na, na.action = 'na.fail', method = 'ML')

```

This is our *global model* from which I will perform model selection. To do this, I'll compute *AIC* scores for every possible combination of covariates and interactions and select the model with the lowest *AIC* score. I'll automate this process with the *dredge* function in the *MuMIn* library. The best fit model is reported below.

```
## model selection
dredge.x1 <- dredge(mod.x1, trace = 2)
best.mod.x1 <- get.models(dredge.x1, 1)[[1]]
anova(best.mod.x1)
summary(best.mod.x1)

```

### Ordination *X*-axis
![Fig. 5](figures/x1_anova.png =x150)

Ordination scores on t-SNE axis 1 are correlated with *sample type* and *mean shoot density*. There are also *sample type* by *longest leaf*, *mean shoot density* and *mesograzer biomass* interactions, indicating that these were correlated with different *sample types* in different ways. Let's explore the model results a bit closer.

![Fig. 6](figures/x1_summary.png =x260)

It looks like these variables are only correlated with environmental microbiomes, and not seagrass-associated ones. We'll keep that in mind as we continue to explore the data.

### Ordination *Y*-axis
I'll repeat the same thing for the 2nd t-SNE axis.

```
## 2nd axis
mod.x2 <- lme(X2 ~
                as.factor(sample.type) * lat +
                as.factor(sample.type) * longest.leaf.cm +
                as.factor(sample.type) * zmarina.above.biomass +
                as.factor(sample.type) * zmarina.below.biomass +
                as.factor(sample.type) * mean.zmarina.shoots.m2 +
                as.factor(sample.type) * mean.mesograzer.b +
                as.factor(sample.type) * std.crustacean.b,
              random = ~1|site, data = tsne.mat.na, na.action = 'na.fail', method = 'ML')

## model selection
dredge.x2 <- dredge(mod.x2, trace = 2)
best.mod.x2 <- get.models(dredge.x2, 1)[[1]]
anova(best.mod.x2)
summary(best.mod.x2)
```

