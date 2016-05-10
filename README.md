# <img src="figures/plant_vec.jpg" width=50, height=61)> The ZEN Microbiome Project

> May 10, 2016

### 1. Introduction 
The first step in understanding the symbiotic relationship between seagrass microbes and their host is to characterize the baseline microbiota and the differences that are associated with host genotype and environment. In the Moore Foundation proposal, we proposed to address the following set of questions:

- How much is microbial community composition influenced by genetic and ecological variation in the host, and is this different for different host tissues?
- Is there a ‘core microbiome’ or set of commonly occurring microorganisms that appear in all assemblages associated with seagrasses?
- Are there significant co-occurrence or co-exclusion relationships between pairs of microbial taxa within the seagrass microbiome?

### 2. Analyses
These analyses depend heavily on Python scripts bundled with [*macQIIME* v1.9](http://www.wernerlab.org/software/macqiime) and the *biom*, *ggplot2*, *vegan*, *igraph* and *ccrepe* packages in the statistical programming environment `R`.

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

Finally, I'll normalize and `cube-root` transform my BIOM table for analysis downstream. I found that this transformation yielded the most reliable ordination.

```r
## add relative abundances to new biom table
dat.rel <- t(t(dat)/rowSums(t(dat)))

## transform data
dat.rel.trans <- dat.rel^(1/3)
```

#### 2a. Community Ordination
I want to get an initial sense of bacterial community composition across our different sample types (Leaf, Root, Sediment and Water). To do this, I'm first going to compute Canberra *distances* of `cube-root` transformed relative abundances between community pairs. I'll then visualize these in 2-dimensions using nonmetric multidimensional scaling (left-hand Fig. 1a). I'll repeat this for nonweighted Unifrac distances computed in *QIIME* (right-hand Fig. 1b).

#### Figure 1

<img src="figures/fig-1a.jpg" width=400, height=275)> <img src="figures/fig-1b.jpg" width=400, height=275)>

The points in Fig. 1 represent microbial communities, colored by sample type. Points that are closer together in NMDS space have more similar bacterial communities, or lower beta-diversity. The ellipses represent 95% confidence intervals. The the `cube-root` transofmration yielded an ordination stress of `0.196`, indicating that an acceptable solution was found. 

There is a clear pattern of differentiation among microbial communities found on different parts of the plant and environment. Roughly, the *x*-axis seperates above- and belowground microbiomes, while the *y*-axis seperates plant-associated communities from environmental ones (Fig. 3). Moreover, root and sediment microbiomes have higher beta-diversity than do leaves and water; there is clear seperation between root and soil microbiomes in NMDS space, and high overlap between leaf and water communities. Assuming that water and sediment are the major sources of colonists for leaves and roots respectively, this would be worth exploring further.

The results in Fig. 1b look qualitatively similar to the results of the Canberra ordination, which reveals about the drivers of community composition in the seagrass microbiome. Namely, it suggests that the OTUs that really differ in presence and relative abundance between our four sample types also exhibit a low level of shared phylogenetic branch length on average. If the differences between the four groups were driven primarily by differences in the relative abundance of close relatives, then we would not expect to see such comparable (to Canberra) beta-diversity patterns in NMDS space when the unweighted unifrac metric was applied.

#### Summary of section 2a
I analyzed beta-diversity in community composition using two different metrics: `Canberra`, which computes community dissimilarities using differences in the presence and relative abundance of OTUs and `unweighted Unifrac`, which weighs community dissimilarity based on the phylogenetic distance between OTUs that are present in the system. Because they yield similar beta-diversity patterns, we can conclude that the major differences between microbiomes sampled from different habitats (i.e., leaf, root, water or sediment) were due to differences in the presence and abundance of phylogenetically distinct OTUs. Moreover, these differences were large between soil and root communities but not between water and leaf communities. If we assume that soil and water are the major colonist sources for roots and leaves, then it is conceivable that different assembly mechanisms operate on the *Z. marina* leaf and root microbiomes.



