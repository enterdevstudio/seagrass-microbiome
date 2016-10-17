## functions for later.
## load libraries. will use 'em all eventually
library(biom)
library(ggplot2)
library(mppa)
library(reshape2)
library(vegan)
# library(fields)
library(e1071)
# library(biotools)
library(FNN)
library(plyr)
library(lda)
library(dplyr)
library(minerva)
library(gplots)
library(ellipse)
library(corpcor)
library(RColorBrewer)
library(ccrepe)
library(igraph)
library(SpiecEasi)
library(Rtsne)
library(MuMIn)
library(stats)
library(indicspecies)
library(edgeR)
library(ape)
library(phyloseq)
library(DESeq2)

# detach('package:fields', unload=TRUE)
# detach('package:spam', unload=TRUE)

## function to take genome and load filef from directory
load.edge <- function(x){edge.list <- read.csv(paste('/Users/Ashkaan/Dropbox/ZEN-microbiome/data/metabolic-models/', x, '.csv', sep = ''))
                         edge.list
}

## geometric mean
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

## function to subsample comp matrix for leaf, sediment, root or water
subsample.tax.matrix <- function(taxon.list, meta.matrix){
  site.tax <- taxon.list
  new.adj <- matrix(nrow = length(site.tax), ncol = length(site.tax), 0)
  colnames(new.adj) <- site.tax
  rownames(new.adj) <- site.tax
  for(i in 1:length(site.tax)){
    for(j in 1:length(site.tax)){
      tryCatch({
      row = site.tax[i]
      col = site.tax[j]
      new.adj[row, col] <- meta.matrix[row, col]
      }, error = function(e){})
    }
    cat('Matrix construction is', round(i/length(site.tax)*100, 2), '% complete\n')
  }
  new.adj
}

## function to subsample comp matrix
comp.to.samp <- function(comp, site.id){
  curr.site <- sub.dat.rel[, which(colnames(sub.dat.rel) == site.id)]
  site.tax <- names(curr.site[curr.site > 0])
  new.adj <- matrix(nrow = length(site.tax), ncol = length(site.tax), 0)
  colnames(new.adj) <- site.tax
  rownames(new.adj) <- site.tax
  for(i in 1:length(site.tax)){
    for(j in 1:length(site.tax)){
      row = site.tax[i]
      col = site.tax[j]
      new.adj[row, col] <- competition.mat[row,col]
    }
#     cat('Matrix construction is', round(i/length(site.tax)*100, 2), '% complete\n')
  }
  new.adj
}

## function to subsample jaccard matrix
jacc.to.samp <- function(jacc, site.id){
  curr.site <- sub.dat.rel[, which(colnames(sub.dat.rel) == site.id)]
  site.tax <- names(curr.site[curr.site > 0])
  new.adj <- matrix(nrow = length(site.tax), ncol = length(site.tax), 0)
  colnames(new.adj) <- site.tax
  rownames(new.adj) <- site.tax
  for(i in 1:length(site.tax)){
    for(j in 1:length(site.tax)){
      row = site.tax[i]
      col = site.tax[j]
      new.adj[row, col] <- jacc.mat[row,col]
    }
    #     cat('Matrix construction is', round(i/length(site.tax)*100, 2), '% complete\n')
  }
  new.adj
}


## function to match blasted species IDs to OTU IDs
otu.matcher <- function(cor.table, blast.out, interaction.table){
  final.table <- cor.table
  int <- interaction
  blast <- blast
  final.table$sppA <- rep(0)
  final.table$sppB <- rep(0)
  
  for(g in 1:length(final.table$sppA)){
    final.table$sppA[g] <- as.character(blast$Species_ID[which(blast$Query_Otu_ID == as.character(final.table$otuA[g]))])
    final.table$sppB[g] <- as.character(blast$Species_ID[which(blast$Query_Otu_ID == as.character(final.table$otuB[g]))])
    cat('Assignment is', round(g / length(final.table$sppA)*100, 2), '% done\n')
  }
  
  final.table$GR.A.solo <- rep(0)
  final.table$GR.B.solo <- rep(0)
  final.table$GR.A.full <- rep(0)
  final.table$GR.B.full <- rep(0)
  final.table$interaction <- rep(0)
  
  for(h in 1:length(final.table$GR.A.solo)){
    final.table$GR.A.solo[h] <- int$GRASolo[which(int$GenomeIDSpeciesA == final.table$sppA[h] & int$GenomeIDSpeciesB == final.table$sppB[h])]
    final.table$GR.B.solo[h] <- int$GRBSolo[which(int$GenomeIDSpeciesA == final.table$sppA[h] & int$GenomeIDSpeciesB == final.table$sppB[h])]
    final.table$GR.A.full[h] <- int$GRSpeciesAFull[which(int$GenomeIDSpeciesA == final.table$sppA[h] & int$GenomeIDSpeciesB == final.table$sppB[h])]
    final.table$GR.B.full[h] <- int$GRSpeciesBFull[which(int$GenomeIDSpeciesA == final.table$sppA[h] & int$GenomeIDSpeciesB == final.table$sppB[h])]
    final.table$interaction[h] <- as.character(int$Type.of.interaction[which(int$GenomeIDSpeciesA == final.table$sppA[h] & int$GenomeIDSpeciesB == final.table$sppB[h])])
    cat('Matching growth rates is', round(h / length(final.table$sppA)*100, 2), '% done\n')
  }
  final.table
}


## autocurve function
autocurve.edges2 <- function (graph, start = 0.5)
{cm <- count.multiple(graph)
 mut <- is.mutual(graph)  #are connections mutual?
 el <- apply(get.edgelist(graph, names = FALSE), 1, paste, collapse = ":")
 ord <- order(el)
 res <- numeric(length(ord))
 p <- 1
 while (p <= length(res)) {
   m <- cm[ord[p]]
   mut.obs <-mut[ord[p]] #are the connections mutual for this point?
   idx <- p:(p + m - 1)
   if (m == 1 & mut.obs == FALSE) { #no mutual conn = no curve
     r <- 0}
   else {r <- seq(-start, start, length = m)}
   res[ord[idx]] <- r
   p <- p + m}
 res
}

## function to merge OTUs by family
merge.by.fam <- function(mat){
  ## replace OTU labels with families
  m <- mat
  otus <- rownames(m)
  these.otus <- tax[match(otus, rownames(tax))]
  tax.dat <- as.data.frame(matrix(ncol = 9, nrow = length(these.otus)))
  names(tax.dat) <- c('otu', 'highest', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'sp')
  tax.dat$otu <- otus
  tax.dat$highest <- rep(999)
  ## seperate string by ';' and unlist
  for(j in 1:length(these.otus)){ 
    tax.dat$kingdom[j] <- unlist(strsplit(as.character(these.otus[[j]][1]), '__'))[2]
    tax.dat$phylum[j] <- unlist(strsplit(as.character(these.otus[[j]][2]), '__'))[2]
    tax.dat$class[j] <- unlist(strsplit(as.character(these.otus[[j]][3]), '__'))[2]
    tax.dat$order[j] <- unlist(strsplit(as.character(these.otus[[j]][4]), '__'))[2]
    tax.dat$family[j] <- unlist(strsplit(as.character(these.otus[[j]][5]), '__'))[2]
    tax.dat$genus[j] <- NA
    tax.dat$sp[j] <- NA
    cat('Loop is', round(j/length(these.otus)*100, 2), '% done\n')
  }
  ## function to find highest taxonomic resolution
  for(y in 1:length(tax.dat$otu)){tax.dat$highest[y] <- tax.dat[y, min(which(is.na(tax.dat[y, ]))) - 1]}
  ## add highest tax to m
  m$highest <- tax.dat$highest
  m.merged <- ddply(m, 'highest', numcolwise(sum))
  rownames(m.merged) <- m.merged$highest
  m.merged <- m.merged[, -1]
  m.merged
}

merge.by.ord <- function(mat){
  ## replace OTU labels with families
  m <- mat
  otus <- rownames(m)
  these.otus <- tax[match(otus, rownames(tax))]
  tax.dat <- as.data.frame(matrix(ncol = 9, nrow = length(these.otus)))
  names(tax.dat) <- c('otu', 'highest', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'sp')
  tax.dat$otu <- otus
  tax.dat$highest <- rep(999)
  ## seperate string by ';' and unlist
  for(j in 1:length(these.otus)){ 
    tax.dat$kingdom[j] <- unlist(strsplit(as.character(these.otus[[j]][1]), '__'))[2]
    tax.dat$phylum[j] <- unlist(strsplit(as.character(these.otus[[j]][2]), '__'))[2]
    tax.dat$class[j] <- unlist(strsplit(as.character(these.otus[[j]][3]), '__'))[2]
    tax.dat$order[j] <- unlist(strsplit(as.character(these.otus[[j]][4]), '__'))[2]
    tax.dat$family[j] <- NA
    tax.dat$genus[j] <- NA
    tax.dat$sp[j] <- NA
    cat('Loop is', round(j/length(these.otus)*100, 2), '% done\n')
  }
  ## function to find highest taxonomic resolution
  for(y in 1:length(tax.dat$otu)){tax.dat$highest[y] <- tax.dat[y, min(which(is.na(tax.dat[y, ]))) - 1]}
  ## add highest tax to m
  m$highest <- tax.dat$highest
  m.merged <- ddply(m, 'highest', numcolwise(sum))
  rownames(m.merged) <- m.merged$highest
  m.merged <- m.merged[, -1]
  m.merged
}

merge.by.class <- function(mat){
  ## replace OTU labels with families
  m <- mat
  otus <- rownames(m)
  these.otus <- tax[match(otus, rownames(tax))]
  tax.dat <- as.data.frame(matrix(ncol = 9, nrow = length(these.otus)))
  names(tax.dat) <- c('otu', 'highest', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'sp')
  tax.dat$otu <- otus
  tax.dat$highest <- rep(999)
  ## seperate string by ';' and unlist
  for(j in 1:length(these.otus)){ 
    tax.dat$kingdom[j] <- unlist(strsplit(as.character(these.otus[[j]][1]), '__'))[2]
    tax.dat$phylum[j] <- unlist(strsplit(as.character(these.otus[[j]][2]), '__'))[2]
    tax.dat$class[j] <- unlist(strsplit(as.character(these.otus[[j]][3]), '__'))[2]
    tax.dat$order[j] <- NA
    tax.dat$family[j] <- NA
    tax.dat$genus[j] <- NA
    tax.dat$sp[j] <- NA
    cat('Loop is', round(j/length(these.otus)*100, 2), '% done\n')
  }
  ## function to find highest taxonomic resolution
  for(y in 1:length(tax.dat$otu)){tax.dat$highest[y] <- tax.dat[y, min(which(is.na(tax.dat[y, ]))) - 1]}
  ## add highest tax to m
  m$highest <- tax.dat$highest
  m.merged <- ddply(m, 'highest', numcolwise(sum))
  rownames(m.merged) <- m.merged$highest
  m.merged <- m.merged[, -1]
  m.merged
}

merge.by.phy <- function(mat){
  ## replace OTU labels with families
  m <- mat
  otus <- rownames(m)
  these.otus <- tax[match(otus, rownames(tax))]
  tax.dat <- as.data.frame(matrix(ncol = 9, nrow = length(these.otus)))
  names(tax.dat) <- c('otu', 'highest', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'sp')
  tax.dat$otu <- otus
  tax.dat$highest <- rep(999)
  ## seperate string by ';' and unlist
  for(j in 1:length(these.otus)){ 
    tax.dat$kingdom[j] <- unlist(strsplit(as.character(these.otus[[j]][1]), '__'))[2]
    tax.dat$phylum[j] <- unlist(strsplit(as.character(these.otus[[j]][2]), '__'))[2]
    tax.dat$class[j] <- NA
    tax.dat$order[j] <- NA
    tax.dat$family[j] <- NA
    tax.dat$genus[j] <- NA
    tax.dat$sp[j] <- NA
    cat('Loop is', round(j/length(these.otus)*100, 2), '% done\n')
  }
  ## function to find highest taxonomic resolution
  for(y in 1:length(tax.dat$otu)){tax.dat$highest[y] <- tax.dat[y, min(which(is.na(tax.dat[y, ]))) - 1]}
  ## add highest tax to m
  m$highest <- tax.dat$highest
  m.merged <- ddply(m, 'highest', numcolwise(sum))
  rownames(m.merged) <- m.merged$highest
  m.merged <- m.merged[, -1]
  m.merged
}

## for ccrepe:
bray.sim <- function(x, y = NA){
  if(is.vector(x) && is.vector(y)) return(as.matrix(vegdist(t(data.frame('x' = x, 'y' = y)), method = 'bray'))[2])
  if(is.matrix(x) && is.na(y)) return(as.matrix(vegdist(t(x), method = 'bray')))
  if(is.data.frame(x) && is.na(y)) return(as.matrix(vegdist(t(x), method = 'bray')))
  else stop('ERROR')
}

## function to group nodes based on community membership
layout.modular <- function(G,c){
  nm <- length(levels(as.factor(c)))
  gr <- 2
  while(gr^2 < nm){
    gr <- gr + 1
  }
  i <- j <- 0
  for(cc in levels(as.factor(c))){
    F <- delete.vertices(G, c != cc)
    F$layout <- layout.fruchterman.reingold(F)
    F$layout <- layout.norm(F$layout, i, i + 0.5, j, j + 0.5)
    G$layout[c == cc,] <- F$layout
    if(i == gr){
      i <- 0
      if(j == gr){
        j <- 0
      } else{
        j <- j + 1
      }
    }else{
      i <- i + 1
    }
  }
  return(G$layout)
}

## function to specify sample ID and generate adjacency matrix from meta-network
meta.to.samp <- function(meta, site.id){
  curr.site <- sub.dat.rel[, which(colnames(sub.dat.rel) == site.id)]
  site.tax <- names(curr.site[curr.site > 0])
  new.adj <- matrix(nrow = length(site.tax), ncol = length(site.tax), 0)
  colnames(new.adj) <- site.tax
  rownames(new.adj) <- site.tax
  for(i in 1:length(site.tax)){
    for(j in 1:length(site.tax)){
      row = site.tax[i]
      col = site.tax[j]
      tryCatch({
      new.adj[row, col] <- merged.dat[row,col]
      }, error = function(e){})
    }
#     cat('Matrix construction is', round(i/length(site.tax)*100, 2), '% complete\n')
  }
  new.adj
}

## function to take tissue type and generate adjacency matrix from meta-network
meta.to.type <- function(meta, type){
  curr.type <- sub.dat.rel[, which(colnames(sub.dat.rel) %in% rownames(subset(tsne.mat.na, sample.type == type)))]
  type.tax <- rownames(curr.type[rowSums(curr.type) > 0, ])
  new.adj <- matrix(nrow = length(type.tax), ncol = length(type.tax), 0)
  colnames(new.adj) <- type.tax
  rownames(new.adj) <- type.tax
  for(i in 1:length(type.tax)){
    for(j in 1:length(type.tax)){
      row = type.tax[i]
      col = type.tax[j]
#       tryCatch({
      new.adj[row, col] <- meta[row, col]
#       }, error = function(e){})
    }
        cat('Matrix construction is', round(i / length(type.tax)*100, 2), '% complete\n')
  }
  new.adj
}


## function to subset metanetwork in edge list format
meta.edge.to.samp <- function(edge.list, site.id){
  v.list <- rownames(sub.dat.rel)[which(sub.dat.rel[, site.id] > 0)]
  contained <- data.frame('a' = as.numeric(edge.list$OtuIDA %in% v.list), 'b' = as.numeric(edge.list$OtuIDB %in% v.list))
  keepers <- c(which(rowSums(contained) == 2))
  new.edge.list <- edge.list[keepers, ]
  new.edge.list
}

# meta$SubsiteNumber[which(meta$SubsiteNumber == 1)] <- 'A'
# meta$SubsiteNumber[which(meta$SubsiteNumber == 2)] <- 'B'
# meta$sub.code <- paste(meta$ZenSite, meta$SubsiteNumber, sep = '.')

rel.abund.nodes <- function(site.igraph.obj, site.id){
  rel.abund.size.vec <-c()
  site.dat <- sub.dat.rel[, site.id]
  curr.nodes <- names(V(site.igraph.obj))
  for(i in 1:length(curr.nodes)){
    rel.abund.size.vec[i] <- site.dat[which(names(site.dat) == curr.nodes[i])]  
  }
  rel.abund.size.vec
}

### tax table to long form function
tax.to.long <- function(mat){
  ## replace OTU labels with families
  m <- mat
  otus <- rownames(m)
  these.otus <- tax[match(otus, rownames(tax))]
  
  tax.dat <- as.data.frame(matrix(ncol = 8, nrow = length(these.otus)))
  names(tax.dat) <- c('otu', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'sp')
  tax.dat$otu <- otus
  tax.dat$highest <- rep(999)
  
  ## seperate string by ';' and unlist
  for(j in 1:length(these.otus)){ 
    tax.dat$domain[j] <- unlist(strsplit(as.character(these.otus[[j]][1]), '__'))[2]
    tax.dat$phylum[j] <- unlist(strsplit(as.character(these.otus[[j]][2]), '__'))[2]
    tax.dat$class[j] <- unlist(strsplit(as.character(these.otus[[j]][3]), '__'))[2]
    tax.dat$order[j] <- unlist(strsplit(as.character(these.otus[[j]][4]), '__'))[2]
    tax.dat$family[j] <- unlist(strsplit(as.character(these.otus[[j]][5]), '__'))[2]
    tax.dat$genus[j] <- unlist(strsplit(as.character(these.otus[[j]][6]), '__'))[2]
    tax.dat$sp[j] <- unlist(strsplit(as.character(these.otus[[j]][7]), '__'))[2]
    if(j == 1){bar.vec <- c(na.omit(seq(1:length(these.otus))[1:length(these.otus) * round(length(these.otus) / 10)]))
               cat('|')}
    if(j %in% bar.vec == TRUE){cat('=====|')}
  }
  tax.dat
}

#################################################################################################################
# Function to test the significance of the species-level and community-level change in deviance.		#
# See Baeten et al. 2013 Methods in Ecology and Evolution for more details on the approach.			#
#														#
# Usage														#
# dDEV(data, sites, survey, surv.old, surv.new, permutations)							#
#														#
# Arguments of the function bDEV.test:										#
# data		the function expects a single site x species data matrix, i.e., with data from two time periods	#
# sites		vector specifying unique site names (the same name for the two time periods)			#
# survey	vector specifying the time period								#
# surv.old	name of the initial time period									#
# surv.new	name of the final time period									#
# permutations	number of permutations used to test the hypothesis of no temporal change (default = 2000)	#	
#														#
# Value														#
# The function returns a list with the first element summarizing the multiple-site results (number of plots,	#
# number of species, delta_D and significance). The second element shows species-level results (delta_Di and	#
# significance).  												#
#################################################################################################################

dDEV<-function(data,sites,survey,surv.old,surv.new,permutations=2000){
  calc.Di <- function(x){
    p <- sum(x)/nsites
    nsites * -2*( p*log(p+1.e-8)+(1-p)*log(1-p+1.e-8) )
  }	
  smpl <- function(){
    out <- matrix(NA,2,nsites*nspecies)
    time <- matrix(1:2,2,nsites*nspecies)
    rtime <- sample(1:2,nsites*nspecies,replace=T)
    out[,rtime==1] <- c(1,2)
    out[,rtime==2] <- c(2,1)
    matrix(out,2*nsites,nspecies)
  }
  t1<-surv.old;t2<-surv.new
  fsites <- factor(sites);fsurv <- factor(survey, levels=c(surv.old, surv.new))
  nsites <- nlevels(fsites);nspecies <- ncol(data)
  m <- as.matrix(data);om <- m[order(fsites, fsurv),]
  # observed deviance
  Di.t1 <- sapply(1:nspecies, function(x) calc.Di(m[fsurv==t1,x]))
  Di.t2 <- sapply(1:nspecies, function(x) calc.Di(m[fsurv==t2,x]))
  dDi <- Di.t2 - Di.t1
  dD <- sum(dDi)
  # random deviance
  rdDi <- matrix(NA,nrow=permutations,ncol=nspecies)
  for(i in 1:permutations){
    smpl.i <- smpl()
    rDi.t1 <- sapply(1:nspecies, function(x) calc.Di(om[smpl.i[,x]==1,x]))
    rDi.t2 <- sapply(1:nspecies, function(x) calc.Di(om[smpl.i[,x]==2,x]))
    rdDi[i,] <- rDi.t2 - rDi.t1
  }
  rdD <- apply(rdDi, 1, sum) 
  # significance
  eps <- 1.e-8
  rS.Pspecies <- rowSums( abs(t(rdDi)) > (abs(dDi) - eps) )
  Pspecies <- (rS.Pspecies + 1) / (permutations + 1)
  
  s.Pstudy <- sum( abs(rdD) > (abs(dD) - eps) )
  Pstudy <- (s.Pstudy + 1) / (permutations + 1)
  # output
  list(dD=data.frame(nsites, nspecies, dD, Pstudy), 
       dDi=data.frame(dDi, Pspecies))
  
}

##----------
## MC bootstrap functions
diffMeans <- function(data, hasTrt){
  # computes our test statistics: the difference of means
  # hasTrt: boolean vector, TRUE if has treatment
  test_stat <- mean(data[hasTrt]) - mean(data[!hasTrt])
  return(test_stat)
}

simPermDsn <- function(data, hasTrt, testStat, k=1000){
  # Simulates the permutation distribution for our data
  # hasTrt: boolean indicating whether group is treatment or control
  # testStat: function defining the test statistic to be computed
  # k: # of samples to build dsn.      
  sim_data   <- replicate(k, sample(data))
  test_stats <- apply(sim_data, 2,
                      function(x) {testStat(x, hasTrt)})
  return( test_stats)
}

permutationTest <- function(data, hasTrt, testStat, k=1000){
  # Performs permutation test for our data, given some pre-selected
  # test statistic, testStat
  currentStat    <- testStat(data, hasTrt)
  simulatedStats <- simPermDsn(data, hasTrt,testStat, k)
  
  # 2-way P-value
  pValue         <- sum(abs(simulatedStats) >= currentStat)  / k
  
  return(pValue)
}



