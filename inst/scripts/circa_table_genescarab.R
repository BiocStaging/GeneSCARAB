# Script for generation of inst/extdata/org_Otaurireduced_eg_db.sqlite.gz
# dataset

# This dataset was generated from the data of tha paper "Multiomics integration
# unveils photoperiodic plasticity in the molecular rhythms of marine 
# phytoplankton" following the script provided by the authors in 
# https://github.com/fran-romero-campero/SANDAL/blob/master/seasonal_and_diurnal_cycles_in_ostreococcus.Rmd,
# which we will summarize here

# First load expression data (downloaded from 
# https://github.com/fran-romero-campero/SANDAL/blob/master/gene_expression.tsv)
# and determine rhythmic genes in both conditions using the RAIN package
gene.expression <- read.table(file ="tables/gene_expression.tsv",sep="\t",header=T,as.is=T)
head(gene.expression[,1:7])

gene.ids <- gene.expression$X

gene.expression <- as.matrix(gene.expression[,2:ncol(gene.expression)])
rownames(gene.expression) <- gene.ids
head(gene.expression[,1:6])

ostta.genes <- rownames(gene.expression)
number.genes <- length(ostta.genes)
number.genes

length(which(apply(X = gene.expression,MARGIN = 1,FUN = max) == 0))
length(which(apply(X = gene.expression,MARGIN = 1,FUN = max) < 1))
length(which(apply(X = gene.expression,MARGIN = 1,FUN = max) < 5))
length(which(apply(X = gene.expression,MARGIN = 1,FUN = max) < 10))

ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
ld.zt.i <- sapply(X = ld.zt,FUN = function(x){ paste(x,1:3,sep="_")})
sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
sd.zt.i <- sapply(X = sd.zt,FUN = function(x){ paste(x,1:3,sep="_")})

ld.sd.gene.expression <- gene.expression[,c(ld.zt.i,sd.zt.i)]
head(ld.sd.gene.expression[,1:6])
dim(ld.sd.gene.expression)

ld.gene.expression <- gene.expression[,ld.zt.i]
sd.gene.expression <- gene.expression[,sd.zt.i]

# LD
library(rain)

results.ld <- rain(t(ld.gene.expression), deltat=4, period=24, verbose=FALSE, nr.series=3)
sum(results.ld$pVal < 0.05)/number.genes

rhythmic.genes.ld <- rownames(subset(results.ld, pVal < 0.05))
length(rhythmic.genes.ld)

results.ld.12 <- rain(t(ld.gene.expression), deltat=4, period=12, verbose=FALSE, nr.series=3)
sum(results.ld.12$pVal < 0.05)/number.genes
rhythmic.genes.ld.12 <- rownames(subset(results.ld.12, pVal < 0.05))
length(rhythmic.genes.ld.12)

complete.ld.rhythmic.genes <- unique(c(rhythmic.genes.ld,rhythmic.genes.ld.12))

# SD
results.sd <- rain(t(sd.gene.expression), deltat=4, period=24, verbose=FALSE, nr.series=3)
sum(results.sd$pVal < 0.05)/number.genes

rhythmic.genes.sd <- rownames(subset(results.sd, pVal < 0.05))
length(rhythmic.genes.sd)

results.sd.12 <- rain(t(sd.gene.expression), deltat=4, period=12, verbose=FALSE, nr.series=3)
sum(results.sd.12$pVal < 0.05)/number.genes
rhythmic.genes.sd.12 <- rownames(subset(results.sd.12, pVal < 0.05))
length(rhythmic.genes.sd.12)

complete.sd.rhythmic.genes <- unique(c(rhythmic.genes.sd,rhythmic.genes.sd.12))

complete.ld.sd.rhythmic.genes <- intersect(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes)

# Now infer phases for those genes using CircaCompare
library(circacompare)

ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")

circacompare.ld.sd <- matrix(nrow=length(complete.ld.sd.rhythmic.genes),ncol=15)
rownames(circacompare.ld.sd) <- complete.ld.sd.rhythmic.genes

for(i in 1:length(complete.ld.sd.rhythmic.genes))
{
  gene.i <- complete.ld.sd.rhythmic.genes[i]
  
  ld.expression.i <- gene.expression[gene.i,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"),paste(ld.zt,3,sep="_"))]
  sd.expression.i <- gene.expression[gene.i,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"),paste(sd.zt,3,sep="_"))]
  
  time.points <- seq(from=0,by=4,length.out = 18)
  
  ld.sd.df <- data.frame(time=c(time.points,time.points),
                         measure=c(ld.expression.i, sd.expression.i),
                         group=c(rep("ld",18),rep("sd",18)))
  
  out.i <- circacompare(x = ld.sd.df, col_time = "time", col_group = "group", col_outcome = "measure",alpha_threshold = 1)
  circacompare.ld.sd[i,] <- out.i[[2]][,2]
}

colnames(circacompare.ld.sd) <- out.i[[2]][,1]

# After that, create final df
circa_table <- data.frame(
  ld.peak.time.hours = circacompare.ld.sd["ld peak time hours"],
  sd.peak.time.hours = circacompare.ld.sd["sd peak time hours"])

# More accessory details can be found in the original script at
# https://github.com/fran-romero-campero/SANDAL/blob/master/seasonal_and_diurnal_cycles_in_ostreococcus.Rmd
