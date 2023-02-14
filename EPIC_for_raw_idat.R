
#EPIC calculations for Illumina Human Methylation EPIC/850K by using raw .idat files

#adjust session directory, do not forget to put input file (.csv) which will be used by R as a guide

#put .idat files for each sample into the same directory

#install and load necessary libraries

# BiocManager::install("xxxxxxx")

#####################################################
library(minfi)
library(ggplot2)
library(DNAmArray)
library(minfiData)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(dplyr)
library(Gviz)
library(GenomicRanges)
library(ComplexHeatmap)

baseDir <- system.file("extdata", package = "minfiData")
targets <- read.metharray.sheet(baseDir)

RGSet <- read.metharray.exp(targets = targets,extended=TRUE)

#Normalization quantile
RGSet.norm <- preprocessQuantile(RGSet)

#Filtering probes
RGSet_filtered=probeFiltering(RGSet)

#Normalization Functional Normalization

RGSet.funnorm <- preprocessFunnorm(RGSet,nPCs=4)
mybeta <- reduce(RGSet.funnorm,RGSet_filtered,what="beta")


######################################################
##annotation of CpGs to genome

data(Locations)
positions=data.frame(row.names=Locations@rownames,
                     chr=Locations@listData$chr,
                     start=Locations@listData$pos,
                     end=Locations@listData$pos+1,
                     strand=Locations@listData$strand)

positions.gr=makeGRangesFromDataFrame(positions)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

positions.annot=annotatePeak(positions.gr,tssRegion=c(-2500, 2500),
                             TxDb=txdb, annoDb="org.Hs.eg.db")
positions.annot.df=as.data.frame(positions.annot,row.names=positions.annot@anno@ranges@NAMES)

######select PCDHA and promoter

PD <- positions.annot.df[grep("PCDHA", positions.annot.df$SYMBOL), ]
PD2 <- PD[,c('annotation', 'SYMBOL', 'start', 'end')]
PD2 <- PD2[grep("Promoter", PD$annotation), ]
PD2_unique <- distinct(PD2, SYMBOL, .keep_all= TRUE)

#combine individual beta values (beta) with pd2, from UCSC script
combo <- merge(PD2, mybeta, by = 0)

#make average beta value for each pcdha gene promoter
sorted_combo <- combo[order(combo$SYMBOL),]
j <- 1
average_table <- data.frame("1","2","3","4","5","6","7","8","9","10")

for (i in 1:nrow(sorted_combo)){
  s <- sorted_combo[i,3]
  print(s)
  x <- grep(s,sorted_combo$SYMBOL) 
  print(x)
  average <- c(0,0,0,0,0,0,0,0,0,0)
  for (k in x){
    for (m in 6:10){
      average[m]=average[m]+sorted_combo[k,m]
    }
  }
  average = average/length(x)
  average[3] = s
  average_table[j,]=average
  j = j + 1
}

average_table <- unique(average_table)
average_table <- average_table[,-c(1,2,4,5)]
names(average_table)[names(average_table) == 'X.3.'] <- 'SYMBOL' #rename column
data <- merge(PD2_unique,average_table, by = "SYMBOL")
names(data)[names(data) == 'X.6.'] <- 'Healthy 3' #rename column
names(data)[names(data) == 'X.7.'] <- 'Healthy 4' #rename column
names(data)[names(data) == 'X.8.'] <- 'Healthy 6' #rename column
names(data)[names(data) == 'X.9.'] <- 'ICF3 7' #rename column
names(data)[names(data) == 'X.10.'] <- 'ICF3 8' #rename column
##################################


##########
#hg19
genome <- "hg19"
chrom <- "chr"
from <- 140153729
to <- 140374011

axTrack <- GenomeAxisTrack()
idxTrack <- IdeogramTrack(genome = "hg19", chromosome = "chr5")

#plotTracks(list(idxTrack, axTrack, knownGenes), from = from, to = to, showTitle = FALSE)

refGenes <- UcscTrack(genome = "hg19", chromosome = "chr5",
                      track = "xenoRefGene", from = from, to = to,
                      trackType = "GeneRegionTrack", rstarts = "exonStarts",
                      rends = "exonEnds", gene = "name", symbol = "name2",
                      transcript = "name", strand = "strand", fill = "black",
                      stacking = "dense", name = "Pcdha")

cpgIslands <- UcscTrack(genome = "hg19", chromosome = "chr5", 
                        track = "cpgIslandExt", from = from, to = to,
                        trackType = "AnnotationTrack", 
                        start = "chromStart", end = "chromEnd", 
                        id = "name", shape = "box", 
                        name = "CpGs")
plotTracks(cpgIslands)

## make data track
dTrack <- DataTrack(data,chromosome = "chr5", genome = "hg19", name = "Methylation", groups = rep(c("Healthy 3", "Healthy 4", "Healthy 6","ICF3 7","ICF3 8")),type = c ("p", "a"))

plotTracks(list(idxTrack, axTrack, refGenes, cpgIslands, dTrack), from = from, to = to, showTitle = TRUE)
plotTracks(list(idxTrack, axTrack, cpgIslands, dTrack), from = from, to = to, showTitle = TRUE)
#########################

###### making Heatmaps with beta values##################
PD <- positions.annot.df[grep("PCDHA", positions.annot.df$SYMBOL), ]
combo <- merge(PD, mybeta, by = 0)

names(combo)[names(combo) == '203740920087_R04C01'] <- 'Healthy 4' #rename column
names(combo)[names(combo) == '203740920087_R05C01'] <- 'Healthy 5' #rename column
names(combo)[names(combo) == '202764140117_R01C01'] <- 'Healthy 1' #rename column
names(combo)[names(combo) == '202764140117_R02C01'] <- 'ICF3 2' #rename column
names(combo)[names(combo) == '202764140117_R05C01'] <- 'ICF3 5' #rename column

#select the beta values, Healthy + ICF3 cpgs,  and gene name from df
combo2 <- combo[,c('Healthy 4','Healthy 5','Healthy 1', 'ICF3 2', 'ICF3 5')]
combo2_matrix <- data.matrix(combo2)
Heatmap(combo2_matrix, name = "beta", column_title = "Pcdha",show_row_names = FALSE)
