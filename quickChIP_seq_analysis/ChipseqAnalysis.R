
library(clusterProfiler)
library(UpSetR)
library(ChIPseeker)
library(ggpubr)
library(vroom)
library(dplyr)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(ReactomePA)
#PARAMETERS
bedfiles = list.files("/Volumes/LaCie/InstitutNeuroMyoGene_DATA/autre_data_caroline/Chipseq/GSE44756", pattern= ".bedGraph", full.names=T)
bedfiles[3] = list.files("/Volumes/LaCie/InstitutNeuroMyoGene_DATA/autre_data_caroline/Chipseq/GSE11062",pattern=".bed",full.names=T)
txdb = TxDb.Mmusculus.UCSC.mm9.knownGene
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
#READING FILES
peak = lapply(bedfiles, function(x) vroom(x,skip=1,col_names=c("chr","start","end","score")))
names(peak) = c("GSE44756_IP","GSE44756_CT","GSE11062")

#DATA MANAGEMENT vroom > GRanges
peakdf = lapply(peak, makeGRangesFromDataFrame)
peakdf[[1]]$score = peak[[1]]$score #25817258
peakdf[[2]]$score = peak[[2]]$score #25816277
peakdf[[3]]$score = NA #5273
#FILTERING
summary(peakdf[[1]]$score) 
summary(peakdf[[2]]$score)

peakdf_filtered = list()
peakdf_filtered[[1]] = peakdf[[1]][peakdf[[1]]$score > 2*mean(peakdf[[1]]$score),] # 2154260
peakdf_filtered[[2]] = peakdf[[2]][peakdf[[2]]$score > 2*mean(peakdf[[2]]$score),] #2144286
peakdf_filtered[[3]] = peakdf[[3]]
names(peakdf_filtered) = c("GSE44756_IP","GSE44756_CT","GSE11062")
#cat "/Volumes/LaCie/InstitutNeuroMyoGene_DATA/autre_data_caroline/Chipseq/GSE44756/GSE44756_IP.bed" | sed '1,1d' | awk '{print $2 "\t" $3 "\t" $4 "\t" $7};' | bedtools -i > GSE44756_IP_merged.bed
peakAnnoList <- lapply(peakdf_filtered, annotatePeak, TxDb=txdb,tssRegion=c(-3000, 3000), verbose=FALSE)
names(peakAnnoList) = c("GSE44756_IP","GSE44756_CT","GSE11062")
bar = plotAnnoBar(peakAnnoList)
dist = plotDistToTSS(peakAnnoList)
tagMatrixList <- lapply(peakdf_filtered, getTagMatrix, windows=promoter)
avg = plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
vennplot1 = upsetplot(peakAnnoList[[1]], vennpie=TRUE)
vennplot2 = upsetplot(peakAnnoList[[2]], vennpie=TRUE)

#FILTERING PROMOTER
peakAnnoList_filtering = list()
peakAnnoList_filtering[[1]] = peakAnnoList$GSE44756_IP@anno[stringr::str_detect(string=peakAnnoList$GSE44756_IP@anno$annotation,pattern="3' UTR|5' UTR|Promoter"),]
peakAnnoList_filtering[[2]] = peakAnnoList$GSE44756_CT@anno[stringr::str_detect(string=peakAnnoList$GSE44756_CT@anno$annotation,pattern="3' UTR|5' UTR|Promoter"),]
peakAnnoList_filtering[[3]] = peakAnnoList$GSE11062@anno[stringr::str_detect(string=peakAnnoList$GSE11062@anno$annotation,pattern="3' UTR|5' UTR|Promoter"),]
names(peakAnnoList_filtering) = c("GSE44756_IP","GSE44756_CT","GSE11062")


genes = lapply(peakAnnoList_filtering, function(i) unique(as.data.frame(i)$geneId))
names(genes) = sub("_", "\n", names(genes))
compKEGG1 <- compareCluster(geneCluster   = genes,
                           fun           = "enrichKEGG",
                           organism = "mmu",
                           keyType = "kegg",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")
compKEGG1_filtered = filter(compKEGG1@compareClusterResult,!stringr::str_detect(string=Description,pattern="disease|ataxia|cancer|infection|osis|carcinoma|Hepatitis|Glioma|leukemia|diabetic"))
#######################
compKEGG1_filtered$tmp = stringr::str_split(string=compKEGG1_filtered$GeneRatio,pattern="/")
pathway1 = dotplot(compKEGG1_filtered, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")

plot1 = ggarrange(dist,bar,avg,pathway1,nrow=2,ncol=2,labels=c("A","B","C","D"))
##################

genes = lapply(peakAnnoList, function(i) unique(as.data.frame(i)$geneId))
names(genes) = sub("_", "\n", names(genes))
compKEGG2 <- compareCluster(geneCluster   = genes,
                           fun           = "enrichKEGG",
                           organism = "mmu",
                           keyType = "kegg",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")
pathway2 = dotplot(compKEGG2, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")
plot2 = ggarrange(pathway1,pathway2,ncol=2,labels=c("A","B"),common.legend = T)




