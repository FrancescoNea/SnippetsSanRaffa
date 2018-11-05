#!/bin/sh



 #************************************************************DATABASE***************************************************************************
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
library("org.Hs.eg.db")
gtf <- rtracklayer::import('gencode.v26.chr_patch_hapl_scaff.annotation.gtf')
#subsetto gtf a seconda dei geni che troviamo nell'rnaseq
toptags_RNA=read.table('TopTags_RNA_DiMicco_Analisi_NuovoGTF.txt',head=T)
gencodegtf_rnaseq=gtf[(elementMetadata(gtf)[,"gene_id"] %in% toptags_RNA$genes)]
gencode_rnaseqtxdb=makeTxDbFromGRanges(gencodegtf_rnaseq)
saveRDS(gencode_rnaseqtxdb,'../../Documents/Pilot_DiMicco/Gencode_RNAseqTXDB.rds')


gencode_rnaseqtxdb=readRDS('Gencode_RNAseqTXDB.rds')

#************************************************************RNA-seq***************************************************************************

newcounts=featureCounts(files=c(
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_1.bam',
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_2.bam',
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_3.bam',
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_4.bam',
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_5.bam',
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_6.bam',
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_7.bam',
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_8.bam',
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_9.bam',
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_10.bam',
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_11.bam',
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_12.bam',
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_13.bam',
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_14.bam'),
annot.ext='/lustre1/genomes/hg38/annotation/gencode.v26.chr_patch_hapl_scaff.annotation.gtf',
isGTFAnnotationFile = TRUE,GTF.attrType='gene_name',chrAliases='Aliases_GTF.txt')



colnames(newcounts$counts)=c("UPN3_1","UPN3_2","UPN3_3","UPN3_4","UPN3_5","UPN3_6","UPN3_7","UPN3_8","UPN3_9","UPN3_10","UPN3_11","UPN3_12","UPN3_13","UPN3_14")

rna_newcounts = DGEList(newcounts$counts[,c(1,2,3,4,8,9,10,11)], genes = newcounts$annotation$GeneID)
rna_newcounts = calcNormFactors(rna_newcounts,method='TMM')
rna_newcounts$samples$group = as.factor(group[row.names(rna_newcounts$samples)])
rna_keep_newcounts = rowSums(cpm(rna_newcounts) > 1) >= 3
rna_newcounts = rna_newcounts[rna_keep_newcounts, ]


library(qvalue)
df = targets[row.names(rna_newcounts$samples),]
D = as.factor(df$disease)
P = as.factor(df$passage) 
#	m_design = model.matrix(~D * P)
design = model.matrix(~D+P)
 design
#  (Intercept) DR P1
#1           1  0  0
#2           1  0  1
#3           1  0  1
#4           1  0  1
#5           1  1  0
#6           1  1  1
#7           1  1  1
#8           1  1  1

rna_newcounts = estimateDisp(rna_newcounts, design)
fit = glmFit(rna_newcounts, design)
lrt <- glmLRT(fit, coef=2)
toptags_RNA=as.data.frame(topTags(lrt,sort.by="PValue",n="all"))
 q = qvalue(toptags_RNA$PValue)
toptags_RNA$Qvalue2 = q$qvalues
toptags_RNA_DEG=toptags_RNA[toptags_RNA$Qvalue<=0.01,]
nrow(toptags_RNA_DEG)
rna_newcounts$samples$group=as.factor(design[,2])
plotSmear(rna_newcounts,smooth.scatter=T)


rna_newcounts_cpm=cbind(as.data.frame(rna_newcounts$genes),cpm(rna_newcounts,log=T)))
rna_cpm
#MDS plot su tutte conte DGE dopo normalizzazione
mds=plotMDS(rna_newcounts,dim.plot=c(1,2),gene.selection="common")
mds2=data.frame(x=mds$x,y=mds$y)
samples=rownames(mds2)
infos = targets[row.names(rna_newcounts$samples),]
mds2=cbind(mds2,as.data.frame(samples),as.data.frame(infos$disease),as.data.frame(infos$passage))
colnames(mds2)=c("x","y","samples","condition","passage")
pdf("PCA_AllSamples_DGE_RNA_Common_NewGTF.pdf")
ggplot(mds2, aes(x, y)) + geom_point(aes(colour = condition,shape=as.character(passage)))+ geom_text(aes(label=samples),hjust=0, vjust=0,size=2) + coord_cartesian(xlim = c(-5, 5))
dev.off()





#************************************************************ATAC***************************************************************************
#Dopo aver calcolato i picchi, li unisco tutti insieme
##in /lustre2/scratch/fsantaniello/DiMicco/ATAC_seq/Peaks

cat ATAC_UPN3_*_Macs_peaks.narrowPeak | cut -f 1,2,3 | bedSort stdin stdout | mergeBed -i - >  All_ATAC_Peaks.bed

#Calcolo il coverage per i bam files
 multiBamCov -bams /lustre1/workspace/DiMicco/prj_ALL_pilot/491_ATACseq/BAM/ATAC_UPN3_1.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/491_ATACseq/BAM/ATAC_UPN3_2.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/491_ATACseq/BAM/ATAC_UPN3_3.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/491_ATACseq/BAM/ATAC_UPN3_4.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/491_ATACseq/BAM/ATAC_UPN3_5.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/491_ATACseq/BAM/ATAC_UPN3_6.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/491_ATACseq/BAM/ATAC_UPN3_7.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/491_ATACseq/BAM/ATAC_UPN3_8.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/491_ATACseq/BAM/ATAC_UPN3_9.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/491_ATACseq/BAM/ATAC_UPN3_10.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/491_ATACseq/BAM/ATAC_UPN3_11.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/491_ATACseq/BAM/ATAC_UPN3_13.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/491_ATACseq/BAM/ATAC_UPN3_14.bam -bed ../Peaks/All_AtacPeaks_NewRuns.bed -q 1 | intersectBed -v -a stdin -b /lustre2/scratch/dcittaro/Ref/hg38.mask.bed | grep -v chrY  > All_Atac_Peaks_Counts_MultiCov_NoReps.txt



atac_all=readPeakFile("All_Atac_Peaks_Counts_MultiCov_NoReps.txt",head=F)
atac_all_anno <- annotatePeak(peak = atac_all, tssRegion=c(-2000, 2000),TxDb=gencode_rnaseqtxdb) #,annoDb="org.Hs.eg.db")
atac_all_anno_df=as.data.frame(atac_all_anno)
colnames(atac_all_anno_df)=c(colnames(atac_all_anno_df)[1:5],"UPN3_1","UPN3_2","UPN3_3","UPN3_4","UPN3_5","UPN3_6","UPN3_7","UPN3_8","UPN3_9","UPN3_10","UPN3_11","UPN3_13","UPN3_14",colnames(atac_all_anno_df)[19:27])

#ANALISI STANDARD EDGER normalizzazione RLE e con le lib.sizes sotto i picchi
atac_all_anno_df$annotation = sub(" .*", "",atac_all_anno_df$annotation)
atac_counts_new=atac_all_anno_df[,c(1,2,3,4,25,19,27,6:18)]
atac_y_new=DGEList(counts=atac_counts_new[,c(8:11,15:18)],genes=atac_counts_new[,c(1:7)],group=c("D_0","D_1","D_1","D_1","R_0","R_1","R_1","R_1"))
#atac_y_new = calcNormFactors(atac_y_new,method='RLE')#,Acutoff=-18.5)
#atac_new_keep = rowSums(cpm(atac_y_new) > 5) >= 3
#atac_y_new = atac_y_new[atac_new_keep, ]
keep <- rowSums(cpm(atac_y_new)>5) >= 3 #PRENDO SOLO I GENI CON CPM>1 IN ALMENO 8 CAMPIONI
atac_y_new<- atac_y_new[keep, , keep.lib.sizes=FALSE]
atac_y_new = calcNormFactors(atac_y_new,method='RLE')#,Acutoff=-18.5)
targets=read.table('Targets_for_DEAnalysis.txt')
targets
df = targets[row.names(atac_y_new$samples),]
df
D = as.factor(df$disease)
P = as.factor(df$passage)
design = model.matrix(~D+P)
design
atac_y_new = estimateDisp(atac_y_new, design)
fit = glmFit(atac_y_new, design) 
lrt <- glmLRT(fit, coef=2)
toptags_ATAC=as.data.frame(topTags(lrt,sort.by="PValue",n="all"))
q = qvalue(toptags_ATAC$PValue)
toptags_ATAC$Qvalue = q$qvalues
atac_y_new$samples$group=as.factor(design[,2])
plotSmear(lrt,smooth.scatter=T)
lrt_info=cbind(lrt$genes,lrt$table)
q_lrt=qvalue(lrt_info$PValue)
lrt_info$Qvalue=q_lrt$qvalues

cloud1=lrt_info[lrt_info$Qvalue<=0.05&lrt_info$logFC< -2& lrt_info$cpm>=2& lrt_info$logCPM<6 ,]
cloud2=lrt_info[lrt_info$Qvalue<=0.05&lrt_info$logFC>2 ,]
cloud3=lrt_info[lrt_info$Qvalue<=0.05&lrt_info$logFC>1 & lrt_info$logCPM>5,]
cloud4=lrt_info[lrt_info$Qvalue<=0.05&lrt_info$logCPM>8,]
cloud_atac=unique(rbind(cloud1,cloud2,cloud3,cloud4))

toptags_ATAC_DEGs=toptags_ATAC[toptags_ATAC$Qvalue<=0.01,]
ATAC_up=toptags_ATAC[toptags_ATAC$Qvalue<=0.01&toptags_ATAC$logFC>0,]
ATAC_down=toptags_ATAC[toptags_ATAC$Qvalue<=0.01&toptags_ATAC$logFC<0,]

nrow(ATAC_up)
[1] 3481

nrow(ATAC_down)
[1] 2836


#PIEPLOT ANNOTAZIONI

picchiDown=round(prop.table(table(ATAC_down$annotation))*100,1)
legend=names(picchiDown)
labels=round(prop.table(table(ATAC_down$annotation))*100,1)
picchiDown=as.vector(as.numeric(t(as.data.frame(picchiDown))[2,]))
 # Plot the chart.
# Plot the chart.
pie(picchiDown, labels = paste(picchiDown,"%",sep=''), main = "Genomic Annotation of ATAC DO peaks",col = brewer.pal(n = 7, name = "Set2"),cex=1)
legend("topleft", legend, cex = 0.8,fill = brewer.pal(n = 7, name = "Set2"))
 
 picchiUp=round(prop.table(table(ATAC_up$annotation))*100,1)
legend=names(picchiUp)
labels=round(prop.table(table(ATAC_up$annotation))*100,1)
picchiUp=as.vector(as.numeric(t(as.data.frame(picchiUp))[2,]))
 # Plot the chart.
# Plot the chart.
pie(picchiUp, labels = paste(picchiUp,"%",sep=''), main = "Genomic Annotation of ATAC DO peaks",col = brewer.pal(n = 7, name = "Set2"),cex=1)
legend("topleft", legend, cex = 0.8,fill = brewer.pal(n = 7, name = "Set2"))
 

#GO CON CLUSTERPROFILER

ATAC_up_prom=ATAC_down[ATAC_up$annotation=='Promoter',]
ATAC_down_prom=ATAC_down[ATAC_down$annotation=='Promoter',]
ATAC_DEGs_prom=toptags_ATAC_DEGs[toptags_ATAC_DEGs$annotation=='Promoter',]


ATACseq_genes_down_prom<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=ATAC_down_prom$geneId,mart=mart)
ATACseq_genes_up_prom<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=ATAC_up_prom$geneId,mart=mart)
ATACseq_deg_genes_prom<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=ATAC_DEGs_prom$geneId,mart=mart)
ATACseq_all_genes=getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=toptags_ATAC$geneId,mart=mart)
atacDown_ego_prom <- enrichGO(gene          = ATACseq_genes_down_prom$entrezgene,
universe      = as.character(ATACseq_all_genes$entrezgene),
OrgDb         = org.Hs.eg.db,
ont           = "BP",
pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = T)
dotplot(atacDown_ego_prom)

atacUp_ego_prom <- enrichGO(gene          = ATACseq_genes_up_prom$entrezgene,
universe      = as.character(ATACseq_all_genes$entrezgene),
OrgDb         = org.Hs.eg.db,
ont           = "BP",
pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = T)
dotplot(atacUp_ego_prom)

atacDeg_ego_prom <- enrichGO(gene          = ATACseq_deg_genes_prom$entrezgene,
universe      = as.character(ATACseq_all_genes$entrezgene),
OrgDb         = org.Hs.eg.db,
ont           = "BP",
pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = T)
dotplot(atacDeg_ego_prom)




#SMEARPLOT

pdf('PlotSmear_ATAC_RLE_Normalized.pdf')
plotSmear(lrt,smooth.scatter=T)
dev.off()

#MDS ORTOGONALI
pdf('OrtogonalPlot_MDS_ATAC_RLE_Normalized.pdf')
#MDS plot su tutte conte DGE dopo normalizzazione
mds=plotMDS(atac_y_new,dim.plot=c(1,2),gene.selection="common")
mds2=data.frame(x=mds$x,y=mds$y)
samples=rownames(mds2)
infos = targets[row.names(atac_y_new$samples),]
mds2=cbind(mds2,as.data.frame(samples),as.data.frame(infos$disease),as.data.frame(infos$passage))
colnames(mds2)=c("MD1","MD2","samples","condition","passage")
#pdf("PCA_AllSamples_DGE_RNA_Common.pdf")
ggplot(mds2, aes(MD1, MD2)) + geom_point(aes(colour = condition,shape=as.character(passage)))+ geom_text(aes(label=samples),hjust=0, vjust=0,size=2) + scale_x_continuous(limits = c(-4, 4)) 


#MDS plot su tutte conte DGE dopo normalizzazione
mds=plotMDS(atac_y_new,dim.plot=c(2,3),gene.selection="common")
mds2=data.frame(x=mds$x,y=mds$y)
samples=rownames(mds2)
infos = targets[row.names(atac_y_new$samples),]
mds2=cbind(mds2,as.data.frame(samples),as.data.frame(infos$disease),as.data.frame(infos$passage))
colnames(mds2)=c("MD2","MD3","samples","condition","passage")
#pdf("PCA_AllSamples_DGE_RNA_Common.pdf")
ggplot(mds2, aes(MD3, MD2)) + geom_point(aes(colour = condition,shape=as.character(passage)))+ geom_text(aes(label=samples),hjust=0, vjust=0,size=2) + scale_x_continuous(limits = c(-4, 4)) 

mds=plotMDS(atac_y_new,dim.plot=c(1,3),gene.selection="common")
mds2=data.frame(x=mds$x,y=mds$y)
samples=rownames(mds2)
infos = targets[row.names(atac_y_new$samples),]
mds2=cbind(mds2,as.data.frame(samples),as.data.frame(infos$disease),as.data.frame(infos$passage))
colnames(mds2)=c("MD1","MD3","samples","condition","passage")
#pdf("PCA_AllSamples_DGE_RNA_Common.pdf")
ggplot(mds2, aes(MD1, MD3)) + geom_point(aes(colour = condition,shape=as.character(passage)))+ geom_text(aes(label=samples),hjust=0, vjust=0,size=2)  + scale_x_continuous(limits = c(-4, 4)) 

dev.off()






toptags_ATAC_prom=toptags_ATAC[toptags_ATAC$annotation=='Promoter',]
colnames(toptags_ATAC_prom)=c(colnames(toptags_ATAC_prom)[1:4],'genes',colnames(toptags_ATAC_prom)[6:13])

rna_atac_merge=merge(toptags_RNA,toptags_ATAC_prom,by='genes')
nrow(rna_atac_merge[rna_atac_merge$Qvalue.x<=0.01&rna_atac_merge$logFC.x<0&rna_atac_merge$Qvalue.y<=0.01&rna_atac_merge$logFC.y<0,])
[1] 204
nrow(rna_atac_merge[rna_atac_merge$Qvalue.x<=0.01&rna_atac_merge$logFC.x<0&rna_atac_merge$Qvalue.y<=0.01&rna_atac_merge$logFC.y>0,])
[1] 32
nrow(rna_atac_merge[rna_atac_merge$Qvalue.x<=0.01&rna_atac_merge$logFC.x>0&rna_atac_merge$Qvalue.y<=0.01&rna_atac_merge$logFC.y>0,])
[1] 132
nrow(rna_atac_merge[rna_atac_merge$Qvalue.x<=0.01&rna_atac_merge$logFC.x>0&rna_atac_merge$Qvalue.y<=0.01&rna_atac_merge$logFC.y<0,])
[1] 62

rna_atac=matrix(c(78,196,15,519),nrow=2,ncol=2)
colnames(rna_atac)=c('UP_rna','DOWN_rna')
rownames(rna_atac)=c('DOWN_atac','UP_atac')
```
          UP_rna DOWN_rna
DOWN_atac     78       15
UP_atac      196      519
```



 
 
#************************************************************K4me3***************************************************************************

#Dopo aver calcolato i picchi, li unisco tutti insieme
##in /lustre2/scratch/fsantaniello/DiMicco/K4me3/Peaks

cat *_Macs_peaks.narrowPeak | cut -f 1,2,3 | bedSort stdin stdout | mergeBed -i - >  All_K4me3_Peaks.bed


#Calcolo il coverage per i bam files
multiBamCov -bams /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K4me3_UPN3_1.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K4me3_UPN3_2.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K4me3_UPN3_3.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K4me3_UPN3_4.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K4me3_UPN3_5.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K4me3_UPN3_6.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K4me3_UPN3_7.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K4me3_UPN3_8.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K4me3_UPN3_9.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K4me3_UPN3_10.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K4me3_UPN3_11.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K4me3_UPN3_12.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K4me3_UPN3_13.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K4me3_UPN3_14.bam -bed /lustre2/scratch/fsantaniello/DiMicco/K4me3/Peaks/All_K4me3_Peaks.bed -q 1 |  intersectBed -v -a stdin -b /lustre2/scratch/dcittaro/Ref/hg38.mask.bed | grep -v chrY > /lustre2/scratch/fsantaniello/DiMicco/K4me3/Peaks/All_K4me3_Peaks_Counts_MultiCov_New_NoReps.txt


#ANALISI EDGER CON METODO RLE e con lib.sizes sotto picchi

k4me3_all=readPeakFile("All_K4me3_Peaks_Counts_MultiCov_New_NoReps_New.txt",head=F)
k4me3_all_anno <- annotatePeak(peak = k4me3_all, tssRegion=c(-2000, 2000),TxDb=gencode_rnaseqtxdb)#,annoDb="org.Hs.eg.db")
 #plotAnnoPie(k4me3_all_anno)


k4me3_all_anno_df=as.data.frame(k4me3_all_anno)
k4me3_all_anno_df$annotation = sub(" .*", "", k4me3_all_anno_df$annotation)
colnames(k4me3_all_anno_df)=c(colnames(k4me3_all_anno_df)[1:5],"UPN3_1","UPN3_2","UPN3_3","UPN3_4","UPN3_5","UPN3_6","UPN3_7","UPN3_8","UPN3_9","UPN3_10","UPN3_11","UPN3_12","UPN3_13","UPN3_14",colnames(k4me3_all_anno_df)[20:28])

#head(k4me3_all_anno_df)
k4me3_counts_new=k4me3_all_anno_df[,c(1,2,3,26,20,28,6:19)]
k4me3_y_new=DGEList(counts=k4me3_counts_new[,c(7:10,14:17)],genes=k4me3_counts_new[,c(1:6)],group=c("D_0","D_1","D_1","D_1","R_0","R_1","R_1","R_1"))#,lib.size = c(22089014,23497699,30146073,21596634,24792817,22820342,26900163,45421712))
#k4me3_y_new = calcNormFactors(k4me3_y_new,method='RLE')#,Acutoff=-18.5)
#k4me3_new_keep = rowSums(cpm(k4me3_y_new) > 5) >= 3
#k4me3_y_new = k4me3_y_new[k4me3_new_keep, ]
keep <- rowSums(cpm(k4me3_y_new)>5) >= 3 #PRENDO SOLO I GENI CON CPM>1 IN ALMENO 8 CAMPIONI
k4me3_y_new<- k4me3_y_new[keep, , keep.lib.sizes=FALSE]
k4me3_y_new = calcNormFactors(k4me3_y_new,method='RLE')#,Acutoff=-18.5)
library(qvalue)
 df = targets[row.names(k4me3_y_new$samples),]
 D = as.factor(df$disease)
 P = as.factor(df$passage)
 design = model.matrix(~D+P)
 design
 k4me3_y_new= estimateDisp(k4me3_y_new, design)
 fit = glmFit(k4me3_y_new, design)
 lrt <- glmLRT(fit, coef=2)
 toptags_k4me3=as.data.frame(topTags(lrt,sort.by="none",n="all"))
 q = qvalue(toptags_k4me3$PValue)
  toptags_k4me3$Qvalue = q$qvalues
 toptags_k4me3_DEG=toptags_k4me3[toptags_k4me3$Qvalue<0.05,]
k4me3_y_new$samples$group=as.factor(design[,2])

lrt_info=cbind(lrt$genes,lrt$table)
q_lrt=qvalue(lrt_info$PValue)
lrt_info$Qvalue=q_lrt$qvalues
cloud1=lrt_info[lrt_info$Qvalue<=0.05&lrt_info$logFC< -4 ,]
cloud2=lrt_info[lrt_info$Qvalue<=0.05&lrt_info$logFC< -2& lrt_info$logCPM>4 ,]
cloud3=lrt_info[lrt_info$Qvalue<=0.05&lrt_info$logFC< -1.4& lrt_info$logCPM>6 ,]
cloud4=lrt_info[lrt_info$Qvalue<=0.05&lrt_info$logFC>2& lrt_info$logCPM>3 ,]
cloud5=lrt_info[lrt_info$Qvalue<=0.05&lrt_info$logFC02& lrt_info$logCPM>8 ,]
cloud_k4me3=unique(rbind(cloud1,cloud2,cloud3,cloud4,cloud5))


k4me3_y_new_cpm_all=cbind(as.data.frame(k4me3_y_new$genes),cpm(k4me3_y_new,log=T))

toptags_k4me3_down=toptags_k4me3[toptags_k4me3$Qvalue<=0.05&toptags_k4me3$logFC<0,]
toptags_k4me3_up=toptags_k4me3[toptags_k4me3$Qvalue<=0.05&toptags_k4me3$logFC>0,]

nrow(toptags_k4me3_down)
[1] 1927
nrow(toptags_k4me3_up)
[1] 1573
pdf('PlotSmear_K4me3_RLE_Normalized.pdf')
plotSmear(lrt,smooth.scatter=T)
dev.off()
pdf('OrtogonalPlot_MDS_K4me3_RLE_Normalized.pdf')
#MDS plot su tutte conte DGE dopo normalizzazione
mds=plotMDS(k4me3_y_new,dim.plot=c(1,2))
mds2=data.frame(x=mds$x,y=mds$y)
samples=rownames(mds2)
infos = targets[row.names(k4me3_y_new$samples),]
mds2=cbind(mds2,as.data.frame(samples),as.data.frame(infos$disease),as.data.frame(infos$passage))
colnames(mds2)=c("MD1","MD2","samples","condition","passage")
#pdf("PCA_AllSamples_DGE_RNA_Common.pdf")
ggplot(mds2, aes(MD1, MD2)) + geom_point(aes(colour = condition,shape=as.character(passage)))+ geom_text(aes(label=samples),hjust=0, vjust=0,size=2) + scale_x_continuous(limits = c(-3.5, 3.5)) 


#MDS plot su tutte conte DGE dopo normalizzazione
mds=plotMDS(k4me3_y_new,dim.plot=c(2,3))
mds2=data.frame(x=mds$x,y=mds$y)
samples=rownames(mds2)
infos = targets[row.names(k4me3_y_new$samples),]
mds2=cbind(mds2,as.data.frame(samples),as.data.frame(infos$disease),as.data.frame(infos$passage))
colnames(mds2)=c("MD2","MD3","samples","condition","passage")
#pdf("PCA_AllSamples_DGE_RNA_Common.pdf")
ggplot(mds2, aes(MD3, MD2)) + geom_point(aes(colour = condition,shape=as.character(passage)))+ geom_text(aes(label=samples),hjust=0, vjust=0,size=2) + scale_x_continuous(limits = c(-3.5, 3.5)) 

mds=plotMDS(k4me3_y_new,dim.plot=c(1,3))
mds2=data.frame(x=mds$x,y=mds$y)
samples=rownames(mds2)
infos = targets[row.names(k4me3_y_new$samples),]
mds2=cbind(mds2,as.data.frame(samples),as.data.frame(infos$disease),as.data.frame(infos$passage))
colnames(mds2)=c("MD1","MD3","samples","condition","passage")
#pdf("PCA_AllSamples_DGE_RNA_Common.pdf")
ggplot(mds2, aes(MD1, MD3)) + geom_point(aes(colour = condition,shape=as.character(passage)))+ geom_text(aes(label=samples),hjust=0, vjust=0,size=2)  + scale_x_continuous(limits = c(-3.5, 3.5)) 
dev.off()


#************************************************************K27ac***************************************************************************

#Dopo aver calcolato i picchi, li unisco tutti insieme
##in /lustre2/scratch/fsantaniello/DiMicco/K27ac/Peaks

cat *_Macs_peaks.narrowPeak | cut -f 1,2,3 | bedSort stdin stdout | mergeBed -i - >  All_K4me3_Peaks.bed

#PBS -l select=1:app=java:ncpus=6:mem=60gb
#PBS -P 161
#PBS -q workq
#PBS -m ae

multiBamCov -bams /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K27ac_UPN3_1.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K27ac_UPN3_2.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K27ac_UPN3_3.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K27ac_UPN3_4.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K27ac_UPN3_5.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K27ac_UPN3_6.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K27ac_UPN3_7.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K27ac_UPN3_8.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K27ac_UPN3_9.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K27ac_UPN3_10.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K27ac_UPN3_11.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K27ac_UPN3_12.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K27ac_UPN3_13.bam /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/H3K27ac_UPN3_14.bam -bed /lustre2/scratch/fsantaniello/DiMicco/K27Ac/Peaks/All_K27ac_Peaks.bed  -q 1 | bedtools intersect -v -a stdin -b /lustre1/genomes/hg38/annotation/hg38.blacklist.bed > /lustre2/scratch/fsantaniello/DiMicco/K27Ac/Peaks/All_K27ac_Peaks_Counts_MultiCov.txt



targets=read.table('Targets_for_DEAnalysis.txt')
k27ac_all=readPeakFile("All_K27ac_Peaks_Counts_MultiCov_NoReps.txt",head=F)
k27ac_all_anno <- annotatePeak(peak = k27ac_all, tssRegion=c(-2000, 2000),TxDb=gencode_rnaseqtxdb) #,annoDb="org.Hs.eg.db")
k27ac_all_anno_df=as.data.frame(k27ac_all_anno)
k27ac_all_anno_df$annotation = sub(" .*", "", k27ac_all_anno_df$annotation)
colnames(k27ac_all_anno_df)=c(colnames(k27ac_all_anno_df)[1:5],"UPN3_1","UPN3_2","UPN3_3","UPN3_4","UPN3_5","UPN3_6","UPN3_7","UPN3_8","UPN3_9","UPN3_10","UPN3_11","UPN3_12","UPN3_13","UPN3_14",colnames(k27ac_all_anno_df)[20:28])
k27ac_counts_new=k27ac_all_anno_df[,c(1,2,3,4,26,20,28,6:19)]
k27ac_y_new=DGEList(counts=k27ac_counts_new[,c(8:11,15:18)],genes=k27ac_counts_new[,c(1:7)],group=c("D_0","D_1","D_1","D_1","R_0","R_1","R_1","R_1"))
#k27ac_y_new = calcNormFactors(k27ac_y_new,method='RLE')#,Acutoff=-18.5)
#k27ac_new_keep = rowSums(cpm(k27ac_y_new) > 5) >= 3
#k27ac_y_new = k27ac_y_new[k27ac_new_keep, ]
keep <- rowSums(cpm(k27ac_y_new)>5) >= 3 #PRENDO SOLO I GENI CON CPM>1 IN ALMENO 8 CAMPIONI
k27ac_y_new<- k27ac_y_new[keep, , keep.lib.sizes=FALSE]
k27ac_y_new = calcNormFactors(k27ac_y_new,method='RLE')#,Acutoff=-18.5)
library(qvalue)
df = targets[row.names(k27ac_y_new$samples),]
D = as.factor(df$disease)
P = as.factor(df$passage)
design = model.matrix(~D+P)
design
k27ac_y_new= estimateDisp(k27ac_y_new, design)
fit = glmFit(k27ac_y_new, design)
lrt <- glmLRT(fit, coef=2)
toptags_k27ac=as.data.frame(topTags(lrt,sort.by="PValue",n="all"))
q = qvalue(toptags_k27ac$PValue)
toptags_k27ac$Qvalue = q$qvalues
 k27ac_y_new$samples$group=as.factor(design[,2])
 plotSmear(lrt,smooth.scatter=T)


cloud1=lrt_info[lrt_info$Qvalue<=0.05&lrt_info$logFC<= -3,]
cloud2=lrt_info[lrt_info$Qvalue<=0.05&lrt_info$logFC< 0 & lrt_info$logCPM >8 ,]
cloud3=lrt_info[lrt_info$Qvalue<=0.05&lrt_info$logFC< -1.5 & lrt_info$logCPM >6 & lrt_info$logCPM<8 ,]
cloud4=lrt_info[lrt_info$Qvalue<=0.05&lrt_info$logFC > 2.5,]
cloud5=lrt_info[lrt_info$Qvalue<=0.05&lrt_info$logFC > 1.7 & lrt_info$logCPM>6.3&lrt_info$logCPM<8,]
cloud6=lrt_info[lrt_info$Qvalue<=0.05&lrt_info$logFC >= 0 & lrt_info$logCPM>8,]
cloud_k27ac=unique(rbind(cloud1,cloud2,cloud3,cloud4,cloud5,cloud6))


toptags_k27ac_DEG=toptags_k27ac[toptags_k27ac$Qvalue<=0.05,]
nrow(toptags_k27ac_DEG[toptags_k27ac_DEG$logFC<0,])
[1] 12840
nrow(toptags_k27ac_DEG[toptags_k27ac_DEG$logFC>0,])
[1] 11322



pdf('PlotSmear_K27_RLE_Normalized.pdf')
plotSmear(lrt,smooth.scatter=T)
dev.off()
pdf('OrtogonalPlot_MDS_K27_RLE_Normalized.pdf')
#MDS plot su tutte conte DGE dopo normalizzazione
mds=plotMDS(k27ac_y_new,dim.plot=c(1,2))
mds2=data.frame(x=mds$x,y=mds$y)
samples=rownames(mds2)
infos = targets[row.names(k27ac_y_new$samples),]
mds2=cbind(mds2,as.data.frame(samples),as.data.frame(infos$disease),as.data.frame(infos$passage))
colnames(mds2)=c("MD1","MD2","samples","condition","passage")
#pdf("PCA_AllSamples_DGE_RNA_Common.pdf")
ggplot(mds2, aes(MD1, MD2)) + geom_point(aes(colour = condition,shape=as.character(passage)))+ geom_text(aes(label=samples),hjust=0, vjust=0,size=2) + scale_x_continuous(limits = c(-4, 4)) 


#MDS plot su tutte conte DGE dopo normalizzazione
mds=plotMDS(k27ac_y_new,dim.plot=c(2,3))
mds2=data.frame(x=mds$x,y=mds$y)
samples=rownames(mds2)
infos = targets[row.names(k27ac_y_new$samples),]
mds2=cbind(mds2,as.data.frame(samples),as.data.frame(infos$disease),as.data.frame(infos$passage))
colnames(mds2)=c("MD2","MD3","samples","condition","passage")
#pdf("PCA_AllSamples_DGE_RNA_Common.pdf")
ggplot(mds2, aes(MD3, MD2)) + geom_point(aes(colour = condition,shape=as.character(passage)))+ geom_text(aes(label=samples),hjust=0, vjust=0,size=2) + scale_x_continuous(limits = c(-4, 4)) 

mds=plotMDS(k27ac_y_new,dim.plot=c(1,3))
mds2=data.frame(x=mds$x,y=mds$y)
samples=rownames(mds2)
infos = targets[row.names(k27ac_y_new$samples),]
mds2=cbind(mds2,as.data.frame(samples),as.data.frame(infos$disease),as.data.frame(infos$passage))
colnames(mds2)=c("MD1","MD3","samples","condition","passage")
#pdf("PCA_AllSamples_DGE_RNA_Common.pdf")
ggplot(mds2, aes(MD1, MD3)) + geom_point(aes(colour = condition,shape=as.character(passage)))+ geom_text(aes(label=samples),hjust=0, vjust=0,size=2)  + scale_x_continuous(limits = c(-4, 4)) 

dev.off()





#*******************************************************Overlap Picchi con GenomicRanges************************************************************************


#toptags_ATAC_DEG_gr=makeGRangesFromDataFrame(toptags_ATAC_DEG,keep.extra.columns = T)
toptags_ATAC_gr=makeGRangesFromDataFrame(toptags_ATAC,keep.extra.columns = T)
toptags_k27ac_gr=makeGRangesFromDataFrame(toptags_k27ac,keep.extra.columns = T)
toptags_k4me3_gr=makeGRangesFromDataFrame(toptags_k4me3,keep.extra.columns = T)

#ATAC + K4me3
df_atac_k4me3=mergeByOverlaps(toptags_ATAC_gr,toptags_k4me3_gr)
Merged_Atac_K4me3=cbind(as.data.frame(attr(df_atac_k4me3[1:10],"listData")),as.data.frame(attr(df_atac_k4me3[11:20],"listData")))
Merged_Atac_K4me3=Merged_Atac_K4me3[,c(1:15,24:37)]
colnames(Merged_Atac_K4me3)=c('seqnames','start','end',colnames(Merged_Atac_K4me3)[4:29])

#ATAC + K27ac
df_atac_k27ac=mergeByOverlaps(toptags_ATAC_gr,toptags_k27ac_gr)
Merged_Atac_K27ac=cbind(as.data.frame(attr(df_atac_k27ac[1:10],"listData")),as.data.frame(attr(df_atac_k27ac[11:20],"listData")))
Merged_Atac_K27ac=Merged_Atac_K27ac[,c(1:15,24:37)]
colnames(Merged_Atac_K27ac)=c('seqnames','start','end',colnames(Merged_Atac_K27ac)[4:29])


#Creo il tabellone
colnames(toptags_ATAC)=c(colnames(toptags_ATAC)[1:4],'genes',colnames(toptags_ATAC)[6:13])
ATAC_RNA=merge(toptags_ATAC,toptags_RNA,by='genes')
#ATAC_K4me3=merge(toptags_ATAC,Merged_Atac_K4me3,by=c('seqnames','start','end'),all.x=T)
#ATAC_K4me3=ATAC_K4me3[,c(1:13,26:39)]
#ATAC_K27Ac=merge(toptags_ATAC,Merged_Atac_K27ac,by=c('seqnames','start','end'),all.x=T)
#ATAC_K27Ac=ATAC_K27Ac[,c(1:13,26:39)]
#ATAC_chip=merge(ATAC_K4me3,ATAC_K27Ac,by=c('seqnames','start','end'),all.x=T,all.y=T)

Atac_K4me3_K27ac=merge(Merged_Atac_K4me3,Merged_Atac_K27ac,by=c('seqnames','start','end'),all.x=T,all.y=T)

Def=merge(ATAC_RNA,Atac_K4me3_K27ac,all.x=T,by=c('seqnames','start','end'),all.y=T)

Def$K4me3Overlap=ifelse(is.na(Def$toptags_k4me3_gr.seqnames),'-','+')
Def$K27acOverlap=ifelse(is.na(Def$toptags_k27ac_gr.seqnames),'-','+')
colnames=c("ATAC_chr","ATAC_start","ATAC_end","ATAC_logFC","ATAC_PValue","ATAC_Qvalue","genes","RNA_logFC","RNA_PValue","RNA_Qvalue","ATAC_annotation","ATAC_distanceToTSS","K4me3Overlap","K27acOverlap","K4me3_chr","K4me3_start","K4me3_end","K4me3_annotation","K4me3_distanceToTSS","K4me3_logFC","K4me3_PValue","K4me3_Qvalue"	,"K27ac_chr","K27ac_start","K27ac_end","K27ac_annotation","K27ac_distanceToTSS","K27ac_logFC","K27ac_Pvalue","K27ac_Qvalue")
Def_towrite=Def[,c(1,2,3,8,11,13,4,14,17,19,6,7,72,73,32,33,34,38,39,40,43,45,58,59,60,64,65,66,69,71)]

write.table(Def_towrite,'Integrated_Annotation_ATAC_RNA_K4me3_K27Ac.txt',sep='\t',col.names=T,row.names=F,quote=F)




enhancers=Def_towrite[Def_towrite$K4me3Overlap=='-'&Def_towrite$K27acOverlap=='+',1:3]
enhancers_deg=Def_towrite[Def_towrite$K4me3Overlap=='-'&Def_towrite$K27acOverlap=='+'&Def_towrite$ATAC_Qvalue<=0.25,1:3]
promoters=Def_towrite[Def_towrite$ATAC_annotation=='Promoter',1:3]



df_atac_degs_k4me3=cbind(as.data.frame(attr(df_atac_degs_k4me3[1:10],"listData")),as.data.frame(attr(df_atac_degs_k4me3[11:20],"listData")))

Merged_Atac_K4me3_gr=makeGRangesFromDataFrame(Merged_Atac_K4me3,keep.extra.columns = T)

 cor.test(Merged_Atac_K4me3$toptags_ATAC_gr.logFC,Merged_Atac_K4me3$toptags_k4me3_gr.logFC)

	Pearsons product-moment correlation

data:  Merged_Atac_K4me3$toptags_ATAC_gr.logFC and Merged_Atac_K4me3$toptags_k4me3_gr.logFC
t = -0.23142, df = 12728, p-value = 0.817
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.01942219  0.01532099
sample estimates:
         cor 
-0.002051222



cor.test(Merged_Atac_K27ac$toptags_ATAC_gr.logFC,Merged_Atac_K27ac$toptags_k27ac_gr.logFC)

	Pearsons product-moment correlation

data:  Merged_Atac_K27ac$toptags_ATAC_gr.logFC and Merged_Atac_K27ac$toptags_k27ac_gr.logFC
t = 26.079, df = 26407, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.1466737 0.1701897
sample estimates:
      cor 
0.1584541



df_atac_k4me3_k27ac=mergeByOverlaps(Merged_Atac_K4me3_gr,Merged_Atac_K27ac_gr)



#************************************************************Ho trovato i siti di binding di Erg***************************************************************************

Def_gr==makeGRangesFromDataFrame(Def,keep.extra.columns = T)
Erg_gr=makeGRangesFromDataFrame(erg,keep.extra.columns = T)

#ATAC + K4me3
df_Def_erg=mergeByOverlaps(Def_gr,Erg_gr)
Erg_data_frame=as.data.frame(attr(df_atac_erg[,10:11],"listData"))
Erg_data_frame=Erg_data_frame[,2:4]

Merged_Atac_Erg=cbind(as.data.frame(attr(df_atac_erg[1:10],"listData")),Erg_data_frame)
Merged_Atac_K4me3=Merged_Atac_K4me3[,c(1:15,24:37)]
colnames(Merged_Atac_K4me3)=c('seqnames','start','end',colnames(Merged_Atac_K4me3)[4:29])





#************************************************************Altre metanalisi da sistemare***************************************************************************





#Faccio gene ontology sui picchi significativi ATAC che overlappano con K4me3

Merged_Atac_K4me3_sigatac=Merged_Atac_K4me3[Merged_Atac_K4me3$toptags_ATAC_gr.Qvalue<=0.05,]

Atac_k4me3_sig<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=Merged_Atac_K4me3_sigatac$toptags_ATAC_gr.SYMBOL,mart=mart)
Atac_k4me3_allgenes<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=Merged_Atac_K4me3$toptags_ATAC_gr.SYMBOL,mart=mart)

atack4me3_ego <- enrichGO(gene = Atac_k4me3_sig$entrezgene,universe = as.character(Atac_k4me3_allgenes$entrezgene),OrgDb= org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = T)
dotplot(atack4me3_ego)


#Faccio gene ontology sui picchi significativi ATAC che overlappano con K4me3 e che hanno logFC>0 per ATAc e logFC<0 per K4me3
Merged_Atac_K4me3_ControCorrelati=Merged_Atac_K4me3[Merged_Atac_K4me3$toptags_ATAC_gr.logFC > 0.5 &Merged_Atac_K4me3$toptags_k4me3_gr.logFC< -2,]

Atac_k4me3_controCorr<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=Merged_Atac_K4me3_ControCorrelati$toptags_ATAC_gr.SYMBOL,mart=mart)
Atac_k4me3_allgenes<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=Merged_Atac_K4me3$toptags_ATAC_gr.SYMBOL,mart=mart)

atack4me3_contrcor_ego <- enrichGO(gene = Atac_k4me3_controCorr$entrezgene,universe = as.character(Atac_k4me3_allgenes$entrezgene),OrgDb= org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = T)
dotplot(atack4me3_contrcor_ego)

#************************************************************RNA-seq + ATAC***************************************************************************

#Top100 Up e Down Genes
toptags_RNA_up_top100=toptags_RNA[toptags_RNA$Qvalue<=0.01&toptags_RNA$logFC>0,]
toptags_RNA_up_top100=toptags_RNA_up_top100[order(-toptags_RNA_up_top100$logFC),]
toptags_RNA_up_top100=toptags_RNA_up_top100[order(-toptags_RNA_up_top100$logFC),]
toptags_RNA_up_top100=toptags_RNA_up_top100[1:100,]

toptags_RNA_down_top100=toptags_RNA[toptags_RNA$Qvalue<=0.01&toptags_RNA$logFC<0,]
toptags_RNA_down_top100=toptags_RNA_down_top100[order(toptags_RNA_down_top100$logFC),]
toptags_RNA_down_top100=toptags_RNA_down_top100[1:100,]
write.table(toptags_RNA_down_top100,"RNAseq_DiMicco_NewAnalysis_Top100DownGenes.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(toptags_RNA_up_top100,"RNAseq_DiMicco_NewAnalysis_Top100UpGenes.txt",sep="\t",col.names=T,row.names=F,quote=F)


#************************************************************RNA-seq + ATAC Promoter***********************************************************************

toptags_ATAC_prom=toptags_ATAC[toptags_ATAC$annotation=="Promoter",]
colnames(toptags_ATAC_prom)=c(colnames(toptags_ATAC_prom)[1:3],"genes",colnames(toptags_ATAC_prom)[5:12])
toptags_RNA_ATAC_prom=merge(toptags_RNA,toptags_ATAC_prom,by="genes")
toptags_RNA_ATAC_prom_top100down=toptags_RNA_ATAC_prom[toptags_RNA_ATAC_prom$Qvalue.x<=0.01&toptags_RNA_ATAC_prom$logFC.x<0,]
toptags_RNA_ATAC_prom_top100down=toptags_RNA_ATAC_prom_top100down[order(toptags_RNA_ATAC_prom_top100down$logFC.x),]
toptags_RNA_ATAC_prom_top100down=toptags_RNA_ATAC_prom_top100down[1:100,]

toptags_RNA_ATAC_prom_top100up=toptags_RNA_ATAC_prom[toptags_RNA_ATAC_prom$Qvalue.x<=0.01&toptags_RNA_ATAC_prom$logFC.x>0,]
toptags_RNA_ATAC_prom_top100up=toptags_RNA_ATAC_prom_top100up[order(-toptags_RNA_ATAC_prom_top100up$logFC.x),]
toptags_RNA_ATAC_prom_top100up=toptags_RNA_ATAC_prom_top100up[1:100,]

write.table(toptags_RNA_ATAC_prom_top100up,"RNAseq_DiMicco_NewAnalysis_Top100UpGenes_AtacPromoter.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(toptags_RNA_ATAC_prom_top100down,"RNAseq_DiMicco_NewAnalysis_Top100DownGenes_AtacPromoter.txt",sep="\t",col.names=T,row.names=F,quote=F)


#************************************************************RNA-seq + ATAC Enhancers***********************************************************************

#Prendo solo i picchi di K27 che NON overlappano con K4
toptags_k27ac_gr_NoK4me3=toptags_k27ac_gr[toptags_k27ac_gr %outside% toptags_k4me3_gr,] 
 
 
df_atac_enhancers=mergeByOverlaps(toptags_ATAC_gr,toptags_k27ac_gr_NoK4me3)
Merged_Atac_Enhancers=cbind(as.data.frame(attr(df_atac_enhancers[1:10],"listData")),as.data.frame(attr(df_atac_enhancers[11:20],"listData")))

plot(Merged_Atac_Enhancers$toptags_ATAC_gr.logFC,Merged_Atac_Enhancers$toptags_k27ac_gr_NoK4me3.logFC)

Atac_enhancers_down<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=Merged_Atac_Enhancers_down$toptags_ATAC_gr.genes,mart=mart)
Atac_enhancers_up<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=Merged_Atac_Enhancers_up$toptags_ATAC_gr.genes,mart=mart)

AtacDown_enhancers_ego <- enrichGO(gene          = Atac_enhancers_down$entrezgene,
universe      = as.character(RNAseq_all_genes$entrezgene),
OrgDb         = org.Hs.eg.db,
ont           = "BP",
pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = T)
dotplot(rnaDown_ego)


 colnames(Merged_Atac_Enhancers)=c("chr","start","end","width","strand","genes","anotation","distanceToTss","logFC","logCPM","LR","PValue","FDR","Qvalue")
 
RNAseq_merged_Atac_Enhancers=merge(toptags_RNA,Merged_Atac_Enhancers,by="genes")
RNAseq_merged_Atac_Enhancers=RNAseq_merged_Atac_Enhancers[RNAseq_merged_Atac_Enhancers$Qvalue.x<=0.01,]
RNAseq_merged_Atac_Enhancers_Top100Down=RNAseq_merged_Atac_Enhancers[RNAseq_merged_Atac_Enhancers$logFC.x<0,]
RNAseq_merged_Atac_Enhancers_Top100Down=RNAseq_merged_Atac_Enhancers_Top100Down[order(RNAseq_merged_Atac_Enhancers_Top100Down$logFC.x),]
 RNAseq_merged_Atac_Enhancers_Top100Down_genesOnly=RNAseq_merged_Atac_Enhancers_Top100Down[,c(1:7,13)]
 RNAseq_merged_Atac_Enhancers_Top100Down_genesOnly=unique(RNAseq_merged_Atac_Enhancers_Top100Down_genesOnly)
 RNAseq_merged_Atac_Enhancers_Top100Down_genesOnly=RNAseq_merged_Atac_Enhancers_Top100Down_genesOnly[1:100,]
 
 
 RNAseq_merged_Atac_Enhancers=merge(toptags_RNA,Merged_Atac_Enhancers,by="genes")
RNAseq_merged_Atac_Enhancers=RNAseq_merged_Atac_Enhancers[RNAseq_merged_Atac_Enhancers$Qvalue.x<=0.01,]
RNAseq_merged_Atac_Enhancers_Top100Up=RNAseq_merged_Atac_Enhancers[RNAseq_merged_Atac_Enhancers$logFC.x>0,]
RNAseq_merged_Atac_Enhancers_Top100Up=RNAseq_merged_Atac_Enhancers_Top100Up[order(-RNAseq_merged_Atac_Enhancers_Top100Up$logFC.x),]
 RNAseq_merged_Atac_Enhancers_Top100Up_genesOnly=RNAseq_merged_Atac_Enhancers_Top100Up[,1:7]
 RNAseq_merged_Atac_Enhancers_Top100Up_genesOnly=unique(RNAseq_merged_Atac_Enhancers_Top100Up_genesOnly)
 RNAseq_merged_Atac_Enhancers_Top100Up_genesOnly=RNAseq_merged_Atac_Enhancers_Top100Up_genesOnly[1:100,]
 

#******************************************************Gene Ontologies delle diverse tecniche***************************************************


Merged_Atac_Enhancers=as.data.frame(attr(df_atac_enhancers[1:10],"listData"))




#******************************************************Gene Ontologies delle diverse tecniche***************************************************
library(DOSE)
library(GO.db)
library(org.Hs.eg.db)
library(topGO)
library(GSEABase)
library(clusterProfiler)
library(biomaRt)
mart=useMart("ENSEMBL_MART_ENSEMBL",host="uswest.ensembl.org",dataset="hsapiens_gene_ensembl")

#RNASEQ
RNAseq_genes_down<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=RNA_down$genes,mart=mart)
RNAseq_all_genes<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=toptags_RNA$genes,mart=mart)
RNAseq_deg_genes<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=toptags_RNA_DEG$genes,mart=mart)
RNAseq_genes_up<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=RNA_up$genes,mart=mart)


#RNASEQ
ATACseq_genes_down<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=ATAC_down$genes,mart=mart)
ATACseq_all_genes<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=toptags_ATAC$genes,mart=mart)
ATACseq_deg_genes<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=toptags_ATAC_DEG$geneId,mart=mart)
ATACseq_genes_up<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=ATAC_up$genes,mart=mart)



rnaDown_ego <- enrichGO(gene          = RNAseq_deg_genes$entrezgene,
universe      = as.character(RNAseq_all_genes$entrezgene),
OrgDb         = org.Hs.eg.db,
ont           = "BP",
pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = T)
dotplot(rnaDown_ego)


rnaUp_ego <- enrichGO(gene          = RNAseq_genes_up$entrezgene,
                         universe      = as.character(RNAseq_all_genes$entrezgene),
                         OrgDb         = org.Hs.eg.db,
                         ont           = "BP",
                         pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = T)

dotplot(rnaUp_ego)



#ATAC-seq
ATAC_up=toptags_ATAC[toptags_ATAC$Qvalue<=0.01&toptags_ATAC$logFC>0,]
ATAC_down=toptags_ATAC[toptags_ATAC$Qvalue<=0.01&toptags_ATAC$logFC<0,]
ATAC_degs=toptags_ATAC_DEG
ATACseq_genes_down<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=ATAC_down$genes,mart=mart)
ATACseq_all_genes<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=toptags_ATAC$genes,mart=mart)
ATACseq_genes_up<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=ATAC_up$SYMBOL,mart=mart)
ATAC_degs<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=toptags_ATAC_DEG$genes,mart=mart)

atacDegs_ego <- enrichGO(gene = ATACseq_genes_down$entrezgene,
universe      = as.character(ATACseq_all_genes$entrezgene),
OrgDb         = org.Hs.eg.db,
ont           = "BP",
pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = T)
dotplot(atacDegs_ego)

atacUp_ego <- enrichGO(gene          = ATACseq_genes_up$entrezgene,
universe      = as.character(ATACseq_all_genes$entrezgene),
OrgDb         = org.Hs.eg.db,
ont           = "BP",
pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = T)
dotplot(atacUp_ego)


atacDown_ego <- enrichGO(gene          = ATACseq_genes_up$entrezgene,
universe      = as.character(ATACseq_all_genes$entrezgene),
OrgDb         = org.Hs.eg.db,
ont           = "BP",
pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = T)
dotplot(atacDown_ego)

ATAC_up_prom=ATAC_up[ATAC_up$annotation=="Promoter",]
ATAC_down_prom=ATAC_up[ATAC_down$annotation=="Promoter",]
ATACseq_genes_down_prom<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=ATAC_down_prom$SYMBOL,mart=mart)
ATACseq_genes_up_prom<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=ATAC_up_prom$SYMBOL,mart=mart)

atacDown_ego_prom <- enrichGO(gene          = ATACseq_genes_down_prom$entrezgene,
universe      = as.character(ATACseq_all_genes$entrezgene),
OrgDb         = org.Hs.eg.db,
ont           = "BP",
pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = T)
dotplot(atacDown_ego_prom)

atacUp_ego_prom <- enrichGO(gene          = ATACseq_genes_up_prom$entrezgene,
universe      = as.character(ATACseq_all_genes$entrezgene),
OrgDb         = org.Hs.eg.db,
ont           = "BP",
pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = T)
dotplot(atacUp_ego_prom)



#ATAC significativi + K4me3

View(Merged_Atac_K4me3)
Merged_Atac_K4me3_sig=Merged_Atac_K4me3[Merged_Atac_K4me3$toptags_ATAC_gr.Qvalue<=0.05,]
Merged_Atac_K4me3_sig_down=Merged_Atac_K4me3_sig[Merged_Atac_K4me3_sig$toptags_ATAC_gr.logFC<0,]
Merged_Atac_K4me3_sig_up=Merged_Atac_K4me3_sig[Merged_Atac_K4me3_sig$toptags_ATAC_gr.logFC>0,]




ATACseq_sig_K4me3_down<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=Merged_Atac_K4me3_sig_down$toptags_ATAC_gr.SYMBOL,mart=mart)
ATACseq_sig_K4me3_up<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=Merged_Atac_K4me3_sig_up$toptags_ATAC_gr.SYMBOL,mart=mart)
ATACseq_K4me3_all<-  getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=Merged_Atac_K4me3$toptags_ATAC_gr.SYMBOL,mart=mart)

atac_K4me3_sig<-getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=Merged_Atac_K4me3_sig_up$toptags_ATAC_gr.SYMBOL,mart=mart)

atac_sig_down_K4me3 <- enrichGO(gene          = ATACseq_sig_K4me3_down$entrezgene,
universe      = as.character(ATACseq_K4me3_all$entrezgene),
OrgDb         = org.Hs.eg.db,
ont           = "BP",
pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = T)
dotplot(atac_sig_down_K4me3)

atac_sig_up_K4me3 <- enrichGO(gene          = ATACseq_sig_K4me3_up$entrezgene,
universe      = as.character(ATACseq_K4me3_all$entrezgene),
OrgDb         = org.Hs.eg.db,
ont           = "BP",
pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = T)
dotplot(atac_sig_up_K4me3)


atac_sig_all_K4me3 <- enrichGO(gene          = ATACseq_sig_K4me3_sig$entrezgene,
universe      = as.character(ATACseq_K4me3_all$entrezgene),
OrgDb         = org.Hs.eg.db,
ont           = "BP",
pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = T)









#Mettere su dropbox bigwigs normalizzati (epi)
#Mettere su dropbox tabelle CPM delle diverse tecniche
#identificare superenhancers
#Rose - algorithm for superenancers
#ANALYSI CON HOMER



#PBS -l select=1:app=java:ncpus=4:mem=20gb
#PBS -P 161 
#PBS -q workq
#PBS -m ae

while read bam peaks;
do  prefix=$(echo $peaks| sed 's/.narrowPeak/.coverage/g');
echo $bam $peaks $prefix;
coverageBed -a $peaks  -b $bam > $prefix ;
done < tt



while read bams peaks; do
echo /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/${file}.bam; 
samtools view -F 0xD40 -u  -q 15  /lustre1/workspace/DiMicco/prj_ALL_pilot/ChIPseq/BAM/${file}.bam |
bamToBed |slopBed -i stdin -l 0 -r 150 -g /lustre2/scratch/gbucci/Genome/hs38DH.tab |
genomeCoverageBed -i stdin -bg -g /lustre2/scratch/gbucci/Genome/hs38DH.tab -scale ${normfactor} |
intersectBed -v -a stdin -b /lustre2/scratch/dcittaro/Ref/hg38.mask.bed  | grep -v chrY |
wigToBigWig stdin /lustre2/scratch/gbucci/Genome/hs38DH.tab /lustre2/scratch/fsantaniello/DiMicco/K4me3/BAM/${file}.bw;
done < /lustre2/scratch/fsantaniello/DiMicco/K4me3/BAM/NormFactors_K4me3seq.txt





ATAC_prom_genes<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=unique(atac_prom),mart=mart)



ATAC_RNA_sig=ATAC_RNA[ATAC_RNA$Qvalue.x<=0.01&ATAC_RNA$Qvalue.y<=0.01,]
ATAC_RNA_sig_down=ATAC_RNA[ATAC_RNA$Qvalue.x<=0.01&ATAC_RNA$Qvalue.y<=0.01&ATAC_RNA$logFC.x<0&ATAC_RNA$logFC.y<0,]
ATAC_RNA_sig_up=ATAC_RNA[ATAC_RNA$Qvalue.x<=0.01&ATAC_RNA$Qvalue.y<=0.01&ATAC_RNA$logFC.x>0&ATAC_RNA$logFC.y>0,]
RNAseq_all_genes<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=toptags_RNA$genes,mart=mart)
ATAC_RNAseq_genes_down<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=ATAC_RNA_sig_down$genes,mart=mart)
ATAC_RNAseq_genes_up<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=ATAC_RNA_sig_up$genes,mart=mart)


atac_rna_sig_down <- enrichGO(gene = ATAC_RNAseq_genes_down$entrezgene,
universe      = as.character(RNAseq_all_genes$entrezgene),
OrgDb         = org.Hs.eg.db,
ont           = "BP",
pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = T)
dotplot(atac_rna_sig_down)

atac_rna_sig_up <- enrichGO(gene = ATAC_RNAseq_genes_up$entrezgene,
universe      = as.character(RNAseq_all_genes$entrezgene),
OrgDb         = org.Hs.eg.db,
ont           = "BP",
pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = T)

dotplot(atac_rna_sig_up)





#PBS -l select=1:app=java:ncpus=4:mem=20gb
#PBS -P 161 
#PBS -q workq
#PBS -m ae

while read file  libsize; do
echo /lustre1/workspace/DiMicco/prj_ALL_pilot/491_ATACseq/BAM/${file}.bam; 
samtools view -F 0xD40 -u  -q 15  /lustre1/workspace/DiMicco/prj_ALL_pilot/491_ATACseq/BAM/${file}.bam |
bamToBed |slopBed -i stdin -l 0 -r 150 -g /lustre2/scratch/gbucci/Genome/hs38DH.tab |
genomeCoverageBed -i stdin -bg -g /lustre2/scratch/gbucci/Genome/hs38DH.tab  |
awk '{OFS="\t"; print $1,$2,$3,($4/libsize)*1000000}' libsize="$libsize" - | 
intersectBed -v -a stdin -b /lustre2/scratch/dcittaro/Ref/hg38.mask.bed  | grep -v chrY |
wigToBigWig stdin /lustre2/scratch/gbucci/Genome/hs38DH.tab /lustre2/scratch/fsantaniello/DiMicco/ATAC_seq/NormedBigWigs/${file}_NormedFullLibSize.bw;
done < /lustre2/scratch/fsantaniello/DiMicco/ATAC_seq/BAM/LibSizes_ATAC.txt






"ATAC_chr","ATAC_start","ATAC_end","ATAC_logFC","ATAC_PValue","ATAC_Qvalue","genes","RNA_logFC","RNA_PValue","RNA_Qvalue","ATAC_annotation","ATAC_distanceToTSS","K4me3Overlap","K27acOverlap","K4me3_chr","K4me3_start","K4me3_end","K4me3_annotation","K4me3_distanceToTSS","K4me3_logFC","K4me3_PValue","K4me3_Qvalue","K27ac_chr","K27ac_start","K27ac_end","K27ac_annotation","K27ac_distanceToTSS","K27ac_logFC","K27ac_Pvalue","K27ac_Qvalue"



breakList=seq(-5, 5, by = 0.1) 
pheatmap( immuno_genes [,2:9], cluster_rows=T,cluster_cols=F,color=colorRampPalette( c("green", "black", "red"), space="rgb")(length(breakList)))



immuno_genes=c('DYNC1I2',
'TREM2',
'ACTR10',
'HLA-DPB1',
'HLA-DQA2',
'HLA-DRB1',
'HLA-DQB1',
'SEC31A',
'AP2A2',
'DYNLL1',
'LGMN',
'AP1S1',
'DCTN6',
'HLA-DRB5',
'SEC24C',
'MARCH1',
'KIF2A',
'CLTA',
'THBS1',
'RILP',
'DCTN5',
'CTSE',
'KIF3B',
'ACTR1A',
'HLA-DPA1',
'FCER1G',
'DYNC1H1',
'DCTN1',
'CENPE',
'DYNC1I1',
'KIF5A',
'OSBPL1A',
'LAG3',
'CLTC',
'SEC24D',
'KIF2B',
'HLA-DOB',
'HLA-DRB4',
'HLA-DQA1',
'HLA-DRA',
'IFI30',
'DCTN4',
'CTSF',
'CANX',
'CTSD',
'DCTN3',
'KIF26A',
'AP1B1',
'HLA-DOA',
'CAPZA1',
'RACGAP1',
'KLC2',
'AP2B1',
'KIFAP3',
'TRAF6',
'KIF2C',
'KIF15',
'SEC23A',
'AP1M2',
'SH3GL2',
'AP2M1',
'AP2S1',
'FCGR2B',
'HLA-DMB',
'HLA-DMA',
'KIF4B',
'SEC13',
'KLC1',
'KIF23',
'DNM2',
'DYNLL2',
'KIF3A',
'DYNC1LI2',
'KIF22',
'CAPZB',
'CAPZA2',
'PYCARD',
'ARF1',
'CTSL',
'CTSS',
'KIF4A',
'HLA-DQB2',
'SEC24B',
'SEC24A',
'AP1G1',
'KIF18A',
'KIF11',
'SAR1B',
'AP2A1',
'DYNC1LI1',
'CAPZA3',
'ACTR1B',
'E7ENX8',
'HLA-DRB3',
'DCTN2',
'RAB7A',
'AP1M1',
'SPTBN2',
'AP1S3',
'KIF3C',
'CD74',
'CTSV',
'AP1S2')


geni_cristina=c('CDKN2C|STATH|SIX5|PALM',
'C120d57',
'PPP1R14A',
'KCTD15',
'RPS2',
'AOP11',
'LOC650698',
'PTHLH',
'FBXL6',
'RAPGEF3',
'GPR172A',
'LOC649447',
'PPP1R14B',
'RPL34',
'LOC642755',
'RBM9',
'GPC4',
'FKBP9L',
'SLC35D2',
'RMND5B',
'WIBG',
'C1orf97',
'CSNK2B',
'FLJ39632',
'LOC440160',
'C11orf9',
'CD24',
'ATP7B',
'PDE9A',
'A4GALT',
'C6orf59',
'BRDG1',
'CD276',
'C0320',
'C10orf125',
'DDAH1',
'DLK1',
'NRBP2',
'RXRA',
'F5',
'SCGB3A1',
'SSR4',
'IGFBP2',
'MMP15',
'NAT14',
'NT5DC2',
'DBN1',
'AVEN',
'BATF3',
'DOCK6',
'OFOTM1',
'TRPC1',
'C19orf51',
'C9orf61',
'GZMA',
'PVRL2',
'LARP6',
'TRIP6',
'PSEN2',
'C21orf63',
'LOC644390',
'PPT2',
'STK33',
'MEIS1',
'IRX5',
'LASS1',
'LOC643911',
'ADAM28',
'ZIK1',
'ZNF544',
'F8X021',
'ZNF135',
'LOC644133',
'SCRN1',
'TOB1',
'EREG',
'F13A1',
'EIF4E3',
'PROK2',
'ANGPT1',
'ITGAL',
'LAT2',
'TNFSF4',
'HIST2H2BE',
'NPTX2',
'BEX2',
'MAP7',
'NLRP3',
'PHF11',
'UBA3',
'SBF2',
'TYROBP',
'HLA-DMA',
'HLA-DRA',
'CIITA',
'HLA-DOA',
'HLA-DMB',
'LOC652479',
'CPNE8',
'C10or54',
'C1orf24',
'PRAME',
'MAN1C1',
'SLC22A15',
'PRKCZ',
'MTX3',
'PLEKHA1')





#Confronto analisi K4me3 DESeq2 vs EdgeR
#ANALISI STANDARD EDGER (TMM, ACUTOFF)

k4me3_all=readPeakFile("All_K4me3_Peaks_Counts_MultiCov_New_NoReps_New.txt",head=F)
k4me3_all_anno <- annotatePeak(peak = k4me3_all, tssRegion=c(-2000, 2000),TxDb=gencode_rnaseqtxdb)#,annoDb="org.Hs.eg.db")
 #plotAnnoPie(k4me3_all_anno)


k4me3_all_anno_df=as.data.frame(k4me3_all_anno)
k4me3_all_anno_df$annotation = sub(" .*", "", k4me3_all_anno_df$annotation)
colnames(k4me3_all_anno_df)=c(colnames(k4me3_all_anno_df)[1:5],"UPN3_1","UPN3_2","UPN3_3","UPN3_4","UPN3_5","UPN3_6","UPN3_7","UPN3_8","UPN3_9","UPN3_10","UPN3_11","UPN3_12","UPN3_13","UPN3_14",colnames(k4me3_all_anno_df)[20:28])

#head(k4me3_all_anno_df)
k4me3_counts_new=k4me3_all_anno_df[,c(1,2,3,26,20,28,6:19)]
k4me3_y_new=DGEList(counts=k4me3_counts_new[,c(7:10,14:17)],genes=k4me3_counts_new[,c(1:6)],group=c("D_0","D_1","D_1","D_1","R_0","R_1","R_1","R_1"))#,lib.size = c(22089014,23497699,30146073,21596634,24792817,22820342,26900163,45421712))
k4me3_y_new = calcNormFactors(k4me3_y_new,Acutoff=-15)
k4me3_new_keep = rowSums(cpm(k4me3_y_new) > 5) >= 3
k4me3_y_new= k4me3_y_new[k4me3_new_keep, ]
library(qvalue)
 df = targets[row.names(k4me3_y_new$samples),]
 D = as.factor(df$disease)
 P = as.factor(df$passage)
 design = model.matrix(~D+P)
 design
 k4me3_y_new= estimateDisp(k4me3_y_new, design)
 fit = glmFit(k4me3_y_new, design)
 lrt <- glmLRT(fit, coef=2)
 toptags_k4me3=as.data.frame(topTags(lrt,sort.by="none",n="all"))
 q = qvalue(toptags_k4me3$PValue)
  toptags_k4me3$Qvalue = q$qvalues
 toptags_k4me3_DEG=toptags_k4me3[toptags_k4me3$Qvalue<0.05,]
k4me3_y_new$samples$group=as.factor(design[,2])


toptags_k4me3_down=toptags_k4me3[toptags_k4me3$Qvalue<=0.05&toptags_k4me3$logFC<0,]
toptags_k4me3_up=toptags_k4me3[toptags_k4me3$Qvalue<=0.05&toptags_k4me3$logFC>0,]

nrow(toptags_k4me3_down)
[1] 4352
nrow(toptags_k4me3_up)
[1] 481
plotSmear(lrt,smooth.scatter=T)




#ANALISI DESEQ2
df_deseq=df
df_deseq$passage=as.factor(df_deseq$passage)
dds=DESeqDataSetFromMatrix(k4me3_counts_new[,c(7:10,14:17)],df_deseq,~passage+ disease)
featureData <- data.frame(gene=k4me3_counts_new[,1:6])
mcols(dds) <- DataFrame(mcols(dds), featureData)
keep = rowSums(counts(dds) > 5) >= 3
dds <- dds[keep,]
dds$disease <- relevel(dds$disease, ref = "D")

dds<- DESeq(dds)#,betaPrior = F,test = "LRT",reduced=model.matrix(~passage,data=df),full = model.matrix(~passage+disease,data=df))
res=results(object = dds,name='disease_R_vs_D')
summary(res)
plotMA(res)

```
out of 71505 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)     : 3780, 5.3% 
LFC < 0 (down)   : 2205, 3.1% 
outliers [1]     : 142, 0.2% 
low counts [2]   : 16634, 23% 
(mean count < 6)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```
res_complete=cbind(mcols(dds),res)
colnames(res_complete)=gsub("gene\\.",'',colnames(res_complete))
res_complete=res_complete[,c(1:6,33,37)]

#Confronto con EdgeR
merge_deseq_edger=merge(res_complete,toptags_k4me3,by=colnames(res_complete)[1:6])
merge_deseq_edger_prom=merge_deseq_edger[merge_deseq_edger$annotation=='Promoter',]

fit<- lm(merge_deseq_edger$log2FoldChange~merge_deseq_edger$logFC)
summary(fit)
```
Call:
lm(formula = merge_deseq_edger$log2FoldChange ~ merge_deseq_edger$logFC)

Residuals:
     Min       1Q   Median       3Q      Max 
-0.69894 -0.03589  0.00148  0.03696  0.98367 

Coefficients:
                         Estimate Std. Error t value Pr(>|t|)    
(Intercept)             0.5973416  0.0003459    1727   <2e-16 ***
merge_deseq_edger$logFC 1.0234761  0.0003006    3405   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.0644 on 40480 degrees of freedom
Multiple R-squared:  0.9965,	Adjusted R-squared:  0.9965 
F-statistic: 1.159e+07 on 1 and 40480 DF,  p-value: < 2.2e-16
```


fit<- lm(merge_deseq_edger_prom$log2FoldChange~merge_deseq_edger_prom$logFC)
summary(fit)
```
Call:
lm(formula = merge_deseq_edger_prom$log2FoldChange ~ merge_deseq_edger_prom$logFC)

Residuals:
     Min       1Q   Median       3Q      Max 
-0.44470 -0.03379 -0.00240  0.03189  0.92749 

Coefficients:
                             Estimate Std. Error t value Pr(>|t|)    
(Intercept)                  0.592486   0.000497    1192   <2e-16 ***
merge_deseq_edger_prom$logFC 1.033363   0.000501    2063   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.05813 on 14103 degrees of freedom
Multiple R-squared:  0.9967,	Adjusted R-squared:  0.9967 
F-statistic: 4.254e+06 on 1 and 14103 DF,  p-value: < 2.2e-16
```


#ANALISI STANDARD EDGER SENZA ACUTOFF (TMM only)

k4me3_all=readPeakFile("All_K4me3_Peaks_Counts_MultiCov_New_NoReps_New.txt",head=F)
k4me3_all_anno <- annotatePeak(peak = k4me3_all, tssRegion=c(-2000, 2000),TxDb=gencode_rnaseqtxdb)#,annoDb="org.Hs.eg.db")
 #plotAnnoPie(k4me3_all_anno)


k4me3_all_anno_df=as.data.frame(k4me3_all_anno)
k4me3_all_anno_df$annotation = sub(" .*", "", k4me3_all_anno_df$annotation)
colnames(k4me3_all_anno_df)=c(colnames(k4me3_all_anno_df)[1:5],"UPN3_1","UPN3_2","UPN3_3","UPN3_4","UPN3_5","UPN3_6","UPN3_7","UPN3_8","UPN3_9","UPN3_10","UPN3_11","UPN3_12","UPN3_13","UPN3_14",colnames(k4me3_all_anno_df)[20:28])

#head(k4me3_all_anno_df)
k4me3_counts_new=k4me3_all_anno_df[,c(1,2,3,26,20,28,6:19)]
k4me3_y_new=DGEList(counts=k4me3_counts_new[,c(7:10,14:17)],genes=k4me3_counts_new[,c(1:6)],group=c("D_0","D_1","D_1","D_1","R_0","R_1","R_1","R_1"))#,lib.size = c(22089014,23497699,30146073,21596634,24792817,22820342,26900163,45421712))
k4me3_y_new = calcNormFactors(k4me3_y_new)#,Acutoff=-15)
k4me3_new_keep = rowSums(cpm(k4me3_y_new) > 5) >= 3
k4me3_y_new= k4me3_y_new[k4me3_new_keep, ]
library(qvalue)
 df = targets[row.names(k4me3_y_new$samples),]
 D = as.factor(df$disease)
 P = as.factor(df$passage)
 design = model.matrix(~D+P)
 design
 k4me3_y_new= estimateDisp(k4me3_y_new, design)
 fit = glmFit(k4me3_y_new, design)
 lrt <- glmLRT(fit, coef=2)
 toptags_k4me3=as.data.frame(topTags(lrt,sort.by="none",n="all"))
 q = qvalue(toptags_k4me3$PValue)
  toptags_k4me3$Qvalue = q$qvalues
 toptags_k4me3_DEG=toptags_k4me3[toptags_k4me3$Qvalue<0.05,]
k4me3_y_new$samples$group=as.factor(design[,2])


toptags_k4me3_down=toptags_k4me3[toptags_k4me3$Qvalue<=0.05&toptags_k4me3$logFC<0,]
toptags_k4me3_up=toptags_k4me3[toptags_k4me3$Qvalue<=0.05&toptags_k4me3$logFC>0,]

nrow(toptags_k4me3_down)
[1] 1266
nrow(toptags_k4me3_up)
[1] 4518

plotSmear(lrt,smooth.scatter=T)



#Confronto con EdgeR
merge_deseq_edger=merge(res_complete,toptags_k4me3,by=colnames(res_complete)[1:6])
merge_deseq_edger_prom=merge_deseq_edger[merge_deseq_edger$annotation=='Promoter',]

fit<- lm(merge_deseq_edger$log2FoldChange~merge_deseq_edger$logFC)
summary(fit)
```
Call:
lm(formula = merge_deseq_edger$log2FoldChange ~ merge_deseq_edger$logFC)

Residuals:
     Min       1Q   Median       3Q      Max 
-0.47824 -0.01410 -0.00342  0.00802  1.33247 

Coefficients:
                         Estimate Std. Error t value Pr(>|t|)    
(Intercept)             -0.116211   0.000263  -441.9   <2e-16 ***
merge_deseq_edger$logFC  1.017006   0.000240  4237.6   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.04819 on 37709 degrees of freedom
Multiple R-squared:  0.9979,	Adjusted R-squared:  0.9979 
F-statistic: 1.796e+07 on 1 and 37709 DF,  p-value: < 2.2e-16
```



fit<- lm(merge_deseq_edger_prom$log2FoldChange~merge_deseq_edger_prom$logFC)
summary(fit)

```
Call:
lm(formula = merge_deseq_edger_prom$log2FoldChange ~ merge_deseq_edger_prom$logFC)

Residuals:
     Min       1Q   Median       3Q      Max 
-0.37765 -0.01217 -0.00272  0.00816  1.28327 

Coefficients:
                               Estimate Std. Error t value Pr(>|t|)    
(Intercept)                  -0.1210818  0.0003843    -315   <2e-16 ***
merge_deseq_edger_prom$logFC  1.0162876  0.0003427    2965   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.03922 on 13884 degrees of freedom
Multiple R-squared:  0.9984,	Adjusted R-squared:  0.9984 
F-statistic: 8.793e+06 on 1 and 13884 DF,  p-value: < 2.2e-16
```



#ANALISI EDGER CON METODO RLE e  ACUTOFF

k4me3_all=readPeakFile("All_K4me3_Peaks_Counts_MultiCov_New_NoReps_New.txt",head=F)
k4me3_all_anno <- annotatePeak(peak = k4me3_all, tssRegion=c(-2000, 2000),TxDb=gencode_rnaseqtxdb)#,annoDb="org.Hs.eg.db")
 #plotAnnoPie(k4me3_all_anno)


k4me3_all_anno_df=as.data.frame(k4me3_all_anno)
k4me3_all_anno_df$annotation = sub(" .*", "", k4me3_all_anno_df$annotation)
colnames(k4me3_all_anno_df)=c(colnames(k4me3_all_anno_df)[1:5],"UPN3_1","UPN3_2","UPN3_3","UPN3_4","UPN3_5","UPN3_6","UPN3_7","UPN3_8","UPN3_9","UPN3_10","UPN3_11","UPN3_12","UPN3_13","UPN3_14",colnames(k4me3_all_anno_df)[20:28])

#head(k4me3_all_anno_df)
k4me3_counts_new=k4me3_all_anno_df[,c(1,2,3,26,20,28,6:19)]
k4me3_y_new=DGEList(counts=k4me3_counts_new[,c(7:10,14:17)],genes=k4me3_counts_new[,c(1:6)],group=c("D_0","D_1","D_1","D_1","R_0","R_1","R_1","R_1"))#,lib.size = c(22089014,23497699,30146073,21596634,24792817,22820342,26900163,45421712))
k4me3_y_new = calcNormFactors(k4me3_y_new,method='RLE',Acutoff=-15)
k4me3_new_keep = rowSums(cpm(k4me3_y_new) > 5) >= 3
k4me3_y_new= k4me3_y_new[k4me3_new_keep, ]
library(qvalue)
 df = targets[row.names(k4me3_y_new$samples),]
 D = as.factor(df$disease)
 P = as.factor(df$passage)
 design = model.matrix(~D+P)
 design
 k4me3_y_new= estimateDisp(k4me3_y_new, design)
 fit = glmFit(k4me3_y_new, design)
 lrt <- glmLRT(fit, coef=2)
 toptags_k4me3=as.data.frame(topTags(lrt,sort.by="none",n="all"))
 q = qvalue(toptags_k4me3$PValue)
  toptags_k4me3$Qvalue = q$qvalues
 toptags_k4me3_DEG=toptags_k4me3[toptags_k4me3$Qvalue<0.05,]
k4me3_y_new$samples$group=as.factor(design[,2])


toptags_k4me3_down=toptags_k4me3[toptags_k4me3$Qvalue<=0.05&toptags_k4me3$logFC<0,]
toptags_k4me3_up=toptags_k4me3[toptags_k4me3$Qvalue<=0.05&toptags_k4me3$logFC>0,]

nrow(toptags_k4me3_down)
[1] 1211
nrow(toptags_k4me3_up)
[1] 4254

plotSmear(lrt,smooth.scatter=T)


#Confronto con EdgeR
merge_deseq_edger=merge(res_complete,toptags_k4me3,by=colnames(res_complete)[1:6])
merge_deseq_edger_prom=merge_deseq_edger[merge_deseq_edger$annotation=='Promoter',]

fit<- lm(merge_deseq_edger$log2FoldChange~merge_deseq_edger$logFC)
summary(fit)
```
Call:
lm(formula = merge_deseq_edger$log2FoldChange ~ merge_deseq_edger$logFC)

Residuals:
     Min       1Q   Median       3Q      Max 
-0.48233 -0.01262 -0.00368  0.00572  1.36827 

Coefficients:
                          Estimate Std. Error t value Pr(>|t|)    
(Intercept)             -0.0979372  0.0002690    -364   <2e-16 ***
merge_deseq_edger$logFC  1.0160249  0.0002466    4120   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.04904 on 37250 degrees of freedom
Multiple R-squared:  0.9978,	Adjusted R-squared:  0.9978 
F-statistic: 1.698e+07 on 1 and 37250 DF,  p-value: < 2.2e-16
```




fit<- lm(merge_deseq_edger_prom$log2FoldChange~merge_deseq_edger_prom$logFC)
summary(fit)

```
Call:
lm(formula = merge_deseq_edger_prom$log2FoldChange ~ merge_deseq_edger_prom$logFC)

Residuals:
     Min       1Q   Median       3Q      Max 
-0.38122 -0.01048 -0.00317  0.00603  1.31038 

Coefficients:
                               Estimate Std. Error t value Pr(>|t|)    
(Intercept)                  -0.1027046  0.0003793  -270.7   <2e-16 ***
merge_deseq_edger_prom$logFC  1.0144109  0.0003403  2981.4   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.0389 on 13861 degrees of freedom
Multiple R-squared:  0.9984,	Adjusted R-squared:  0.9984 
F-statistic: 8.889e+06 on 1 and 13861 DF,  p-value: < 2.2e-16
```




#ANALISI EDGER CON METODO RLE e  ACUTOFF

k4me3_all=readPeakFile("All_K4me3_Peaks_Counts_MultiCov_New_NoReps_New.txt",head=F)
k4me3_all_anno <- annotatePeak(peak = k4me3_all, tssRegion=c(-2000, 2000),TxDb=gencode_rnaseqtxdb)#,annoDb="org.Hs.eg.db")
 #plotAnnoPie(k4me3_all_anno)


k4me3_all_anno_df=as.data.frame(k4me3_all_anno)
k4me3_all_anno_df$annotation = sub(" .*", "", k4me3_all_anno_df$annotation)
colnames(k4me3_all_anno_df)=c(colnames(k4me3_all_anno_df)[1:5],"UPN3_1","UPN3_2","UPN3_3","UPN3_4","UPN3_5","UPN3_6","UPN3_7","UPN3_8","UPN3_9","UPN3_10","UPN3_11","UPN3_12","UPN3_13","UPN3_14",colnames(k4me3_all_anno_df)[20:28])

#head(k4me3_all_anno_df)
k4me3_counts_new=k4me3_all_anno_df[,c(1,2,3,26,20,28,6:19)]
k4me3_y_new=DGEList(counts=k4me3_counts_new[,c(7:10,14:17)],genes=k4me3_counts_new[,c(1:6)],group=c("D_0","D_1","D_1","D_1","R_0","R_1","R_1","R_1"))#,lib.size = c(22089014,23497699,30146073,21596634,24792817,22820342,26900163,45421712))
k4me3_y_new = calcNormFactors(k4me3_y_new,method='RLE')#,Acutoff=-15)
k4me3_new_keep = rowSums(cpm(k4me3_y_new) > 5) >= 3
k4me3_y_new= k4me3_y_new[k4me3_new_keep, ]
library(qvalue)
 df = targets[row.names(k4me3_y_new$samples),]
 D = as.factor(df$disease)
 P = as.factor(df$passage)
 design = model.matrix(~D+P)
 design
 k4me3_y_new= estimateDisp(k4me3_y_new, design)
 fit = glmFit(k4me3_y_new, design)
 lrt <- glmLRT(fit, coef=2)
 toptags_k4me3=as.data.frame(topTags(lrt,sort.by="none",n="all"))
 q = qvalue(toptags_k4me3$PValue)
  toptags_k4me3$Qvalue = q$qvalues
 toptags_k4me3_DEG=toptags_k4me3[toptags_k4me3$Qvalue<0.05,]
k4me3_y_new$samples$group=as.factor(design[,2])


toptags_k4me3_down=toptags_k4me3[toptags_k4me3$Qvalue<=0.05&toptags_k4me3$logFC<0,]
toptags_k4me3_up=toptags_k4me3[toptags_k4me3$Qvalue<=0.05&toptags_k4me3$logFC>0,]

nrow(toptags_k4me3_down)
[1] 1211
nrow(toptags_k4me3_up)
[1] 4254

plotSmear(lrt,smooth.scatter=T)


#Confronto con EdgeR
merge_deseq_edger=merge(res_complete,toptags_k4me3,by=colnames(res_complete)[1:6])
merge_deseq_edger_prom=merge_deseq_edger[merge_deseq_edger$annotation=='Promoter',]

fit<- lm(merge_deseq_edger$log2FoldChange~merge_deseq_edger$logFC)
summary(fit)
```
Call:
lm(formula = merge_deseq_edger$log2FoldChange ~ merge_deseq_edger$logFC)

Residuals:
     Min       1Q   Median       3Q      Max 
-0.48233 -0.01262 -0.00368  0.00572  1.36827 

Coefficients:
                          Estimate Std. Error t value Pr(>|t|)    
(Intercept)             -0.0979372  0.0002690    -364   <2e-16 ***
merge_deseq_edger$logFC  1.0160249  0.0002466    4120   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.04904 on 37250 degrees of freedom
Multiple R-squared:  0.9978,	Adjusted R-squared:  0.9978 
F-statistic: 1.698e+07 on 1 and 37250 DF,  p-value: < 2.2e-16
```


fit<- lm(merge_deseq_edger_prom$log2FoldChange~merge_deseq_edger_prom$logFC)
summary(fit)

```
Call:
lm(formula = merge_deseq_edger_prom$log2FoldChange ~ merge_deseq_edger_prom$logFC)

Residuals:
     Min       1Q   Median       3Q      Max 
-0.38122 -0.01048 -0.00317  0.00603  1.31038 

Coefficients:
                               Estimate Std. Error t value Pr(>|t|)    
(Intercept)                  -0.1027046  0.0003793  -270.7   <2e-16 ***
merge_deseq_edger_prom$logFC  1.0144109  0.0003403  2981.4   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.0389 on 13861 degrees of freedom
Multiple R-squared:  0.9984,	Adjusted R-squared:  0.9984 
F-statistic: 8.889e+06 on 1 and 13861 DF,  p-value: < 2.2e-16
```





#Confronto analisi ATAC DESeq2 vs EdgeR
#ANALISI STANDARD EDGER (TMM, ACUTOFF)




atac_all=readPeakFile("All_Atac_Peaks_Counts_MultiCov_NoReps.txt",head=F)
atac_all_anno <- annotatePeak(peak = atac_all, tssRegion=c(-2000, 2000),TxDb=gencode_rnaseqtxdb) #,annoDb="org.Hs.eg.db")
atac_all_anno_df=as.data.frame(atac_all_anno)
colnames(atac_all_anno_df)=c(colnames(atac_all_anno_df)[1:5],"UPN3_1","UPN3_2","UPN3_3","UPN3_4","UPN3_5","UPN3_6","UPN3_7","UPN3_8","UPN3_9","UPN3_10","UPN3_11","UPN3_13","UPN3_14",colnames(atac_all_anno_df)[19:27])

atac_all_anno_df$annotation = sub(" .*", "",atac_all_anno_df$annotation)
atac_counts_new=atac_all_anno_df[,c(1,2,3,4,25,19,27,6:18)]
atac_y_new=DGEList(counts=atac_counts_new[,c(8:11,15:18)],genes=atac_counts_new[,c(1:7)],group=c("D_0","D_1","D_1","D_1","R_0","R_1","R_1","R_1"))
atac_y_new = calcNormFactors(atac_y_new,Acutoff=-18.5)
atac_new_keep = rowSums(cpm(atac_y_new) > 5) >= 3
atac_y_new = atac_y_new[atac_new_keep, ]
targets=read.table('Targets_for_DEAnalysis.txt')
targets
df = targets[row.names(atac_y_new$samples),]
df
D = as.factor(df$disease)
P = as.factor(df$passage)
design = model.matrix(~D+P)
design
atac_y_new = estimateDisp(atac_y_new, design)
fit = glmFit(atac_y_new, design) 
lrt <- glmLRT(fit, coef=2)
toptags_ATAC=as.data.frame(topTags(lrt,sort.by="PValue",n="all"))
q = qvalue(toptags_ATAC$PValue)
toptags_ATAC$Qvalue = q$qvalues
atac_y_new$samples$group=as.factor(design[,2])
plotSmear(atac_y_new,smooth.scatter=T)
ATAC_up=toptags_ATAC[toptags_ATAC$Qvalue<=0.01&toptags_ATAC$logFC>0,]
ATAC_down=toptags_ATAC[toptags_ATAC$Qvalue<=0.01&toptags_ATAC$logFC<0,]

nrow(ATAC_up)
[1] 2880
nrow(ATAC_down)
[1] 3853





#ANALISI DESEQ2
df_deseq=df
df_deseq$passage=as.factor(df_deseq$passage)
dds=DESeqDataSetFromMatrix(atac_counts_new[,c(8:11,15:18)],df_deseq,~passage+ disease)
featureData <- data.frame(gene=atac_counts_new[,1:7])
mcols(dds) <- DataFrame(mcols(dds), featureData)
keep = rowSums(counts(dds) > 5) >= 3
dds <- dds[keep,]
dds$disease <- relevel(dds$disease, ref = "D")

dds<- DESeq(dds)#,betaPrior = F,test = "LRT",reduced=model.matrix(~passage,data=df),full = model.matrix(~passage+disease,data=df))
res=results(object = dds,name='disease_R_vs_D')
summary(res)
```
out of 124650 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)     : 9859, 7.9% 
LFC < 0 (down)   : 11227, 9% 
outliers [1]     : 25, 0.02% 
low counts [2]   : 0, 0% 
(mean count < 2)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```
plotMA(res)

res_complete=cbind(mcols(dds),res)
colnames(res_complete)=gsub("gene\\.",'',colnames(res_complete))
res_complete=res_complete[,c(1:7,34,38)]




#Confronto con EdgeR
merge_deseq_edger=merge(res_complete,toptags_ATAC,by=colnames(res_complete)[1:7])
merge_deseq_edger_prom=merge_deseq_edger[merge_deseq_edger$annotation=='Promoter',]

fit<- lm(merge_deseq_edger$log2FoldChange~merge_deseq_edger$logFC)
summary(fit)
```
Call:
lm(formula = merge_deseq_edger$log2FoldChange ~ merge_deseq_edger$logFC)

Residuals:
    Min      1Q  Median      3Q     Max 
-3.4679 -0.0042  0.0009  0.0054  0.5980 

Coefficients:
                         Estimate Std. Error t value Pr(>|t|)    
(Intercept)             -0.156320   0.000104   -1504   <2e-16 ***
merge_deseq_edger$logFC  1.007286   0.000126    7994   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.025 on 58226 degrees of freedom
Multiple R-squared:  0.9991,	Adjusted R-squared:  0.9991 
F-statistic: 6.391e+07 on 1 and 58226 DF,  p-value: < 2.2e-16
```


fit<- lm(merge_deseq_edger_prom$log2FoldChange~merge_deseq_edger_prom$logFC)
summary(fit)

```
Call:
lm(formula = merge_deseq_edger_prom$log2FoldChange ~ merge_deseq_edger_prom$logFC)

Residuals:
    Min      1Q  Median      3Q     Max 
-3.4541 -0.0033  0.0011  0.0054  0.1543 

Coefficients:
                               Estimate Std. Error t value Pr(>|t|)    
(Intercept)                  -0.1558662  0.0003133  -497.5   <2e-16 ***
merge_deseq_edger_prom$logFC  1.0103708  0.0004958  2038.0   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.0349 on 13374 degrees of freedom
Multiple R-squared:  0.9968,	Adjusted R-squared:  0.9968 
F-statistic: 4.153e+06 on 1 and 13374 DF,  p-value: < 2.2e-16
```



#ANALISI STANDARD EDGER normalizzazione RLE e con le lib.sizes native
atac_all_anno_df$annotation = sub(" .*", "",atac_all_anno_df$annotation)
atac_counts_new=atac_all_anno_df[,c(1,2,3,4,25,19,27,6:18)]
atac_y_new=DGEList(counts=atac_counts_new[,c(8:11,15:18)],genes=atac_counts_new[,c(1:7)],group=c("D_0","D_1","D_1","D_1","R_0","R_1","R_1","R_1"),lib.size = c(139987133,28373982,43162692,27291027,17232760,40830216,45646788,41369221))
atac_y_new = calcNormFactors(atac_y_new,method='RLE')#,Acutoff=-18.5)
atac_new_keep = rowSums(cpm(atac_y_new) > 5) >= 3
atac_y_new = atac_y_new[atac_new_keep, ]
targets=read.table('Targets_for_DEAnalysis.txt')
targets
df = targets[row.names(atac_y_new$samples),]
df
D = as.factor(df$disease)
P = as.factor(df$passage)
design = model.matrix(~D+P)
design
atac_y_new = estimateDisp(atac_y_new, design)
fit = glmFit(atac_y_new, design) 
lrt <- glmLRT(fit, coef=2)
toptags_ATAC=as.data.frame(topTags(lrt,sort.by="PValue",n="all"))
q = qvalue(toptags_ATAC$PValue)
toptags_ATAC$Qvalue = q$qvalues
atac_y_new$samples$group=as.factor(design[,2])
plotSmear(atac_y_new,smooth.scatter=T)
ATAC_up=toptags_ATAC[toptags_ATAC$Qvalue<=0.01&toptags_ATAC$logFC>0,]
ATAC_down=toptags_ATAC[toptags_ATAC$Qvalue<=0.01&toptags_ATAC$logFC<0,]

nrow(ATAC_up)
[1] 18
nrow(ATAC_down)
[1] 1255











#
#toptags_k4me3_prom=toptags_k4me3[toptags_k4me3$annotation=='Promoter',]
#colnames(toptags_k4me3_prom)=c(colnames(toptags_k4me3_prom)[1:3],'genes',colnames(toptags_k4me3_prom)[5:12])
#
#rna_k4me3_merge=merge(toptags_RNA,toptags_k4me3_prom,by='genes')
#nrow(rna_k4me3_merge[rna_k4me3_merge$Qvalue.x<=0.01&rna_k4me3_merge$logFC.x<0&rna_k4me3_merge$Qvalue.y<=0.05&rna_k4me3_merge$logFC.y<0,])
#[1] 249
#nrow(rna_k4me3_merge[rna_k4me3_merge$Qvalue.x<=0.01&rna_k4me3_merge$logFC.x<0&rna_k4me3_merge$Qvalue.y<=0.05&rna_k4me3_merge$logFC.y>0,])
#[1] 115
#nrow(rna_k4me3_merge[rna_k4me3_merge$Qvalue.x<=0.01&rna_k4me3_merge$logFC.x>0&rna_k4me3_merge$Qvalue.y<=0.05&rna_k4me3_merge$logFC.y>0,])
#[1] 328
#nrow(rna_k4me3_merge[rna_k4me3_merge$Qvalue.x<=0.01&rna_k4me3_merge$logFC.x>0&rna_k4me3_merge$Qvalue.y<=0.05&rna_k4me3_merge$logFC.y<0,])
#[1] 42
#
#rna_k4me3=matrix(c(711,21,330,193),nrow=2,ncol=2)
#colnames(rna_k4me3)=c('UP_rna','DOWN_rna')
#rownames(rna_k4me3)=c('UP_K4me3','DOWN_K4me3')
#```
#          UP_rna DOWN_rna
#UP_K4me3      711      330
#DOWN_K4me3     21      193
#```
#
#
#atac_k4me3_merge=merge(toptags_ATAC_prom,toptags_k4me3_prom,by='genes')
#nrow(atac_k4me3_merge[atac_k4me3_merge$Qvalue.x<=0.01&atac_k4me3_merge$logFC.x<0&atac_k4me3_merge$Qvalue.y<=0.05&atac_k4me3_merge$logFC.y<0,])
#[1] 171
#nrow(atac_k4me3_merge[atac_k4me3_merge$Qvalue.x<=0.01&atac_k4me3_merge$logFC.x<0&atac_k4me3_merge$Qvalue.y<=0.05&atac_k4me3_merge$logFC.y>0,])
#[1] 357
#nrow(atac_k4me3_merge[atac_k4me3_merge$Qvalue.x<=0.01&atac_k4me3_merge$logFC.x>0&atac_k4me3_merge$Qvalue.y<=0.05&atac_k4me3_merge$logFC.y>0,])
#[1] 69
#nrow(atac_k4me3_merge[atac_k4me3_merge$Qvalue.x<=0.01&atac_k4me3_merge$logFC.x>0&atac_k4me3_merge$Qvalue.y<=0.05&atac_k4me3_merge$logFC.y<0,])
#[1] 16
#
#
#atac_k4me3=matrix(c(69,16,357,171),nrow=2,ncol=2)
#colnames(rna_k4me3)=c('UP_atac','DOWN_atac')
#rownames(rna_k4me3)=c('UP_K4me3','DOWN_K4me3')
#
#
#
#fsantaniello@login2:/lustre1/workspace/DiMicco/prj_ALL_pilot/491_ATACseq/BAM$ for i in *UPN3_*bam; do echo -n "$i "; samtools idxstats $i  | awk '{sum += $3} END {print sum}'; done
#ATAC_UPN3_1.bam 139987133
#ATAC_UPN3_10.bam 45646788
#ATAC_UPN3_11.bam 41369221
#ATAC_UPN3_12.bam 62310702
#ATAC_UPN3_13.bam 39057377
#ATAC_UPN3_14.bam 41091232
#ATAC_UPN3_1_Resampled.bam [bam_idxstats] fail to load the index.
#
#ATAC_UPN3_2.bam 28373982
#ATAC_UPN3_3.bam 43162692
#ATAC_UPN3_4.bam 27291027
#ATAC_UPN3_5.bam 26611964
#ATAC_UPN3_6.bam 27633416
#ATAC_UPN3_7.bam 31284686
#ATAC_UPN3_8.bam 17232760
#ATAC_UPN3_9.bam 40830216
#



#Metto insieme tutte le nuvole

#toptags_ATAC_DEG_gr=makeGRangesFromDataFrame(toptags_ATAC_DEG,keep.extra.columns = T)
cloud_atac_gr=makeGRangesFromDataFrame(cloud_atac,keep.extra.columns = T)
cloud_k27ac_gr=makeGRangesFromDataFrame(cloud_k27ac,keep.extra.columns = T)
cloud_k4me3_gr=makeGRangesFromDataFrame(cloud_k4me3,keep.extra.columns = T)

#ATAC + K4me3
cloud_atac_k4me3=mergeByOverlaps(cloud_atac_gr,cloud_k4me3_gr)
Merged_cloud_atac_k4me3=as.data.frame(attr(cloud_atac_k4me3[1:13],"listData"))
Merged_cloud_atac_k4me3=Merged_cloud_atac_k4me3[,c(1:13,22:34)]
colnames(Merged_cloud_atac_k4me3)=c('seqnames','start','end',colnames(Merged_cloud_atac_k4me3)[4:26])

#ATAC + K27ac
cloud_atac_k27ac=mergeByOverlaps(cloud_atac_gr,cloud_k27ac_gr)
Merged_cloud_atac_k27ac=as.data.frame(attr(cloud_atac_k27ac[1:13],"listData"))
Merged_cloud_atac_k27ac=Merged_cloud_atac_k27ac[,c(1:13,22:34)]
colnames(Merged_cloud_atac_k27ac)=c('seqnames','start','end',colnames(Merged_cloud_atac_k27ac)[4:26])


#ATAC_K4me3_K27ac
Cloud_Atac_K27ac_K4me3=merge(Merged_cloud_atac_k4me3,Merged_cloud_atac_k27ac,by=c('seqnames','start','end'),all=T)
Cloud_Atac_K27ac_K4me3_def=Cloud_Atac_K27ac_K4me3[,c(1,2,3,19,42,7,8,9,10,13,20,21,22,23,26,43,44,45,46)]

colnames(Cloud_Atac_K27ac_K4me3_def)=c("seqnames" ,"start", "end","gene_k4me3"        
,"gene_k27ac"        , "annotation_atac"   
,"distanceToTSS_atac", "logFC_atac"        
,"logCPM_atac"       , "Qvalue_atac"       
,"annotation_k4me3"    , "distanceToTSS_k4me3" 
,"logFC_k4me3"         , "logCPM_k4me3"        
,"Qvalue_k4me3"        , "annotation_k27ac"    
,"distanceToTSS_k27ac" , "logFC_k27ac"         
,"logCPM_k27ac")

##K27 annotata su genehancer
#k27ac_all=readPeakFile("../../Desktop/Scratch/K27Ac/Peaks/genehancer.countsk27.upn3.bed",head=F)
#k27ac_all_anno <- annotatePeak(peak = k27ac_all, tssRegion=c(-2000, 2000),TxDb=gencode_rnaseqtxdb) #,annoDb="org.Hs.eg.db")
#k27ac_all_anno_df=as.data.frame(k27ac_all_anno)
#k27ac_all_anno_df$annotation = sub(" .*", "", k27ac_all_anno_df$annotation)
#colnames(k27ac_all_anno_df)=c(colnames(k27ac_all_anno_df)[1:5],'GeneHancer1','GeneHancer2','Genehancer3',"UPN3_10","UPN3_11","UPN3_12","UPN3_13","UPN3_14","UPN3_1","UPN3_2","UPN3_3","UPN3_4","UPN3_5","UPN3_6","UPN3_7","UPN3_8","UPN3_9",colnames(k27ac_all_anno_df)[23:31])
#k27ac_counts_new=k27ac_all_anno_df[,c(1,2,3,4,29,23,31,9:22)]
#k27ac_y_new=DGEList(counts=k27ac_counts_new[,c(13:16,20,21,8,9)],genes=k27ac_counts_new[,c(1:7)],group=c("D_0","D_1","D_1","D_1","R_0","R_1","R_1","R_1"))
##k27ac_y_new = calcNormFactors(k27ac_y_new,method='RLE')#,Acutoff=-18.5)
##k27ac_new_keep = rowSums(cpm(k27ac_y_new) > 5) >= 3
##k27ac_y_new = k27ac_y_new[k27ac_new_keep, ]
#keep <- rowSums(cpm(k27ac_y_new)>10) >= 3 #PRENDO SOLO I GENI CON CPM>1 IN ALMENO 8 CAMPIONI
#k27ac_y_new<- k27ac_y_new[keep, , keep.lib.sizes=FALSE]
#k27ac_y_new = calcNormFactors(k27ac_y_new,method='RLE')#,Acutoff=-18.5)
#library(qvalue)
#df = targets[row.names(k27ac_y_new$samples),]
#D = as.factor(df$disease)
#P = as.factor(df$passage)
#design = model.matrix(~D+P)
#design
#k27ac_y_new= estimateDisp(k27ac_y_new, design)
#fit = glmFit(k27ac_y_new, design)
#lrt <- glmLRT(fit, coef=2)
#toptags_k27ac=as.data.frame(topTags(lrt,sort.by="PValue",n="all"))
#q = qvalue(toptags_k27ac$PValue)
#toptags_k27ac$Qvalue = q$qvalues
# k27ac_y_new$samples$group=as.factor(design[,2])
# plotSmear(lrt,smooth.scatter=T)
##
#
#toptags_k27ac_DEG=toptags_k27ac[toptags_k27ac$Qvalue<=0.05,]
#nrow(toptags_k27ac_DEG[toptags_k27ac_DEG$logFC<0,])
#[1] 10449
#nrow(toptags_k27ac_DEG[toptags_k27ac_DEG$logFC>0,])
#[1] 12568
#
#
#
#
#pdf('PlotSmear_ATAC_RLE_Normalized.pdf')
#plotSmear(lrt,smooth.scatter=T)
#dev.off()
#pdf('OrtogonalPlot_MDS_ATAC_RLE_Normalized.pdf')
##MDS plot su tutte conte DGE dopo normalizzazione
#mds=plotMDS(k27ac_y_new,dim.plot=c(1,2),gene.selection="common")
#mds2=data.frame(x=mds$x,y=mds$y)
#samples=rownames(mds2)
#infos = targets[row.names(k27ac_y_new$samples),]
#mds2=cbind(mds2,as.data.frame(samples),as.data.frame(infos$disease),as.data.frame(infos$passage))
#colnames(mds2)=c("MD1","MD2","samples","condition","passage")
##pdf("PCA_AllSamples_DGE_RNA_Common.pdf")
#ggplot(mds2, aes(MD1, MD2)) + geom_point(aes(colour = condition,shape=as.character(passage)))+ geom_text(aes(label=samples),hjust=0, vjust=0,size=2) + scale_x_continuous(limits = c(-4, 4)) 
#
#
##MDS plot su tutte conte DGE dopo normalizzazione
#mds=plotMDS(k27ac_y_new,dim.plot=c(2,3),gene.selection="common")
#mds2=data.frame(x=mds$x,y=mds$y)
#samples=rownames(mds2)
#infos = targets[row.names(k27ac_y_new$samples),]
#mds2=cbind(mds2,as.data.frame(samples),as.data.frame(infos$disease),as.data.frame(infos$passage))
#colnames(mds2)=c("MD2","MD3","samples","condition","passage")
##pdf("PCA_AllSamples_DGE_RNA_Common.pdf")
#ggplot(mds2, aes(MD3, MD2)) + geom_point(aes(colour = condition,shape=as.character(passage)))+ geom_text(aes(label=samples),hjust=0, vjust=0,size=2) + scale_x_continuous(limits = c(-4, 4)) 
#
#mds=plotMDS(k27ac_y_new,dim.plot=c(1,3),gene.selection="common")
#mds2=data.frame(x=mds$x,y=mds$y)
#samples=rownames(mds2)
#infos = targets[row.names(k27ac_y_new$samples),]
#mds2=cbind(mds2,as.data.frame(samples),as.data.frame(infos$disease),as.data.frame(infos$passage))
#colnames(mds2)=c("MD1","MD3","samples","condition","passage")
##pdf("PCA_AllSamples_DGE_RNA_Common.pdf")
#ggplot(mds2, aes(MD1, MD3)) + geom_point(aes(colour = condition,shape=as.character(passage)))+ geom_text(aes(label=samples),hjust=0, vjust=0,size=2)  + scale_x_continuous(limits = c(-4, 4)) 
#
#dev.off()

toptags_RNA


RNA_DEGs=toptags_RNA[toptags_RNA$Qvalue<=0.01,]


RNAseq_genes_down_prom<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=ATAC_down_prom$geneId,mart=mart)
ATACseq_genes_up_prom<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=ATAC_up_prom$geneId,mart=mart)
RNA_deg_genes<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=RNA_DEGs$genes,mart=mart)
RNAseq_all_genes=getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=toptags_RNA$genes,mart=mart)
atacDown_ego_prom <- enrichGO(gene          = ATACseq_genes_down_prom$entrezgene,
universe      = as.character(ATACseq_all_genes$entrezgene),
OrgDb         = org.Hs.eg.db,
ont           = "BP",
pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = T)
dotplot(atacDown_ego_prom)

atacUp_ego_prom <- enrichGO(gene          = ATACseq_genes_up_prom$entrezgene,
universe      = as.character(ATACseq_all_genes$entrezgene),
OrgDb         = org.Hs.eg.db,
ont           = "BP",
pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = T)
dotplot(atacUp_ego_prom)

rnaDeg_ego<- enrichGO(gene          = RNA_deg_genes$entrezgene,
universe      = as.character(RNAseq_all_genes$entrezgene),
OrgDb         = org.Hs.eg.db,
ont           = "BP",
pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = T)
dotplot(rnaDeg_ego)
rnaDeg_ego_df=rnaDeg_ego

#***********************************************************METILAZIONE*************************************************************************************


library(methylKit)
rrbs_bams='/lustre1/workspace/DiMicco/prj_ALL_pilot/RRBS/BAM/'
files=list(list.files(path=rrbs_bams,pattern='bam$'))
files.ls= Filter(function(x) grepl("UPN3", x), files[[1]])
files.ls=paste(rrbs_bams,files.ls[c(1,7,8,9,13,14,2,3)],sep='')

objs=processBismarkAln(location=files.ls[1:8]              ,sample.id=list('D_1','D_2'),assembly="hg38",             save.folder=NULL,save.context=NULL,read.context="CpG",             nolap=FALSE,mincov=10,minqual=20,phred64=FALSE,             treatment=c(0,1))
             
             


objs=processBismarkAln(location=files.ls[1:8] ,sample.id=list('D_0','D_1','D_2','D_3','R_0','R_1','R_2','R_3'),assembly="hg38", save.folder=getwd(),save.context= "CpG",read.context="CpG", nolap=FALSE,mincov=10,minqual=20,phred64=FALSE, treatment=c(0,0,0,0,1,1,1,1))

cat D_0_CpG.txt D_1_CpG.txt D_2_CpG.txt D_3_CpG.txt R_0_CpG.txt R_1_CpG.txt R_2_CpG.txt R_3_CpG.txt | tail -n +2 | awk -v OFS='\t' '{print $2,$3-1,$3}' > Tutte_Regioni_Metilazione.txt
bedSort Tutte_Regioni_Metilazione.txt Tutte_Regioni_Metilazione.txt

 for i in  *CpG.txt; do prefix=$(echo $i | sed 's/.txt/_for_R.txt/g'); tail -n +2 $i | awk -v OFS='\t' '{print $2,$3-1,$3,$6}' - >   $prefix; done

D0_met=read.table('D_0_CpG_for_R.txt',head=F,col.names = c('pos','start','end','Met_C'))
D1_met=read.table('D_1_CpG_for_R.txt',head=F,col.names = c('pos','start','end','Met_C'))
D2_met=read.table('D_2_CpG_for_R.txt',head=F,col.names = c('pos','start','end','Met_C'))
D3_met=read.table('D_3_CpG_for_R.txt',head=F,col.names = c('pos','start','end','Met_C'))
R0_met=read.table('R_0_CpG_for_R.txt',head=F,col.names = c('pos','start','end','Met_C'))
R1_met=read.table('R_1_CpG_for_R.txt',head=F,col.names = c('pos','start','end','Met_C'))
R2_met=read.table('R_2_CpG_for_R.txt',head=F,col.names = c('pos','start','end','Met_C'))
R3_met=read.table('R_3_CpG_for_R.txt',head=F,col.names = c('pos','start','end','Met_C'))

All_Regions_met=read.table('Tutte_Regioni_Metilazione.txt',head=F,col.names = c('pos','start','end'))


Mets=list(All_Regions_met,D0_met,D1_met,D2_met,D3_met,R0_met,R1_met,R2_met,R3_met)

Metilazione=join_all(Mets,by=c('pos','start','end'))
colnames(Metilazione_filt)=c('chr','start','end','UPN3_1','UPN3_2','UPN3_3','UPN3_4','UPN3_8','UPN3_9','UPN3_10','UPN3_11')
Metilazione=Metilazione[Metilazione$start!=-1,]
saveRDS(Metilazione,'Tabella_Metilazione.RDS')
Metilazione=
Metilazione[,4:11][is.na(Metilazione[,4:11])]<-0
write.table(Metilazione,'Metilazione_Campioni.txt',sep='\t',col.names=T,row.names=F,quote=F)
Metilazione_peaks=readPeakFile("Metilazione_Campioni.txt",head=T)
Metilazione_anno <- annotatePeak(peak = Metilazione_peaks, tssRegion=c(-200, 200),TxDb=gencode_rnaseqtxdb) #,annoDb="org.Hs.eg.db")
Metilazione_anno_df=as.data.frame(Metilazione_anno)
#Metilazione_anno_df=Metilazione_anno_df[,c(1,2,3,6:13,20)]
rows=rowSums(Metilazione_anno_df[,6:13]>20)>2
Metilazione_anno_df=Metilazione_anno_df[rows,]
Metilazione_anno_df=Metilazione_anno_df[!duplicated(Metilazione_anno_df),]
rownames(Metilazione_anno_df)=make.names(Metilazione_anno_df$geneId,unique = T)
saveRDS(Metilazione_anno_df,'Metilazione_anno_df.RDS')
meth=Metilazione_anno_df[,6:13]
meth=meth/100
idx <- which(apply(meth, 1, function(x) all(x == 1)) == T)
Meth=meth[-idx,]

rownames(atac_y_new_cpm)=make.names(atac_y_new_cpm$geneId,unique = T)
rownames(k27ac_y_new_cpm)=make.names(k27ac_y_new_cpm$geneId,unique = T)
rownames(k4me3_y_new_cpm)=make.names(k4me3_y_new_cpm$geneId,unique = T)

k4me3_y_new_cpm$x <- do.call(paste, c(k4me3_y_new_cpm[1:3], sep="-"))
rownames(k4me3_y_new_cpm)=make.names(k4me3_y_new_cpm$x)
k27ac_y_new_cpm$x <- do.call(paste, c(k27ac_y_new_cpm[1:3], sep="-"))
rownames(k27ac_y_new_cpm)=make.names(k27ac_y_new_cpm$x)
atac_y_new_cpm$x <- do.call(paste, c(atac_y_new_cpm[1:3], sep="-"))
rownames(atac_y_new_cpm)=make.names(atac_y_new_cpm$x)


ll=list()
ll$acetil=k27ac_y_new_cpm [,8:15]
ll$atac=atac_y_new_cpm[,8:15]
ll$rna=rna_cpm
ll$trimetil=k4me3_y_new_cpm[,7:14]
ll$meth=methyl[,4:11]
targets=read.table('Targets_for_DEAnalysis.txt')



Saving model in /var/folders/7z/12qs_17x17gcgy1qs3lxc4srzw4_z5/T//RtmpX3Oxct/file698e6beccbb...

[1] "Error loading the training data..."
Error in `factorNames<-`(`*tmp*`, value = c("intercept", "LF1", "LF0")) : 
  Length of factor names does not match the dimensionality of the latent variable matrix
  
  
  
pilot_targets=targets[c(1,2,3,4,8,9,10,11),]
maepilot <- MultiAssayExperiment(experiments = ll, colData = pilot_targets )


MOFAobject_meth <- createMOFAobject(ll)
plotTilesData(MOFAobject_meth)
DataOptions <- getDefaultDataOptions()
DataOptions
ModelOptions <- getDefaultModelOptions(MOFAobject_meth)
TrainOptions <- getDefaultTrainOptions()
TrainOptions$DropFactorThreshold <- 0.02 
TrainOptions


MOFAobject_meth <- prepareMOFA(
MOFAobject_meth,
DataOptions = DataOptions,
ModelOptions = ModelOptions,
TrainOptions = TrainOptions)
MOFAobject_meth <- runMOFA(MOFAobject_meth, outfile=tempfile())

```
Setting random seed 300581...
```


MOFAobject
```
Trained MOFA model with the following characteristics: 
 Number of views: 4 
 View names: k27 rna k4 atac 
 Number of features per view: 39080 13711 38168 57711 
 Number of samples: 8 
 Number of factors: 4 
```


r2 <- calculateVarianceExplained(MOFAobject)
r2
```
r2
$R2Total
      k27       rna        k4      atac 
0.5356315 0.6254190 0.9115040 0.5915935 

$R2PerFactor
           k27         rna          k4       atac
LF1 0.11632767 0.256672143 0.600677856 0.15744056
LF2 0.20569116 0.292741859 0.109059708 0.32646445
LF3 0.16488193 0.072091931 0.068515642 0.05446214
LF4 0.01328265 0.000293540 0.155112306 0.01438974
LF5 0.04578521 0.003170524 0.003246636 0.03910484
```

RNA_Factors1_2_nometh<- getWeights(MOFAobject,      views = "rna", factors = c(1,2) ,as.data.frame = TRUE)
RNA_Factors1_2_nometh_LFT1=RNA_Factors1_2_nometh[RNA_Factors1_2_nometh$factor=='LF1',]
RNA_Factors1_2_nometh_LFT2=RNA_Factors1_2_nometh[RNA_Factors1_2_nometh$factor=='LF2',]
RNA_Factors1_2_nometh_Factors=cbind(RNA_Factors1_2_nometh_LFT1,RNA_Factors1_2_nometh_LFT2)
RNA_Factors1_2_nometh_Factors=RNA_Factors1_2_nometh_Factors[,c(1,3,7)]

colnames(RNA_Factors1_2_nometh_Factors)=c('Genes','LF1','LF2')
RNA_Factors1_2_nometh_Factors_LF1=RNA_Factors1_2_nometh_Factors[,1:2]
RNA_Factors1_2_nometh_Factors_LF2=RNA_Factors1_2_nometh_Factors[,c(1,3)]

RNA_Factors3_nometh<- getWeights(MOFAobject,      views = "rna", factors = 3 ,as.data.frame = TRUE)


write.table(RNA_Factors1_2_nometh_Factors_LF1[order(-RNA_Factors1_2_nometh_Factors_LF1$LF1),],'MOFA_RNA_Factors_NoMeth_Values_LF1.rnk',sep='\t',col.names = F,row.names=F,quote=F)
write.table(RNA_Factors1_2_nometh_Factors_LF2[order(-RNA_Factors1_2_nometh_Factors_LF2$LF2),],'MOFA_RNA_Factors_NoMeth_Values_LF2.rnk',sep='\t',col.names = F,row.names=F,quote=F)
write.table(RNA_Factors3_nometh[order(-RNA_Factors3_nometh$factor),c(1,3)],'MOFA_RNA_Factors_NoMeth_Values_LF3.rnk',sep='\t',col.names = F,row.names=F,quote=F))

K4_Factors1_2_nometh <- getWeights(MOFAobject,      views = "k27", factors = c(1,2),  as.data.frame = TRUE)
K4_Factors1_2_nometh_LFT1=K4_Factors1_2_nometh[K4_Factors1_2_nometh$factor=='LF1',]
K4_Factors1_2_nometh_LFT2=K4_Factors1_2_nometh[K4_Factors1_2_nometh$factor=='LF2',]
K4_Factors1_2_nometh_Factors=cbind(K4_Factors1_2_nometh_LFT1,K4_Factors1_2_nometh_LFT2)
K4_Factors1_2_nometh_Factors=K4_Factors1_2_nometh_Factors[,c(1,3,7)]


K4_Factors1_2_nometh_Factors_meanLF1=K4_Factors1_2_nometh_Factors[,c(2,4)]
K4_Factors1_2_nometh_Factors_meanLF1=aggregate(K4_Factors1_2_nometh_Factors_meanLF1[,1], list(Genes=K4_Factors1_2_nometh_Factors_meanLF1[,2]), FUN=mean)
K4_Factors1_2_nometh_Factors_meanLF2=K4_Factors1_2_nometh_Factors[,c(3,4)]
K4_Factors1_2_nometh_Factors_meanLF2=aggregate(K4_Factors1_2_nometh_Factors_meanLF2[,1], list(Genes=K4_Factors1_2_nometh_Factors_meanLF2[,2]), FUN=mean)

colnames(K4_Factors1_2_nometh_Factors_meanLF1)=c('Genes','Mean_LF1')
colnames(K4_Factors1_2_nometh_Factors_meanLF2)=c('Genes','Mean_LF2')

K4_Factors1_2_nometh_Factors_meanValues=merge(K4_Factors1_2_nometh_Factors_meanLF1,K4_Factors1_2_nometh_Factors_meanLF2,by='Genes')

write.table(K4_Factors1_2_nometh_Factors_meanLF1[order(-K4_Factors1_2_nometh_Factors_meanLF1$Mean_LF1),],'MOFA_K4_Factors_NoMeth_MeanValues_LF1.rnk',sep='\t',col.names = F,row.names=F,quote=F)
write.table(K4_Factors1_2_nometh_Factors_meanLF2[order(-K4_Factors1_2_nometh_Factors_meanLF2$Mean_LF2),],'MOFA_K4_Factors_NoMeth_MeanValues_LF2.rnk',sep='\t',col.names = F,row.names=F,quote=F)


K27_Factors1_2_nometh <- getWeights(MOFAobject,      views = "k4", factors = c(1,2),  as.data.frame = TRUE)
K27_Factors1_2_nometh_LFT1=K27_Factors1_2_nometh[K27_Factors1_2_nometh$factor=='LF1',]
K27_Factors1_2_nometh_LFT2=K27_Factors1_2_nometh[K27_Factors1_2_nometh$factor=='LF2',]
K27_Factors1_2_nometh_Factors=cbind(K27_Factors1_2_nometh_LFT1,K27_Factors1_2_nometh_LFT2)
K27_Factors1_2_nometh_Factors=K27_Factors1_2_nometh_Factors[,c(1,3,7)]

K27_Factors1_2_nometh_Factors_meanLF1=K27_Factors1_2_nometh_Factors[,c(2,4)]
K27_Factors1_2_nometh_Factors_meanLF1=aggregate(K27_Factors1_2_nometh_Factors_meanLF1[,1], list(Genes=K27_Factors1_2_nometh_Factors_meanLF1[,2]), FUN=mean)
K27_Factors1_2_nometh_Factors_meanLF2=K27_Factors1_2_nometh_Factors[,c(3,4)]
K27_Factors1_2_nometh_Factors_meanLF2=aggregate(K27_Factors1_2_nometh_Factors_meanLF2[,1], list(Genes=K27_Factors1_2_nometh_Factors_meanLF2[,2]), FUN=mean)

colnames(K27_Factors1_2_nometh_Factors_meanLF1)=c('Genes','Mean_LF1')
colnames(K27_Factors1_2_nometh_Factors_meanLF2)=c('Genes','Mean_LF2')

K27_Factors1_2_nometh_Factors_meanValues=merge(K27_Factors1_2_nometh_Factors_meanLF1,K27_Factors1_2_nometh_Factors_meanLF2,by='Genes')

write.table(K27_Factors1_2_nometh_Factors_meanLF1[order(-K27_Factors1_2_nometh_Factors_meanLF1$Mean_LF1),],'MOFA_K27_Factors_NoMeth_MeanValues_LF1.rnk',sep='\t',col.names = F,row.names=F,quote=F)
write.table(K27_Factors1_2_nometh_Factors_meanLF2[order(-K27_Factors1_2_nometh_Factors_meanLF2$Mean_LF2),],'MOFA_K27_Factors_NoMeth_MeanValues_LF2.rnk',sep='\t',col.names = F,row.names=F,quote=F)




ATAC_Factors1_2_nometh <- getWeights(MOFAobject,      views = "atac", factors = c(1,2),  as.data.frame = TRUE)
ATAC_Factors1_2_nometh_LFT1=ATAC_Factors1_2_nometh[ATAC_Factors1_2_nometh$factor=='LF1',]
ATAC_Factors1_2_nometh_LFT2=ATAC_Factors1_2_nometh[ATAC_Factors1_2_nometh$factor=='LF2',]
ATAC_Factors1_2_nometh_Factors=cbind(ATAC_Factors1_2_nometh_LFT1,ATAC_Factors1_2_nometh_LFT2)
ATAC_Factors1_2_nometh_Factors=ATAC_Factors1_2_nometh_Factors[,c(1,3,7)]

ATAC_Factors1_2_nometh_Factors_meanLF1=ATAC_Factors1_2_nometh_Factors[,c(2,4)]
ATAC_Factors1_2_nometh_Factors_meanLF1=aggregate(ATAC_Factors1_2_nometh_Factors_meanLF1[,1], list(Genes=ATAC_Factors1_2_nometh_Factors_meanLF1[,2]), FUN=mean)
ATAC_Factors1_2_nometh_Factors_meanLF2=ATAC_Factors1_2_nometh_Factors[,c(3,4)]
ATAC_Factors1_2_nometh_Factors_meanLF2=aggregate(ATAC_Factors1_2_nometh_Factors_meanLF2[,1], list(Genes=ATAC_Factors1_2_nometh_Factors_meanLF2[,2]), FUN=mean)

colnames(ATAC_Factors1_2_nometh_Factors_meanLF1)=c('Genes','Mean_LF1')
colnames(ATAC_Factors1_2_nometh_Factors_meanLF2)=c('Genes','Mean_LF2')

ATAC_Factors1_2_nometh_Factors_meanValues=merge(ATAC_Factors1_2_nometh_Factors_meanLF1,ATAC_Factors1_2_nometh_Factors_meanLF2,by='Genes')

write.table(ATAC_Factors1_2_nometh_Factors_meanLF1[order(-ATAC_Factors1_2_nometh_Factors_meanLF1$Mean_LF1),],'MOFA_ATAC_Factors_NoMeth_MeanValues_LF1.rnk',sep='\t',col.names = F,row.names=F,quote=F)
write.table(ATAC_Factors1_2_nometh_Factors_meanLF2[order(-ATAC_Factors1_2_nometh_Factors_meanLF2$Mean_LF2),],'MOFA_ATAC_Factors_NoMeth_MeanValues_LF2.rnk',sep='\t',col.names = F,row.names=F,quote=F)



tabellone=read.table('Tabella_GSEA_MOFA_AllTechniques_LF1LF2.txt',sep='\t',head=T)
tabellone_nes=tabellone[,c(1,2,3,6,8,10,12,14,16)]
tabellone_fdr=tabellone[,c(1,4,5,7,9,11,13,15,17)]
tabellone_nes=melt(tabellone_nes)
tabellone_fdr=melt(tabellone_fdr)
tabellone_heatmap=cbind(tabellone_nes,tabellone_fdr)
tabellone_heatmap=tabellone_heatmap[,c(1,2,3,6)]
colnames(tabellone_heatmap)=c('Hallmark','Technique_LFs', 'NES','FDR')

p=ggplot(tabellone_heatmap, aes(Technique_LFs,Hallmark)) +     geom_point(aes(size = -log10(FDR), colour=NES))   +  theme(panel.background = element_rect(fill = "white",colour = "black",size = 0.5, linetype = "solid"),panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "black"),panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour ="black")) + theme(legend.position="right") +theme_grey(base_size = 10,) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_gradient2(low = "midnightblue", mid = "white",high = "orange4", midpoint = 0, space = "Lab",na.value = "grey50", guide = "colourbar", aesthetics = "colour") + coord_flip()
pdf('Plot_MOFA_FActors.pdf')
print(p, vp=viewport(angle=-90))
dev.off()

Tabella_GSEA_MOFA_RNA_LF1LF2LF3=read.table('Tabella_GSEA_MOFA_RNA_LF1LF2LF3.txt',head=T,sep='\t')
Tabella_GSEA_MOFA_RNA_LF1LF2LF3=Tabella_GSEA_MOFA_RNA_LF1LF2LF3[,1:7]
Tabella_GSEA_MOFA_RNA_LF1LF2LF3_NES=Tabella_GSEA_MOFA_RNA_LF1LF2LF3[,c(1,2,4,6)]
Tabella_GSEA_MOFA_RNA_LF1LF2LF3_FDR=Tabella_GSEA_MOFA_RNA_LF1LF2LF3[,c(1,3,5,7)]
Tabella_GSEA_MOFA_RNA_LF1LF2LF3_NES=melt(Tabella_GSEA_MOFA_RNA_LF1LF2LF3_NES)
Tabella_GSEA_MOFA_RNA_LF1LF2LF3_FDR=melt(Tabella_GSEA_MOFA_RNA_LF1LF2LF3_FDR)
Tabella_GSEA_MOFA_RNA_LF1LF2LF3_heatmap=cbind(Tabella_GSEA_MOFA_RNA_LF1LF2LF3_NES,Tabella_GSEA_MOFA_RNA_LF1LF2LF3_FDR)
Tabella_GSEA_MOFA_RNA_LF1LF2LF3_heatmap=Tabella_GSEA_MOFA_RNA_LF1LF2LF3_heatmap[,c(1,2,3,6)]
colnames(Tabella_GSEA_MOFA_RNA_LF1LF2LF3_heatmap)=c('Hallmark','Technique_LFs', 'NES','FDR') 

p=ggplot(Tabella_GSEA_MOFA_RNA_LF1LF2LF3_heatmap, aes(Technique_LFs,Hallmark)) +     geom_point(aes(size = -log10(FDR), colour=NES))   +  theme(panel.background = element_rect(fill = "white",colour = "black",size = 0.5, linetype = "solid"),panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "black"),panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour ="black")) + theme(legend.position="right") +theme_grey(base_size = 10,) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_gradient2(low = "midnightblue", mid = "white",high = "orange4", midpoint = 0, space = "Lab",na.value = "grey50", guide = "colourbar", aesthetics = "colour") + coord_flip()
pdf('Plot_MOFA_RNAseq_FActors123.pdf')
print(p, vp=viewport(angle=-90))
dev.off()




ATAC_Factors1_2_nometh_LFT1_moredetails=cbind(atac_y_new_cpm[1:5],ATAC_Factors1_2_nometh_LFT1)
colnames(ATAC_Factors1_2_nometh_LFT1_moredetails)=c("ATAC_chr"     ,"ATAC_start" ,  "ATAC_end","width",    "geneId",   "genes",    "factor",  "value"    ,"view")
ATAC_LatentFactors1_TabelloneAnnotazione=merge(ATAC_Factors1_2_nometh_LFT1_moredetails,Def_towrite_cols,by=c("ATAC_chr"     ,"ATAC_start" ,  "ATAC_end"))
ATAC_LatentFactors1_TabelloneAnnotazione=ATAC_LatentFactors1_TabelloneAnnotazione[c(1,2,3,4,5,7,8,9,10,11,12,13)]
ATAC_LatentFactors1_TabelloneAnnotazione=unique(ATAC_LatentFactors1_TabelloneAnnotazione)
ATAC_LatentFactors1_TabelloneAnnotazione$color=ifelse(ATAC_LatentFactors1_TabelloneAnnotazione$K27acOverlap=='+'&ATAC_LatentFactors1_TabelloneAnnotazione$K4me3Overlap=='-','red',ifelse(ATAC_LatentFactors1_TabelloneAnnotazione$K27acOverlap=='+'&ATAC_LatentFactors1_TabelloneAnnotazione$K4me3Overlap=='+','green',ifelse(ATAC_LatentFactors1_TabelloneAnnotazione$K27acOverlap=='-'&ATAC_LatentFactors1_TabelloneAnnotazione$K4me3Overlap=='+','green','black')))
ATAC_LatentFactors1_TabelloneAnnotazione=ATAC_LatentFactors1_TabelloneAnnotazione[order(-ATAC_LatentFactors1_TabelloneAnnotazione$value),]

ATAC_LatentFactors1_TabelloneAnnotazione_sign=ATAC_LatentFactors1_TabelloneAnnotazione[ATAC_LatentFactors1_TabelloneAnnotazione$ATAC_Qvalue<=0.05,]
ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCposBlack=ATAC_LatentFactors1_TabelloneAnnotazione_sign[ATAC_LatentFactors1_TabelloneAnnotazione_sign$color=='black'&ATAC_LatentFactors1_TabelloneAnnotazione_sign$value>0&ATAC_LatentFactors1_TabelloneAnnotazione_sign$ATAC_logFC>0,]
nrow(ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCposBlack)
ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCposBlack=ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCposBlack[1:100,]
ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCposBlack=ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCposBlack[,1:3]
write.table(ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCposBlack,"../../Desktop/Home3/ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCposBlack.bed",col.names = F,row.names = F,quote=F,sep='\t')

ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCnegBlack=ATAC_LatentFactors1_TabelloneAnnotazione_sign[ATAC_LatentFactors1_TabelloneAnnotazione_sign$color=='black'&ATAC_LatentFactors1_TabelloneAnnotazione_sign$value>0&ATAC_LatentFactors1_TabelloneAnnotazione_sign$ATAC_logFC<0,]
nrow(ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCnegBlack)
ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCnegBlack=ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCnegBlack[1:100,]
ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCnegBlack=ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCnegBlack[,1:3]
write.table(ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCnegBlack,"../../Desktop/Home3/ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCnegBlack.bed",col.names = F,row.names = F,quote=F,sep='\t')

ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCposBlack=ATAC_LatentFactors1_TabelloneAnnotazione_sign[ATAC_LatentFactors1_TabelloneAnnotazione_sign$color=='black'&ATAC_LatentFactors1_TabelloneAnnotazione_sign$value<0&ATAC_LatentFactors1_TabelloneAnnotazione_sign$ATAC_logFC>0,]
nrow(ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCposBlack)
ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCposBlack=ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCposBlack[1:100,]
ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCposBlack=ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCposBlack[,1:3]
write.table(ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCposBlack,"../../Desktop/Home3/ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCposBlack.bed",col.names = F,row.names = F,quote=F,sep='\t')

ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCnegBlack=ATAC_LatentFactors1_TabelloneAnnotazione_sign[ATAC_LatentFactors1_TabelloneAnnotazione_sign$color=='black'&ATAC_LatentFactors1_TabelloneAnnotazione_sign$value<0&ATAC_LatentFactors1_TabelloneAnnotazione_sign$ATAC_logFC<0,]
nrow(ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCnegBlack)
ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCnegBlack=ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCnegBlack[1:100,]
ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCnegBlack=ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCnegBlack[,1:3]
write.table(ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCnegBlack,"../../Desktop/Home3/ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCnegBlack.bed",col.names = F,row.names = F,quote=F,sep='\t')


ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCposAcetil=ATAC_LatentFactors1_TabelloneAnnotazione_sign[ATAC_LatentFactors1_TabelloneAnnotazione_sign$color=='red'&ATAC_LatentFactors1_TabelloneAnnotazione_sign$value>0&ATAC_LatentFactors1_TabelloneAnnotazione_sign$ATAC_logFC>0,]
nrow(ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCposAcetil)
ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCposAcetil=ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCposAcetil[1:100,]
ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCposAcetil=ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCposAcetil[,1:3]
write.table(ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCposAcetil,"../../Desktop/Home3/ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCposAcetil.bed",col.names = F,row.names = F,quote=F,sep='\t')

ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCnegAcetil=ATAC_LatentFactors1_TabelloneAnnotazione_sign[ATAC_LatentFactors1_TabelloneAnnotazione_sign$color=='red'&ATAC_LatentFactors1_TabelloneAnnotazione_sign$value>0&ATAC_LatentFactors1_TabelloneAnnotazione_sign$ATAC_logFC<0,]
nrow(ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCnegAcetil)
ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCnegAcetil=ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCnegAcetil[1:100,]
ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCnegAcetil=ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCnegAcetil[,1:3]
write.table(ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCnegAcetil,"../../Desktop/Home3/ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCnegAcetil.bed",col.names = F,row.names = F,quote=F,sep='\t')

ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCposAcetil=ATAC_LatentFactors1_TabelloneAnnotazione_sign[ATAC_LatentFactors1_TabelloneAnnotazione_sign$color=='red'&ATAC_LatentFactors1_TabelloneAnnotazione_sign$value<0&ATAC_LatentFactors1_TabelloneAnnotazione_sign$ATAC_logFC>0,]
nrow(ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCposAcetil)
ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCposAcetil=ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCposAcetil[1:100,]
ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCposAcetil=ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCposAcetil[,1:3]
write.table(ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCposAcetil,"../../Desktop/Home3/ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCposAcetil.bed",col.names = F,row.names = F,quote=F,sep='\t')

ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCnegAcetil=ATAC_LatentFactors1_TabelloneAnnotazione_sign[ATAC_LatentFactors1_TabelloneAnnotazione_sign$color=='red'&ATAC_LatentFactors1_TabelloneAnnotazione_sign$value<0&ATAC_LatentFactors1_TabelloneAnnotazione_sign$ATAC_logFC<0,]
nrow(ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCnegAcetil)
ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCnegAcetil=ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCnegAcetil[1:100,]
ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCnegAcetil=ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCnegAcetil[,1:3]
write.table(ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCnegAcetil,"../../Desktop/Home3/ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCnegAcetil.bed",col.names = F,row.names = F,quote=F,sep='\t')


ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCposTrimetil=ATAC_LatentFactors1_TabelloneAnnotazione_sign[ATAC_LatentFactors1_TabelloneAnnotazione_sign$color=='green'&ATAC_LatentFactors1_TabelloneAnnotazione_sign$value>0&ATAC_LatentFactors1_TabelloneAnnotazione_sign$ATAC_logFC>0,]
nrow(ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCposTrimetil)
ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCposTrimetil=ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCposTrimetil[,1:3]
write.table(ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCposTrimetil,"../../Desktop/Home3/ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCposTrimetil.bed",col.names = F,row.names = F,quote=F,sep='\t')

ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCnegTrimetil=ATAC_LatentFactors1_TabelloneAnnotazione_sign[ATAC_LatentFactors1_TabelloneAnnotazione_sign$color=='green'&ATAC_LatentFactors1_TabelloneAnnotazione_sign$value>0&ATAC_LatentFactors1_TabelloneAnnotazione_sign$ATAC_logFC<0,]
nrow(ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCnegTrimetil)
ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCnegTrimetil=ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCnegTrimetil[1:100,]
ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCnegTrimetil=ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCnegTrimetil[,1:3]
write.table(ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCnegTrimetil,"../../Desktop/Home3/ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFposFCnegTrimetil.bed",col.names = F,row.names = F,quote=F,sep='\t')

ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCposTrimetil=ATAC_LatentFactors1_TabelloneAnnotazione_sign[ATAC_LatentFactors1_TabelloneAnnotazione_sign$color=='green'&ATAC_LatentFactors1_TabelloneAnnotazione_sign$value<0&ATAC_LatentFactors1_TabelloneAnnotazione_sign$ATAC_logFC>0,]
nrow(ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCposTrimetil)
ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCposTrimetil=ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCposTrimetil[1:100,]
ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCposTrimetil=ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCposTrimetil[,1:3]
write.table(ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCposTrimetil,"../../Desktop/Home3/ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCposTrimetil.bed",col.names = F,row.names = F,quote=F,sep='\t')

ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCnegTrimetil=ATAC_LatentFactors1_TabelloneAnnotazione_sign[ATAC_LatentFactors1_TabelloneAnnotazione_sign$color=='green'&ATAC_LatentFactors1_TabelloneAnnotazione_sign$value<0&ATAC_LatentFactors1_TabelloneAnnotazione_sign$ATAC_logFC<0,]
nrow(ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCnegTrimetil)
ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCnegTrimetil=ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCnegTrimetil[,1:3]
write.table(ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCnegTrimetil,"../../Desktop/Home3/ATAC_LatentFactors1_TabelloneAnnotazione_sign_LFnegFCnegTrimetil.bed",col.names = F,row.names = F,quote=F,sep='\t')


#
#
#
##Rifaccio senza K4
#
#MOFAobject_NoK4 <- createMOFAobject(ll)
#plotTilesData(MOFAobject_NoK4)
#DataOptions <- getDefaultDataOptions()
#DataOptions
#ModelOptions <- getDefaultModelOptions(MOFAobject_NoK4)
#TrainOptions <- getDefaultTrainOptions()
#TrainOptions$DropFactorThreshold <- 0.02
#TrainOptions
#
#
#MOFAobject_NoK4 <- prepareMOFA(
#MOFAobject_NoK4,
#DataOptions = DataOptions,
#ModelOptions = ModelOptions,
#TrainOptions = TrainOptions)
#MOFAobject_NoK4 <- runMOFA(MOFAobject_NoK4, outfile=tempfile())
#
#r2 <- calculateVarianceExplained(MOFAobject_NoK4)
#r2
#
#```
#$R2Total
#      k27       rna        k4      atac 
#0.5088427 0.6252040 0.9325992 0.5686192 
#
#$R2PerFactor
#           k27          rna         k4       atac
#LF1 0.11150386 0.2547107204 0.60046364 0.15873568
#LF2 0.20638771 0.2927571472 0.11066199 0.32497330
#LF3 0.16718935 0.0708103624 0.06975823 0.05307114
#LF4 0.00941138 0.0007203005 0.16207743 0.01311091
#LF5 0.02098114 0.0058172822 0.01986787 0.01859151
#```
#
#
#
#RNA_Factors1_2_nometh_NoK4<- getWeights(MOFAobject_NoK4,      views = "rna", factors = c(1,2) ,as.data.frame = TRUE)
#RNA_Factors1_2_nometh_NoK4_LFT1=RNA_Factors1_2_nometh_NoK4[RNA_Factors1_2_nometh_NoK4$factor=='LF1',]
#RNA_Factors1_2_nometh_NoK4_LFT2=RNA_Factors1_2_nometh_NoK4[RNA_Factors1_2_nometh_NoK4$factor=='LF2',]
#RNA_Factors1_2_nometh_NoK4_Factors=cbind(RNA_Factors1_2_nometh_NoK4_LFT1,RNA_Factors1_2_nometh_NoK4_LFT2)
#RNA_Factors1_2_nometh_NoK4_Factors=RNA_Factors1_2_nometh_NoK4_Factors[,c(1,3,7)]
#
#colnames(RNA_Factors1_2_nometh_NoK4_Factors)=c('Genes','LF1','LF2')
#RNA_Factors1_2_nometh_NoK4_Factors_LF1=RNA_Factors1_2_nometh_NoK4_Factors[,1:2]
#RNA_Factors1_2_nometh_NoK4_Factors_LF2=RNA_Factors1_2_nometh_NoK4_Factors[,c(1,3)]
#
#
#write.table(RNA_Factors1_2_nometh_NoK4_Factors_LF1[order(-RNA_Factors1_2_nometh_NoK4_Factors_LF1$LF1),],'MOFA_RNA_Factors_NoMeth_NoK4_Values_LF1.rnk',sep='\t',col.names = F,row.names=F,quote=F)
#write.table(RNA_Factors1_2_nometh_NoK4_Factors_LF2[order(-RNA_Factors1_2_nometh_NoK4_Factors_LF2$LF2),],'MOFA_RNA_Factors_NoMeth_NoK4_Values_LF2.rnk',sep='\t',col.names = F,row.names=F,quote=F)
#
#
#rownames(atac_y_new_cpm)=make.names(atac_y_new_cpm$geneId,unique = T)
#rownames(k27ac_y_new_cpm)=make.names(k27ac_y_new_cpm$geneId,unique = T)
#rownames(k4me3_y_new_cpm)=make.names(k4me3_y_new_cpm$geneId,unique = T)
#
#
#K27_Factors1_2_nometh_NoK4<- getWeights(MOFAobject_NoK4,      views = "acetil", factors = c(1,2) ,as.data.frame = TRUE)
#K27_Factors1_2_nometh_NoK4$feature=k27ac_y_new_cpm$geneId
#K27_Factors1_2_nometh_NoK4_LFT1=K27_Factors1_2_nometh_NoK4[K27_Factors1_2_nometh_NoK4$factor=='LF1',]
#K27_Factors1_2_nometh_NoK4_LFT2=K27_Factors1_2_nometh_NoK4[K27_Factors1_2_nometh_NoK4$factor=='LF2',]
#K27_Factors1_2_nometh_NoK4_Factors=cbind(K27_Factors1_2_nometh_NoK4_LFT1,K27_Factors1_2_nometh_NoK4_LFT2)
#K27_Factors1_2_nometh_NoK4_Factors=K27_Factors1_2_nometh_NoK4_Factors[,c(1,3,7)]
#
#
#K27_Factors1_2_nometh_NoK4_Factors_meanLF1=K27_Factors1_2_nometh_NoK4_Factors[,c(1,2)]
#K27_Factors1_2_nometh_NoK4_Factors_meanLF1=aggregate(K27_Factors1_2_nometh_NoK4_Factors_meanLF1[,2], list(Genes=K27_Factors1_2_nometh_NoK4_Factors_meanLF1[,1]), FUN=mean)
#K27_Factors1_2_nometh_NoK4_Factors_meanLF2=K27_Factors1_2_nometh_NoK4_Factors[,c(1,3)]
#K27_Factors1_2_nometh_NoK4_Factors_meanLF2=aggregate(K27_Factors1_2_nometh_NoK4_Factors_meanLF2[,2], list(Genes=K27_Factors1_2_nometh_NoK4_Factors_meanLF2[,1]), FUN=mean)
#
#
#
#
#write.table(K27_Factors1_2_nometh_NoK4_Factors_meanLF1[order(-K27_Factors1_2_nometh_NoK4_Factors_meanLF1$x),],'MOFA_K27_Factors_NoMeth_NoK4_Values_LF1.rnk',sep='\t',col.names = F,row.names=F,quote=F)
#write.table(K27_Factors1_2_nometh_NoK4_Factors_meanLF2[order(-K27_Factors1_2_nometh_NoK4_Factors_meanLF2$x),],'MOFA_K27_Factors_NoMeth_NoK4_Values_LF2.rnk',sep='\t',col.names = F,row.names=F,quote=F)
#
#
#
#write.table(K27_Factors1_2_nometh_NoK4_Factors_LF1[order(-K27_Factors1_2_nometh_NoK4_Factors_LF1$LF1),],'MOFA_K27_Factors_NoMeth_NoK4_Values_LF1.rnk',sep='\t',col.names = F,row.names=F,quote=F)
#write.table(K27_Factors1_2_nometh_NoK4_Factors_LF2[order(-K27_Factors1_2_nometh_NoK4_Factors_LF2$LF2),],'MOFA_K27_Factors_NoMeth_NoK4_Values_LF2.rnk',sep='\t',col.names = F,row.names=F,quote=F)
#
#
#ATAC_Factors1_2_nometh_NoK4<- getWeights(MOFAobject_NoK4,      views = "atac", factors = c(1,2) ,as.data.frame = TRUE)
#ATAC_Factors1_2_nometh_NoK4$feature=atac_y_new_cpm$geneId
#ATAC_Factors1_2_nometh_NoK4_LFT1=ATAC_Factors1_2_nometh_NoK4[ATAC_Factors1_2_nometh_NoK4$factor=='LF1',]
#ATAC_Factors1_2_nometh_NoK4_LFT2=ATAC_Factors1_2_nometh_NoK4[ATAC_Factors1_2_nometh_NoK4$factor=='LF2',]
#ATAC_Factors1_2_nometh_NoK4_Factors=cbind(ATAC_Factors1_2_nometh_NoK4_LFT1,ATAC_Factors1_2_nometh_NoK4_LFT2)
#ATAC_Factors1_2_nometh_NoK4_Factors=ATAC_Factors1_2_nometh_NoK4_Factors[,c(1,3,7)]
#
#
#ATAC_Factors1_2_nometh_NoK4_Factors_meanLF1=ATAC_Factors1_2_nometh_NoK4_Factors[,c(1,2)]
#ATAC_Factors1_2_nometh_NoK4_Factors_meanLF1=aggregate(ATAC_Factors1_2_nometh_NoK4_Factors_meanLF1[,2], list(Genes=ATAC_Factors1_2_nometh_NoK4_Factors_meanLF1[,1]), FUN=mean)
#ATAC_Factors1_2_nometh_NoK4_Factors_meanLF2=ATAC_Factors1_2_nometh_NoK4_Factors[,c(1,3)]
#ATAC_Factors1_2_nometh_NoK4_Factors_meanLF2=aggregate(ATAC_Factors1_2_nometh_NoK4_Factors_meanLF2[,2], list(Genes=ATAC_Factors1_2_nometh_NoK4_Factors_meanLF2[,1]), FUN=mean)
#
#
#

#write.table(ATAC_Factors1_2_nometh_NoK4_Factors_meanLF1[order(-ATAC_Factors1_2_nometh_NoK4_Factors_meanLF1$x),],'MOFA_ATAC_Factors_NoMeth_NoK4_Values_LF1.rnk',sep='\t',col.names = F,row.names=F,quote=F)
#write.table(ATAC_Factors1_2_nometh_NoK4_Factors_meanLF2[order(-ATAC_Factors1_2_nometh_NoK4_Factors_meanLF2$x),],'MOFA_ATAC_Factors_NoMeth_NoK4_Values_LF2.rnk',sep='\t',col.names = F,row.names=F,quote=F)

#
#
#Dat_Factor1_pos=Dat_Factor1[order(-Dat_Factor1$value),]
#Dat_Factor1_pos=Dat_Factor1_pos[1:50,]
#Dat_Factor1_neg=Dat_Factor1[order(Dat_Factor1$value),]
#Dat_Factor1_neg=Dat_Factor1_neg[1:50,]
#Dat_Factor1_Top100=rbind(Dat_Factor1_pos,Dat_Factor1_neg)
#
#
#Dat_Factor1_pos_GO=Dat_Factor1[order(-Dat_Factor1$value),]
#Dat_Factor1_pos_GO=Dat_Factor1_pos_GO[1:50,]
#Dat_Factor1_neg_GO=Dat_Factor1[order(Dat_Factor1$value),]
#Dat_Factor1_neg_GO=Dat_Factor1_neg_GO[1:50,]
#Dat_Factor1_GO_Top100=rbind(Dat_Factor1_pos_GO,Dat_Factor1_neg_GO)
#
#
#library(biomaRt)
#mart=useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="hsapiens_gene_ensembl")
##RNASEQ
#
#RNAseq_all_genes<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=toptags_RNA$genes,mart=mart)
#TopFeatures_RNA_Factor1<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=Dat_Factor1_GO_Top100$feature,mart=mart)
#
#rna_topFeatures_Factor1 <- enrichGO(gene          = TopFeatures_RNA_Factor1$entrezgene,
#universe      = as.character(RNAseq_all_genes$entrezgene),
#OrgDb         = org.Hs.eg.db,
#ont           = "BP",
#pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = T)
#dotplot(rna_topFeatures_Factor1)
#
#
#
#Dat_Factor2_pos=Dat_Factor2[order(-Dat_Factor2$value),]
#Dat_Factor2_pos=Dat_Factor2_pos[1:50,]
#Dat_Factor2_neg=Dat_Factor2[order(Dat_Factor2$value),]
#Dat_Factor2_neg=Dat_Factor2_neg[1:50,]
#Dat_Factor2_Top100=rbind(Dat_Factor2_pos,Dat_Factor2_neg)
#
#
#TopFeatures_RNA_Factor2<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=Dat_Factor2_Top100$feature,mart=mart)
#
#rna_topFeatures_Factor2 <- enrichGO(gene          = TopFeatures_RNA_Factor2$entrezgene,
#universe      = as.character(RNAseq_all_genes$entrezgene),
#OrgDb         = org.Hs.eg.db,
#ont           = "BP",
#pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = T)
#dotplot(rna_topFeatures_Factor2)
#
#
#
#
#Dat_Factor3_pos=Dat_Factor3[order(-Dat_Factor3$value),]
#Dat_Factor3_pos=Dat_Factor3_pos[1:50,]
#Dat_Factor3_neg=Dat_Factor3[order(Dat_Factor3$value),]
#Dat_Factor3_neg=Dat_Factor3_neg[1:50,]
#Dat_Factor3_Top100=rbind(Dat_Factor3_pos,Dat_Factor3_neg)
#
#TopFeatures_RNA_Factor3<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=Dat_Factor3_Top100$feature,mart=mart)
#
#rna_topFeatures_Factor3 <- enrichGO(gene          = TopFeatures_RNA_Factor3$entrezgene,
#universe      = as.character(RNAseq_all_genes$entrezgene),
#OrgDb         = org.Hs.eg.db,
#ont           = "BP",
#pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = T)
#dotplot(rna_topFeatures_Factor3)
#
#
#
#Dat_Factor4_pos=Dat_Factor4[order(-Dat_Factor4$value),]
#Dat_Factor4_pos=Dat_Factor4_pos[1:50,]
#Dat_Factor4_neg=Dat_Factor4[order(Dat_Factor4$value),]
#Dat_Factor4_neg=Dat_Factor4_neg[1:50,]
#Dat_Factor4_Top100=rbind(Dat_Factor4_pos,Dat_Factor4_neg)
#
#
#TopFeatures_RNA_Factor4<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=Dat_Factor4_Top100$feature,mart=mart)
#
#rna_topFeatures_Factor4 <- enrichGO(gene          = TopFeatures_RNA_Factor4$entrezgene,
#universe      = as.character(RNAseq_all_genes$entrezgene),
#OrgDb         = org.Hs.eg.db,
#ont           = "BP",
#pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = T)
#dotplot(rna_topFeatures_Factor4)
#
#
#
#
#write.table(as.data.frame(Dat_Factor1_Top100[,c(1,3)]),'RNA_Factor1_MOFA_NoMeth.rnk',sep='\t',col.names=F,row.names=F,quote=F)
#write.table(as.data.frame(Dat_Factor2_Top100[,c(1,3)]),'RNA_Factor2_MOFA_NoMeth.rnk',sep='\t',col.names=F,row.names=F,quote=F)
#write.table(as.data.frame(Dat_Factor3_Top100[,c(1,3)]),'RNA_Factor3_MOFA_NoMeth.rnk',sep='\t',col.names=F,row.names=F,quote=F)
#write.table(as.data.frame(Dat_Factor4_Top100[,c(1,3)]),'RNA_Factor4_MOFA_NoMeth.rnk',sep='\t',col.names=F,row.names=F,quote=F)
#
#plotWeights(  MOFAobject, view = "rna",      factor = 1,      nfeatures = 50 )
#plotWeights(    MOFAobject_nometh,      view = "rna",      factor = 1,     nfeatures = 50
#
#
#
#
#'''
#Faccio di nuovo analisi levando il constraint del piccolo numero di numfactors
#'''
#
#
#MOFAobject_noconstrain <- createMOFAobject(ll)
#plotTilesData(MOFAobject_noconstrain)
#DataOptions <- getDefaultDataOptions()
#DataOptions
#ModelOptions <- getDefaultModelOptions(MOFAobject_noconstrain)
#TrainOptions <- getDefaultTrainOptions()
#TrainOptions$DropFactorThreshold <- 0.02
#TrainOptions
#
#
#MOFAobject_noconstrain <- prepareMOFA(
#MOFAobject_noconstrain,
#DataOptions = DataOptions,
#ModelOptions = ModelOptions,
#TrainOptions = TrainOptions)
#MOFAobject_noconstrain <- runMOFA(MOFAobject_noconstrain, outfile=tempfile())
#
#```
#```
#
#MOFAweights <- getWeights(MOFAobject_nometh,      views = "rna", factors = 1,  as.data.frame = TRUE)
#
# MOFAweights=MOFAweights[order(MOFAweights$value),]
# MOFAweights_bottom=MOFAweights[1:25,]
# MOFAweights=MOFAweights[order(-MOFAweights$value),]
# MOFAweights_top=MOFAweights[1:25,]
# MOFAweights_topfeatures=rbind(MOFAweights_top,MOFAweights_bottom)
# 
#plotFactorScatter(MOFAobject_nometh,
#                  factors = 1,
#                  color_by = "CD74")
#
#
#MOFAweights <- getWeights(MOFAobject_nometh,      views = "rna", factors = 1,  as.data.frame = TRUE)
#
# MOFAweights=MOFAweights[order(MOFAweights$value),]
# MOFAweights_bottom=MOFAweights[1:25,]
# MOFAweights=MOFAweights[order(-MOFAweights$value),]
# MOFAweights_top=MOFAweights[1:25,]
# MOFAweights_topfeatures=rbind(MOFAweights_top,MOFAweights_bottom)
# 
# 
# 
#
#RNAseq_all_genes<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=toptags_RNA$genes,mart=mart)
#TopFeatures_RNA_Factor2<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=MOFAweights_topfeatures$feature,mart=mart)
#
#rna_topFeatures_Factor2 <- enrichGO(gene          = TopFeatures_RNA_Factor2$entrezgene,
#universe      = as.character(RNAseq_all_genes$entrezgene),
#OrgDb         = org.Hs.eg.db,
#ont           = "BP",
#pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = T)
#dotplot(rna_topFeatures_Factor2)
#
#




#RIFACCIO TUTTO CON I SECONDI PASSAGGI
#RNA

newcounts=featureCounts(files=c(
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_1.bam',
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_2.bam',
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_3.bam',
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_4.bam',
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_5.bam',
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_6.bam',
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_7.bam',
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_8.bam',
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_9.bam',
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_10.bam',
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_11.bam',
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_12.bam',
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_13.bam',
'/lustre1/workspace/DiMicco/prj_ALL_pilot/497_RNAseq/BAM/RNA_UPN3_14.bam'),
annot.ext='/lustre1/genomes/hg38/annotation/gencode.v26.chr_patch_hapl_scaff.annotation.gtf',
isGTFAnnotationFile = TRUE,GTF.attrType='gene_name',chrAliases='Aliases_GTF.txt')



colnames(newcounts$counts)=c("UPN3_1","UPN3_2","UPN3_3","UPN3_4","UPN3_5","UPN3_6","UPN3_7","UPN3_8","UPN3_9","UPN3_10","UPN3_11","UPN3_12","UPN3_13","UPN3_14")

rna_newcounts = DGEList(newcounts$counts[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)], genes = newcounts$annotation$GeneID)
rna_newcounts = calcNormFactors(rna_newcounts,method='TMM')
rna_newcounts$samples$group = as.factor(group[row.names(rna_newcounts$samples)])
rna_keep_newcounts = rowSums(cpm(rna_newcounts) > 1) >= 3
rna_newcounts = rna_newcounts[rna_keep_newcounts, ]
rna_newcounts_cpm=cbind(as.data.frame(rna_newcounts$genes),cpm(rna_newcounts,log=T))
df = targets[row.names(rna_newcounts$samples),]
D = as.factor(df$disease)
P = as.factor(df$passage)
design = model.matrix(~D+P)
design
rna_newcounts= estimateDisp(rna_newcounts, design)
fit = glmFit(rna_newcounts, design)
lrt <- glmLRT(fit, coef=2)
toptags_atac=as.data.frame(topTags(lrt,sort.by="none",n="all"))
write.table(toptags_atac,'TopTags_RNA_AllSamples_ToQvalue.txt',sep='\t',col.names=T,row.names=F,quote=F)
 





#ATAC
 atac_all=readPeakFile("All_Atac_Peaks_Counts_MultiCov_NoReps.txt",head=F)
 atac_all_anno <- annotatePeak(peak = atac_all, tssRegion=c(-2000, 2000),TxDb=gencode_rnaseqtxdb) #,annoDb="org.Hs.eg.db")
 atac_all_anno_df=as.data.frame(atac_all_anno)
 colnames(atac_all_anno_df)=c(colnames(atac_all_anno_df)[1:5],"UPN3_1","UPN3_2","UPN3_3","UPN3_4","UPN3_5","UPN3_6","UPN3_7","UPN3_8","UPN3_9","UPN3_10","UPN3_11","UPN3_13","UPN3_14",colnames(atac_all_anno_df)[19:27])

#ANALISI STANDARD EDGER normalizzazione RLE e con le lib.sizes sotto i picchi
atac_all_anno_df$annotation = sub(" .*", "",atac_all_anno_df$annotation)
atac_counts_new=atac_all_anno_df[,c(1,2,3,4,25,19,27,6:18)]
atac_y_new=DGEList(counts=atac_counts_new[,8:20],genes=atac_counts_new[,c(1:7)],group=c("D_0","D_1","D_1","D_1","D_2","D_2","D_2","R_0","R_1","R_1","R_1","R_2","R_2"))
#atac_y_new = calcNormFactors(atac_y_new,method='RLE')#,Acutoff=-18.5)
#atac_new_keep = rowSums(cpm(atac_y_new) > 5) >= 3
#atac_y_new = atac_y_new[atac_new_keep, ]
keep <- rowSums(cpm(atac_y_new)>5) >= 3 #PRENDO SOLO I GENI CON CPM>1 IN ALMENO 8 CAMPIONI
atac_y_new<- atac_y_new[keep, , keep.lib.sizes=FALSE]
atac_y_new = calcNormFactors(atac_y_new,method='RLE')#,Acutoff=-18.5)
atac_y_new_cpm_allsamples=cbind(atac_y_new$genes,cpm(atac_y_new,log = T))
UPN3_12=rep(NA,nrow(atac_y_new_cpm_allsamples))
atac_y_new_cpm_allsamples=cbind(atac_y_new_cpm_allsamples[,1:18],as.data.frame(UPN3_12),atac_y_new_cpm_allsamples[,19:20])
df = targets[row.names(atac_y_new$samples),]
D = as.factor(df$disease)
P = as.factor(df$passage)
design = model.matrix(~D+P)
design
atac_y_new= estimateDisp(atac_y_new, design)
fit = glmFit(atac_y_new, design)
lrt <- glmLRT(fit, coef=2)
toptags_atac=as.data.frame(topTags(lrt,sort.by="none",n="all"))
write.table(toptags_atac,'TopTags_Atac_AllSamples_ToQvalue.txt',sep='\t',col.names=T,row.names=F,quote=F)
 



#K4me3
k4me3_all=readPeakFile("All_K4me3_Peaks_Counts_MultiCov_New_NoReps_New.txt",head=F)
k4me3_all_anno <- annotatePeak(peak = k4me3_all, tssRegion=c(-2000, 2000),TxDb=gencode_rnaseqtxdb)#,annoDb="org.Hs.eg.db")
 #plotAnnoPie(k4me3_all_anno)


k4me3_all_anno_df=as.data.frame(k4me3_all_anno)
k4me3_all_anno_df$annotation = sub(" .*", "", k4me3_all_anno_df$annotation)
colnames(k4me3_all_anno_df)=c(colnames(k4me3_all_anno_df)[1:5],"UPN3_1","UPN3_2","UPN3_3","UPN3_4","UPN3_5","UPN3_6","UPN3_7","UPN3_8","UPN3_9","UPN3_10","UPN3_11","UPN3_12","UPN3_13","UPN3_14",colnames(k4me3_all_anno_df)[20:28])

#head(k4me3_all_anno_df)
k4me3_counts_new=k4me3_all_anno_df[,c(1,2,3,26,20,28,6:19)]
k4me3_y_new=DGEList(counts=k4me3_counts_new[,c(7:20)],genes=k4me3_counts_new[,c(1:6)],group=c("D_0","D_1","D_1","D_1","D_2","D_2","D_2","R_0","R_1","R_1","R_1","R_2","R_2","R_2"))#,lib.size = c(22089014,23497699,30146073,21596634,24792817,22820342,26900163,45421712))
#k4me3_y_new = calcNormFactors(k4me3_y_new,method='RLE')#,Acutoff=-18.5)
#k4me3_new_keep = rowSums(cpm(k4me3_y_new) > 5) >= 3
#k4me3_y_new = k4me3_y_new[k4me3_new_keep, ]
keep <- rowSums(cpm(k4me3_y_new)>5) >= 3 #PRENDO SOLO I GENI CON CPM>1 IN ALMENO 8 CAMPIONI
k4me3_y_new<- k4me3_y_new[keep, , keep.lib.sizes=FALSE]
k4me3_y_new = calcNormFactors(k4me3_y_new,method='RLE')#,Acutoff=-18.5)
k4me3_y_new_cpm_allsamples=cbind(k4me3_y_new$genes,cpm(k4me3_y_new,log = T))
library(qvalue)
 df = targets[row.names(k4me3_y_new$samples),]
 D = as.factor(df$disease)
 P = as.factor(df$passage)
 design = model.matrix(~D+P)
 design
 k4me3_y_new= estimateDisp(k4me3_y_new, design)
 fit = glmFit(k4me3_y_new, design)
 lrt <- glmLRT(fit, coef=2)
 toptags_k4me3=as.data.frame(topTags(lrt,sort.by="none",n="all"))
 write.table(toptags_k4me3,'TopTags_K4me3_AllSamples_ToQvalue.txt',sep='\t',col.names=T,row.names=F)
 
 q = qvalue(toptags_k4me3$PValue)
  toptags_k4me3$Qvalue = q$qvalues
 toptags_k4me3_DEG=toptags_k4me3[toptags_k4me3$Qvalue<0.05,]
k4me3_y_new$samples$group=as.factor(design[,2])




#K27ac

k27ac_all=readPeakFile("All_K27ac_Peaks_Counts_MultiCov_NoReps.txt",head=F)
k27ac_all_anno <- annotatePeak(peak = k27ac_all, tssRegion=c(-2000, 2000),TxDb=gencode_rnaseqtxdb) #,annoDb="org.Hs.eg.db")
k27ac_all_anno_df=as.data.frame(k27ac_all_anno)
k27ac_all_anno_df$annotation = sub(" .*", "", k27ac_all_anno_df$annotation)
colnames(k27ac_all_anno_df)=c(colnames(k27ac_all_anno_df)[1:5],"UPN3_1","UPN3_2","UPN3_3","UPN3_4","UPN3_5","UPN3_6","UPN3_7","UPN3_8","UPN3_9","UPN3_10","UPN3_11","UPN3_12","UPN3_13","UPN3_14",colnames(k27ac_all_anno_df)[20:28])
k27ac_counts_new=k27ac_all_anno_df[,c(1,2,3,4,26,20,28,6:19)]
k27ac_y_new=DGEList(counts=k27ac_counts_new[,c(8:21)],genes=k27ac_counts_new[,c(1:7)],group=c("D_0","D_1","D_1","D_1","D_2","D_2","D_2","R_0","R_1","R_1","R_1","R_2","R_2","R_2"))
#k27ac_y_new = calcNormFactors(k27ac_y_new,method='RLE')#,Acutoff=-18.5)
#k27ac_new_keep = rowSums(cpm(k27ac_y_new) > 5) >= 3
#k27ac_y_new = k27ac_y_new[k27ac_new_keep, ]
keep <- rowSums(cpm(k27ac_y_new)>5) >= 3 #PRENDO SOLO I GENI CON CPM>1 IN ALMENO 8 CAMPIONI
k27ac_y_new<- k27ac_y_new[keep, , keep.lib.sizes=FALSE]
k27ac_y_new = calcNormFactors(k27ac_y_new,method='RLE')#,Acutoff=-18.5)
k27ac_y_new_cpm_allsamples=cbind(k27ac_y_new$genes,cpm(k27ac_y_new,log = T))
df = targets[row.names(k27ac_y_new$samples),]
D = as.factor(df$disease)
P = as.factor(df$passage)
design = model.matrix(~D+P)
design
k27ac_y_new= estimateDisp(k27ac_y_new, design)
fit = glmFit(k27ac_y_new, design)
lrt <- glmLRT(fit, coef=2)
toptags_k27ac=as.data.frame(topTags(lrt,sort.by="none",n="all"))
write.table(toptags_k27ac,'TopTags_K27ac_AllSamples_ToQvalue.txt',sep='\t',col.names=T,row.names=F,quote=F)
 




ll_allsamples=list()
ll_allsamples$acetil=k27ac_y_new_cpm_allsamples[,8:21]
ll_allsamples$atac=atac_y_new_cpm_allsamples[,8:21]
ll_allsamples$rna=rna_newcounts_cpm[2:15]
ll_allsamples$trimetil=k4me3_y_new_cpm_allsamples[,7:20]


MOFAobject_AllSamples <- MOFAobject_AllSamples(ll_allsamples)
plotTilesData(MOFAobject_AllSamples)
DataOptions <- getDefaultDataOptions()
DataOptions
ModelOptions <- getDefaultModelOptions(MOFAobject_AllSamples)
TrainOptions <- getDefaultTrainOptions()
TrainOptions$DropFactorThreshold <- 0.02
TrainOptions


MOFAobject_AllSamples <- prepareMOFA(
MOFAobject_AllSamples,
DataOptions = DataOptions,
ModelOptions = ModelOptions,
TrainOptions = TrainOptions)
MOFAobject_AllSamples <- runMOFA(MOFAobject_AllSamples, outfile=tempfile())


r2_allsamples <- calculateVarianceExplained(MOFAobject_AllSamples)
r2

plotVarianceExplained(MOFAobject_AllSamples)

```
$R2Total
   acetil      atac       rna  trimetil 
0.4705220 0.8309163 0.4405562 0.5022958 

$R2PerFactor
          acetil         atac          rna     trimetil
LF1 1.104740e-01 5.898239e-01 6.317011e-02 8.721174e-02
LF2 6.015932e-05 2.467785e-02 2.560043e-01 1.435461e-01
LF3 3.148389e-02 6.910538e-02 5.872213e-02 1.950202e-01
LF4 2.222388e-01 2.272495e-05 3.292180e-05 5.562800e-02
LF5 1.230093e-01 2.288807e-05 2.777648e-02 2.379787e-02
LF6 3.391632e-05 9.915800e-02 3.048016e-05 1.087269e-05
LF7 2.649425e-05 5.833288e-02 5.327847e-05 1.249004e-05
LF8 1.518583e-03 3.811892e-05 3.703321e-02 6.745522e-03
```


RNA_Factors1_2_nometh_allSamples<- getWeights(MOFAobject_AllSamples,      views = "trimetil", factors = c(1,2), as.data.frame = TRUE)
RNA_Factors1_2_nometh_allSamples_LFT1=RNA_Factors1_2_nometh_allSamples[RNA_Factors1_2_nometh_allSamples$factor=='LF1',]
RNA_Factors1_2_nometh_allSamples_LFT2=RNA_Factors1_2_nometh_allSamples[RNA_Factors1_2_nometh_allSamples$factor=='LF2',]
RNA_Factors1_2_nometh_allSamples_Factors=cbind(RNA_Factors1_2_nometh_allSamples_LFT1,RNA_Factors1_2_nometh_allSamples_LFT2)
RNA_Factors1_2_nometh_allSamples_Factors$genes=rna_newcounts_cpm$genes
RNA_Factors1_2_nometh_allSamples_Factors=RNA_Factors1_2_nometh_allSamples_Factors[,c(9,3,7)]

colnames(RNA_Factors1_2_nometh_allSamples_Factors)=c('Genes','LF1','LF2')
RNA_Factors1_2_nometh_allSamples_Factors_LF1=RNA_Factors1_2_nometh_allSamples_Factors[,1:2]
RNA_Factors1_2_nometh_allSamples_Factors_LF2=RNA_Factors1_2_nometh_allSamples_Factors[,c(1,3)]

write.table(RNA_Factors1_2_nometh_allSamples_Factors_LF1[order(-RNA_Factors1_2_nometh_allSamples_Factors_LF1$LF1),],'MOFA_RNA_Factors_NoMeth_Values_LF1_AllSamples.rnk',sep='\t',col.names = F,row.names=F,quote=F)
write.table(RNA_Factors1_2_nometh_allSamples_Factors_LF2[order(-RNA_Factors1_2_nometh_allSamples_Factors_LF2$LF2),],'MOFA_RNA_Factors_NoMeth_Values_LF2_AllSamples.rnk',sep='\t',col.names = F,row.names=F,quote=F)

RNA_Factors3_4_nometh_allSamples<- getWeights(MOFAobject_AllSamples,      views = "trimetil", factors = c(3,4), as.data.frame = TRUE)
RNA_Factors3_4_nometh_allSamples_LFT3=RNA_Factors3_4_nometh_allSamples[RNA_Factors3_4_nometh_allSamples$factor=='LF3',]
RNA_Factors3_4_nometh_allSamples_LFT4=RNA_Factors3_4_nometh_allSamples[RNA_Factors3_4_nometh_allSamples$factor=='LF4',]
RNA_Factors3_4_nometh_allSamples_Factors=cbind(RNA_Factors3_4_nometh_allSamples_LFT3,RNA_Factors3_4_nometh_allSamples_LFT4)
RNA_Factors3_4_nometh_allSamples_Factors$genes=rna_newcounts_cpm$genes
RNA_Factors3_4_nometh_allSamples_Factors=RNA_Factors3_4_nometh_allSamples_Factors[,c(9,3,7)]

colnames(RNA_Factors3_4_nometh_allSamples_Factors)=c('Genes','LF3','LF4')
RNA_Factors3_4_nometh_allSamples_Factors_LF3=RNA_Factors3_4_nometh_allSamples_Factors[,1:2]
RNA_Factors3_4_nometh_allSamples_Factors_LF4=RNA_Factors3_4_nometh_allSamples_Factors[,c(1,3)]

write.table(RNA_Factors3_4_nometh_allSamples_Factors_LF3[order(-RNA_Factors3_4_nometh_allSamples_Factors_LF3$LF3),],'MOFA_RNA_Factors_NoMeth_Values_LF3_AllSamples.rnk',sep='\t',col.names = F,row.names=F,quote=F)
write.table(RNA_Factors3_4_nometh_allSamples_Factors_LF4[order(-RNA_Factors3_4_nometh_allSamples_Factors_LF4$LF4),],'MOFA_RNA_Factors_NoMeth_Values_LF4_AllSamples.rnk',sep='\t',col.names = F,row.names=F,quote=F)



Tabella_GSEA_MOFA_RNA_AllSamples_LF1LF2LF3F4=read.table('Tabella_MOFA_RnaSeq_AllSamples_Factors1234.txt',head=T,sep='\t')
Tabella_GSEA_MOFA_RNA_AllSamples_LF1LF2LF3F4_NES=Tabella_GSEA_MOFA_RNA_AllSamples_LF1LF2LF3F4[,c(1,2,4,6,8)]
Tabella_GSEA_MOFA_RNA_AllSamples_LF1LF2LF3F4_FDR=Tabella_GSEA_MOFA_RNA_AllSamples_LF1LF2LF3F4[,c(1,3,5,7,9)]
Tabella_GSEA_MOFA_RNA_AllSamples_LF1LF2LF3F4_NES=melt(Tabella_GSEA_MOFA_RNA_AllSamples_LF1LF2LF3F4_NES)
Tabella_GSEA_MOFA_RNA_AllSamples_LF1LF2LF3F4_FDR=melt(Tabella_GSEA_MOFA_RNA_AllSamples_LF1LF2LF3F4_FDR)
Tabella_GSEA_MOFA_RNA_AllSamples_LF1LF2LF3F4_heatmap=cbind(Tabella_GSEA_MOFA_RNA_AllSamples_LF1LF2LF3F4_NES,Tabella_GSEA_MOFA_RNA_AllSamples_LF1LF2LF3F4_FDR)
Tabella_GSEA_MOFA_RNA_AllSamples_LF1LF2LF3F4_heatmap=Tabella_GSEA_MOFA_RNA_AllSamples_LF1LF2LF3F4_heatmap[,c(1,2,3,6)]
colnames(Tabella_GSEA_MOFA_RNA_AllSamples_LF1LF2LF3F4_heatmap)=c('Hallmark','Technique_LFs', 'NES','FDR') 


p=ggplot(Tabella_GSEA_MOFA_RNA_AllSamples_LF1LF2LF3F4_heatmap, aes(Technique_LFs,Hallmark)) +     geom_point(aes(size = -log10(FDR), colour=NES))   +  theme(panel.background = element_rect(fill = "white",colour = "black",size = 0.5, linetype = "solid"),panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "black"),panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour ="black")) + theme(legend.position="right") +theme_grey(base_size = 10,) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_gradient2(low = "midnightblue", mid = "white",high = "orange4", midpoint = 0, space = "Lab",na.value = "grey50", guide = "colourbar", aesthetics = "colour") + coord_flip()
pdf('Plot_MOFA_RNAseq_AllSamples_FActors1234.pdf')
print(p, vp=viewport(angle=-90))
dev.off()




ATAC_Factors1_2_nometh_allSamples<- getWeights(MOFAobject_AllSamples,      views = "rna", factors = c(1,2), as.data.frame = TRUE)
ATAC_Factors1_2_nometh_allSamples_LFT1=ATAC_Factors1_2_nometh_allSamples[ATAC_Factors1_2_nometh_allSamples$factor=='LF1',]
ATAC_Factors1_2_nometh_allSamples_LFT2=ATAC_Factors1_2_nometh_allSamples[ATAC_Factors1_2_nometh_allSamples$factor=='LF2',]
ATAC_Factors1_2_nometh_allSamples_Factors=cbind(ATAC_Factors1_2_nometh_allSamples_LFT1,ATAC_Factors1_2_nometh_allSamples_LFT2)
ATAC_Factors1_2_nometh_allSamples_Factors$genes=atac_y_new_cpm_allsamples$geneId
ATAC_Factors1_2_nometh_allSamples_Factors=ATAC_Factors1_2_nometh_allSamples_Factors[,c(9,3,7)]

colnames(ATAC_Factors1_2_nometh_allSamples_Factors)=c('Genes','LF1','LF2')
ATAC_Factors1_2_nometh_allSamples_Factors_LF1=ATAC_Factors1_2_nometh_allSamples_Factors[,1:2]

ATAC_Factors1_2_nometh_allSamples_Factors_LF1=aggregate(ATAC_Factors1_2_nometh_allSamples_Factors_LF1[,2], list(Genes=ATAC_Factors1_2_nometh_allSamples_Factors_LF1[,1]), FUN=mean)

ATAC_Factors1_2_nometh_allSamples_Factors_LF2=ATAC_Factors1_2_nometh_allSamples_Factors[,c(1,3)]

ATAC_Factors1_2_nometh_allSamples_Factors_LF2=aggregate(ATAC_Factors1_2_nometh_allSamples_Factors_LF2[,2], list(Genes=ATAC_Factors1_2_nometh_allSamples_Factors_LF2[,1]), FUN=mean)


write.table(ATAC_Factors1_2_nometh_allSamples_Factors_LF1[order(-ATAC_Factors1_2_nometh_allSamples_Factors_LF1$x),],'MOFA_ATAC_Factors_NoMeth_Values_LF1_AllSamples.rnk',sep='\t',col.names = F,row.names=F,quote=F)
write.table(ATAC_Factors1_2_nometh_allSamples_Factors_LF2[order(-ATAC_Factors1_2_nometh_allSamples_Factors_LF2$x),],'MOFA_ATAC_Factors_NoMeth_Values_LF2_AllSamples.rnk',sep='\t',col.names = F,row.names=F,quote=F)


K4_Factors1_2_nometh_allSamples<- getWeights(MOFAobject_AllSamples,      views = "acetil", factors = c(1,2), as.data.frame = TRUE)
K4_Factors1_2_nometh_allSamples_LFT1=K4_Factors1_2_nometh_allSamples[K4_Factors1_2_nometh_allSamples$factor=='LF1',]
K4_Factors1_2_nometh_allSamples_LFT2=K4_Factors1_2_nometh_allSamples[K4_Factors1_2_nometh_allSamples$factor=='LF2',]
K4_Factors1_2_nometh_allSamples_Factors=cbind(K4_Factors1_2_nometh_allSamples_LFT1,K4_Factors1_2_nometh_allSamples_LFT2)
K4_Factors1_2_nometh_allSamples_Factors$genes=k4me3_y_new_cpm_allsamples$geneId
K4_Factors1_2_nometh_allSamples_Factors=K4_Factors1_2_nometh_allSamples_Factors[,c(9,3,7)]

colnames(K4_Factors1_2_nometh_allSamples_Factors)=c('Genes','LF1','LF2')
K4_Factors1_2_nometh_allSamples_Factors_LF1=K4_Factors1_2_nometh_allSamples_Factors[,1:2]

K4_Factors1_2_nometh_allSamples_Factors_LF1=aggregate(K4_Factors1_2_nometh_allSamples_Factors_LF1[,2], list(Genes=K4_Factors1_2_nometh_allSamples_Factors_LF1[,1]), FUN=mean)

K4_Factors1_2_nometh_allSamples_Factors_LF2=K4_Factors1_2_nometh_allSamples_Factors[,c(1,3)]

K4_Factors1_2_nometh_allSamples_Factors_LF2=aggregate(K4_Factors1_2_nometh_allSamples_Factors_LF2[,2], list(Genes=K4_Factors1_2_nometh_allSamples_Factors_LF2[,1]), FUN=mean)

write.table(K4_Factors1_2_nometh_allSamples_Factors_LF1[order(-K4_Factors1_2_nometh_allSamples_Factors_LF1$x),],'MOFA_K4_Factors_NoMeth_Values_LF1_AllSamples.rnk',sep='\t',col.names = F,row.names=F,quote=F)
write.table(K4_Factors1_2_nometh_allSamples_Factors_LF2[order(-K4_Factors1_2_nometh_allSamples_Factors_LF2$x),],'MOFA_K4_Factors_NoMeth_Values_LF2_AllSamples.rnk',sep='\t',col.names = F,row.names=F,quote=F)


K27_Factors1_2_nometh_allSamples<- getWeights(MOFAobject_AllSamples,      views = "atac", factors = c(1,2), as.data.frame = TRUE)
K27_Factors1_2_nometh_allSamples_LFT1=K27_Factors1_2_nometh_allSamples[K27_Factors1_2_nometh_allSamples$factor=='LF1',]
K27_Factors1_2_nometh_allSamples_LFT2=K27_Factors1_2_nometh_allSamples[K27_Factors1_2_nometh_allSamples$factor=='LF2',]
K27_Factors1_2_nometh_allSamples_Factors=cbind(K27_Factors1_2_nometh_allSamples_LFT1,K27_Factors1_2_nometh_allSamples_LFT2)
K27_Factors1_2_nometh_allSamples_Factors$genes=k27ac_y_new_cpm_allsamples$geneId
K27_Factors1_2_nometh_allSamples_Factors=K27_Factors1_2_nometh_allSamples_Factors[,c(9,3,7)]

colnames(K27_Factors1_2_nometh_allSamples_Factors)=c('Genes','LF1','LF2')
K27_Factors1_2_nometh_allSamples_Factors_LF1=K27_Factors1_2_nometh_allSamples_Factors[,1:2]

K27_Factors1_2_nometh_allSamples_Factors_LF1=aggregate(K27_Factors1_2_nometh_allSamples_Factors_LF1[,2], list(Genes=K27_Factors1_2_nometh_allSamples_Factors_LF1[,1]), FUN=mean)

K27_Factors1_2_nometh_allSamples_Factors_LF2=K27_Factors1_2_nometh_allSamples_Factors[,c(1,3)]

K27_Factors1_2_nometh_allSamples_Factors_LF2=aggregate(K27_Factors1_2_nometh_allSamples_Factors_LF2[,2], list(Genes=K27_Factors1_2_nometh_allSamples_Factors_LF2[,1]), FUN=mean)

write.table(K27_Factors1_2_nometh_allSamples_Factors_LF1[order(-K27_Factors1_2_nometh_allSamples_Factors_LF1$x),],'MOFA_K27_Factors_NoMeth_Values_LF1_AllSamples.rnk',sep='\t',col.names = F,row.names=F,quote=F)
write.table(K27_Factors1_2_nometh_allSamples_Factors_LF2[order(-K27_Factors1_2_nometh_allSamples_Factors_LF2$x),],'MOFA_K27_Factors_NoMeth_Values_LF2_AllSamples.rnk',sep='\t',col.names = F,row.names=F,quote=F)


