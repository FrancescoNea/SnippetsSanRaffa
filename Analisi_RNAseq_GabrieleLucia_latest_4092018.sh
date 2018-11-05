#!/bin/sh
# ho fatto un RData pulito che si trova in /home/fsantaiello/Analisi_RNaseqGabriLucia_Clean.RData
library(edgeR)
library(ggplot2)
library(qvalue)
library(pheatmap)
library(Rsubread)
library(tximport)
library(readr)
library(fpc)
samples <- read.table(file.path("./", "samples_wDR1.tsv"),sep="\t", header = FALSE)
files <- file.path(".", paste(samples$V1,"kallisto",sep="_"), "abundance.tsv")
names(files) <- samples$sample
all(file.exists(files))
tx2gene <- read.csv(file.path("./", "Tx2Gene.txt"),sep="\t",header=FALSE)
txi_length <- tximport(files, type = "kallisto", tx2gene = tx2gene, reader = read_tsv,txOut = FALSE,lengthCol=2)
countdata=txi_length$counts
samplenames=c()
for (i in 1:length(files)){samplenames=c(samplenames,unlist(strsplit((files)[i], "/", fixed=TRUE))[2])}
colnames(countdata)=samplenames
condition=rep(c("Diagnosis","Relapse"),15)
countdata_filt <- countdata[rowSums(countdata[, 1:30]) > 0, ]
countdata_filtDGE=DGEList(countdata_filt[,1:30],genes=rownames(countdata_filt),group=condition)
keep <- rowSums(cpm(countdata_filtDGE)>1) >= 8 #PRENDO SOLO I GENI CON CPM>1 IN ALMENO 8 CAMPIONI
countdata_filtDGE_keep<- countdata_filtDGE[keep, , keep.lib.sizes=FALSE]
countdata_filtDGE_keep=calcNormFactors(countdata_filtDGE_keep,method="TMM")
mds=plotMDS(countdata_filtDGE_keep,dim.plot=c(1,2))
mds2=data.frame(x=mds$x,y=mds$y)
samples=rownames(mds2)
samples2=unlist(strsplit(samples,split="_kallisto"))
samples2=samples2[samples2!=""]
mds2=cbind(mds2,as.data.frame(samples2),as.data.frame(condition))
patient=rep(1:15,each=2)
mds2=cbind(mds2,as.data.frame(patient))
colnames(mds2)=c("x","y","samples","condition","patient")
design = model.matrix (~condition + x + y + as.factor(patient),data=mds2)
countdata_filtDGE_keep <- estimateDisp(countdata_filtDGE_keep,design)
fit <- glmFit(countdata_filtDGE_keep, design)
lrt <- glmLRT(fit,coef=2)
#Correggo per qvalue e seleziono solo i significativi per qvalue<=0.25
library(qvalue)
toptags_xypat_Kallisto=as.data.frame(topTags(lrt,sort.by="PValue",n="all"))
q = qvalue(toptags_xypat_Kallisto$PValue)
toptags_xypat_Kallisto$Qvalue = q$qvalues
toptags_xypat_Kallisto_DEG=toptags_xypat_Kallisto[toptags_xypat_Kallisto$Qvalue<=0.25,]
nrow(toptags_xypat_Kallisto_DEG)

toptags_xypat_Kallisto_cpm=cbind(cpm(countdata_filtDGE_keep,log=T),as.data.frame(countdata_filtDGE_keep$genes))
toptags_xypat_Kallisto_cpm=merge(toptags_xypat_Kallisto,toptags_xypat_Kallisto_cpm,by="genes")
toptags_xypat_Kallisto_cpm$AlfeFC=toptags_xypat_Kallisto_cpm$ALFE2_kallisto-toptags_xypat_Kallisto_cpm$ALFE1_kallisto
toptags_xypat_Kallisto_cpm$BesuFC=toptags_xypat_Kallisto_cpm$BESU2_kallisto-toptags_xypat_Kallisto_cpm$BESU1_kallisto
toptags_xypat_Kallisto_cpm$CaluFC=toptags_xypat_Kallisto_cpm$CALU2_kallisto-toptags_xypat_Kallisto_cpm$CALU1_kallisto
toptags_xypat_Kallisto_cpm$DeivFC=toptags_xypat_Kallisto_cpm$DEIV2b_kallisto-toptags_xypat_Kallisto_cpm$DEIV2_kallisto
toptags_xypat_Kallisto_cpm$DestFC=toptags_xypat_Kallisto_cpm$DEST2_kallisto-toptags_xypat_Kallisto_cpm$DEST1_kallisto
toptags_xypat_Kallisto_cpm$Dr1FC=toptags_xypat_Kallisto_cpm$DR1_2_kallisto-toptags_xypat_Kallisto_cpm$DR1_1_kallisto
toptags_xypat_Kallisto_cpm$Dr4FC=toptags_xypat_Kallisto_cpm$DR42_kallisto-toptags_xypat_Kallisto_cpm$DR41_kallisto
toptags_xypat_Kallisto_cpm$Dr5FC=toptags_xypat_Kallisto_cpm$DR52_kallisto-toptags_xypat_Kallisto_cpm$DR51_kallisto
toptags_xypat_Kallisto_cpm$FocaFC=toptags_xypat_Kallisto_cpm$FOCA2_kallisto-toptags_xypat_Kallisto_cpm$FOCA1_kallisto
toptags_xypat_Kallisto_cpm$GagraFC=toptags_xypat_Kallisto_cpm$GAGRA2_kallisto-toptags_xypat_Kallisto_cpm$GAGRA1_kallisto
toptags_xypat_Kallisto_cpm$LuanFC=toptags_xypat_Kallisto_cpm$LUAN2_kallisto-toptags_xypat_Kallisto_cpm$LUAN1_kallisto
toptags_xypat_Kallisto_cpm$MabiFC=toptags_xypat_Kallisto_cpm$MABI2_kallisto-toptags_xypat_Kallisto_cpm$MABI1_kallisto
toptags_xypat_Kallisto_cpm$MogeFC=toptags_xypat_Kallisto_cpm$MOGE2_kallisto-toptags_xypat_Kallisto_cpm$MOGE1_kallisto
toptags_xypat_Kallisto_cpm$PiagFC=toptags_xypat_Kallisto_cpm$PIAG2_kallisto-toptags_xypat_Kallisto_cpm$PIAG1_kallisto
toptags_xypat_Kallisto_cpm$PreluFC=toptags_xypat_Kallisto_cpm$PRELU2_kallisto-toptags_xypat_Kallisto_cpm$PRELU1_kallisto




toptags_xypat_Kallisto_cpm_DEG=merge(toptags_xypat_Kallisto_cpm,toptags_xypat_Kallisto_DEG,by="genes")

#
#pdf("FoldChanges_Kallisto_DEGs_15PAzienti.pdf")
#breakList=seq(-5, 5, by = 0.1) 
#pheatmap( toptags_xypat_Kallisto_cpm_DEG [,38:52], cluster_rows=T,cluster_cols=T,color=colorRampPalette( c("green", "black", "red"), space="rgb")(length(breakList)))
#dev.off()




#library(biomaRt)
#mart=useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="hsapiens_gene_ensembl")
#G_list_DegGenes_wDR1<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=toptags_xypat_Kallisto_DEG$genes,mart=mart)
#write.table(as.data.frame(G_list_DegGenes_wDR1$hgnc_symbol),"/home/fsantaniello/GeniDeregolati_AnalisiConDr1.txt",sep="\t",col.names=F,row.names=F,quote=F)
#
#G_list_AllGenes_wDR1<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=toptags_xypat_Kallisto$genes,mart=mart)
#write.table(as.data.frame(G_list_AllGenes_wDR1$entrezgene),"/home/fsantaniello/Background_Geni_AnalisiConDr1.txt",sep="\t",col.names=F,row.names=F,quote=F)


library(ggcorrplot)
pdf("/home/fsantaniello/Correlazioni_AllGenes_Pazienti_NewAnalysis.pdf")
corr <- round(cor(toptags_xypat_Kallisto_cpm_DEG[,38:52]), 1)
p.mat <- cor_pmat(toptags_xypat_Kallisto_cpm[,38:52])
col<- colorRampPalette(c("blue", "white", "red"))(20)
ggcorrplot(corr, hc.order = TRUE,    outline.col = "white",   p.mat = p.mat,insig = "blank")
dev.off()
library(Hmisc)


#Network Analysis su logFCs

res2allgenes <- rcorr(as.matrix(t(toptags_xypat_Kallisto_cpm[,38:52])))
cor.matrix=res2allgenes$r
colnames(cor.matrix)=toptags_xypat_Kallisto_cpm$genes
rownames(cor.matrix)=toptags_xypat_Kallisto_cpm$genes

cor.matrix_rows=apply(cor.matrix, 1, percent_rank)
cor.matrix_cols=apply(cor.matrix, 2, percent_rank)
cor.matrix_sqrt=sqrt(cor.matrix_rows*cor.matrix_cols)
library(reshape2)
my_cor_df <- melt(cor.matrix_sqrt)
my_cor_df= my_cor_df %>% filter(Var1 != Var2)
my_adj_list <- my_cor_df %>% filter(value > 0.995)
names(my_adj_list) <- c('from', 'to', 'weight') 
g=graph.data.frame(my_adj_list, directed=F)
communities=cluster_louvain(g, weights = NULL)
comms_Dataframe=as.data.frame(as.matrix(membership(communities)))
comms_Dataframe$genes=rownames(comms_Dataframe)
colnames(comms_Dataframe)=c('communities','genes')   
kallisto_comms=merge(toptags_xypat_Kallisto_cpm,comms_Dataframe,by="genes")

#piu stringente
res2allgenes <- rcorr(as.matrix(t(toptags_xypat_Kallisto_cpm[,38:52])))
cor.matrix=res2allgenes$r
colnames(cor.matrix)=toptags_xypat_Kallisto_cpm$genes
rownames(cor.matrix)=toptags_xypat_Kallisto_cpm$genes

cor.matrix_rows=apply(cor.matrix, 1, percent_rank)
cor.matrix_cols=apply(cor.matrix, 2, percent_rank)
cor.matrix_sqrt=sqrt(cor.matrix_rows*cor.matrix_cols)
library(reshape2)
my_cor_df <- melt(cor.matrix_sqrt)
my_cor_df= my_cor_df %>% filter(Var1 != Var2)
my_adj_list <- my_cor_df %>% filter(value > 0.995)
names(my_adj_list) <- c('from', 'to', 'weight') 
g_new=graph.data.frame(my_adj_list, directed=F)
communities_new=cluster_louvain(g, weights = NULL)
comms_Dataframe_new=as.data.frame(as.matrix(membership(communities_new)))
comms_Dataframe_new$genes=rownames(comms_Dataframe_new)
colnames(comms_Dataframe_new)=c('communities','genes')   
kallisto_comms_new=merge(toptags_xypat_Kallisto_cpm,comms_Dataframe_new,by="genes")



breakList=seq(-7, 7, by = 0.1)
for (i in 1:11){pdf(paste('Heatmap_FoldChanges_Comunita_AllGenes_',i,'.pdf',sep=''));a=kallisto_comms[kallisto_comms$communities==i,];rownames(a)=a$genes;pheatmap( a[,38:52], cluster_rows=T,cluster_cols=T,color=colorRampPalette( c("cyan", "white", "yellow"), space="rgb")(length(breakList)),breaks = breakList,main = paste('Community',i,sep=' '),border_color = NA,fontsize_row = 6);dev.off() }


breakList=seq(-7, 7, by = 0.1)
for (i in 1:13){pdf(paste('Heatmap_FoldChanges_Comunita_AllGenes_MoreStringent_',i,'.pdf',sep=''));a=kallisto_comms_new[kallisto_comms_new$communities==i,];rownames(a)=a$genes;pheatmap( a[,38:52], cluster_rows=T,cluster_cols=T,color=colorRampPalette( c("cyan", "black", "yellow"), space="rgb")(length(breakList)),breaks = breakList,main = paste('Community ',i,sep=' '),border_color = NA,fontsize_row = 6);dev.off() }


#Calcolo le sottocommunities 
for (i in 1:11){             nam=kallisto_comms_new[kallisto_comms$communities==i,];
                             print(head(nam));
                             res2allgenes <- rcorr(as.matrix(t(nam[,38:52])));
                             cor.matrix=res2allgenes$r;
                             colnames(cor.matrix)=nam$genes;
                             rownames(cor.matrix)=nam$genes;
                             cor.matrix_rows=apply(cor.matrix, 1, percent_rank);
                             cor.matrix_cols=apply(cor.matrix, 2, percent_rank);
                             cor.matrix_sqrt=sqrt(cor.matrix_rows*cor.matrix_cols);
                             library(reshape2);
                             my_cor_df <- melt(cor.matrix_sqrt);
                             my_cor_df= my_cor_df %>% filter(Var1 != Var2);
                             my_adj_list <- my_cor_df %>% filter(value > 0.999);
                             names(my_adj_list) <- c('from', 'to', 'weight') ;
                             g=graph.data.frame(my_adj_list, directed=F);
                             communities=cluster_louvain(g, weights = NULL);
                             comms_Dataframe=as.data.frame(as.matrix(membership(communities)));
                             comms_Dataframe$genes=rownames(comms_Dataframe);
                             colnames(comms_Dataframe)=c('communities','genes');
                             nam=merge(nam,comms_Dataframe,by="genes");
                             print(head(nam));
                             kalcom <- paste("Kal_CoDeregulated_Community_MoreStringent_", i, sep = "")
                             assign(kalcom,nam);
                             print(head(kalcom));
                             print(paste('finished for Kal_CoDeregulated_Community_MoreStringent_',i,sep=''))}
                              





#Calcolo le sottocommunities utilizzando una maggior stringenza
for (i in 1:11){             nam=kallisto_comms[kallisto_comms$communities==i,];
                             print(head(nam));
                             res2allgenes <- rcorr(as.matrix(t(nam[,38:52])));
                             cor.matrix=res2allgenes$r;
                             colnames(cor.matrix)=nam$genes;
                             rownames(cor.matrix)=nam$genes;
                             cor.matrix_rows=apply(cor.matrix, 1, percent_rank);
                             cor.matrix_cols=apply(cor.matrix, 2, percent_rank);
                             cor.matrix_sqrt=sqrt(cor.matrix_rows*cor.matrix_cols);
                             library(reshape2);
                             my_cor_df <- melt(cor.matrix_sqrt);
                             my_cor_df= my_cor_df %>% filter(Var1 != Var2);
                             my_adj_list <- my_cor_df %>% filter(value > 0.999);
                             names(my_adj_list) <- c('from', 'to', 'weight') ;
                             g=graph.data.frame(my_adj_list, directed=F);
                             communities=cluster_louvain(g, weights = NULL);
                             comms_Dataframe=as.data.frame(as.matrix(membership(communities)));
                             comms_Dataframe$genes=rownames(comms_Dataframe);
                             colnames(comms_Dataframe)=c('communities','genes');
                             nam=merge(nam,comms_Dataframe,by="genes");
                             print(head(nam));
                             kalcom <- paste("Kal_CoDeregulated_Community_", i, sep = "")
                             assign(kalcom,nam);
                             print(head(kalcom));
                             print(paste('finished for Kal_CoDeregulated_Community_',i,sep=''))}
                              

 
 
 
 #Ora scrivo i files RNK per ogni sottocommunity
 
 for (i in 1:15){kal=Kal_CoDeregulated_Community_1[Kal_CoDeregulated_Community_1$communities.y==i,1:2];kal=kal[order(-kal$logFC),];write.table(kal[,1:2],paste('TopTags_Kallysto_15Patients_Community1_SubCommunity_',i,'.rnk',sep=''),sep='\t',col.names=F,row.names=F,quote=F)}

for (i in 1:13){kal=Kal_CoDeregulated_Community_2[Kal_CoDeregulated_Community_2$communities.y==i,1:2];kal=kal[order(-kal$logFC),];write.table(kal[,1:2],paste('TopTags_Kallysto_15Patients_Community2_SubCommunity_',i,'.rnk',sep=''),sep='\t',col.names=F,row.names=F,quote=F)}

for (i in 1:13){kal=Kal_CoDeregulated_Community_3[Kal_CoDeregulated_Community_3$communities.y==i,1:2];kal=kal[order(-kal$logFC),];write.table(kal[,1:2],paste('TopTags_Kallysto_15Patients_Community3_SubCommunity_',i,'.rnk',sep=''),sep='\t',col.names=F,row.names=F,quote=F)}


#La community 4 e' gia' molto eterogenea. Non c'e' bisogno di dividerla in piu' community

for (i in 1:19){kal=Kal_CoDeregulated_Community_5[Kal_CoDeregulated_Community_5$communities.y==i,1:2];kal=kal[order(-kal$logFC),];write.table(kal[,1:2],paste('TopTags_Kallysto_15Patients_Community5_SubCommunity_',i,'.rnk',sep=''),sep='\t',col.names=F,row.names=F,quote=F)}

for (i in 1:11){kal=Kal_CoDeregulated_Community_6[Kal_CoDeregulated_Community_6$communities.y==i,1:2];kal=kal[order(-kal$logFC),];write.table(kal[,1:2],paste('TopTags_Kallysto_15Patients_Community6_SubCommunity_',i,'.rnk',sep=''),sep='\t',col.names=F,row.names=F,quote=F)}
 
for (i in 1:12){kal=Kal_CoDeregulated_Community_7[Kal_CoDeregulated_Community_7$communities.y==i,1:2];kal=kal[order(-kal$logFC),];write.table(kal[,1:2],paste('TopTags_Kallysto_15Patients_Community7_SubCommunity_',i,'.rnk',sep=''),sep='\t',col.names=F,row.names=F,quote=F)}

for (i in 1:10){kal=Kal_CoDeregulated_Community_8[Kal_CoDeregulated_Community_8$communities.y==i,1:2];kal=kal[order(-kal$logFC),];write.table(kal[,1:2],paste('TopTags_Kallysto_15Patients_Community8_SubCommunity_',i,'.rnk',sep=''),sep='\t',col.names=F,row.names=F,quote=F)}

for (i in 1:16){kal=Kal_CoDeregulated_Community_9[Kal_CoDeregulated_Community_9$communities.y==i,1:2];kal=kal[order(-kal$logFC),];write.table(kal[,1:2],paste('TopTags_Kallysto_15Patients_Community9_SubCommunity_',i,'.rnk',sep=''),sep='\t',col.names=F,row.names=F,quote=F)}

for (i in 1:12){kal=Kal_CoDeregulated_Community_10[Kal_CoDeregulated_Community_10$communities.y==i,1:2];kal=kal[order(-kal$logFC),];write.table(kal[,1:2],paste('TopTags_Kallysto_15Patients_Community10_SubCommunity_',i,'.rnk',sep=''),sep='\t',col.names=F,row.names=F,quote=F)}

for (i in 1:13){kal=Kal_CoDeregulated_Community_11[Kal_CoDeregulated_Community_11$communities.y==i,1:2];kal=kal[order(-kal$logFC),];write.table(kal[,1:2],paste('TopTags_Kallysto_15Patients_Community11_SubCommunity_',i,'.rnk',sep=''),sep='\t',col.names=F,row.names=F,quote=F)}




nam=kallisto_comms[kallisto_comms$communities==10,];
print(head(nam));
res2allgenes <- rcorr(as.matrix(t(nam[,38:52])))
cor.matrix=res2allgenes$r
colnames(cor.matrix)=nam$genes
rownames(cor.matrix)=nam$genes
cor.matrix_rows=apply(cor.matrix, 1, percent_rank)
cor.matrix_cols=apply(cor.matrix, 2, percent_rank)
cor.matrix_sqrt=sqrt(cor.matrix_rows*cor.matrix_cols)
library(reshape2)
my_cor_df <- melt(cor.matrix_sqrt)
my_cor_df= my_cor_df %>% filter(Var1 != Var2)
my_adj_list <- my_cor_df %>% filter(value > 0.999)
names(my_adj_list) <- c('from', 'to', 'weight')
g=graph.data.frame(my_adj_list, directed=F)
communities=cluster_louvain(g, weights = NULL)
comms_Dataframe=as.data.frame(as.matrix(membership(communities)))
comms_Dataframe$genes=rownames(comms_Dataframe)
colnames(comms_Dataframe)=c('communities','genes')
nam=merge(nam,comms_Dataframe,by="genes")
print(head(nam));
kalcom <- paste("Kal_CoDeregulated_Community_", i, sep = "")
assign(kalcom,nam);
print(head(kalcom));
print(paste('finished for Kal_CoDeregulated_Community_',i,sep=''))

library(biomaRt)
library(org.Hs.eg.db)
library(GO.db)
library(DOSE)
mart=useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="hsapiens_gene_ensembl")



RNAseq_all_genes<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=toptags_xypat_Kallisto$genes,mart=mart)


RNAseq_HLA_genes<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=Kal_CoDeregulated_Community_10[Kal_CoDeregulated_Community_10$communities.y==1,1],mart=mart)



#Gene Ontology Comunita 1
#Seleziono solo quelle sottocommunities>20
kom_table=as.data.frame(table(Kal_CoDeregulated_Community_1$communities.y))
colnames(kom_table)=c('communities.y','com_counts')
Kal_CoDeregulated_Community_1_filt=plyr::join(Kal_CoDeregulated_Community_1,kom_table,by='communities.y')
Kal_CoDeregulated_Community_1_filt=Kal_CoDeregulated_Community_1_filt[Kal_CoDeregulated_Community_1_filt$com_counts>=20,]
for (i in unique(Kal_CoDeregulated_Community_1_filt$communities.y)) {
                RNA_Subcommunity_Genes<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=Kal_CoDeregulated_Community_1_filt[Kal_CoDeregulated_Community_1_filt$communities.y==i,1],mart=mart);
                
                Kallisto_rna_subcommunity <- enrichGO(gene = RNA_Subcommunity_Genes$entrezgene,
                universe      = as.character(RNAseq_all_genes$entrezgene),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 0.05,readable = T)
                GO_kal <- paste("EnrichGO_GObp_Community1_SubCommunity_", i, sep = "")
                print(GO_kal)
                assign(GO_kal,Kallisto_rna_subcommunity)}
                 
for (i in unique(Kal_CoDeregulated_Community_1_filt$communities.y)) {pdf(paste('Pheatmap_Kal_CoDeregulated_Community_1_filt_',i,'.pdf',sep=''));a=Kal_CoDeregulated_Community_1_filt[Kal_CoDeregulated_Community_1_filt$communities.y==i,];rownames(a)=a$genes;pheatmap( a[,38:52], cluster_rows=T,cluster_cols=T,color=colorRampPalette( c("cyan", "black", "yellow"), space="rgb")(length(breakList)),breaks = breakList,main = paste('Community 1 subCom',i,sep=' '),border_color = NA,fontsize_row = 6);dev.off() }


[1] "EnrichGO_GObp_Community1_SubCommunity_262"
[1] "EnrichGO_GObp_Community1_SubCommunity_244"
write.table(as.data.frame(EnrichGO_GObp_Community1_SubCommunity_244),'EnrichGO_GObp_Community1_SubCommunity_244.txt',sep='\t',col.names=T,row.names=F,quote=F)
[1] "EnrichGO_GObp_Community1_SubCommunity_219"
write.table(as.data.frame(EnrichGO_GObp_Community1_SubCommunity_219),'EnrichGO_GObp_Community1_SubCommunity_219.txt',sep='\t',col.names=T,row.names=F,quote=F)
[1] "EnrichGO_GObp_Community1_SubCommunity_171"
[1] "EnrichGO_GObp_Community1_SubCommunity_248"
[1] "EnrichGO_GObp_Community1_SubCommunity_233"
[1] "EnrichGO_GObp_Community1_SubCommunity_215"
[1] "EnrichGO_GObp_Community1_SubCommunity_209"


for (i in unique(Kal_CoDeregulated_Community_1_filt$communities.y)) {pdf(paste('Pheatmap_Kal_CoDeregulated_Community_1_filt_',i,'.pdf',sep=''));a=Kal_CoDeregulated_Community_1_filt[Kal_CoDeregulated_Community_1_filt$communities.y==i,];rownames(a)=a$genes;pheatmap( a[,38:52], cluster_rows=T,cluster_cols=T,color=colorRampPalette( c("cyan", "black", "yellow"), space="rgb")(length(breakList)),breaks = breakList,main = paste('Community_1 subCom',i,sep=' '),border_color = NA,fontsize_row = 6);dev.off() }


#Gene Ontology Comunita 2
kom_table=as.data.frame(table(Kal_CoDeregulated_Community_2$communities.y))
colnames(kom_table)=c('communities.y','com_counts')
Kal_CoDeregulated_Community_2_filt=plyr::join(Kal_CoDeregulated_Community_2,kom_table,by='communities.y')
Kal_CoDeregulated_Community_2_filt=Kal_CoDeregulated_Community_2_filt[Kal_CoDeregulated_Community_2_filt$com_counts>=20,]
for (i in unique(Kal_CoDeregulated_Community_2_filt$communities.y)) {
                RNA_Subcommunity_Genes<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=Kal_CoDeregulated_Community_2_filt[Kal_CoDeregulated_Community_2_filt$communities.y==i,1],mart=mart);
                
                Kallisto_rna_subcommunity <- enrichGO(gene = RNA_Subcommunity_Genes$entrezgene,
                universe      = as.character(RNAseq_all_genes$entrezgene),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 0.05,readable = T)
                GO_kal <- paste("EnrichGO_GObp_Community2_SubCommunity_", i, sep = "")
                print(GO_kal)
                assign(GO_kal,Kallisto_rna_subcommunity)}

[1] "EnrichGO_GObp_Community2_SubCommunity_286"
write.table(as.data.frame(EnrichGO_GObp_Community2_SubCommunity_286),'EnrichGO_GObp_Community2_SubCommunity_286.txt',sep='\t',col.names=T,row.names=F,quote=F)
[1] "EnrichGO_GObp_Community2_SubCommunity_235"
[1] "EnrichGO_GObp_Community2_SubCommunity_258"
[1] "EnrichGO_GObp_Community2_SubCommunity_263"
[1] "EnrichGO_GObp_Community2_SubCommunity_249"
[1] "EnrichGO_GObp_Community2_SubCommunity_238"
[1] "EnrichGO_GObp_Community2_SubCommunity_265"
[1] "EnrichGO_GObp_Community2_SubCommunity_261"
[1] "EnrichGO_GObp_Community2_SubCommunity_267"
[1] "EnrichGO_GObp_Community2_SubCommunity_280"

for (i in unique(Kal_CoDeregulated_Community_2_filt$communities.y)) {pdf(paste('Pheatmap_Kal_CoDeregulated_Community_2_filt_',i,'.pdf',sep=''));a=Kal_CoDeregulated_Community_2_filt[Kal_CoDeregulated_Community_2_filt$communities.y==i,];rownames(a)=a$genes;pheatmap( a[,38:52], cluster_rows=T,cluster_cols=T,color=colorRampPalette( c("cyan", "black", "yellow"), space="rgb")(length(breakList)),breaks = breakList,main = paste('Community_2 subCom',i,sep=' '),border_color = NA,fontsize_row = 6);dev.off() }



#Gene Ontology Comunita 3
kom_table=as.data.frame(table(Kal_CoDeregulated_Community_3$communities.y))
colnames(kom_table)=c('communities.y','com_counts')
Kal_CoDeregulated_Community_3_filt=plyr::join(Kal_CoDeregulated_Community_3,kom_table,by='communities.y')
Kal_CoDeregulated_Community_3_filt=Kal_CoDeregulated_Community_3_filt[Kal_CoDeregulated_Community_3_filt$com_counts>=10,]
for (i in unique(Kal_CoDeregulated_Community_3_filt$communities.y)) {
                RNA_Subcommunity_Genes<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=Kal_CoDeregulated_Community_3_filt[Kal_CoDeregulated_Community_3_filt$communities.y==i,1],mart=mart);
                
                Kallisto_rna_subcommunity <- enrichGO(gene = RNA_Subcommunity_Genes$entrezgene,
                universe      = as.character(RNAseq_all_genes$entrezgene),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 0.05,readable = T)
                GO_kal <- paste("EnrichGO_GObp_Community3_SubCommunity_", i, sep = "")
                print(GO_kal)
                assign(GO_kal,Kallisto_rna_subcommunity)}



head(EnrichGO_GObp_Community3_SubCommunity_191)
head(EnrichGO_GObp_Community3_SubCommunity_240)
head(EnrichGO_GObp_Community3_SubCommunity_301)
head(EnrichGO_GObp_Community3_SubCommunity_252)
head(EnrichGO_GObp_Community3_SubCommunity_309)
head(EnrichGO_GObp_Community3_SubCommunity_313)
head(EnrichGO_GObp_Community3_SubCommunity_320)
head(EnrichGO_GObp_Community3_SubCommunity_262)
head(EnrichGO_GObp_Community3_SubCommunity_271)
head(EnrichGO_GObp_Community3_SubCommunity_330)
head(EnrichGO_GObp_Community3_SubCommunity_291)
head(EnrichGO_GObp_Community3_SubCommunity_331)
head(EnrichGO_GObp_Community3_SubCommunity_279)
head(EnrichGO_GObp_Community3_SubCommunity_306)
head(EnrichGO_GObp_Community3_SubCommunity_212)
head(EnrichGO_GObp_Community3_SubCommunity_207)
head(EnrichGO_GObp_Community3_SubCommunity_286)
head(EnrichGO_GObp_Community3_SubCommunity_296)
head(EnrichGO_GObp_Community3_SubCommunity_316)
head(EnrichGO_GObp_Community3_SubCommunity_273)
head(EnrichGO_GObp_Community3_SubCommunity_326)
head(EnrichGO_GObp_Community3_SubCommunity_246)
head(EnrichGO_GObp_Community3_SubCommunity_205)
head(EnrichGO_GObp_Community3_SubCommunity_278)
head(EnrichGO_GObp_Community3_SubCommunity_281)
head(EnrichGO_GObp_Community3_SubCommunity_299)
head(EnrichGO_GObp_Community3_SubCommunity_241)
head(EnrichGO_GObp_Community3_SubCommunity_314)
head(EnrichGO_GObp_Community3_SubCommunity_239)
head(EnrichGO_GObp_Community3_SubCommunity_315)
head(EnrichGO_GObp_Community3_SubCommunity_263)
head(EnrichGO_GObp_Community3_SubCommunity_277)
head(EnrichGO_GObp_Community3_SubCommunity_268)
head(EnrichGO_GObp_Community3_SubCommunity_226)
head(EnrichGO_GObp_Community3_SubCommunity_334)


write.table(as.data.frame(EnrichGO_GObp_Community3_SubCommunity_240),'EnrichGO_GObp_Community3_SubCommunity_240.txt',sep='\t',col.names=T,row.names=F,quote=F)
write.table(as.data.frame(EnrichGO_GObp_Community3_SubCommunity_301),'EnrichGO_GObp_Community3_SubCommunity_301.txt',sep='\t',col.names=T,row.names=F,quote=F)
write.table(as.data.frame(EnrichGO_GObp_Community3_SubCommunity_252),'EnrichGO_GObp_Community3_SubCommunity_252.txt',sep='\t',col.names=T,row.names=F,quote=F)
write.table(as.data.frame(EnrichGO_GObp_Community3_SubCommunity_331),'EnrichGO_GObp_Community3_SubCommunity_331.txt',sep='\t',col.names=T,row.names=F,quote=F)
write.table(as.data.frame(EnrichGO_GObp_Community3_SubCommunity_301),'EnrichGO_GObp_Community3_SubCommunity_301.txt',sep='\t',col.names=T,row.names=F,quote=F)
write.table(as.data.frame(EnrichGO_GObp_Community3_SubCommunity_301),'EnrichGO_GObp_Community3_SubCommunity_301.txt',sep='\t',col.names=T,row.names=F,quote=F)
write.table(as.data.frame(EnrichGO_GObp_Community3_SubCommunity_301),'EnrichGO_GObp_Community3_SubCommunity_301.txt',sep='\t',col.names=T,row.names=F,quote=F)
write.table(as.data.frame(EnrichGO_GObp_Community3_SubCommunity_301),'EnrichGO_GObp_Community3_SubCommunity_301.txt',sep='\t',col.names=T,row.names=F,quote=F)


pdf('Community3_SubCommunities_0999Corr_10MinNumberGenes.pdf')
breakList=seq(-8, 8, by = 0.1)
for (i in unique(Kal_CoDeregulated_Community_3_filt$communities.y)) {a=Kal_CoDeregulated_Community_3_filt[Kal_CoDeregulated_Community_3_filt$communities.y==i,];rownames(a)=a$genes;pheatmap( a[,38:52], cluster_rows=T,cluster_cols=T,color=colorRampPalette( c("cyan", "black", "yellow"), space="rgb")(length(breakList)),breaks = breakList,main = paste('Community 3 subCom',i,sep=' '),border_color = NA,fontsize_row = 6)}
dev.off()



for (i in unique(Kal_CoDeregulated_Community_3_filt$communities.y)) {pdf(paste('Pheatmap_Kal_CoDeregulated_Community_3_filt_',i,'.pdf',sep=''));a=Kal_CoDeregulated_Community_3_filt[Kal_CoDeregulated_Community_3_filt$communities.y==i,];rownames(a)=a$genes;pheatmap( a[,38:52], cluster_rows=T,cluster_cols=T,color=colorRampPalette( c("cyan", "black", "yellow"), space="rgb")(length(breakList)),breaks = breakList,main = paste('Community_3 subCom',i,sep=' '),border_color = NA,fontsize_row = 6);dev.off() }

#Gene Ontology Comunita 4
#La community 4 non ha sottocommunities

#Gene Ontology Comunita 3
kom_table=as.data.frame(table(Kal_CoDeregulated_Community_3$communities.y))
colnames(kom_table)=c('communities.y','com_counts')
Kal_CoDeregulated_Community_4_filt=plyr::join(Kal_CoDeregulated_Community_4,kom_table,by='communities.y')
Kal_CoDeregulated_Community_4_filt=Kal_CoDeregulated_Community_4_filt[Kal_CoDeregulated_Community_4_filt$com_counts>=10,]
for (i in unique(Kal_CoDeregulated_Community_4_filt$communities.y)) {
                RNA_Subcommunity_Genes<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=Kal_CoDeregulated_Community_4_filt[Kal_CoDeregulated_Community_4_filt$communities.y==i,1],mart=mart);
                
                Kallisto_rna_subcommunity <- enrichGO(gene = RNA_Subcommunity_Genes$entrezgene,
                universe      = as.character(RNAseq_all_genes$entrezgene),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 0.05,readable = T)
                GO_kal <- paste("EnrichGO_GObp_Community4_SubCommunity_", i, sep = "")
                print(GO_kal)
                assign(GO_kal,Kallisto_rna_subcommunity)}

                
#Gene Ontology Comunita 5
#La community 5 non ha sottocommunities 

#Gene Ontology Comunita 6
kom_table=as.data.frame(table(Kal_CoDeregulated_Community_6$communities.y))
colnames(kom_table)=c('communities.y','com_counts')
Kal_CoDeregulated_Community_6_filt=plyr::join(Kal_CoDeregulated_Community_6,kom_table,by='communities.y')
Kal_CoDeregulated_Community_6_filt=Kal_CoDeregulated_Community_6_filt[Kal_CoDeregulated_Community_6_filt$com_counts>=10,]
for (i in unique(Kal_CoDeregulated_Community_6_filt$communities.y)) {
                RNA_Subcommunity_Genes<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=Kal_CoDeregulated_Community_6_filt[Kal_CoDeregulated_Community_6_filt$communities.y==i,1],mart=mart);
                
                Kallisto_rna_subcommunity <- enrichGO(gene = RNA_Subcommunity_Genes$entrezgene,
                universe      = as.character(RNAseq_all_genes$entrezgene),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 0.05,readable = T)
                GO_kal <- paste("EnrichGO_GObp_Community6_SubCommunity_", i, sep = "")
                print(GO_kal)
                assign(GO_kal,Kallisto_rna_subcommunity)}
[1] "EnrichGO_GObp_Community6_SubCommunity_394"

pdf('Community6_SubCommunities_0999Corr_10MinNumberGenes.pdf')
breakList=seq(-8, 8, by = 0.1)
for (i in unique(Kal_CoDeregulated_Community_6_filt$communities.y)) {a=Kal_CoDeregulated_Community_6_filt[Kal_CoDeregulated_Community_6_filt$communities.y==i,];rownames(a)=a$genes;pheatmap( a[,38:52], cluster_rows=T,cluster_cols=T,color=colorRampPalette( c("cyan", "black", "yellow"), space="rgb")(length(breakList)),breaks = breakList,main = paste('Community 10 subCom',i,sep=' '),border_color = NA,fontsize_row = 6)}
dev.off()


for (i in unique(Kal_CoDeregulated_Community_6_filt$communities.y)) {pdf(paste('Pheatmap_Kal_CoDeregulated_Community_6_filt_',i,'.pdf',sep=''));a=Kal_CoDeregulated_Community_6_filt[Kal_CoDeregulated_Community_6_filt$communities.y==i,];rownames(a)=a$genes;pheatmap( a[,38:52], cluster_rows=T,cluster_cols=T,color=colorRampPalette( c("cyan", "black", "yellow"), space="rgb")(length(breakList)),breaks = breakList,main = paste('Community_6 subCom',i,sep=' '),border_color = NA,fontsize_row = 6);dev.off() }


#Gene Ontology Comunita 7
kom_table=as.data.frame(table(Kal_CoDeregulated_Community_7$communities.y))
colnames(kom_table)=c('communities.y','com_counts')
Kal_CoDeregulated_Community_7_filt=plyr::join(Kal_CoDeregulated_Community_7,kom_table,by='communities.y')
Kal_CoDeregulated_Community_7_filt=Kal_CoDeregulated_Community_7_filt[Kal_CoDeregulated_Community_7_filt$com_counts>=10,]
for (i in unique(Kal_CoDeregulated_Community_7_filt$communities.y)) {
                RNA_Subcommunity_Genes<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=Kal_CoDeregulated_Community_7_filt[Kal_CoDeregulated_Community_7_filt$communities.y==i,1],mart=mart);
                
                Kallisto_rna_subcommunity <- enrichGO(gene = RNA_Subcommunity_Genes$entrezgene,
                universe      = as.character(RNAseq_all_genes$entrezgene),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 0.05,readable = T)
                GO_kal <- paste("EnrichGO_GObp_Community7_SubCommunity_", i, sep = "")
                print(GO_kal)
                assign(GO_kal,Kallisto_rna_subcommunity)}

for (i in unique(Kal_CoDeregulated_Community_7_filt$communities.y)) {pdf(paste('Pheatmap_Kal_CoDeregulated_Community_7_filt_',i,'.pdf',sep=''));a=Kal_CoDeregulated_Community_7_filt[Kal_CoDeregulated_Community_7_filt$communities.y==i,];rownames(a)=a$genes;pheatmap( a[,38:52], cluster_rows=T,cluster_cols=T,color=colorRampPalette( c("cyan", "black", "yellow"), space="rgb")(length(breakList)),breaks = breakList,main = paste('Community 11 subCom',i,sep=' '),border_color = NA,fontsize_row = 6);dev.off() }



#Gene Ontology Comunita 8
kom_table=as.data.frame(table(Kal_CoDeregulated_Community_8$communities.y))
colnames(kom_table)=c('communities.y','com_counts')
Kal_CoDeregulated_Community_8_filt=plyr::join(Kal_CoDeregulated_Community_8,kom_table,by='communities.y')
Kal_CoDeregulated_Community_8_filt=Kal_CoDeregulated_Community_8_filt[Kal_CoDeregulated_Community_8_filt$com_counts>=20,]
for (i in unique(Kal_CoDeregulated_Community_8_filt$communities.y)) {
                RNA_Subcommunity_Genes<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=Kal_CoDeregulated_Community_8_filt[Kal_CoDeregulated_Community_8_filt$communities.y==i,1],mart=mart);
                
                Kallisto_rna_subcommunity <- enrichGO(gene = RNA_Subcommunity_Genes$entrezgene,
                universe      = as.character(RNAseq_all_genes$entrezgene),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 0.05,readable = T)
                GO_kal <- paste("EnrichGO_GObp_Community8_SubCommunity_", i, sep = "")
                print(GO_kal)
                assign(GO_kal,Kallisto_rna_subcommunity)}


[1] "EnrichGO_GObp_Community8_SubCommunity_258"
[1] "EnrichGO_GObp_Community8_SubCommunity_321"
[1] "EnrichGO_GObp_Community8_SubCommunity_330"
write.table(as.data.frame(EnrichGO_GObp_Community8_SubCommunity_330),'EnrichGO_GObp_Community8_SubCommunity_330.txt',sep='\t',col.names=T,row.names=F,quote=F)

[1] "EnrichGO_GObp_Community8_SubCommunity_333"
[1] "EnrichGO_GObp_Community8_SubCommunity_319"
[1] "EnrichGO_GObp_Community7_SubCommunity_278"
write.table(as.data.frame(EnrichGO_GObp_Community7_SubCommunity_278),'EnrichGO_GObp_Community7_SubCommunity_278.txt',sep='\t',col.names=T,row.names=F,quote=F)

[1] "EnrichGO_GObp_Community7_SubCommunity_293"
[1] "EnrichGO_GObp_Community7_SubCommunity_298"
 
for (i in unique(Kal_CoDeregulated_Community_8_filt$communities.y)) {pdf(paste('Pheatmap_Kal_CoDeregulated_Community_8_filt_',i,'.pdf',sep=''));a=Kal_CoDeregulated_Community_8_filt[Kal_CoDeregulated_Community_8_filt$communities.y==i,];rownames(a)=a$genes;pheatmap( a[,38:52], cluster_rows=T,cluster_cols=T,color=colorRampPalette( c("cyan", "black", "yellow"), space="rgb")(length(breakList)),breaks = breakList,main = paste('Community_8 subCom',i,sep=' '),border_color = NA,fontsize_row = 6);dev.off() }


#Gene Ontology Comunita 9
kom_table=as.data.frame(table(Kal_CoDeregulated_Community_9$communities.y))
colnames(kom_table)=c('communities.y','com_counts')
Kal_CoDeregulated_Community_9_filt=plyr::join(Kal_CoDeregulated_Community_9,kom_table,by='communities.y')
Kal_CoDeregulated_Community_9_filt=Kal_CoDeregulated_Community_9_filt[Kal_CoDeregulated_Community_9_filt$com_counts>=20,]
for (i in unique(Kal_CoDeregulated_Community_9_filt$communities.y)) {
                RNA_Subcommunity_Genes<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=Kal_CoDeregulated_Community_9_filt[Kal_CoDeregulated_Community_9_filt$communities.y==i,1],mart=mart);
                
                Kallisto_rna_subcommunity <- enrichGO(gene = RNA_Subcommunity_Genes$entrezgene,
                universe      = as.character(RNAseq_all_genes$entrezgene),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 0.05,readable = T)
                GO_kal <- paste("EnrichGO_GObp_Community9_SubCommunity_", i, sep = "")
                print(GO_kal)
                assign(GO_kal,Kallisto_rna_subcommunity)}




#Gene Ontology Comunita 10
kom_table=as.data.frame(table(Kal_CoDeregulated_Community_10$communities.y))
colnames(kom_table)=c('communities.y','com_counts')
Kal_CoDeregulated_Community_10_filt=plyr::join(Kal_CoDeregulated_Community_10,kom_table,by='communities.y')
Kal_CoDeregulated_Community_10_filt=Kal_CoDeregulated_Community_10_filt[Kal_CoDeregulated_Community_10_filt$com_counts>=10,]
for (i in unique(Kal_CoDeregulated_Community_10_filt$communities.y)) {
                RNA_Subcommunity_Genes<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=Kal_CoDeregulated_Community_10_filt[Kal_CoDeregulated_Community_10_filt$communities.y==i,1],mart=mart)
                Kallisto_rna_subcommunity <- enrichGO(gene = RNA_Subcommunity_Genes$entrezgene,
                universe      = as.character(RNAseq_all_genes$entrezgene),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 0.05,readable = T)
                GO_kal <- paste("EnrichGO_GObp_Community10_SubCommunity_", i, sep = "")
                print(GO_kal)
                assign(GO_kal,Kallisto_rna_subcommunity)}

pdf('Community10_SubCommunities_0999Corr_10MinNumberGenes.pdf')
breakList=seq(-8, 8, by = 0.1)
for (i in unique(Kal_CoDeregulated_Community_10_filt$communities.y)) {a=Kal_CoDeregulated_Community_10_filt[Kal_CoDeregulated_Community_10_filt$communities.y==i,];rownames(a)=a$genes;pheatmap( a[,38:52], cluster_rows=T,cluster_cols=T,color=colorRampPalette( c("cyan", "black", "yellow"), space="rgb")(length(breakList)),breaks = breakList,main = paste('Community 10 subCom',i,sep=' '),border_color = NA,fontsize_row = 6)}
dev.off()

head(EnrichGO_GObp_Community10_SubCommunity_26)
head(EnrichGO_GObp_Community10_SubCommunity_15)
head(EnrichGO_GObp_Community10_SubCommunity_21)
head(EnrichGO_GObp_Community10_SubCommunity_12)
head(EnrichGO_GObp_Community10_SubCommunity_17)
head(EnrichGO_GObp_Community10_SubCommunity_11)
head(EnrichGO_GObp_Community10_SubCommunity_10)
head(EnrichGO_GObp_Community10_SubCommunity_22)
head(EnrichGO_GObp_Community10_SubCommunity_9)
head(EnrichGO_GObp_Community10_SubCommunity_14)
head(EnrichGO_GObp_Community10_SubCommunity_16)
head(EnrichGO_GObp_Community10_SubCommunity_2)
head(EnrichGO_GObp_Community10_SubCommunity_31)
head(EnrichGO_GObp_Community10_SubCommunity_19)
head(EnrichGO_GObp_Community10_SubCommunity_7)
head(EnrichGO_GObp_Community10_SubCommunity_37)
head(EnrichGO_GObp_Community10_SubCommunity_13)
head(EnrichGO_GObp_Community10_SubCommunity_24)
head(EnrichGO_GObp_Community10_SubCommunity_41)
head(EnrichGO_GObp_Community10_SubCommunity_8)
head(EnrichGO_GObp_Community10_SubCommunity_34)
head(EnrichGO_GObp_Community10_SubCommunity_20)
head(EnrichGO_GObp_Community10_SubCommunity_36)
head(EnrichGO_GObp_Community10_SubCommunity_33)
head(EnrichGO_GObp_Community10_SubCommunity_5)
head(EnrichGO_GObp_Community10_SubCommunity_42)
head(EnrichGO_GObp_Community10_SubCommunity_28)
head(EnrichGO_GObp_Community10_SubCommunity_32)


write.table(as.data.frame(EnrichGO_GObp_Community10_SubCommunity_26),'EnrichGO_GObp_Community10_SubCommunity_26.txt',sep='\t',col.names=T,row.names=F,quote=F)
write.table(as.data.frame(EnrichGO_GObp_Community10_SubCommunity_9),'EnrichGO_GObp_Community10_SubCommunity_9.txt',sep='\t',col.names=T,row.names=F,quote=F)
write.table(as.data.frame(EnrichGO_GObp_Community10_SubCommunity_2),'EnrichGO_GObp_Community10_SubCommunity_2.txt',sep='\t',col.names=T,row.names=F,quote=F)
write.table(as.data.frame(EnrichGO_GObp_Community10_SubCommunity_7),'EnrichGO_GObp_Community10_SubCommunity_7.txt',sep='\t',col.names=T,row.names=F,quote=F)
write.table(as.data.frame(EnrichGO_GObp_Community10_SubCommunity_42),'EnrichGO_GObp_Community10_SubCommunity_42.txt',sep='\t',col.names=T,row.names=F,quote=F)
write.table(as.data.frame(EnrichGO_GObp_Community10_SubCommunity_28),'EnrichGO_GObp_Community10_SubCommunity_28.txt',sep='\t',col.names=T,row.names=F,quote=F)


pdf('Community10_SubCommunities_0999Corr_10MinNumberGenes.pdf')
for (i in unique(Kal_CoDeregulated_Community_10_filt$communities.y)) {a=Kal_CoDeregulated_Community_10_filt[Kal_CoDeregulated_Community_10_filt$communities.y==i,];rownames(a)=a$genes;pheatmap( a[,38:52], cluster_rows=T,cluster_cols=T,color=colorRampPalette( c("cyan", "black", "yellow"), space="rgb")(length(breakList)),breaks = breakList,main = paste('Community 10 subCom',i,sep=' '),border_color = NA,fontsize_row = 6)}


#Gene Ontology Comunita 11
kom_table=as.data.frame(table(Kal_CoDeregulated_Community_11$communities.y))
colnames(kom_table)=c('communities.y','com_counts')
Kal_CoDeregulated_Community_11_filt=plyr::join(Kal_CoDeregulated_Community_11,kom_table,by='communities.y')
Kal_CoDeregulated_Community_11_filt=Kal_CoDeregulated_Community_11_filt[Kal_CoDeregulated_Community_11_filt$com_counts>=10,]
for (i in unique(Kal_CoDeregulated_Community_11_filt$communities.y)) {
                RNA_Subcommunity_Genes<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=Kal_CoDeregulated_Community_11_filt[Kal_CoDeregulated_Community_11_filt$communities.y==i,1],mart=mart);
                
                Kallisto_rna_subcommunity <- enrichGO(gene = RNA_Subcommunity_Genes$entrezgene,
                universe      = as.character(RNAseq_all_genes$entrezgene),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 0.05,readable = T)
                GO_kal <- paste("EnrichGO_GObp_Community11_SubCommunity_", i, sep = "")
                print(GO_kal)
                assign(GO_kal,Kallisto_rna_subcommunity)}




head(EnrichGO_GObp_Community11_SubCommunity_299)
head(EnrichGO_GObp_Community11_SubCommunity_200)
head(EnrichGO_GObp_Community11_SubCommunity_255)
head(EnrichGO_GObp_Community11_SubCommunity_286)
head(EnrichGO_GObp_Community11_SubCommunity_301)
head(EnrichGO_GObp_Community11_SubCommunity_230)
head(EnrichGO_GObp_Community11_SubCommunity_190)
head(EnrichGO_GObp_Community11_SubCommunity_249)
head(EnrichGO_GObp_Community11_SubCommunity_307)
head(EnrichGO_GObp_Community11_SubCommunity_178)
head(EnrichGO_GObp_Community11_SubCommunity_245)
head(EnrichGO_GObp_Community11_SubCommunity_293)
head(EnrichGO_GObp_Community11_SubCommunity_279)
head(EnrichGO_GObp_Community11_SubCommunity_213)
head(EnrichGO_GObp_Community11_SubCommunity_277)
head(EnrichGO_GObp_Community11_SubCommunity_217)
head(EnrichGO_GObp_Community11_SubCommunity_259)
head(EnrichGO_GObp_Community11_SubCommunity_288)
head(EnrichGO_GObp_Community11_SubCommunity_219)
head(EnrichGO_GObp_Community11_SubCommunity_202)
head(EnrichGO_GObp_Community11_SubCommunity_274)
head(EnrichGO_GObp_Community11_SubCommunity_258)
head(EnrichGO_GObp_Community11_SubCommunity_256)
head(EnrichGO_GObp_Community11_SubCommunity_246)
head(EnrichGO_GObp_Community11_SubCommunity_232)
head(EnrichGO_GObp_Community11_SubCommunity_269)
head(EnrichGO_GObp_Community11_SubCommunity_236)
head(EnrichGO_GObp_Community11_SubCommunity_239)
head(EnrichGO_GObp_Community11_SubCommunity_266)


for (i in unique(Kal_CoDeregulated_Community_11_filt$communities.y)) {pdf(paste('Pheatmap_Kal_CoDeregulated_Community_11_filt_',i,'.pdf',sep=''));a=Kal_CoDeregulated_Community_11_filt[Kal_CoDeregulated_Community_11_filt$communities.y==i,];rownames(a)=a$genes;pheatmap( a[,38:52], cluster_rows=T,cluster_cols=T,color=colorRampPalette( c("cyan", "black", "yellow"), space="rgb")(length(breakList)),breaks = breakList,main = paste('Community 11 subCom',i,sep=' '),border_color = NA,fontsize_row = 6);dev.off() }

pdf('Community11_SubCommunities_0999Corr_10MinNumberGenes.pdf')
breakList=seq(-8, 8, by = 0.1)
for (i in unique(Kal_CoDeregulated_Community_11_filt$communities.y)) {a=Kal_CoDeregulated_Community_11_filt[Kal_CoDeregulated_Community_11_filt$communities.y==i,];rownames(a)=a$genes;pheatmap( a[,38:52], cluster_rows=T,cluster_cols=T,color=colorRampPalette( c("cyan", "black", "yellow"), space="rgb")(length(breakList)),breaks = breakList,main = paste('Community 11 subCom',i,sep=' '),border_color = NA,fontsize_row = 6)}
dev.off()









res2allgenes <- rcorr(as.matrix(t(Kal_CoDeregulated_Community_10[,38:52])))
cor.matrix=res2allgenes$r
colnames(cor.matrix)=Kal_CoDeregulated_Community_10$genes
rownames(cor.matrix)=Kal_CoDeregulated_Community_10$genes

cor.matrix_rows=apply(cor.matrix, 1, percent_rank)
cor.matrix_cols=apply(cor.matrix, 2, percent_rank)
cor.matrix_sqrt=sqrt(cor.matrix_rows*cor.matrix_cols)
library(reshape2)
my_cor_df <- melt(cor.matrix_sqrt)
my_cor_df= my_cor_df %>% filter(Var1 != Var2)
my_adj_list <- my_cor_df %>% filter(value > 0.995)
names(my_adj_list) <- c('from', 'to', 'weight') 
g=graph.data.frame(my_adj_list, directed=F)
communities=cluster_louvain(g, weights = NULL)
comms_Dataframe=as.data.frame(as.matrix(membership(communities)))
comms_Dataframe$genes=rownames(comms_Dataframe)
colnames(comms_Dataframe)=c('communities','genes')   
kallisto_comms_community10=merge(Kal_CoDeregulated_Community_10,comms_Dataframe,by="genes")






for (i in 1:12){kal=kallisto_comms_community10[kallisto_comms_community10$communities.y==i,1:2];kal=kal[order(-kal$logFC),];write.table(kal[,1:2],paste('TopTags_Kallysto_15Patients_Community10_SubCommunity_',i,'.rnk',sep=''),sep='\t',col.names=F,row.names=F,quote=F)}
 


#Network Analysis su Diagnosis logCPM

diagnosis_res2allgenes <- rcorr(as.matrix(t(toptags_xypat_Kallisto_cpm[,c(8,10,12,14,16,18,20,22,24,26,28,30,32,34,36)])))
diagnosis_cor.matrix=diagnosis_res2allgenes$r
colnames(diagnosis_cor.matrix)=toptags_xypat_Kallisto_cpm$genes
rownames(diagnosis_cor.matrix)=toptags_xypat_Kallisto_cpm$genes

diagnosis_cor.matrix_rows=apply(diagnosis_cor.matrix, 1, percent_rank)
diagnosis_cor.matrix_cols=apply(diagnosis_cor.matrix, 2, percent_rank)
diagnosis_cor.matrix_sqrt=sqrt(diagnosis_cor.matrix_rows*diagnosis_cor.matrix_cols)
library(reshape2)
diagnosis_my_cor_df <- melt(diagnosis_cor.matrix_sqrt)
diagnosis_my_cor_df = diagnosis_my_cor_df %>% filter(Var1 != Var2)
diagnosis_my_adj_list <- diagnosis_my_cor_df %>% filter(value > 0.995)
names(diagnosis_my_adj_list) <- c('from', 'to', 'weight') 
diagnosis_g=graph.data.frame(diagnosis_my_adj_list, directed=F)
saveRDS(diagnosis_g,'Grafo_Diagnosi_LuciaGabri.RDS')
diagnosis_communities=cluster_louvain(diagnosis_g,weights=NULL)
diagnosis_communities_genes=as.data.frame(as.matrix(membership(diagnosis_communities)))
diagnosis_communities_genes$genes=rownames(diagnosis_communities_genes)
diagnosis_communities_genes=diagnosis_communities_genes[,c(2,1)]
colnames(diagnosis_communities_genes)=c('genes','diagnosis_community')


#Network Analysis su Relapse logCPM

relapse_res2allgenes <- rcorr(as.matrix(t(toptags_xypat_Kallisto_cpm[,c(9,11,13,15,17,19,21,23,25,27,29,31,33,35,37)])))
relapse_cor.matrix=relapse_res2allgenes$r
colnames(relapse_cor.matrix)=toptags_xypat_Kallisto_cpm$genes
rownames(relapse_cor.matrix)=toptags_xypat_Kallisto_cpm$genes

relapse_cor.matrix_rows=apply(relapse_cor.matrix, 1, percent_rank)
relapse_cor.matrix_cols=apply(relapse_cor.matrix, 2, percent_rank)
relapse_cor.matrix_sqrt=sqrt(relapse_cor.matrix_rows*relapse_cor.matrix_cols)
library(reshape2)
relapse_my_cor_df <- melt(relapse_cor.matrix_sqrt)
relapse_my_cor_df = relapse_my_cor_df %>% filter(Var1 != Var2)
relapse_my_adj_list <- relapse_my_cor_df %>% filter(value > 0.995)
names(relapse_my_adj_list) <- c('from', 'to', 'weight') 
relapse_g=graph.data.frame(relapse_my_adj_list, directed=F)
saveRDS(relapse_g,'Grafo_Relapse_LuciaGabri.RDS')
relapse_communities=cluster_louvain(relapse_g,weights=NULL)
relapse_communities_genes=as.data.frame(as.matrix(membership(relapse_communities)))
relapse_communities_genes$genes=rownames(relapse_communities_genes)
relapse_communities_genes=relapse_communities_genes[,c(2,1)]
colnames(relapse_communities_genes)=c('genes','relapse_community')


kallisto_comms=toptags_xypat_Kallisto_cpm
kallisto_comms=merge(kallisto_comms,diagnosis_communities_genes,by='genes')
kallisto_comms=merge(kallisto_comms,relapse_communities_genes,by='genes')

table(kallisto_comms$diagnosis_community)

   1    2    3    4    5    6    7    8    9   10   11   12 
 838 1490 1818 1166 2067  807 1313 2036  205  768  348 3247 

table(kallisto_comms$relapse_community)

   1    2    3    4    5    6    7    8    9   10 
 416 1714 1152 1646 2187 2647 1821 1634 1997  889 

for(i in 1:12) { 
    nam <- paste("Kal_Diagnosis_", i, sep = "")
    assign(nam, kallisto_comms[kallisto_comms$diagnosis_community==i,1])
}

diag_list=list(Kal_Diagnosis_1,Kal_Diagnosis_2,Kal_Diagnosis_3,Kal_Diagnosis_4,Kal_Diagnosis_5,Kal_Diagnosis_6,Kal_Diagnosis_7,Kal_Diagnosis_8,Kal_Diagnosis_9,Kal_Diagnosis_10,Kal_Diagnosis_11,Kal_Diagnosis_12)


for(i in 1:10) { 
    nam <- paste("Kal_Relapse_", i, sep = "")
    assign(nam, kallisto_comms[kallisto_comms$relapse_community==i,1])
}

list=list(Kal_Diagnosis_1,Kal_Diagnosis_2,Kal_Diagnosis_3,Kal_Diagnosis_4,Kal_Diagnosis_5,Kal_Diagnosis_6,Kal_Diagnosis_7,Kal_Diagnosis_8,Kal_Diagnosis_9,Kal_Diagnosis_10,Kal_Diagnosis_11,Kal_Diagnosis_12,Kal_Relapse_1,Kal_Relapse_2,Kal_Relapse_3,Kal_Relapse_4,Kal_Relapse_5,Kal_Relapse_6,Kal_Relapse_7,Kal_Relapse_8,Kal_Relapse_9,Kal_Relapse_10)
    
names(list)=c(paste('Diagnosis_Community_',seq(1,12,by = 1)),paste('Relapse_Community_',seq(1,10,by = 1))) 

 dist <- unlist(lapply(combn(list, 2, simplify = FALSE), function(x) {
+     length(intersect(x[[1]], x[[2]]))/length(union(x[[1]], x[[2]])) }))
 JAccard_Index_Communities=as.data.frame(cbind(t(combn(22,2)), dist))

Unique_R7_notD5=setdiff(Kal_Relapse_7,Kal_Diagnosis_5)
Unique_D5_notR7=setdiff(Kal_Diagnosis_5,Kal_Relapse_7)


g2 <- induced_subgraph(relapse_g, V(relapse_g)[Unique_R7_notD5])

V(g2)$size <- 8
V(g2)$frame.color <- "white"
V(g2)$color <- "orange"
V(g2)$label <- ""
E(g2)$arrow.mode <- 0



rownames(as.matrix(V(relapse_g)))
V(relapse_g)$label=rownames(as.matrix(V(relapse_g)))
V(relapse_g)$label
prova=subgraph.edges(relapse_g, E(g)[label!=Unique_R7_notD5])
g <- make_ring(10)
V(g)
g2 <- make_ring(10) %>%
set_vertex_attr("name", value = letters[1:10])
V(g2)
V(relapse_g)$name



#
#res2allgenes <- rcorr(as.matrix(t(toptags_xypat_Kallisto_cpm[,38:52])))
#cor.matrix=res2allgenes$r
#colnames(cor.matrix)=toptags_xypat_Kallisto_cpm$genes
#rownames(cor.matrix)=toptags_xypat_Kallisto_cpm$genes
#
#cor.matrix_rows=apply(cor.matrix, 1, percent_rank)
#cor.matrix_cols=apply(cor.matrix, 2, percent_rank)
#cor.matrix_sqrt=sqrt(cor.matrix_rows*cor.matrix_cols)
#library(reshape2)
#my_cor_df <- melt(cor.matrix_sqrt)
#my_cor_df= my_cor_df %>% filter(Var1 != Var2)
#my_adj_list <- my_cor_df %>% filter(value > 0.995)
#names(my_adj_list) <- c('from', 'to', 'weight') 
#g=graph.data.frame(my_adj_list, directed=F)
#communities=cluster_louvain(g, weights = NULL)
#
#
#
#res2allgenes_DEG <- rcorr(as.matrix(t(toptags_xypat_Kallisto_cpm_DEG[,38:52])))
#cor.matrix_DEG=res2allgenes_DEG$r
#colnames(cor.matrix_DEG)=toptags_xypat_Kallisto_cpm_DEG$genes
#rownames(cor.matrix_DEG)=toptags_xypat_Kallisto_cpm_DEG$genes
#
#cor.matrix_rows_DEG=apply(cor.matrix_DEG, 1, percent_rank)
#cor.matrix_cols_DEG=apply(cor.matrix_DEG, 2, percent_rank)
#cor.matrix_sqrt_DEG=sqrt(cor.matrix_rows_DEG*cor.matrix_cols_DEG)
#library(reshape2)
#my_cor_df_DEG <- melt(cor.matrix_sqrt_DEG)
#my_cor_df_DEG= my_cor_df_DEG %>% filter(Var1 != Var2)
#my_adj_list_DEG <- my_cor_df_DEG %>% filter(value > 0.995)
#names(my_adj_list_DEG) <- c('from', 'to', 'weight') 
#g_DEG=graph.data.frame(my_adj_list_DEG, directed=F)
#communities_DEG=cluster_louvain(g_DEG, weights = NULL)
#
#


#
#
#
#library(dplyr)
#my_adj_list <- my_cor_df %>% filter(value > 0.995)
#names(my_adj_list) <- c('from', 'to', 'weight')
# g=graph.data.frame(my_adj_list)
# E(g)$weight <- my_adj_list$weight
#
#
#cor_sqr_df=as.data.frame(cor.matrix_sqrt) 
#
#
#res_g <- simplify(g, membership(communities),remove.loops=T) 
# res_communities=cluster_louvain(res_g, weights = NULL)
#
#
#
#vertex.color=rainbow(3, alpha=0.6)[guga$membership]
#
#layout = layout.auto
#
#plot(g, layout=layout,vertex.label=NA,vertex.size=2,vertex.color=rainbow(11, alpha=0.6)[communities$membership])
#
#cricche=cliques(g)



colnames(res2allgenes$P)=toptags_xypat_Kallisto_cpm$genes
rownames(res2allgenes$P)=toptags_xypat_Kallisto_cpm$genes
pval.matrix=res2allgenes$P


m=as.matrix(res2allgenes$r)

df=data.frame(genes1=rownames(m)[row(m)], genes2=colnames(m)[col(m)], corr=c(m))

m2=as.matrix(res2allgenes$P)

df2=data.frame(genes1=rownames(m2)[row(m2)], genes2=colnames(m2)[col(m2)], pval=c(m2))
df_cor=merge(df,df2,by=c('genes1','genes2'))




E(grafo)$weight <- abs(E(grafo)$weight)
grafo <- delete_edges(grafo, E(grafo)[which(E(grafo)$weight<0.8)])

grafo.communities <- edge.betweenness.community(grafo, directed=FALSE)
grafo.clustering <- make_clusters(grafo, membership=grafo.communities$membership)
commSummary <- data.frame(  grafo.communities$names,  grafo.communities$membership)



m=as.matrix(res2allgenes$r)

df=data.frame(row=rownames(m)[row(m)], col=colnames(m)[col(m)], corr=c(m))


g <- graph.adjacency(
  as.matrix(as.dist(cor(t(toptags_xypat_Kallisto_cpm[,38:52]), method="pearson"))),

  mode="undirected",
  weighted=TRUE,
  diag=FALSE
)

g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE)
E(g)$weight <- abs(E(g)$weight)
g <- delete_edges(g, E(g)[which(E(g)$weight<0.8)])
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
vSizes <- (scale01(apply(toptags_xypat_Kallisto_cpm[,38:52], 1, mean)) + 1.0) * 10

#Amplify or decrease the width of the edges
edgeweights <- E(g)$weight * 2.0

#Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
 mst <- mst(g, algorithm="prim")
 mst.communities <- edge.betweenness.community(mst, directed=FALSE)
 mst.clustering <- make_clusters(mst, membership=mst.communities$membership)
 V(mst)$color <- mst.communities$membership + 1

par(mfrow=c(1,2))
plot(
  mst.clustering, mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.size=vSizes,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="My first graph"
)

degree(mst)

#Output information for each community, including vertex-to-community assignments and modularity
commSummary2 <- data.frame(  mst.communities$names,  mst.communities$membership,mst.clustering$membership)
  
colnames(commSummary2) <- c("Gene", "Community","Clustering")
options(scipen=999)
commSummary

lovaino=cluster_louvain(grafo)
Geni_comunita=cbind(lovaino$names,lovaino$membership)
colnames(Geni_comunita)=c("genes","comunita")
toptags_xypat_Kallisto_cpm_comm=merge(toptags_xypat_Kallisto_cpm,Geni_comunita,by="genes")
toptags_xypat_Kallisto_cpm_comm$Significant=ifelse(toptags_xypat_Kallisto_cpm_comm$Qvalue<=0.25,"Significant","NonSignificant")

pdf("/home/fsantaniello/PCA_Colorata_Comunita_16kGenes_new.pdf")
cbPalette <- c('#e6194b','#3cb44b','#ffe119','#0082c8','#f58231','#911eb4','#46f0f0','#f032e6','#d2f53c','#fabebe','#008080','#e6beff',
'#aa6e28','#fffac8','#800000','#aaffc3','#808000')

toptags_xypat_Kallisto_cpm_comm=merge(toptags_xypat_Kallisto_cpm_DEG,comms,by="genes")

pca_Kallisto_prcomp_communities=prcomp(t(toptags_xypat_Kallisto_cpm_comm[,38:52]))
mds=as.data.frame(pca_Kallisto_prcomp_communities$x)
mds2=data.frame(x=mds$PC1,y=mds$PC2)
mds2=cbind(mds2,as.data.frame(toptags_xypat_Kallisto_cpm_comm$community))
colnames(mds2)=c("PC1","PC2","comunita")


ggplot(mds2, aes(PC1, PC2)) + geom_point(aes(colour =factor(comunita)))+scale_color_manual(values=cbPalette)










library(biomaRt)
library(org.Hs.eg.db)
library(GO.db)
library(DOSE)
mart=useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="hsapiens_gene_ensembl")


RNAseq_all_genes<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=toptags_xypat_Kallisto$genes,mart=mart)
RNAseq_DEG_genes<- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=toptags_xypat_Kallisto_DEG$genes,mart=mart)



Kallisto_rna <- enrichGO(gene = RNAseq_DEG_genes$entrezgene,
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





pdf("PCA_Conte_Pazienti_RNAseq_wDR1_NewGTF.pdf")
 mds=plotMDS(countdata_filtDGE_keep,dim.plot=c(1,2))
 mds2=data.frame(x=mds$x,y=mds$y)
samples=rownames(mds2)
 samples2=unlist(strsplit(samples,split="_kallisto"))
 samples2=samples2[samples2!=""]
 mds2=cbind(mds2,as.data.frame(samples2),as.data.frame(condition))
 patient=rep(1:15,each=2)
mds2=cbind(mds2,as.data.frame(patient))
 colnames(mds2)=c("x","y","samples","condition","patient")
 ggplot(mds2,aes(x,y))+geom_point(aes(shape=condition,col=as.character(patient))) + geom_text(label= gsub("_kallisto","",samples),size=3)
 

#ORa calcolo le logCPM per tutti i pazienti
toptags_xypat_Kallisto_cpm=cbind(cpm(countdata_filtDGE_keep,log=T),as.data.frame(countdata_filtDGE_keep$genes))
toptags_xypat_Kallisto_cpm=merge(toptags_xypat_Kallisto,toptags_xypat_Kallisto_cpm,by="genes")
toptags_xypat_Kallisto_cpm$AlfeFC=toptags_xypat_Kallisto_cpm$ALFE2_kallisto-toptags_xypat_Kallisto_cpm$ALFE1_kallisto
toptags_xypat_Kallisto_cpm$BesuFC=toptags_xypat_Kallisto_cpm$BESU2_kallisto-toptags_xypat_Kallisto_cpm$BESU1_kallisto
toptags_xypat_Kallisto_cpm$CaluFC=toptags_xypat_Kallisto_cpm$CALU2_kallisto-toptags_xypat_Kallisto_cpm$CALU1_kallisto
toptags_xypat_Kallisto_cpm$DeivFC=toptags_xypat_Kallisto_cpm$DEIV2b_kallisto-toptags_xypat_Kallisto_cpm$DEIV2_kallisto
toptags_xypat_Kallisto_cpm$DestFC=toptags_xypat_Kallisto_cpm$DEST2_kallisto-toptags_xypat_Kallisto_cpm$DEST1_kallisto
toptags_xypat_Kallisto_cpm$Dr1FC=toptags_xypat_Kallisto_cpm$DR1_2_kallisto-toptags_xypat_Kallisto_cpm$DR1_1_kallisto
toptags_xypat_Kallisto_cpm$Dr4FC=toptags_xypat_Kallisto_cpm$DR42_kallisto-toptags_xypat_Kallisto_cpm$DR41_kallisto
toptags_xypat_Kallisto_cpm$Dr5FC=toptags_xypat_Kallisto_cpm$DR52_kallisto-toptags_xypat_Kallisto_cpm$DR51_kallisto
toptags_xypat_Kallisto_cpm$FocaFC=toptags_xypat_Kallisto_cpm$FOCA2_kallisto-toptags_xypat_Kallisto_cpm$FOCA1_kallisto
toptags_xypat_Kallisto_cpm$GagraFC=toptags_xypat_Kallisto_cpm$GAGRA2_kallisto-toptags_xypat_Kallisto_cpm$GAGRA1_kallisto
toptags_xypat_Kallisto_cpm$LuanFC=toptags_xypat_Kallisto_cpm$LUAN2_kallisto-toptags_xypat_Kallisto_cpm$LUAN1_kallisto
toptags_xypat_Kallisto_cpm$MabiFC=toptags_xypat_Kallisto_cpm$MABI2_kallisto-toptags_xypat_Kallisto_cpm$MABI1_kallisto
toptags_xypat_Kallisto_cpm$MogeFC=toptags_xypat_Kallisto_cpm$MOGE2_kallisto-toptags_xypat_Kallisto_cpm$MOGE1_kallisto
toptags_xypat_Kallisto_cpm$PiagFC=toptags_xypat_Kallisto_cpm$PIAG2_kallisto-toptags_xypat_Kallisto_cpm$PIAG1_kallisto
toptags_xypat_Kallisto_cpm$PreluFC=toptags_xypat_Kallisto_cpm$PRELU2_kallisto-toptags_xypat_Kallisto_cpm$PRELU1_kallisto



#Ora su questi dati creo dei files RNK


toptags_xypat_Kallisto_cpm$AlfeRNK=   ifelse((toptags_xypat_Kallisto_cpm$ALFE2_kallisto  >= 0.58  | toptags_xypat_Kallisto_cpm$ALFE1_kallisto   >= 0.58 ) & ( toptags_xypat_Kallisto_cpm$AlfeFC   <= -0.58 | toptags_xypat_Kallisto_cpm$AlfeFC   >= 0.58),toptags_xypat_Kallisto_cpm$AlfeFC  ,NA) 
toptags_xypat_Kallisto_cpm$BesuRNK=   ifelse((toptags_xypat_Kallisto_cpm$BESU2_kallisto  >= 0.58  | toptags_xypat_Kallisto_cpm$BESU1_kallisto   >= 0.58 ) & ( toptags_xypat_Kallisto_cpm$BesuFC   <= -0.58 | toptags_xypat_Kallisto_cpm$BesuFC   >= 0.58),toptags_xypat_Kallisto_cpm$BesuFC  ,NA) 
toptags_xypat_Kallisto_cpm$CaluRNK=   ifelse((toptags_xypat_Kallisto_cpm$CALU2_kallisto  >= 0.58  | toptags_xypat_Kallisto_cpm$CALU1_kallisto   >= 0.58 ) & ( toptags_xypat_Kallisto_cpm$CaluFC   <= -0.58 | toptags_xypat_Kallisto_cpm$CaluFC   >= 0.58),toptags_xypat_Kallisto_cpm$CaluFC  ,NA) 
toptags_xypat_Kallisto_cpm$DeivRNK=   ifelse((toptags_xypat_Kallisto_cpm$DEIV2b_kallisto >= 0.58  | toptags_xypat_Kallisto_cpm$DEIV2_kallisto  >= 0.58 ) & (  toptags_xypat_Kallisto_cpm$DeivFC   <= -0.58 | toptags_xypat_Kallisto_cpm$DeivFC   >= 0.58),toptags_xypat_Kallisto_cpm$DeivFC  ,NA) 
toptags_xypat_Kallisto_cpm$DestRNK=   ifelse((toptags_xypat_Kallisto_cpm$DEST2_kallisto  >= 0.58  | toptags_xypat_Kallisto_cpm$DEST1_kallisto   >= 0.58 ) & ( toptags_xypat_Kallisto_cpm$DestFC   <= -0.58 | toptags_xypat_Kallisto_cpm$DestFC   >= 0.58),toptags_xypat_Kallisto_cpm$DestFC  ,NA) 
toptags_xypat_Kallisto_cpm$Dr1RNK=    ifelse((toptags_xypat_Kallisto_cpm$DR1_2_kallisto  >= 0.58  | toptags_xypat_Kallisto_cpm$DR1_1_kallisto   >= 0.58 ) & ( toptags_xypat_Kallisto_cpm$Dr1FC    <= -0.58 | toptags_xypat_Kallisto_cpm$Dr1FC    >= 0.58),toptags_xypat_Kallisto_cpm$Dr1FC   ,NA) 
toptags_xypat_Kallisto_cpm$Dr4RNK=    ifelse((toptags_xypat_Kallisto_cpm$DR42_kallisto   >= 0.58  | toptags_xypat_Kallisto_cpm$DR41_kallisto   >= 0.58 ) & (  toptags_xypat_Kallisto_cpm$Dr4FC    <= -0.58 | toptags_xypat_Kallisto_cpm$Dr4FC    >= 0.58),toptags_xypat_Kallisto_cpm$Dr4FC   ,NA) 
toptags_xypat_Kallisto_cpm$Dr5RNK=    ifelse((toptags_xypat_Kallisto_cpm$DR52_kallisto   >= 0.58  | toptags_xypat_Kallisto_cpm$DR51_kallisto   >= 0.58 ) & (  toptags_xypat_Kallisto_cpm$Dr5FC    <= -0.58 | toptags_xypat_Kallisto_cpm$Dr5FC    >= 0.58),toptags_xypat_Kallisto_cpm$Dr5FC   ,NA) 
toptags_xypat_Kallisto_cpm$FocaRNK=   ifelse((toptags_xypat_Kallisto_cpm$FOCA2_kallisto  >= 0.58  | toptags_xypat_Kallisto_cpm$FOCA1_kallisto   >= 0.58 ) & ( toptags_xypat_Kallisto_cpm$FocaFC   <= -0.58 | toptags_xypat_Kallisto_cpm$FocaFC   >= 0.58),toptags_xypat_Kallisto_cpm$FocaFC  ,NA) 
toptags_xypat_Kallisto_cpm$GagraRNK = ifelse((toptags_xypat_Kallisto_cpm$GAGRA2_kallisto >= 0.58  | toptags_xypat_Kallisto_cpm$GAGRA1_kallisto  >= 0.58 ) & ( toptags_xypat_Kallisto_cpm$GagraFC  <= -0.58 | toptags_xypat_Kallisto_cpm$GagraFC  >= 0.58),toptags_xypat_Kallisto_cpm$GagraFC ,NA) 
toptags_xypat_Kallisto_cpm$LuanRNK =  ifelse((toptags_xypat_Kallisto_cpm$LUAN2_kallisto  >= 0.58  | toptags_xypat_Kallisto_cpm$LUAN1_kallisto   >= 0.58 ) & ( toptags_xypat_Kallisto_cpm$LuanC    <= -0.58 | toptags_xypat_Kallisto_cpm$LuanC    >= 0.58),toptags_xypat_Kallisto_cpm$LuanC   ,NA) 
toptags_xypat_Kallisto_cpm$MabiRNK =  ifelse((toptags_xypat_Kallisto_cpm$MABI2_kallisto  >= 0.58  | toptags_xypat_Kallisto_cpm$MABI1_kallisto   >= 0.58 ) & ( toptags_xypat_Kallisto_cpm$MabiFC   <= -0.58 | toptags_xypat_Kallisto_cpm$MabiFC   >= 0.58),toptags_xypat_Kallisto_cpm$MabiFC  ,NA) 
toptags_xypat_Kallisto_cpm$MogeRNK =  ifelse((toptags_xypat_Kallisto_cpm$MOGE2_kallisto  >= 0.58  | toptags_xypat_Kallisto_cpm$MOGE1_kallisto   >= 0.58 ) & ( toptags_xypat_Kallisto_cpm$MogeFC   <= -0.58 | toptags_xypat_Kallisto_cpm$MogeFC   >= 0.58),toptags_xypat_Kallisto_cpm$MogeFC  ,NA) 
toptags_xypat_Kallisto_cpm$PiagRNK =  ifelse((toptags_xypat_Kallisto_cpm$PIAG2_kallisto  >= 0.58  | toptags_xypat_Kallisto_cpm$PIAG1_kallisto   >= 0.58 ) & ( toptags_xypat_Kallisto_cpm$PiagFC   <= -0.58 | toptags_xypat_Kallisto_cpm$PiagFC   >= 0.58),toptags_xypat_Kallisto_cpm$PiagFC  ,NA) 
toptags_xypat_Kallisto_cpm$PreluRNK = ifelse((toptags_xypat_Kallisto_cpm$PRELU2_kallisto >= 0.58  | toptags_xypat_Kallisto_cpm$PRELU1_kallisto  >= 0.58 ) & ( toptags_xypat_Kallisto_cpm$PreluFC  <= -0.58 | toptags_xypat_Kallisto_cpm$PreluFC  >= 0.58),toptags_xypat_Kallisto_cpm$PreluFC ,NA) 




#Ora vedo come sono correlati i diversi pazienti
#Prima per LogFC

library(ggcorrplot)
pdf("/home/fsantaniello/Correlazioni_AllGenes_Pazienti_NewAnalysis_Onlydiagnoses.pdf")
corr <- round(cor(toptags_xypat_Kallisto_cpm[,c(38:52)]), 1)
#p.mat <- cor_pmat(toptags_xypat_Kallisto_cpm[,38:52])
col<- colorRampPalette(c("blue", "white", "red"))(20)
ggcorrplot(corr, hc.order = TRUE,    outline.col = "white",   p.mat = p.mat,insig = "blank")
dev.off()

#poi se ci sono cambiamenti alle diagnosi o alle relapse
#Diagnosi

library(ggcorrplot)
pdf("/home/fsantaniello/Correlazioni_AllGenes_Pazienti_NewAnalysis_Onlydiagnoses.pdf")
corr_dx <- round(cor(toptags_xypat_Kallisto_cpm[,c(8,10,12,14,16,18,20,22,24,26,28,30,32,34,36)]), 1)
#p.mat <- cor_pmat(toptags_xypat_Kallisto_cpm[,38:52])
col<- colorRampPalette(c("blue", "white", "red"))(20)
ggcorrplot(corr_dx, hc.order = TRUE,    outline.col = "white",   p.mat = p.mat,insig = "blank")
dev.off()

#Relapse
library(ggcorrplot)
pdf("/home/fsantaniello/Correlazioni_AllGenes_Pazienti_NewAnalysis_OnlyRelapses.pdf")
corr_rx <- round(cor(toptags_xypat_Kallisto_cpm[,c(9,11,13,15,17,19,21,23,25,27,29,31,33,35,37)]), 1)
#p.mat <- cor_pmat(toptags_xypat_Kallisto_cpm[,38:52])
col<- colorRampPalette(c("blue", "white", "red"))(20)
ggcorrplot(corr_rx, hc.order = TRUE,    outline.col = "white",   p.mat = p.mat,insig = "blank")
dev.off()


#FAccio dendrogrammi delle correlazioni Diagnosi
pdf("/home/fsantaniello/Dendrogramma_Correlazioni_AllGenes_Pazienti_NewAnalysis_OnlyDiagnosis.pdf")
plot(hclust(dist(corr_dx)))
dev.off()


#FAccio dendrogrammi delle correlazioni Relapse
pdf("/home/fsantaniello/Dendrogramma_Correlazioni_AllGenes_Pazienti_NewAnalysis_OnlyRelapses.pdf")
plot(hclust(dist(corr_rx)))
dev.off()


Accio dendrogrammi delle correlazioni Relapse
pdf("/home/fsantaniello/Dendrogramma_Correlazioni_AllGenes_Pazienti_NewAnalysis_OnlyRelapses.pdf")
plot(hclust(dist(corr_rx)))
dev.off()


#FAccio dendrogrammi  delle correlazioni FoldChanges
pdf("/home/fsantaniello/Dendrogramma_AllGenes_Pazienti_NewAnalysis_OnlyFoldChanges.pdf")
plot(hclust(dist(corr)))
dev.off()


#I risultati sono interessanti. Provo a farci sopra un MDS (Uno per Diagnosi, uno per Relapse)

#Diagnosi
pdf("MDS_15PAzienti_DiagnosisOnly.pdf")
mds=plotMDS(toptags_xypat_Kallisto_cpm[,c(8,10,12,14,16,18,20,22,24,26,28,30,32,34,36)],dim.plot=c(1,2))
mds2=data.frame(x=mds$x,y=mds$y)
samples=rownames(mds2)
samples=unlist(strsplit(samples,split="_kallisto"))
samples=samples[samples!=""]
rownames(mds2)=samples
 ggplot(mds2,aes(x,y))+geom_point(aes(col=as.character(rownames(mds2)))) + geom_text(label= samples,size=3) + scale_x_continuous(limits = c(-5, 5))
dev.off()


#Relapse
pdf("MDS_15PAzienti_RelapseOnly.pdf")
mds=plotMDS(toptags_xypat_Kallisto_cpm[,c(9,11,13,15,17,19,21,23,25,27,29,31,33,35,37)],dim.plot=c(1,2))
mds2=data.frame(x=mds$x,y=mds$y)
samples=rownames(mds2)
samples=unlist(strsplit(samples,split="_kallisto"))
samples=samples[samples!=""]
rownames(mds2)=samples
 ggplot(mds2,aes(x,y))+geom_point(aes(col=as.character(rownames(mds2)))) + geom_text(label= samples,size=3) + scale_x_continuous(limits = c(-5, 5))
dev.off()

#logFCs
pdf("MDS_15PAzienti_logFCsOnly.pdf")
mds=plotMDS(toptags_xypat_Kallisto_cpm[,c(38:52)],dim.plot=c(1,2))
mds2=data.frame(x=mds$x,y=mds$y)
samples=rownames(mds2)
samples=unlist(strsplit(samples,split="_kallisto"))
samples=samples[samples!=""]
rownames(mds2)=samples
 ggplot(mds2,aes(x,y))+geom_point(aes(col=as.character(rownames(mds2)))) + geom_text(label= samples,size=3) 
dev.off()




colnames(toptags_xypat_Kallisto_cpm)
[1] "genes"           "logFC"           "logCPM"          "LR"
[5] "PValue"          "FDR"             "Qvalue"          "ALFE1_kallisto"
[9] "ALFE2_kallisto"  "BESU1_kallisto"  "BESU2_kallisto"  "CALU1_kallisto"
[13] "CALU2_kallisto"  "DEIV2_kallisto"  "DEIV2b_kallisto" "DEST1_kallisto"
[17] "DEST2_kallisto"  "DR1_1_kallisto"  "DR1_2_kallisto"  "DR41_kallisto"
[21] "DR42_kallisto"   "DR51_kallisto"   "DR52_kallisto"   "FOCA1_kallisto"
[25] "FOCA2_kallisto"  "GAGRA1_kallisto" "GAGRA2_kallisto" "LUAN1_kallisto"
[29] "LUAN2_kallisto"  "MABI1_kallisto"  "MABI2_kallisto"  "MOGE1_kallisto"
[33] "MOGE2_kallisto"  "PIAG1_kallisto"  "PIAG2_kallisto"  "PRELU1_kallisto"
[37] "PRELU2_kallisto" "AlfeFC"          "BesuFC"          "CaluFC"
[41] "DeivFC"          "DestFC"          "Dr1FC"           "Dr4FC"
[45] "Dr5FC"           "FocaFC"          "GagraFC"         "LuanFC"
[49] "MabiFC"          "MogeFC"          "PiagFC"          "PreluFC"





 library(mixtools)
mixmdl = normalmixEM(t(toptags_xypat_Kallisto_cpm[,c(8,10,12,14,16,18,20,22,24,26,28,30,32,34,36)]))

#Ora ottengo tutte le conte per i pazienti levando il vincolo di almeno 1 cpm in almeno 8 pazienti

countdata_filt <- countdata[rowSums(countdata[, 1:30]) > 0, ]
countdata_filtDGE_Nothreshold=DGEList(countdata_filt[,1:30],genes=rownames(countdata_filt),group=condition)
keep_nothreshold <- rowSums(cpm(countdata_filtDGE_Nothreshold)>1) >= 1 #PRENDO SOLO I GENI CON CPM>1 IN ALMENO 8 CAMPIONI
countdata_filtDGE_keep_Nothreshold <- countdata_filtDGE_Nothreshold[keep_nothreshold, , keep.lib.sizes=FALSE]
countdata_filtDGE_keep_Nothreshold=calcNormFactors(countdata_filtDGE_keep_Nothreshold,method="TMM")
countdata_filtDGE_Nothreshold_CPM=as.data.frame(cpm(countdata_filtDGE_Nothreshold,log=T))
countdata_filtDGE_Nothreshold_CPM$genes=rownames(countdata_filtDGE_Nothreshold_CPM)

#Ora su questi dati creo dei files RNK
#Inizio con i FCs
countdata_filtDGE_Nothreshold_CPM$AlfeFC  =  countdata_filtDGE_Nothreshold_CPM$ALFE2_kallisto   -   countdata_filtDGE_Nothreshold_CPM$ALFE1_kallisto 
countdata_filtDGE_Nothreshold_CPM$BesuFC  =  countdata_filtDGE_Nothreshold_CPM$BESU2_kallisto   -   countdata_filtDGE_Nothreshold_CPM$BESU1_kallisto 
countdata_filtDGE_Nothreshold_CPM$CaluFC  =  countdata_filtDGE_Nothreshold_CPM$CALU2_kallisto   -   countdata_filtDGE_Nothreshold_CPM$CALU1_kallisto 
countdata_filtDGE_Nothreshold_CPM$DeivFC  =  countdata_filtDGE_Nothreshold_CPM$DEIV2b_kallisto  -   countdata_filtDGE_Nothreshold_CPM$DEIV2_kallisto
countdata_filtDGE_Nothreshold_CPM$DestFC  =  countdata_filtDGE_Nothreshold_CPM$DEST2_kallisto   -   countdata_filtDGE_Nothreshold_CPM$DEST1_kallisto 
countdata_filtDGE_Nothreshold_CPM$Dr1FC   =   countdata_filtDGE_Nothreshold_CPM$DR1_2_kallisto    -  countdata_filtDGE_Nothreshold_CPM$DR1_1_kallisto 
countdata_filtDGE_Nothreshold_CPM$Dr4FC   =   countdata_filtDGE_Nothreshold_CPM$DR42_kallisto     -  countdata_filtDGE_Nothreshold_CPM$DR41_kallisto 
countdata_filtDGE_Nothreshold_CPM$Dr5FC   =   countdata_filtDGE_Nothreshold_CPM$DR52_kallisto     -  countdata_filtDGE_Nothreshold_CPM$DR51_kallisto 
countdata_filtDGE_Nothreshold_CPM$FocaFC  =  countdata_filtDGE_Nothreshold_CPM$FOCA2_kallisto   -   countdata_filtDGE_Nothreshold_CPM$FOCA1_kallisto 
countdata_filtDGE_Nothreshold_CPM$GagraFC = countdata_filtDGE_Nothreshold_CPM$GAGRA2_kallisto -    countdata_filtDGE_Nothreshold_CPM$GAGRA1_kallisto
countdata_filtDGE_Nothreshold_CPM$LuanC   =   countdata_filtDGE_Nothreshold_CPM$LUAN2_kallisto    -  countdata_filtDGE_Nothreshold_CPM$LUAN1_kallisto 
countdata_filtDGE_Nothreshold_CPM$MabiFC  =  countdata_filtDGE_Nothreshold_CPM$MABI2_kallisto   -   countdata_filtDGE_Nothreshold_CPM$MABI1_kallisto 
countdata_filtDGE_Nothreshold_CPM$MogeFC  =  countdata_filtDGE_Nothreshold_CPM$MOGE2_kallisto   -   countdata_filtDGE_Nothreshold_CPM$MOGE1_kallisto 
countdata_filtDGE_Nothreshold_CPM$PiagFC  =  countdata_filtDGE_Nothreshold_CPM$PIAG2_kallisto   -   countdata_filtDGE_Nothreshold_CPM$PIAG1_kallisto 
countdata_filtDGE_Nothreshold_CPM$PreluFC = countdata_filtDGE_Nothreshold_CPM$PRELU2_kallisto -    countdata_filtDGE_Nothreshold_CPM$PRELU1_kallisto


countdata_filtDGE_Nothreshold_CPM$AlfeRNK=   ifelse((countdata_filtDGE_Nothreshold_CPM$ALFE2_kallisto  >= 0.58  | countdata_filtDGE_Nothreshold_CPM$ALFE1_kallisto   >= 0.58 ) & ( countdata_filtDGE_Nothreshold_CPM$AlfeFC   <= -0.58 | countdata_filtDGE_Nothreshold_CPM$AlfeFC   >= 0.58),countdata_filtDGE_Nothreshold_CPM$AlfeFC  ,NA) 
countdata_filtDGE_Nothreshold_CPM$BesuRNK=   ifelse((countdata_filtDGE_Nothreshold_CPM$BESU2_kallisto  >= 0.58  | countdata_filtDGE_Nothreshold_CPM$BESU1_kallisto   >= 0.58 ) & ( countdata_filtDGE_Nothreshold_CPM$BesuFC   <= -0.58 | countdata_filtDGE_Nothreshold_CPM$BesuFC   >= 0.58),countdata_filtDGE_Nothreshold_CPM$BesuFC  ,NA) 
countdata_filtDGE_Nothreshold_CPM$CaluRNK=   ifelse((countdata_filtDGE_Nothreshold_CPM$CALU2_kallisto  >= 0.58  | countdata_filtDGE_Nothreshold_CPM$CALU1_kallisto   >= 0.58 ) & ( countdata_filtDGE_Nothreshold_CPM$CaluFC   <= -0.58 | countdata_filtDGE_Nothreshold_CPM$CaluFC   >= 0.58),countdata_filtDGE_Nothreshold_CPM$CaluFC  ,NA) 
countdata_filtDGE_Nothreshold_CPM$DeivRNK=   ifelse((countdata_filtDGE_Nothreshold_CPM$DEIV2b_kallisto >= 0.58  | countdata_filtDGE_Nothreshold_CPM$DEIV2_kallisto  >= 0.58 ) & (  countdata_filtDGE_Nothreshold_CPM$DeivFC   <= -0.58 | countdata_filtDGE_Nothreshold_CPM$DeivFC   >= 0.58),countdata_filtDGE_Nothreshold_CPM$DeivFC  ,NA) 
countdata_filtDGE_Nothreshold_CPM$DestRNK=   ifelse((countdata_filtDGE_Nothreshold_CPM$DEST2_kallisto  >= 0.58  | countdata_filtDGE_Nothreshold_CPM$DEST1_kallisto   >= 0.58 ) & ( countdata_filtDGE_Nothreshold_CPM$DestFC   <= -0.58 | countdata_filtDGE_Nothreshold_CPM$DestFC   >= 0.58),countdata_filtDGE_Nothreshold_CPM$DestFC  ,NA) 
countdata_filtDGE_Nothreshold_CPM$Dr1RNK=    ifelse((countdata_filtDGE_Nothreshold_CPM$DR1_2_kallisto  >= 0.58  | countdata_filtDGE_Nothreshold_CPM$DR1_1_kallisto   >= 0.58 ) & ( countdata_filtDGE_Nothreshold_CPM$Dr1FC    <= -0.58 | countdata_filtDGE_Nothreshold_CPM$Dr1FC    >= 0.58),countdata_filtDGE_Nothreshold_CPM$Dr1FC   ,NA) 
countdata_filtDGE_Nothreshold_CPM$Dr4RNK=    ifelse((countdata_filtDGE_Nothreshold_CPM$DR42_kallisto   >= 0.58  | countdata_filtDGE_Nothreshold_CPM$DR41_kallisto   >= 0.58 ) & (  countdata_filtDGE_Nothreshold_CPM$Dr4FC    <= -0.58 | countdata_filtDGE_Nothreshold_CPM$Dr4FC    >= 0.58),countdata_filtDGE_Nothreshold_CPM$Dr4FC   ,NA) 
countdata_filtDGE_Nothreshold_CPM$Dr5RNK=    ifelse((countdata_filtDGE_Nothreshold_CPM$DR52_kallisto   >= 0.58  | countdata_filtDGE_Nothreshold_CPM$DR51_kallisto   >= 0.58 ) & (  countdata_filtDGE_Nothreshold_CPM$Dr5FC    <= -0.58 | countdata_filtDGE_Nothreshold_CPM$Dr5FC    >= 0.58),countdata_filtDGE_Nothreshold_CPM$Dr5FC   ,NA) 
countdata_filtDGE_Nothreshold_CPM$FocaRNK=   ifelse((countdata_filtDGE_Nothreshold_CPM$FOCA2_kallisto  >= 0.58  | countdata_filtDGE_Nothreshold_CPM$FOCA1_kallisto   >= 0.58 ) & ( countdata_filtDGE_Nothreshold_CPM$FocaFC   <= -0.58 | countdata_filtDGE_Nothreshold_CPM$FocaFC   >= 0.58),countdata_filtDGE_Nothreshold_CPM$FocaFC  ,NA) 
countdata_filtDGE_Nothreshold_CPM$GagraRNK = ifelse((countdata_filtDGE_Nothreshold_CPM$GAGRA2_kallisto >= 0.58  | countdata_filtDGE_Nothreshold_CPM$GAGRA1_kallisto  >= 0.58 ) & ( countdata_filtDGE_Nothreshold_CPM$GagraFC  <= -0.58 | countdata_filtDGE_Nothreshold_CPM$GagraFC  >= 0.58),countdata_filtDGE_Nothreshold_CPM$GagraFC ,NA) 
countdata_filtDGE_Nothreshold_CPM$LuanRNK =  ifelse((countdata_filtDGE_Nothreshold_CPM$LUAN2_kallisto  >= 0.58  | countdata_filtDGE_Nothreshold_CPM$LUAN1_kallisto   >= 0.58 ) & ( countdata_filtDGE_Nothreshold_CPM$LuanC    <= -0.58 | countdata_filtDGE_Nothreshold_CPM$LuanC    >= 0.58),countdata_filtDGE_Nothreshold_CPM$LuanC   ,NA) 
countdata_filtDGE_Nothreshold_CPM$MabiRNK =  ifelse((countdata_filtDGE_Nothreshold_CPM$MABI2_kallisto  >= 0.58  | countdata_filtDGE_Nothreshold_CPM$MABI1_kallisto   >= 0.58 ) & ( countdata_filtDGE_Nothreshold_CPM$MabiFC   <= -0.58 | countdata_filtDGE_Nothreshold_CPM$MabiFC   >= 0.58),countdata_filtDGE_Nothreshold_CPM$MabiFC  ,NA) 
countdata_filtDGE_Nothreshold_CPM$MogeRNK =  ifelse((countdata_filtDGE_Nothreshold_CPM$MOGE2_kallisto  >= 0.58  | countdata_filtDGE_Nothreshold_CPM$MOGE1_kallisto   >= 0.58 ) & ( countdata_filtDGE_Nothreshold_CPM$MogeFC   <= -0.58 | countdata_filtDGE_Nothreshold_CPM$MogeFC   >= 0.58),countdata_filtDGE_Nothreshold_CPM$MogeFC  ,NA) 
countdata_filtDGE_Nothreshold_CPM$PiagRNK =  ifelse((countdata_filtDGE_Nothreshold_CPM$PIAG2_kallisto  >= 0.58  | countdata_filtDGE_Nothreshold_CPM$PIAG1_kallisto   >= 0.58 ) & ( countdata_filtDGE_Nothreshold_CPM$PiagFC   <= -0.58 | countdata_filtDGE_Nothreshold_CPM$PiagFC   >= 0.58),countdata_filtDGE_Nothreshold_CPM$PiagFC  ,NA) 
countdata_filtDGE_Nothreshold_CPM$PreluRNK = ifelse((countdata_filtDGE_Nothreshold_CPM$PRELU2_kallisto >= 0.58  | countdata_filtDGE_Nothreshold_CPM$PRELU1_kallisto  >= 0.58 ) & ( countdata_filtDGE_Nothreshold_CPM$PreluFC  <= -0.58 | countdata_filtDGE_Nothreshold_CPM$PreluFC  >= 0.58),countdata_filtDGE_Nothreshold_CPM$PreluFC ,NA) 



#Analisi GSEA



ALFE_RNK=countdata_filtDGE_Nothreshold_CPM[countdata_filtDGE_Nothreshold_CPM$ALFE2_kallisto>= 0.58 | countdata_filtDGE_Nothreshold_CPM$ALFE1_kallisto>=0.58,c(31,2,1)]
BESU_RNK=countdata_filtDGE_Nothreshold_CPM[countdata_filtDGE_Nothreshold_CPM$BESU2_kallisto>= 0.58 | countdata_filtDGE_Nothreshold_CPM$BESU1_kallisto>=0.58,c(31,4,3)]
CALU_RNK=countdata_filtDGE_Nothreshold_CPM[countdata_filtDGE_Nothreshold_CPM$CALU2_kallisto>= 0.58 | countdata_filtDGE_Nothreshold_CPM$CALU1_kallisto>=0.58,c(31,6,5)]
DEIV_RNK=countdata_filtDGE_Nothreshold_CPM[countdata_filtDGE_Nothreshold_CPM$DEIV2b_kallisto>= 0.58 | countdata_filtDGE_Nothreshold_CPM$DEIV2_kallisto>=0.58,c(31,8,7)]
DEST_RNK=countdata_filtDGE_Nothreshold_CPM[countdata_filtDGE_Nothreshold_CPM$DEST2_kallisto>= 0.58 | countdata_filtDGE_Nothreshold_CPM$DEST1_kallisto>=0.58,c(31,10,9)]
DR1_RNK=countdata_filtDGE_Nothreshold_CPM[countdata_filtDGE_Nothreshold_CPM$DR1_2_kallisto>= 0.58 | countdata_filtDGE_Nothreshold_CPM$DR1_1_kallisto>=0.58,c(31,12,11)]
DR4_RNK=countdata_filtDGE_Nothreshold_CPM[countdata_filtDGE_Nothreshold_CPM$DR42_kallisto>= 0.58 | countdata_filtDGE_Nothreshold_CPM$DR41_kallisto>=0.58,c(31,14,13)]
DR5_RNK=countdata_filtDGE_Nothreshold_CPM[countdata_filtDGE_Nothreshold_CPM$DR52_kallisto>= 0.58 | countdata_filtDGE_Nothreshold_CPM$DR51_kallisto>=0.58,c(31,16,15)]
FOCA_RNK=countdata_filtDGE_Nothreshold_CPM[countdata_filtDGE_Nothreshold_CPM$FOCA2_kallisto>= 0.58 | countdata_filtDGE_Nothreshold_CPM$FOCA1_kallisto>=0.58,c(31,18,17)]
GAGRA_RNK=countdata_filtDGE_Nothreshold_CPM[countdata_filtDGE_Nothreshold_CPM$GAGRA2_kallisto>= 0.58 | countdata_filtDGE_Nothreshold_CPM$GAGRA1_kallisto>=0.58,c(31,20,19)]
LUAN_RNK=countdata_filtDGE_Nothreshold_CPM[countdata_filtDGE_Nothreshold_CPM$LUAN2_kallisto>= 0.58 | countdata_filtDGE_Nothreshold_CPM$LUAN1_kallisto>=0.58,c(31,22,21)]
MABI_RNK=countdata_filtDGE_Nothreshold_CPM[countdata_filtDGE_Nothreshold_CPM$MABI2_kallisto>= 0.58 | countdata_filtDGE_Nothreshold_CPM$MABI1_kallisto>=0.58,c(31,24,23)]
MOGE_RNK=countdata_filtDGE_Nothreshold_CPM[countdata_filtDGE_Nothreshold_CPM$MOGE2_kallisto>= 0.58 | countdata_filtDGE_Nothreshold_CPM$MOGE1_kallisto>=0.58,c(31,26,25)]
PIAG_RNK=countdata_filtDGE_Nothreshold_CPM[countdata_filtDGE_Nothreshold_CPM$PIAG2_kallisto>= 0.58 | countdata_filtDGE_Nothreshold_CPM$PIAG1_kallisto>=0.58,c(31,28,27)]
PRELU_RNK=countdata_filtDGE_Nothreshold_CPM[countdata_filtDGE_Nothreshold_CPM$PRELU2_kallisto>= 0.58 | countdata_filtDGE_Nothreshold_CPM$PRELU1_kallisto>=0.58,c(31,30,29)]

colnames(ALFE_RNK)=c('GeneSymbol','RELAPSE','DIAGNOSIS')
colnames(BESU_RNK)=c('GeneSymbol','RELAPSE','DIAGNOSIS')
colnames(CALU_RNK)=c('GeneSymbol','RELAPSE','DIAGNOSIS')
colnames(DEIV_RNK)=c('GeneSymbol','RELAPSE','DIAGNOSIS')
colnames(DEST_RNK)=c('GeneSymbol','RELAPSE','DIAGNOSIS')
colnames(DR1_RNK)=c('GeneSymbol','RELAPSE','DIAGNOSIS')
colnames(DR4_RNK)=c('GeneSymbol','RELAPSE','DIAGNOSIS')
colnames(DR5_RNK)=c('GeneSymbol','RELAPSE','DIAGNOSIS')
colnames(FOCA_RNK)=c('GeneSymbol','RELAPSE','DIAGNOSIS')
colnames(GAGRA_RNK)=c('GeneSymbol','RELAPSE','DIAGNOSIS')
colnames(LUAN_RNK)=c('GeneSymbol','RELAPSE','DIAGNOSIS')
colnames(MABI_RNK)=c('GeneSymbol','RELAPSE','DIAGNOSIS')
colnames(MOGE_RNK)=c('GeneSymbol','RELAPSE','DIAGNOSIS')
colnames(PIAG_RNK)=c('GeneSymbol','RELAPSE','DIAGNOSIS')
colnames(PRELU_RNK)=c('GeneSymbol','RELAPSE','DIAGNOSIS')


write.table(ALFE_RNK,"ALFE_Expression_Kallisto_NeoEpitopes.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(BESU_RNK,"BESU_Expression_Kallisto_NeoEpitopes.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(CALU_RNK,"CALU_Expression_Kallisto_NeoEpitopes.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(DEIV_RNK,"DEIV_Expression_Kallisto_NeoEpitopes.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(DEST_RNK,"DEST_Expression_Kallisto_NeoEpitopes.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(DR1_RNK,"DR1_Expression_Kallisto_NeoEpitopes.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(DR4_RNK,"DR4_Expression_Kallisto_NeoEpitopes.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(DR5_RNK,"DR5_Expression_Kallisto_NeoEpitopes.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(FOCA_RNK,"FOCA_Expression_Kallisto_NeoEpitopes.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(GAGRA_RNK,"GAGRA_Expression_Kallisto_NeoEpitopes.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(LUAN_RNK,"LUAN_Expression_Kallisto_NeoEpitopes.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(MABI_RNK,"MABI_Expression_Kallisto_NeoEpitopes.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(MOGE_RNK,"MOGE_Expression_Kallisto_NeoEpitopes.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(PIAG_RNK,"PIAG_Expression_Kallisto_NeoEpitopes.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(PRELU_RNK,"PRELU_Expression_Kallisto_NeoEpitopes.txt",sep="\t",col.names=T,row.names=F,quote=F)



head(ALFE_RNK  )
head(BESU_RNK  )
head(CALU_RNK  )
head(DEIV_RNK  )
head(DEST_RNK  )
head(DR1_RNK   )
head(DR4_RNK   )
head(DR5_RNK   )
head(FOCA_RNK  )
head(GAGRA_RNK )
head(LUAN_RNK  )
head(MABI_RNK  )
head(MOGE_RNK  )
head(PIAG_RNK  )
head(PRELU_RNK )



ALFE_RNK  = as.data.frame(ALFE_RNK[order(-ALFE_RNK$AlfeRNK),])
BESU_RNK  = as.data.frame(BESU_RNK[order(-BESU_RNK$BesuRNK),])
CALU_RNK  = as.data.frame(CALU_RNK[order(-CALU_RNK$CaluRNK),])
DEIV_RNK  = as.data.frame(DEIV_RNK[order(-DEIV_RNK$DeivRNK),])
DEST_RNK  = as.data.frame(DEST_RNK[order(-DEST_RNK$DestRNK),])
DR1_RNK   = as.data.frame(DR1_RNK[order(-DR1_RNK$Dr1RNK),])
DR4_RNK   = as.data.frame(DR4_RNK[order(-DR4_RNK$Dr4RNK),])
DR5_RNK   = as.data.frame(DR5_RNK[order(-DR5_RNK$Dr5RNK),])
FOCA_RNK  = as.data.frame(FOCA_RNK[order(-FOCA_RNK$FocaRNK),])
GAGRA_RNK = as.data.frame(GAGRA_RNK[order(-GAGRA_RNK$GagraRNK),])
LUAN_RNK  = as.data.frame(LUAN_RNK[order(-LUAN_RNK$LuanRNK),])
MABI_RNK  = as.data.frame(MABI_RNK[order(-MABI_RNK$MabiRNK),])
MOGE_RNK  = as.data.frame(MOGE_RNK[order(-MOGE_RNK$MogeRNK),])
PIAG_RNK  = as.data.frame(PIAG_RNK[order(-PIAG_RNK$PiagRNK),])
PRELU_RNK = as.data.frame(PRELU_RNK[order(-PRELU_RNK$PreluRNK),])

ALFE_RNK  = ALFE_RNK[ALFE_RNK$AlfeRNK!="NA",]
BESU_RNK  = BESU_RNK[BESU_RNK$BesuRNK!="NA",]
CALU_RNK  = CALU_RNK[CALU_RNK$CaluRNK!="NA",]
DEIV_RNK  = DEIV_RNK[DEIV_RNK$DeivRNK!="NA",]
DEST_RNK  = DEST_RNK[DEST_RNK$DestRNK!="NA",]
DR1_RNK   = DR1_RNK[DR1_RNK$Dr1RNK!="NA",]
DR4_RNK   = DR4_RNK[DR4_RNK$Dr4RNK!="NA",]
DR5_RNK   = DR5_RNK[DR5_RNK$Dr5RNK!="NA",]
FOCA_RNK  = FOCA_RNK[FOCA_RNK$FocaRNK!="NA",]
GAGRA_RNK = GAGRA_RNK[GAGRA_RNK$GagraRNK!="NA",]
LUAN_RNK  = LUAN_RNK[LUAN_RNK$LuanRNK!="NA",]
MABI_RNK  = MABI_RNK[MABI_RNK$MabiRNK!="NA",]
MOGE_RNK  = MOGE_RNK[MOGE_RNK$MogeRNK!="NA",]
PIAG_RNK  = PIAG_RNK[PIAG_RNK$PiagRNK!="NA",]
PRELU_RNK = PRELU_RNK[PRELU_RNK$PreluRNK!="NA",]







write.table(ALFE_RNK,"ALFE_RNK_NoThreshold.rnk",sep="\t",col.names=F,row.names=F,quote=F)
write.table(BESU_RNK,"BESU_RNK_NoThreshold.rnk",sep="\t",col.names=F,row.names=F,quote=F)
write.table(CALU_RNK,"CALU_RNK_NoThreshold.rnk",sep="\t",col.names=F,row.names=F,quote=F)
write.table(DEIV_RNK,"DEIV_RNK_NoThreshold.rnk",sep="\t",col.names=F,row.names=F,quote=F)
write.table(DEST_RNK,"DEST_RNK_NoThreshold.rnk",sep="\t",col.names=F,row.names=F,quote=F)
write.table(DR1_RNK,"DR1_RNK_NoThreshold.rnk",sep="\t",col.names=F,row.names=F,quote=F)
write.table(DR4_RNK,"DR4_RNK_NoThreshold.rnk",sep="\t",col.names=F,row.names=F,quote=F)
write.table(DR5_RNK,"DR5_RNK_NoThreshold.rnk",sep="\t",col.names=F,row.names=F,quote=F)
write.table(FOCA_RNK,"FOCA_RNK_NoThreshold.rnk",sep="\t",col.names=F,row.names=F,quote=F)
write.table(GAGRA_RNK,"GAGRA_RNK_NoThreshold.rnk",sep="\t",col.names=F,row.names=F,quote=F)
write.table(LUAN_RNK,"LUAN_RNK_NoThreshold.rnk",sep="\t",col.names=F,row.names=F,quote=F)
write.table(MABI_RNK,"MABI_RNK_NoThreshold.rnk",sep="\t",col.names=F,row.names=F,quote=F)
write.table(MOGE_RNK,"MOGE_RNK_NoThreshold.rnk",sep="\t",col.names=F,row.names=F,quote=F)
write.table(PIAG_RNK,"PIAG_RNK_NoThreshold.rnk",sep="\t",col.names=F,row.names=F,quote=F)
write.table(PRELU_RNK,"PRELU_RNK_NoThreshold.rnk",sep="\t",col.names=F,row.names=F,quote=F)



#confronto con dati di cristina toffalori

countdata_filtDGE_Nothreshold_CPM_cri=countdata_filtDGE_Nothreshold_CPM_cri[grep("HLA-DR|HLA-DQ|CD274|CD273|CD276|C10orf54|LGALS9|CD200|PVR|PVRL2|CD86|CD80|ICOSL6|TNFSF4|ITGAL",rownames(countdata_filtDGE_Nothreshold_CPM_cri)),]

#salvo e rimaneggio su excel

write.table(countdata_filtDGE_Nothreshold_CPM_cri,'../Progetto_LuciaGabri/Geni_Heatmap_Cristina.txt',sep='\t',col.names=T,row.names=F,quote=F)
geni_cristina_heat=read.table('../Progetto_LuciaGabri/Geni_Heatmap_Cristina.txt',head=T)
head(geni_cristina_heat)
```
      CaluFC      DestFC    FocaFC     GagraFC     MabiFC     MogeFC    PreluFC    genes
1 -0.3605993  0.02429480 -1.699084  0.47230892 0.38002653 -0.5940172 -0.4786808  HLA-DRA
2 -0.6522032 -0.21273213 -2.026820  0.55672355 0.26013928 -0.6629054 -0.5094507 HLA-DRB1
3 -0.5489297 -0.24177024 -2.112525  0.07473945 0.09590079 -0.7681641 -0.4145400 HLA-DRB5
4 -2.2354618  0.07294259 -1.315932 -0.53982933 0.99494504 -0.8202409 -1.3603720 HLA-DQA1
5 -3.6244392  0.57705045 -1.348464  0.11176613 3.02386435 -0.7625783 -1.2083065 HLA-DQA2
6 -0.3150906  0.15821667 -3.050103  0.21787927 0.36190111 -0.6375350 -1.3670686 HLA-DQB1
```

rownames(geni_cristina_heat)=geni_cristina_heat$genes

#creo heatmap
breakList=seq(-8, 8, by = 0.1)

pheatmap( geni_cristina_heat[,1:7], cluster_rows=F,cluster_cols=T,color=colorRampPalette( c("cyan", "gray15", "yellow1"), space="rgb")(length(breakList)),breaks=breakList)


#GSEA



library(dplyr)
AgI=read.table('Ag_I_genes.grp',head=F)
colnames(AgI)="genes"
AgI=merge(countdata_filtDGE_Nothreshold_CPM,AgI,by='genes')
write.table(AgI,'AgI_List_wInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)


AgII=read.table('Ag_II_genes.grp',head=F)
colnames(AgII)="genes"
AgII=merge(countdata_filtDGE_Nothreshold_CPM,AgII,by='genes')
write.table(AgII,'AgII_List_wInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)


ActMols=read.table('ActivatoryMolecules.grp',head=F)
colnames(ActMols)="genes"
ActMols=merge(countdata_filtDGE_Nothreshold_CPM,ActMols,by='genes')
write.table(ActMols,'ActMols_List_wInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)

CoActCoStim=read.table('CoActCoStim_genes.grp',head=F)
colnames(CoActCoStim)="genes"
CoActCoStim=merge(countdata_filtDGE_Nothreshold_CPM,CoActCoStim,by='genes')
write.table(CoActCoStim,'CoActCoStim_List_wInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)

CoInhibMols=read.table('CoInhibMolecules_genes.grp',head=F)
colnames(CoInhibMols)="genes"
CoInhibMols=merge(countdata_filtDGE_Nothreshold_CPM,CoInhibMols,by='genes')
write.table(CoInhibMols,'CoInhibMols_List_wInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)

hla=read.table('HLA_MHCI.grp',head=F)
colnames(hla)="genes"
hla=merge(countdata_filtDGE_Nothreshold_CPM,hla,by='genes')
write.table(hla,'HLA_mhcI_List_wInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)

hla2=read.table('HLA_MHCII.grp',head=F)
colnames(hla2)="genes"
hla2=merge(countdata_filtDGE_Nothreshold_CPM,hla2,by='genes')
write.table(hla2,'HLA_mhcII_List_wInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)

PosRespIfnG=read.table('PositiveResponseIfnG.grp',head=F)
colnames(PosRespIfnG)="genes"
PosRespIfnG=merge(countdata_filtDGE_Nothreshold_CPM,PosRespIfnG,by='genes')
write.table(PosRespIfnG,'PositiveResponseIfnG_List_wInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)

NegRespIfnG=read.table('NegativeResponseIfnG.grp',head=F)
colnames(NegRespIfnG)="genes"
NegRespIfnG=merge(countdata_filtDGE_Nothreshold_CPM,NegRespIfnG,by='genes')
write.table(NegRespIfnG,'NegativeResponseIfnG_List_wInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)

ImmSopr=read.table('ImmunoSopr_Molecules.grp',head=F)
colnames(ImmSopr)="genes"
ImmSopr=merge(countdata_filtDGE_Nothreshold_CPM,ImmSopr,by='genes')
write.table(ImmSopr,'ImmunoSopr_Molecules_List_wInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)

LeuAGs=read.table('LEU_AGs.grp',head=F)
colnames(LeuAGs)="genes"
LeuAGs=merge(countdata_filtDGE_Nothreshold_CPM,LeuAGs,by='genes')
write.table(LeuAGs,'LEU_AGs_List_wInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)

CTA=read.table('CTA_genes.grp',head=F)
colnames(CTA)="genes"
CTA=merge(countdata_filtDGE_Nothreshold_CPM,CTA,by='genes')
write.table(CTA,'CTA_genes_List_wInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)

SharedProcesses=read.table('SharedProcesses.grp',head=F)
colnames(SharedProcesses)="genes"
SharedProcesses=merge(countdata_filtDGE_Nothreshold_CPM,SharedProcesses,by='genes')
write.table(SharedProcesses,'SharedProcesses_List_wInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)

DnaReplication=read.table('DnaReplication.grp',head=F)
colnames(DnaReplication)="genes"
DnaReplication=merge(countdata_filtDGE_Nothreshold_CPM,DnaReplication,by='genes')
write.table(DnaReplication,'DnaReplication_List_wInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)

DnaRepair=read.table('DnaRepair.grp',head=F)
colnames(DnaRepair)="genes"
DnaRepair=merge(countdata_filtDGE_Nothreshold_CPM,DnaRepair,by='genes')
write.table(DnaRepair,'DnaRepair_List_wInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)

DnaG1STransition=read.table('DnaG1STransition.grp',head=F)
colnames(DnaG1STransition)="genes"
DnaG1STransition=merge(countdata_filtDGE_Nothreshold_CPM,DnaG1STransition,by='genes')
write.table(DnaG1STransition,'DnaG1STransition_List_wInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)

DnaRecomb=read.table('DnaRecomb.grp',head=F)
colnames(DnaRecomb)="genes"
DnaRecomb=merge(countdata_filtDGE_Nothreshold_CPM,DnaRecomb,by='genes')
write.table(DnaRecomb,'DnaRecomb_List_wInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)

DnaDamageResponseDetection=read.table('DnaDamageResponseDetection.grp',head=F)
colnames(DnaDamageResponseDetection)="genes"
DnaDamageResponseDetection=merge(countdata_filtDGE_Nothreshold_CPM,DnaDamageResponseDetection,by='genes')
write.table(DnaDamageResponseDetection,'DnaDamageResponseDetection_List_wInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)

IFNg_Response_Genes=read.table('IFNg_Response_Genes.grp',head=F)
colnames(IFNg_Response_Genes)="genes"
IFNg_Response_Genes=merge(countdata_filtDGE_Nothreshold_CPM,IFNg_Response_Genes,by='genes')
write.table(IFNg_Response_Genes,'IFNg_Response_Genes_List_wInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)







head(BESU_RNK_GSEA)
head(CALU_RNK_GSEA)
head(DEIV_RNK_GSEA)
head(DEST_RNK_GSEA)
head(DR1_RNK_GSEA)
head(DR4_RNK_GSEA)
head(DR5_RNK_GSEA)
head(FOCA_RNK_GSEA)
head(GAGRA_RNK_GSEA)
head(LUAN_RNK_GSEA)
head(MABI_RNK_GSEA)
head(MOGE_RNK_GSEA)
head(PIAG_RNK_GSEA)
head(PRELU_RNK_GSEA)





...
...

write.table(as.data.frame(toptags_xypat_Kallisto_cpm_comm_sig[toptags_xypat_Kallisto_cpm_comm$comunita==1,1]),"/home/fsantaniello/Comunita_1_Genes.txt",sep="\t",col.names=F,row.names=F)
write.table(as.data.frame(toptags_xypat_Kallisto_cpm_comm_sig[toptags_xypat_Kallisto_cpm_comm$comunita==2,1]),"/home/fsantaniello/Comunita_2_Genes.txt",sep="\t",col.names=F,row.names=F)
write.table(as.data.frame(toptags_xypat_Kallisto_cpm_comm_sig[toptags_xypat_Kallisto_cpm_comm$comunita==3,1]),"/home/fsantaniello/Comunita_3_Genes.txt",sep="\t",col.names=F,row.names=F)
...

pdf("Dendro_Pazienti.pdf")
plot(hclust(dist(corr),"complete"))
dev.off()


colnames(countdata_filt)
[1] "ALFE1_kallisto"  "ALFE2_kallisto"  "BESU1_kallisto"  "BESU2_kallisto"
[5] "CALU1_kallisto"  "CALU2_kallisto"  "DEIV2_kallisto"  "DEIV2b_kallisto"
[9] "DEST1_kallisto"  "DEST2_kallisto"  "DR1_1_kallisto"  "DR1_2_kallisto"
[13] "DR41_kallisto"   "DR42_kallisto"   "DR51_kallisto"   "DR52_kallisto"
[17] "FOCA1_kallisto"  "FOCA2_kallisto"  "GAGRA1_kallisto" "GAGRA2_kallisto"
[21] "LUAN1_kallisto"  "LUAN2_kallisto"  "MABI1_kallisto"  "MABI2_kallisto"
[25] "MOGE1_kallisto"  "MOGE2_kallisto"  "PIAG1_kallisto"  "PIAG2_kallisto"
[29] "PRELU1_kallisto" "PRELU2_kallisto"


Gruppo1= countdata_filt[,c(1,2,5,6,9,10,19,20,25,26,29,30)]
ALFE,CALU,DEST,GAGRA,MOGE,PRELU

Gruppo2=countdata_filt[,c(3,4,7,8,23,24,27,28)]
BESU,DEIV,MABI,PIAG
Gruppo3=countdata_filt[,c(11,12,13,14,15,16,17,18,21,22)]
DR1,DR4,DR5,FOCA,LUAN



condition_gruppo1=rep(c("Diagnosis","Relapse"),6)
gruppo1_filtDGE=DGEList(Gruppo1,genes=rownames(Gruppo1),group=condition_gruppo1)
keep <- rowSums(cpm(gruppo1_filtDGE)>1) >= 2 #PRENDO SOLO I GENI CON CPM>1 IN ALMENO 8 CAMPIONI
gruppo1_filtDGE_keep<- gruppo1_filtDGE[keep, , keep.lib.sizes=FALSE]
gruppo1_filtDGE_keep=calcNormFactors(gruppo1_filtDGE_keep,method="TMM")
mds_gruppo1=plotMDS(gruppo1_filtDGE_keep,dim.plot=c(1,2))
mds2_gruppo1=data.frame(x=mds_gruppo1$x,y=mds_gruppo1$y)
samples_gruppo1=rownames(mds2_gruppo1)
samples2_gruppo1=unlist(strsplit(samples_gruppo1,split="_kallisto"))
samples2_gruppo1=samples2_gruppo1[samples2_gruppo1!=""]
mds2_gruppo1=cbind(mds2_gruppo1,as.data.frame(samples2_gruppo1),as.data.frame(condition_gruppo1))
patient_gruppo1=rep(1:6,each=2)
mds2_gruppo1=cbind(mds2_gruppo1,as.data.frame(patient_gruppo1))
colnames(mds2_gruppo1)=c("x","y","samples","condition","patient")
design_gruppo1 = model.matrix (~condition + as.factor(patient),data=mds2_gruppo1)
gruppo1_filtDGE_keep <- estimateDisp(gruppo1_filtDGE_keep,design_gruppo1)
fit_gruppo1 <- glmFit(gruppo1_filtDGE_keep, design_gruppo1)
lrt_gruppo1 <- glmLRT(fit_gruppo1,coef=2)
#Correggo per qvalue e seleziono solo i significativi per qvalue<=0.25
library(qvalue)
toptags_xypat_Kallisto_gruppo1=as.data.frame(topTags(lrt_gruppo1,sort.by="PValue",n="all"))
q_gruppo1 = qvalue(toptags_xypat_Kallisto_gruppo1$PValue)
toptags_xypat_Kallisto_gruppo1$Qvalue = q_gruppo1$qvalues
toptags_xypat_Kallisto_DEG_gruppo1=toptags_xypat_Kallisto_gruppo1[toptags_xypat_Kallisto_gruppo1$Qvalue<=0.1,]
nrow(toptags_xypat_Kallisto_DEG_gruppo1)
2017

Gruppo1_cpm=cbind(cpm(gruppo1_filtDGE_keep,log=T),as.data.frame(gruppo1_filtDGE_keep$genes))

colnames(Gruppo1_cpm)
[1] "ALFE1_kallisto"  "ALFE2_kallisto"  "CALU1_kallisto"  "CALU2_kallisto"
[5] "DEST1_kallisto"  "DEST2_kallisto"  "GAGRA1_kallisto" "GAGRA2_kallisto"
[9] "MOGE1_kallisto"  "MOGE2_kallisto"  "PRELU1_kallisto" "PRELU2_kallisto"
[13] "genes"

Gruppo1_cpm$ALFE_FC=Gruppo1_cpm$ALFE2_kallisto - Gruppo1_cpm$ALFE1_kallisto
Gruppo1_cpm$CALU_FC=Gruppo1_cpm$CALU2_kallisto - Gruppo1_cpm$CALU1_kallisto
Gruppo1_cpm$DEST_FC=Gruppo1_cpm$DEST2_kallisto - Gruppo1_cpm$DEST1_kallisto
Gruppo1_cpm$GAGRA_FC=Gruppo1_cpm$GAGRA2_kallisto - Gruppo1_cpm$GAGRA1_kallisto
Gruppo1_cpm$MOGE_FC=Gruppo1_cpm$MOGE2_kallisto - Gruppo1_cpm$MOGE1_kallisto
Gruppo1_cpm$PRELU_FC=Gruppo1_cpm$PRELU2_kallisto - Gruppo1_cpm$PRELU1_kallisto

ALFE_RNK_Grp1=Gruppo1_cpm[(Gruppo1_cpm$ALFE2_kallisto>= 0.58 | Gruppo1_cpm$ALFE1_kallisto>=0.58) & (Gruppo1_cpm$ALFE_FC<= -0.58 | Gruppo1_cpm$ALFE_FC >=0.58),c(1,2,14,13)]
CALU_RNK_Grp1=Gruppo1_cpm[(Gruppo1_cpm$CALU2_kallisto>= 0.58 | Gruppo1_cpm$CALU1_kallisto>=0.58) & (Gruppo1_cpm$CALU_FC<= -0.58 | Gruppo1_cpm$CALU_FC >=0.58),c(3,4,15,13)]
DEST_RNK_Grp1=Gruppo1_cpm[(Gruppo1_cpm$DEST2_kallisto>= 0.58 | Gruppo1_cpm$DEST1_kallisto>=0.58) & (Gruppo1_cpm$DEST_FC<= -0.58 | Gruppo1_cpm$DEST_FC >=0.58),c(5,6,16,13)]
GAGRA_RNK_Grp1=Gruppo1_cpm[(Gruppo1_cpm$GAGRA2_kallisto>= 0.58 | Gruppo1_cpm$GAGRA1_kallisto>=0.58) & (Gruppo1_cpm$GAGRA_FC<= -0.58 | Gruppo1_cpm$GAGRA_FC >=0.58),c(7,8,17,13)]
MOGE_RNK_Grp1=Gruppo1_cpm[(Gruppo1_cpm$MOGE2_kallisto>= 0.58 | Gruppo1_cpm$MOGE1_kallisto>=0.58) & (Gruppo1_cpm$MOGE_FC<= -0.58 | Gruppo1_cpm$MOGE_FC >=0.58),c(9,10,18,13)]
PRELU_RNK_Grp1=Gruppo1_cpm[(Gruppo1_cpm$PRELU2_kallisto>= 0.58 | Gruppo1_cpm$PRELU1_kallisto>=0.58) & (Gruppo1_cpm$PRELU_FC<= -0.58 | Gruppo1_cpm$PRELU_FC >=0.58),c(11,12,19,13)]
Geni_gruppo1=as.data.frame(Gruppo1_cpm$genes)
colnames(Geni_gruppo1)="genes"

library(dplyr)
dfs=list(Geni_gruppo1,ALFE_RNK_Grp1,CALU_RNK_Grp1,DEST_RNK_Grp1,GAGRA_RNK_Grp1,MOGE_RNK_Grp1,PRELU_RNK_Grp1)
Gruppo1_RNKs=join_all(dfs,by="genes")


dfs=list(as.data.frame(Gruppo1_cpm$genes),ALFE_RNK_Grp1,CALU_RNK_Grp1,DEST_RNK_Grp1,GAGRA_RNK_Grp1,MOGE_RNK_Grp1,PRELU_RNK_Grp1)
Gruppo1_cpm_DEG=merge(Gruppo1_cpm,toptags_xypat_Kallisto_DEG_gruppo1,by="genes")
pdf("/home/fsantaniello/Heatmap_DEGs_Gruppo1_FCs.pdf")
breakList=seq(-5, 5, by = 0.1) 
pheatmap( Gruppo1_cpm_DEG [,14:19], cluster_rows=T,cluster_cols=T,color=colorRampPalette( c("green", "black", "red"), space="rgb")(length(breakList)))
dev.off()



condition_gruppo2=rep(c("Diagnosis","Relapse"),4)
gruppo2_filtDGE=DGEList(Gruppo2,genes=rownames(Gruppo2),group=condition_gruppo2)
keep <- rowSums(cpm(gruppo2_filtDGE)>1) >= 2 #PRENDO SOLO I GENI CON CPM>1 IN ALMENO 8 CAMPIONI
gruppo2_filtDGE_keep<- gruppo2_filtDGE[keep, , keep.lib.sizes=FALSE]
gruppo2_filtDGE_keep=calcNormFactors(gruppo2_filtDGE_keep,method="TMM")
mds_gruppo2=plotMDS(gruppo2_filtDGE_keep,dim.plot=c(1,2))
mds2_gruppo2=data.frame(x=mds_gruppo2$x,y=mds_gruppo2$y)
samples_gruppo2=rownames(mds2_gruppo2)
samples2_gruppo2=unlist(strsplit(samples_gruppo2,split="_kallisto"))
samples2_gruppo2=samples2_gruppo2[samples2_gruppo2!=""]
mds2_gruppo2=cbind(mds2_gruppo2,as.data.frame(samples2_gruppo2),as.data.frame(condition_gruppo2))
patient_gruppo2=rep(1:4,each=2)
mds2_gruppo2=cbind(mds2_gruppo2,as.data.frame(patient_gruppo2))
colnames(mds2_gruppo2)=c("x","y","samples","condition","patient")
design_gruppo2 = model.matrix (~condition + x + y + as.factor(patient),data=mds2_gruppo2)
gruppo2_filtDGE_keep <- estimateDisp(gruppo2_filtDGE_keep,design_gruppo2)
fit_gruppo2 <- glmFit(gruppo2_filtDGE_keep, design_gruppo2)
lrt_gruppo2 <- glmLRT(fit_gruppo2,coef=2)
#Correggo per qvalue e seleziono solo i significativi per qvalue<=0.25
library(qvalue)
toptags_xypat_Kallisto_gruppo2=as.data.frame(topTags(lrt_gruppo2,sort.by="PValue",n="all"))
q_gruppo2 = qvalue(toptags_xypat_Kallisto_gruppo2$PValue)
toptags_xypat_Kallisto_gruppo2$Qvalue = q_gruppo2$qvalues
toptags_xypat_Kallisto_DEG_gruppo2=toptags_xypat_Kallisto_gruppo2[toptags_xypat_Kallisto_gruppo2$Qvalue<=0.1,]
nrow(toptags_xypat_Kallisto_DEG_gruppo2)
184

Gruppo2_cpm=cbind(cpm(gruppo2_filtDGE_keep,log=T),as.data.frame(gruppo2_filtDGE_keep$genes))

colnames(Gruppo2_cpm)
[1] "BESU1_kallisto"  "BESU2_kallisto"  "DEIV2_kallisto"  "DEIV2b_kallisto"
[5] "MABI1_kallisto"  "MABI2_kallisto"  "PIAG1_kallisto"  "PIAG2_kallisto"
[9] "genes"

Gruppo2_cpm_DEGs=merge(Gruppo2_cpm,toptags_xypat_Kallisto_DEG_gruppo2,by="genes")


Gruppo2_cpm$BESU_FC=Gruppo2_cpm$BESU2_kallisto - Gruppo2_cpm$BESU1_kallisto
Gruppo2_cpm$DEIV_FC=Gruppo2_cpm$DEIV2b_kallisto - Gruppo2_cpm$DEIV2_kallisto
Gruppo2_cpm$MABI_FC=Gruppo2_cpm$MABI2_kallisto - Gruppo2_cpm$MABI1_kallisto
Gruppo2_cpm$PIAG_FC=Gruppo2_cpm$PIAG2_kallisto - Gruppo2_cpm$PIAG1_kallisto

colnames(Gruppo2_cpm)
[1] "BESU1_kallisto"  "BESU2_kallisto"  "DEIV2_kallisto"  "DEIV2b_kallisto"
[5] "MABI1_kallisto"  "MABI2_kallisto"  "PIAG1_kallisto"  "PIAG2_kallisto"
[9] "genes"           "BESU_FC"         "DEIV_FC"         "MABI_FC"
[13] "PIAG_FC"


BESU_RNK_Grp2=Gruppo2_cpm[(Gruppo2_cpm$BESU2_kallisto>= 0.58 | Gruppo2_cpm$BESU1_kallisto>=0.58) & (Gruppo2_cpm$BESU_FC<= -0.58 | Gruppo2_cpm$BESU_FC >=0.58),c(1,2,10,9)]
DEIV_RNK_Grp2=Gruppo2_cpm[(Gruppo2_cpm$DEIV2b_kallisto>= 0.58 | Gruppo2_cpm$DEIV2_kallisto>=0.58) & (Gruppo2_cpm$DEIV_FC<= -0.58 | Gruppo2_cpm$DEIV_FC >=0.58),c(3,4,11,9)]
MABI_RNK_Grp2=Gruppo2_cpm[(Gruppo2_cpm$MABI2_kallisto>= 0.58 | Gruppo2_cpm$MABI1_kallisto>=0.58) & (Gruppo2_cpm$MABI_FC<= -0.58 | Gruppo2_cpm$MABI_FC >=0.58),c(5,6,12,9)]
PIAG_RNK_Grp2=Gruppo2_cpm[(Gruppo2_cpm$PIAG2_kallisto>= 0.58 | Gruppo2_cpm$PIAG1_kallisto>=0.58) & (Gruppo2_cpm$PIAG_FC<= -0.58 | Gruppo2_cpm$PIAG_FC >=0.58),c(7,8,13,9)]
Geni_gruppo2=as.data.frame(Gruppo2_cpm$genes)
colnames(Geni_gruppo2)="genes"

library(dplyr)
dfs=list(Geni_gruppo2,BESU_RNK_Grp2,DEIV_RNK_Grp2,MABI_RNK_Grp2,PIAG_RNK_Grp2)
Gruppo2_RNKs=join_all(dfs,by="genes")



condition_gruppo3=rep(c("Diagnosis","Relapse"),5)
gruppo3_filtDGE=DGEList(Gruppo3,genes=rownames(Gruppo3),group=condition_gruppo3)
keep <- rowSums(cpm(gruppo3_filtDGE)>1) >= 2 #PRENDO SOLO I GENI CON CPM>1 IN ALMENO 8 CAMPIONI
gruppo3_filtDGE_keep<- gruppo3_filtDGE[keep, , keep.lib.sizes=FALSE]
gruppo3_filtDGE_keep=calcNormFactors(gruppo3_filtDGE_keep,method="TMM")
mds_gruppo3=plotMDS(gruppo3_filtDGE_keep,dim.plot=c(1,2))
mds2_gruppo3=data.frame(x=mds_gruppo3$x,y=mds_gruppo3$y)
samples_gruppo3=rownames(mds2_gruppo3)
samples2_gruppo3=unlist(strsplit(samples_gruppo3,split="_kallisto"))
samples2_gruppo3=samples2_gruppo3[samples2_gruppo3!=""]
mds2_gruppo3=cbind(mds2_gruppo3,as.data.frame(samples2_gruppo3),as.data.frame(condition_gruppo3))
patient_gruppo3=rep(1:5,each=2)
mds2_gruppo3=cbind(mds2_gruppo3,as.data.frame(patient_gruppo3))
colnames(mds2_gruppo3)=c("x","y","samples","condition","patient")
design_gruppo3 = model.matrix (~condition + x + y + as.factor(patient),data=mds2_gruppo3)
gruppo3_filtDGE_keep <- estimateDisp(gruppo3_filtDGE_keep,design_gruppo3)
fit_gruppo3 <- glmFit(gruppo3_filtDGE_keep, design_gruppo3)
lrt_gruppo3 <- glmLRT(fit_gruppo3,coef=2)
#Correggo per qvalue e seleziono solo i significativi per qvalue<=0.25
library(qvalue)
toptags_xypat_Kallisto_gruppo3=as.data.frame(topTags(lrt_gruppo3,sort.by="PValue",n="all"))
q_gruppo3 = qvalue(toptags_xypat_Kallisto_gruppo3$PValue)
toptags_xypat_Kallisto_gruppo3$Qvalue = q_gruppo3$qvalues
toptags_xypat_Kallisto_DEG_gruppo3=toptags_xypat_Kallisto_gruppo3[toptags_xypat_Kallisto_gruppo3$Qvalue<=0.1,]
nrow(toptags_xypat_Kallisto_DEG_gruppo3)
370

Gruppo3_cpm=cbind(cpm(gruppo3_filtDGE_keep,log=T),as.data.frame(gruppo3_filtDGE_keep$genes))
colnames(Gruppo3_cpm)
[1] "DR1_1_kallisto" "DR1_2_kallisto" "DR41_kallisto"  "DR42_kallisto"
[5] "DR51_kallisto"  "DR52_kallisto"  "FOCA1_kallisto" "FOCA2_kallisto"
[9] "LUAN1_kallisto" "LUAN2_kallisto" "genes"

Gruppo3_cpm$DR1_FC=Gruppo3_cpm$DR1_2_kallisto - Gruppo3_cpm$DR1_1_kallisto
Gruppo3_cpm$DR4_FC=Gruppo3_cpm$DR42_kallisto - Gruppo3_cpm$DR41_kallisto
Gruppo3_cpm$DR5_FC=Gruppo3_cpm$DR52_kallisto - Gruppo3_cpm$DR51_kallisto
Gruppo3_cpm$FOCA_FC=Gruppo3_cpm$FOCA2_kallisto - Gruppo3_cpm$FOCA1_kallisto
Gruppo3_cpm$LUAN_FC=Gruppo3_cpm$LUAN2_kallisto - Gruppo3_cpm$LUAN1_kallisto


colnames(Gruppo3_cpm)
[1] "DR1_1_kallisto" "DR1_2_kallisto" "DR41_kallisto"  "DR42_kallisto"
[5] "DR51_kallisto"  "DR52_kallisto"  "FOCA1_kallisto" "FOCA2_kallisto"
[9] "LUAN1_kallisto" "LUAN2_kallisto" "genes"          "DR1_FC"
[13] "DR4_FC"         "DR5_FC"         "FOCA_FC"        "LUAN_FC"

Gruppo3_cpm_DEGs=merge(Gruppo3_cpm,toptags_xypat_Kallisto_DEG_gruppo3,by="genes")

DR1_RNK_Grp3=Gruppo3_cpm[(Gruppo3_cpm$DR1_2_kallisto>= 0.58 | Gruppo3_cpm$DR1_1_kallisto>=0.58) & (Gruppo3_cpm$DR1_FC<= -0.58 | Gruppo3_cpm$DR1_FC >=0.58),c(1,2,12,11)]
DR4_RNK_Grp3=Gruppo3_cpm[(Gruppo3_cpm$DR42_kallisto>= 0.58 | Gruppo3_cpm$DR41_kallisto>=0.58) & (Gruppo3_cpm$DR4_FC<= -0.58 | Gruppo3_cpm$DR4_FC >=0.58),c(3,4,13,11)]
DR5_RNK_Grp3=Gruppo3_cpm[(Gruppo3_cpm$DR52_kallisto>= 0.58 | Gruppo3_cpm$DR51_kallisto>=0.58) & (Gruppo3_cpm$DR5_FC<= -0.58 | Gruppo3_cpm$DR5_FC >=0.58),c(5,6,14,11)]
FOCA_RNK_Grp3=Gruppo3_cpm[(Gruppo3_cpm$FOCA2_kallisto>= 0.58 | Gruppo3_cpm$FOCA1_kallisto>=0.58) & (Gruppo3_cpm$FOCA_FC<= -0.58 | Gruppo3_cpm$FOCA_FC >=0.58),c(7,8,15,11)]
LUAN_RNK_Grp3=Gruppo3_cpm[(Gruppo3_cpm$LUAN2_kallisto>= 0.58 | Gruppo3_cpm$LUAN1_kallisto>=0.58) & (Gruppo3_cpm$LUAN_FC<= -0.58 | Gruppo3_cpm$LUAN_FC >=0.58),c(9,10,16,11)]


Geni_gruppo3=as.data.frame(Gruppo3_cpm$genes)
colnames(Geni_gruppo3)="genes"

library(dplyr)
dfs=list(Geni_gruppo3,DR1_RNK_Grp3,DR4_RNK_Grp3,DR5_RNK_Grp3,FOCA_RNK_Grp3,LUAN_RNK_Grp3)
Gruppo3_RNKs=join_all(dfs,by="genes")


library(biomart)
mart=useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="hsapiens_gene_ensembl")

G_list_gruppo1 <- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=toptags_xypat_Kallisto_DEG_gruppo1$genes,mart=mart)
colnames(G_list_gruppo1)=c("genes","entrezgene")
toptags_xypat_Kallisto_DEG_gruppo1_entrez=merge(G_list_gruppo1,toptags_xypat_Kallisto_DEG_gruppo1,by="genes")
write.table(toptags_xypat_Kallisto_DEG_gruppo1_entrez,'Geni_DEGs_Gruppo1_AllInfos_Entrez.txt',col.names=T,row.names=F,quote=F)


G_list_gruppo2 <- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=toptags_xypat_Kallisto_DEG_gruppo2$genes,mart=mart)
colnames(G_list_gruppo2)=c("genes","entrezgene")
toptags_xypat_Kallisto_DEG_gruppo2_entrez=merge(G_list_gruppo2,toptags_xypat_Kallisto_DEG_gruppo2,by="genes")
write.table(toptags_xypat_Kallisto_DEG_gruppo3_entrez,'Geni_DEGs_Gruppo2_AllInfos_Entrez.txt',col.names=T,row.names=F,quote=F)



G_list_gruppo3 <- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=toptags_xypat_Kallisto_DEG_gruppo3$genes,mart=mart)
colnames(G_list_gruppo3)=c("genes","entrezgene")
toptags_xypat_Kallisto_DEG_gruppo3_entrez=merge(G_list_gruppo3,toptags_xypat_Kallisto_DEG_gruppo3,by="genes")
write.table(toptags_xypat_Kallisto_DEG_gruppo3_entrez,'Geni_DEGs_Gruppo3_AllInfos_Entrez.txt',col.names=T,row.names=F,quote=F)


write.table(as.data.frame(toptags_xypat_Kallisto_DEG_gruppo3[,1]),"/home/fsantaniello/Geni_DEG_Gruppo3.txt",sep="\t",quote=F,col.names=F,row.names=F)
write.table(as.data.frame(toptags_xypat_Kallisto_DEG_gruppo2[,1]),"/home/fsantaniello/Geni_DEG_Gruppo2.txt",sep="\t",quote=F,col.names=F,row.names=F)
write.table(as.data.frame(toptags_xypat_Kallisto_DEG_gruppo1[,1]),"/home/fsantaniello/Geni_DEG_Gruppo1.txt",sep="\t",quote=F,col.names=F,row.names=F)

write.table(as.data.frame(G_list$entrezgene),"/home/fsantaniello/Background_Geni_Gruppo3.txt",sep="\t",col.names=F,row.names=F,quote=F)

Gruppo_3_2_Aggregato=countdata_filt[,c(3,4,7,8,11,12,13,14,15,16,17,18,21,22,23,24,27,28)]
condition_gruppo3_2_aggregato=rep(c("Diagnosis","Relapse"),9)
gruppo3_2_aggregato_filtDGE=DGEList(Gruppo_3_2_Aggregato,genes=rownames(Gruppo_3_2_Aggregato),group=condition_gruppo3_2_aggregato)
keep <- rowSums(cpm(gruppo3_2_aggregato_filtDGE)>1) >= 4 #PRENDO SOLO I GENI CON CPM>1 IN ALMENO 8 CAMPIONI
gruppo3_2_aggregato_filtDGE_keep<- gruppo3_2_aggregato_filtDGE[keep, , keep.lib.sizes=FALSE]
gruppo3_2_aggregato_filtDGE_keep=calcNormFactors(gruppo3_2_aggregato_filtDGE_keep,method="TMM")
mds_gruppo3_2_aggregato=plotMDS(gruppo3_2_aggregato_filtDGE_keep,dim.plot=c(1,2))
mds2_gruppo3_2_aggregato=data.frame(x=mds_gruppo3_2_aggregato$x,y=mds_gruppo3_2_aggregato$y)
samples_gruppo3_2_aggregato=rownames(mds2_gruppo3_2_aggregato)
samples2_gruppo3_2_aggregato=unlist(strsplit(samples_gruppo3_2_aggregato,split="_kallisto"))
samples2_gruppo3_2_aggregato=samples2_gruppo3_2_aggregato[samples2_gruppo3_2_aggregato!=""]
mds2_gruppo3_2_aggregato=cbind(mds2_gruppo3_2_aggregato,as.data.frame(samples2_gruppo3_2_aggregato),as.data.frame(condition_gruppo3_2_aggregato))
patient_gruppo3_2_aggregato=rep(1:9,each=2)
mds2_gruppo3_2_aggregato=cbind(mds2_gruppo3_2_aggregato,as.data.frame(patient_gruppo3_2_aggregato))
colnames(mds2_gruppo3_2_aggregato)=c("x","y","samples","condition","patient")
design_gruppo3_2_aggregato = model.matrix (~condition  + x + y + as.factor(patient),data=mds2_gruppo3_2_aggregato)
gruppo3_2_aggregato_filtDGE_keep <- estimateDisp(gruppo3_2_aggregato_filtDGE_keep,design_gruppo3_2_aggregato)
fit_gruppo3_2_aggregato <- glmFit(gruppo3_2_aggregato_filtDGE_keep, design_gruppo3_2_aggregato)
lrt_gruppo3_2_aggregato <- glmLRT(fit_gruppo3_2_aggregato,coef=2)
#Correggo per qvalue e seleziono solo i significativi per qvalue<=0.25
library(qvalue)
toptags_xypat_Kallisto_gruppo3_2_aggregato=as.data.frame(topTags(lrt_gruppo3_2_aggregato,sort.by="PValue",n="all"))
q_gruppo3_2_aggregato = qvalue(toptags_xypat_Kallisto_gruppo3_2_aggregato$PValue)
toptags_xypat_Kallisto_gruppo3_2_aggregato$Qvalue = q_gruppo3_2_aggregato$qvalues
toptags_xypat_Kallisto_DEG_gruppo3_2_aggregato=toptags_xypat_Kallisto_gruppo3_2_aggregato[toptags_xypat_Kallisto_gruppo3_2_aggregato$Qvalue<=0.25,]

 nrow(toptags_xypat_Kallisto_DEG_gruppo3_2_aggregato)
[1] 308

write.table(as.data.frame(toptags_xypat_Kallisto_DEG_gruppo3_2_aggregato$genes),"Geni_DEGs_Aggregazione_Gruppo3Gruppo2_redone.txt",sep="\t",col.names=F,row.names=F,quote=F)

library(biomart)
mart=useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="hsapiens_gene_ensembl")
G_list_gruppo3 <- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=toptags_xypat_Kallisto_DEG_gruppo3$genes,mart=mart)
write.table(as.data.frame(G_list_gruppo3$entrezgene),"/home/fsantaniello/Background_Geni_Gruppo3.txt",sep="\t",col.names=F,row.names=F,quote=F)

G_list_gruppo2 <- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=toptags_xypat_Kallisto_DEG_gruppo2$genes,mart=mart)
write.table(as.data.frame(G_list_gruppo2$hgnc_symbol),"/home/fsantaniello/DEGs_gruppo2_GabriLucia.txt",sep="\t",col.names=F,row.names=F,quote=F)

G_list_gruppo1 <- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=toptags_xypat_Kallisto_DEG_gruppo1$genes,mart=mart)
write.table(as.data.frame(G_list_gruppo1$hgnc_symbol),"/home/fsantaniello/DEGs_gruppo1_GabriLucia.txt",sep="\t",col.names=F,row.names=F,quote=F)
library(biomart)

G_list_gruppo3_2_aggregato <- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=toptags_xypat_Kallisto_DEG_gruppo3_2_aggregato$genes,mart=mart)
write.table(as.data.frame(G_list_gruppo3_2_aggregato$hgnc_symbol),"/home/fsantaniello/DEGs_Gruppo3_Gruppo2_Aggregato.txt",sep="\t",col.names=F,row.names=F,quote=F)

All_16kGenes_GabriLucia <- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=toptags_xypat_Kallisto$genes,mart=mart)
write.table(as.data.frame(All_16kGenes_GabriLucia$entrezgene),"/home/fsantaniello/Background_All_16kGenes_GabriLucia.txt",sep="\t",col.names=F,row.names=F,quote=F)


Gruppo1_GO=read.delim("GeneOntology_Gruppo1_DEGs.txt",head=T,stringsAsFactors=F,sep="\t",row.names=NULL)
Gruppo1_GO$ID=gsub("\\~.*","",Gruppo1_GO$Term)
Gruppo1_GO$Term2=gsub(".*~","",Gruppo1_GO$Term)
David_Gruppo1=Gruppo1_GO[,c(1,14,15,12,6)]
colnames(David_Gruppo1)=c("category","ID","Term","adj_pval","genes")

Genes_DEG_gruppo1=toptags_xypat_Kallisto_DEG_gruppo1[,c(1,2)]
colnames(Genes_DEG_gruppo1)=c("ID","logFC")
circ <- circle_dat(David_Gruppo1,Genes_DEG_gruppo1)


David_Gruppo1=Gruppo1_GO[,c(1,14,15,12,6)]
colnames(David_Gruppo3)=c("category","ID","Term","adj_pval","genes")
Genes_DEG_gruppo3=toptags_xypat_Kallisto_DEG_gruppo3[,c(1,2)]
colnames(Genes_DEG_gruppo3)=c("ID","logFC")



Gruppo3_Gruppo2_GO=read.delim("GO_Gruppo3_Gruppo2_Aggregati.txt",head=T,stringsAsFactors=F,sep="\t",row.names=NULL)
Gruppo3_Gruppo2_GO$ID=gsub("\\~.*","",Gruppo3_Gruppo2_GO$Term)
Gruppo3_Gruppo2_GO$Term2=gsub(".*~","",Gruppo3_Gruppo2_GO$Term)
David_Gruppo3_Gruppo2=Gruppo3_Gruppo2_GO[,c(1,14,15,12,6)]
colnames(David_Gruppo3_Gruppo2)=c("category","ID","Term","adj_pval","genes")
David_Gruppo3_Gruppo2=David_Gruppo3_Gruppo2[David_Gruppo3_Gruppo2$adj_pval<=0.1,]
Genes_DEG_gruppo3_gruppo2=toptags_xypat_Kallisto_DEG_gruppo3_2_aggregato[,c(1,2)]
colnames(Genes_DEG_gruppo3_gruppo2)=c("ID","logFC")
circ <- circle_dat(David_Gruppo3_Gruppo2,Genes_DEG_gruppo3_gruppo2)

GOBar(circ, display = 'multiple', title = 'Z-score coloured barplot', zsc.col = c('yellow', 'black', 'cyan'))







geni_GOs_Gruppi=c('CLSPN',
'MCM10',
'CDT1',
'CDC45',
'MCM7',
'ORC1',
'RTEL1',
'CDC7',
'CDK1',
'CDC6',
'DTL',
'LIG1',
'MCM2',
'RNASEH2A',
'MCM3',
'MCM4',
'MCM5',
'RAD50',
'MCM6',
'RFC5',
'RFC3',
'RFC4',
'TIMELESS',
'RFC2',
'RRM2',
'TBRG1',
'RRM1',
'DSCC1',
'BLM',
'TIPIN',
'POLA1',
'CHEK1',
'POLA2',
'RPA3',
'DTD1',
'POLE2',
'FEN1',
'EXO1',
'RECQL4',
'GINS2',
'SSRP1',
'GINS3',
'GINS4',
'BRIP1',
'BRCA1',
'CDC25A',
'DNA2',
'POLD1',
'POLD2',
'PCNA',
'CHTF18',
'CHAF1A',
'CHAF1B',
'RBM14-RBM4',
'REPIN1',
'DUT',
'CDC7',
'PRIM1',
'CCNE1',
'PRIM2',
'PRIM1',
'PTTG1',
'DDX11',
'FANCI',
'NSMCE1',
'FANCA',
'GTF2H4',
'EEPD1',
'RAD51',
'NABP1',
'UHRF1',
'FAAP100',
'RFWD3',
'UCHL5',
'RAD18',
'RUVBL2',
'RUVBL1',
'UBE2T',
'PPP5C',
'FOXM1',
'UNG',
'BCCIP',
'CHD1L',
'NPM1',
'FBXO6',
'POLQ',
'APEX1',
'ERCC1',
'MSH6',
'RAD51AP1',
'RAD54L',
'MNAT1',
'PARPBP',
'ALKBH3',
'CENPX',
'HLA-DQB1',
'HLA-DRB1',
'HLA-DRB5',
'HLA-DPA1',
'HLA-DPB1',
'HLA-DMB',
'HLA-DOA',
'HLA-DMA',
'HLA-DRA')

colnames(geni_GOs_Gruppi)='genes'

toptags_xypat_Kallisto_cpm_DEG_goGruppi=merge(geni_GOs_Gruppi,toptags_xypat_Kallisto_cpm_DEG,by='genes',sort=F)

breakList=seq(-5, 5, by = 0.1) 

pheatmap( toptags_xypat_Kallisto_cpm_DEG_goGruppi[,38:52], cluster_rows=F,cluster_cols=T,color=colorRampPalette( c("cyan", "gray15", "yellow1"), space="rgb")(length(breakList)),breaks=breakList)

colnames(Gruppo1_cpm)
[1] "ALFE1_kallisto"  "ALFE2_kallisto"  "CALU1_kallisto"  "CALU2_kallisto"
[5] "DEST1_kallisto"  "DEST2_kallisto"  "GAGRA1_kallisto" "GAGRA2_kallisto"
[9] "MOGE1_kallisto"  "MOGE2_kallisto"  "PRELU1_kallisto" "PRELU2_kallisto"
[13] "genes"

Gruppo1_cpm$ALFE_FC=Gruppo1_cpm$ALFE2_kallisto - Gruppo1_cpm$ALFE1_kallisto
Gruppo1_cpm$CALU_FC=Gruppo1_cpm$CALU2_kallisto - Gruppo1_cpm$CALU1_kallisto
Gruppo1_cpm$DEST_FC=Gruppo1_cpm$DEST2_kallisto - Gruppo1_cpm$DEST1_kallisto
Gruppo1_cpm$GAGRA_FC=Gruppo1_cpm$GAGRA2_kallisto - Gruppo1_cpm$GAGRA1_kallisto
Gruppo1_cpm$MOGE_FC=Gruppo1_cpm$MOGE2_kallisto - Gruppo1_cpm$MOGE1_kallisto
Gruppo1_cpm$PRELU_FC=Gruppo1_cpm$PRELU2_kallisto - Gruppo1_cpm$PRELU1_kallisto

Gruppo1_cpm_DEG=merge(Gruppo1_cpm,toptags_xypat_Kallisto_DEG_gruppo1,by="genes")
pdf("/home/fsantaniello/Heatmap_DEGs_Gruppo1_FCs.pdf")
breakList=seq(-5, 5, by = 0.1) 
pheatmap( Gruppo1_cpm_DEG [,c(2,4,6,8,10,12,3,5,7,9,11,13)], cluster_rows=T,cluster_cols=F,color=colorRampPalette( c('yellow', 'black', 'cyan'), space="rgb")(length(breakList)),breaks=breakList,scale='row')
dev.off()

Gruppo2_cpm=cbind(cpm(gruppo2_filtDGE_keep,log=T),as.data.frame(gruppo2_filtDGE_keep$genes))

colnames(Gruppo2_cpm)
[1] "BESU1_kallisto"  "BESU2_kallisto"  "DEIV2_kallisto"  "DEIV2b_kallisto"
[5] "MABI1_kallisto"  "MABI2_kallisto"  "PIAG1_kallisto"  "PIAG2_kallisto"
[9] "genes"


Gruppo2_cpm$BESU_FC=Gruppo2_cpm$BESU2_kallisto - Gruppo2_cpm$BESU1_kallisto
Gruppo2_cpm$DEIV_FC=Gruppo2_cpm$DEIV2b_kallisto - Gruppo2_cpm$DEIV2_kallisto
Gruppo2_cpm$MABI_FC=Gruppo2_cpm$MABI2_kallisto - Gruppo2_cpm$MABI1_kallisto
Gruppo2_cpm$PIAG_FC=Gruppo2_cpm$PIAG2_kallisto - Gruppo2_cpm$PIAG1_kallisto

Gruppo2_cpm_DEG=merge(Gruppo2_cpm,toptags_xypat_Kallisto_DEG_gruppo2,by="genes")
pdf("/home/fsantaniello/Heatmap_DEGs_Gruppo2_FCs.pdf")
breakList=seq(-5, 5, by = 0.1) 
pheatmap( Gruppo2_cpm_DEG [,c(2,4,6,8,3,5,7,9)], cluster_rows=T,cluster_cols=F,color=colorRampPalette( c('yellow', 'black', 'cyan'), space="rgb")(length(breakList)),breaks=breakList)
dev.off()

Gruppo3_cpm=cbind(cpm(gruppo3_filtDGE_keep,log=T),as.data.frame(gruppo3_filtDGE_keep$genes))
colnames(Gruppo3_cpm)
[1] "DR1_1_kallisto" "DR1_2_kallisto" "DR41_kallisto"  "DR42_kallisto"
[5] "DR51_kallisto"  "DR52_kallisto"  "FOCA1_kallisto" "FOCA2_kallisto"
[9] "LUAN1_kallisto" "LUAN2_kallisto" "genes"

Gruppo3_cpm$DR1_FC=Gruppo3_cpm$DR1_2_kallisto - Gruppo3_cpm$DR1_1_kallisto
Gruppo3_cpm$DR4_FC=Gruppo3_cpm$DR42_kallisto - Gruppo3_cpm$DR41_kallisto
Gruppo3_cpm$DR5_FC=Gruppo3_cpm$DR52_kallisto - Gruppo3_cpm$DR51_kallisto
Gruppo3_cpm$FOCA_FC=Gruppo3_cpm$FOCA2_kallisto - Gruppo3_cpm$FOCA1_kallisto
Gruppo3_cpm$LUAN_FC=Gruppo3_cpm$LUAN2_kallisto - Gruppo3_cpm$LUAN1_kallisto

Gruppo3_cpm_DEG=merge(Gruppo3_cpm,toptags_xypat_Kallisto_DEG_gruppo3,by="genes")
pdf("/home/fsantaniello/Heatmap_DEGs_Gruppo3_FCs.pdf")
breakList=seq(-5, 5, by = 0.1) 
pheatmap( Gruppo3_cpm_DEG [,c(2,4,6,8,10,3,5,7,9,11)], cluster_rows=T,cluster_cols=F,color=colorRampPalette( c('yellow', 'black', 'cyan'), space="rgb")(length(breakList)))
dev.off()

Gruppo_aggregato_23_cpm=cbind(cpm(gruppo3_2_aggregato_filtDGE_keep,log=T),as.data.frame(gruppo3_2_aggregato_filtDGE_keep$genes))
Gruppo_aggregato_23_cpm_DEG=merge(Gruppo3Gruppo_aggregato_23_cpm_cpm,toptags_xypat_Kallisto_DEG_gruppo3_2_aggregato,by="genes")
breakList=seq(-10, 10, by = 0.1) 
pheatmap( Gruppo_aggregato_23_cpm_DEG [,c(2,4,6,8,10,12,14,16,18,3,5,7,9,11,13,15,17,19)], cluster_rows=T,cluster_cols=F,color=colorRampPalette( c('yellow', 'black', 'cyan'), space="rgb")(length(breakList)),breaks=breakList)
dev.off()


library(org.Hs.eg.db)
library(GO.db)
library(GOSemSim)
library(clusterProfiler)
library(data.table)
library(plyr)
hsGO2 <- godata('org.Hs.eg.db', keytype = "SYMBOL", ont="BP", computeIC=FALSE)


ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)



GeneSimilarity_Gruppo1=mgeneSim(Gruppo1_cpm_DEG$genes,semData=hsGO2, measure="Wang", combine="BMA", verbose=FALSE)
Gruppo1_GO=read.delim("GeneOntology_Gruppo1_DEGs.txt",head=T,stringsAsFactors=F,sep="\t",row.names=NULL)
Gruppo1_GO= Gruppo1_GO[Gruppo1_GO$Benjamini<=0.25,]
geni_gruppo1=as.data.frame(unlist(strsplit(Gruppo1_GO$Genes,",")))
colnames(geni_gruppo1)="genes"


a_gruppo1=data.frame()
for (i in 1:nrow(geni_gruppo1)){a_gruppo1=rbind(a_gruppo1,cbind(as.data.frame(geni_gruppo1[i,]),as.data.frame(Gruppo1_GO[grep(geni_gruppo1[i,],Gruppo1_GO$Genes),2])))}
colnames(a_gruppo1)=c("genes","GO_TERM")
a_gruppo1=unique(a_gruppo1)
a_gruppo1$genes=as.character(a_gruppo1$genes)
a_gruppo1$genes=gsub(" ","",a_gruppo1$genes)
b_gruppo1=as.data.frame(as.character(rownames(GeneSimilarity_Gruppo1)))
colnames(b_gruppo1)="genes"
b_gruppo1$genes=as.character(b_gruppo1$genes)
a1_gruppo1=dcast(setDT(a_gruppo1), genes~rowid(genes, prefix="GO"), value.var="GO_TERM")
a1_gruppo1=a1_gruppo1[,1:2]
a1_gruppo1=as.data.frame(a1_gruppo1)
a1_gruppo1$genes=as.character(a1_gruppo1$genes)
a1_gruppo1$GO1=gsub(" ","_",a1_gruppo1$GO1)
a1_gruppo1$GO1=gsub(".*~","",a1_gruppo1$GO1)

#a1$GO1=gsub("\\*-~","",a1$GO1)

prova_gruppo1=merge(a1_gruppo1,b_gruppo1,by="genes",all.y=T)
prova_gruppo1$GO1[is.na(prova_gruppo1$GO1)]<-"NO_GO"
annotation_gruppo1=as.data.frame(prova_gruppo1$GO1)
colnames(annotation_gruppo1)="goterm"
rownames(annotation_gruppo1)=prova_gruppo1$genes
nrow(as.data.frame(unique(annotation_gruppo1)))
write.table(as.data.frame(unique(annotation_gruppo1)),"GO_Univoche_Gruppo1_DEGs.txt",sep="\t",col.names=F,row.names=F)

# sample(colors(), 24)
#  [1] "lightpink2"      "gold4"           "forestgreen"     "red4"
#  [5] "springgreen3"    "grey24"          "lightblue1"      "darksalmon"
#  [9] "olivedrab4"      "lightblue2"      "seashell4"       "indianred3"
# [13] "lightsalmon"     "aliceblue"       "slateblue3"      "orchid4"
# [17] "blue2"           "papayawhip"      "turquoise4"      "grey14"
# [21] "springgreen4"    "darkorchid"      "maroon4"         "gray76"
# [25] "grey34"          "grey32"          "darkviolet"      "grey33"
# [29] "ivory3"          "chartreuse4"     "blue1"           "cadetblue2"
# [33] "grey39"          "gray67"          "grey54"          "midnightblue"
# [37] "chartreuse3"     "burlywood2"      "moccasin"        "snow1"
# [41] "grey1"           "tomato"          "cyan1"           "grey49"
# [45] "gray36"          "mediumpurple"    "gray84"          "gray59"
# [49] "pink2"           "dodgerblue3"     "limegreen"       "orangered4"
# [53] "slateblue2"      "grey56"          "grey76"          "gray7"
# [57] "wheat4"          "yellow"          "lightsteelblue3" "coral"
# [61] "orchid"          "mistyrose"       "darkorchid2"     "mediumpurple2"
# [65] "cornsilk3"       "gray30"          "deeppink4"       "grey10"
# [69] "gray14"          "mediumslateblue" "brown1"          "grey98"
# [73] "slategray2"      "darkgray"        "magenta1"        "cornflowerblue"
# [77] "lightsalmon3"    "tan"             "aquamarine3"     "dimgrey"
# [81] "indianred2"      "gray19"          "hotpink2"        "royalblue4"
# [85] "gray95"          "grey59"          "peachpuff3"      "gray52"
# [89] "gray23"          "grey86"          "darkseagreen3"   "gray48"
# [93] "red2"            "brown"           "mediumturquoise" "gold1"
# [97] "seashell"        "gray74"          "lightyellow3"    "azure2"
#[101] "slateblue"       "orange4"         "turquoise2"      "orangered3"
#[105] "gray40"          "gray47"          "navajowhite3"    "wheat1"
#[109] "gray91"          "peachpuff4"      "turquoise3"      "violet"
#[113] "grey90"          "bisque3"         "grey80"          "mediumblue"
#[117] "grey77"          "grey85"          "cadetblue3"      "gray37"
#[121] "grey35"          "pink3"           "gray94"          "palegoldenrod"
#[125] "gray46"          "lightskyblue4"   "grey65"          "rosybrown"
#[129] "grey72"          "tan2"            "gray"            "gray62"
#[133] "grey78"          "dodgerblue"      "coral4"          "deeppink1"
#[137] "aquamarine"      "khaki3"          "grey63"          "thistle"
#[141] "lightpink3"      "yellow1"         "gray9"           "darkslategray3"
#[145] "grey51"

n=24
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
sample(col_vector, n)
write.table(as.data.frame(sample(col_vector, n)),"/home/fsantaniello/Colori_Pheatmap_Gruppo1.txt",sep="\t",col.names=F,row.names=F,quote=T)


pdf("/home/fsantaniello/GeneSemanticSimilarity_Heatmap_AnnotazioneGO_Gruppo1.pdf",width=60,height=60)

colors_gruppo1=list(goterm = c("NO_GO"="#F781BF"
,"response_to_lipopolysaccharide"="#E5D8BD"
,"G1/S_transition_of_mitotic_cell_cycle"="#FCCDE5"
,"purine_nucleotide_biosynthetic_process"="#FED9A6"
,"tRNA_aminoacylation_for_protein_translation"="#A6761D"
,"nucleobase-containing_compound_metabolic_process"="#7FC97F"
,"DNA_repair"="#66C2A5"
,"cell_division"="#FFFFB3"
,"cell_proliferation"="#386CB0"
,"DNA_duplex_unwinding"="#E6F5C9"
,"mitotic_nuclear_division"="#FFFF99"
,"mitochondrion_organization"="#984EA3"
,"DNA_replication"="#B3E2CD"
,"rRNA_processing"="#B2DF8A"
,"sister_chromatid_cohesion"="#CAB2D6"
,"spliceosomal_snRNP_assembly"="#CCCCCC"
,"DNA_biosynthetic_process"="#CCEBC5"
,"DNA_recombination"="#377EB8"
,"nuclear_import"="#33A02C"
,"RNA_catabolic_process"="#B3CDE3"
,"mismatch_repair"="#FDDAEC"
,"DNA_damage_response,_detection_of_DNA_damage"="#CCEBC5"
,"nucleotide-excision_repair,_DNA_incision,_5'-to_lesion"="#BEAED4"
,"telomere_maintenance_via_recombination"="#E41A1C"
))

pheatmap(GeneSimilarity_Gruppo1,cluster_rows=T,cluster_cols=T,border_color=NA, annotation_row=annotation_gruppo1,annotation_colors = colors_gruppo1[1])
dev.off()






pdf("/home/fsantaniello/GeneSemanticSimilarity_Heatmap_Gruppo1.pdf",width=30,height=30)
pheatmap(GeneSimilarity,cluster_rows=T,cluster_cols=T,border_color=NA)
dev.off()






GeneSimilarity_Gruppo2=mgeneSim(Gruppo2_cpm_DEG$genes,semData=hsGO2, measure="Wang", combine="BMA", verbose=FALSE)
pdf("/home/fsantaniello/GeneSemanticSimilarity_Heatmap_Gruppo2.pdf",width=30,height=30)
pheatmap(GeneSimilarity,cluster_rows=T,cluster_cols=T,border_color=NA)
dev.off()



GeneSimilarity_Gruppo3=mgeneSim(Gruppo3_cpm_DEG$genes,semData=hsGO2, measure="Wang", combine="BMA", verbose=FALSE)

Gruppo3_GO=read.table('GeneOntology_Gruppo3_DEGs.txt',head=T,stringsAsFactors=F,sep="\t")
Gruppo3_GO=Gruppo3_GO[Gruppo3_GO$Benjamini<=0.25,]
#Gruppo3_GO_Pval= Gruppo3_GO[Gruppo3_GO$Benjamini<=0.05,]
geni=as.data.frame(unlist(strsplit(Gruppo3_GO $Genes,",")))
colnames(geni)="genes"




a=data.frame()
for (i in 1:nrow(geni)){a=rbind(a,cbind(as.data.frame(geni[i,]),as.data.frame(Gruppo3_GO[grep(geni[i,],Gruppo3_GO$Genes),2])))}
colnames(a)=c("genes","GO_TERM")
a=unique(a)
a$genes=as.character(a$genes)
a$genes=gsub(" ","",a$genes)
b=as.data.frame(as.character(rownames(GeneSimilarity_Gruppo3)))
colnames(b)="genes"
b$genes=as.character(b$genes)
a1=dcast(setDT(a), genes~rowid(genes, prefix="GO"), value.var="GO_TERM")
a1=a1[,1:2]
a1=as.data.frame(a1)
a1$genes=as.character(a1$genes)
a1$GO1=gsub(" ","_",a1$GO1)
a1$GO1=gsub(".*~","",a1$GO1)

#a1$GO1=gsub("\\*-~","",a1$GO1)

prova=merge(a1,b,by="genes",all.y=T)
prova$GO1[is.na(prova$GO1)]<-"NO_GO"
annotation_gruppo3=as.data.frame(prova$GO1)
colnames(annotation_gruppo3)="goterm"
rownames(annotation_gruppo3)=prova$genes
nrow(as.data.frame(unique(annotation_gruppo3)))
write.table(as.data.frame(unique(annotation_gruppo3)),"GO_Univoche_Gruppo3_DEGs.txt",sep="\t",col.names=T)

library(RColorBrewer)
n=26
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
sample(col_vector, n)
1] "#BF5B17" "#8DD3C7" "#E78AC3" "#E7298A" "#FC8D62" "#666666" "#FB9A99"
8] "#E6AB02" "#CCEBC5" "#B2DF8A" "#A65628" "#FFFF99" "#B3E2CD" "#7570B3"
5] "#1B9E77" "#E41A1C" "#FFED6F" "#E6F5C9" "#CCEBC5" "#BEAED4" "#D9D9D9"
2] "#E5D8BD" "#666666" "#FED9A6" "#FF7F00" "#E31A1C" "#FFFF99" "#F1E2CC"
9] "#A6CEE3" "#80B1D3" "#4DAF4A" "#B3DE69" "#A6761D" "#F4CAE4" "#FDC086"
6] "#FFFF33" "#D95F02" "#7FC97F" "#8DA0CB" "#B15928" "#F2F2F2" "#F781BF"
3] "#FBB4AE" "#1F78B4" "#FDCDAC" "#E5C494" "#DECBE4" "#FF7F00" "#984EA3"
0] "#999999" "#377EB8" "#CAB2D6" "#6A3D9A" "#CCCCCC" "#FFFFCC" "#FDBF6F"
7] "#FFFFB3"
write.table(as.data.frame(sample(col_vector, n)),"/home/fsantaniello/Colori_Pheatmap_Gruppo3_GO.txt",sep="\t",col.names=F,row.names=F,quote=T)


pdf("/home/fsantaniello/GeneSemanticSimilarity_Heatmap_AnnotazioneGO_Gruppo3.pdf",width=60,height=60)

colors_gruppo3=list(goterm = c('platelet_degranulation'="#FDB462"
,'NO_GO'="#377EB8"
,'immune_response'="#666666"
,'receptor-mediated_endocytosis'="#BC80BD"
,'positive_regulation_of_blood_vessel_endothelial_cell_migration'="#FFD92F"
,'angiogenesis'="#BEBADA"
,'extracellular_matrix_organization'="#B3CDE3"
,'G-protein_coupled_receptor_signaling_pathway'="#E6F5C9"
,'cell_adhesion'="#7570B3"
,'positive_regulation_of_cell_proliferation'="#FF7F00"
,'ventricular_septum_morphogenesis'="#CBD5E8"
,'inflammatory_response'="#8DA0CB"
,'regulation_of_cell_proliferation'="#FDDAEC"
,'extracellular_matrix_disassembly'="#F4CAE4"
,'antigen_processing_and_presentation'="#7FC97F"
,'artery_morphogenesis'="#66C2A5"
,'blood_coagulation'="#CCEBC5"
,'oxygen_transport'="#FDC086"
,'antigen_processing_and_presentation_of_peptide_or_polysaccharide_antigen_via_MHC_class_II'="#E78AC3"
,'response_to_wounding'="#A6D854"
,'neutrophil_chemotaxis'="#FF7F00"
,'collagen_catabolic_process'="#A65628"
,'positive_regulation_of_mesenchymal_cell_proliferation'="#666666"
,'chemotaxis'="#F2F2F2"
,'positive_regulation_of_phosphatidylinositol_3-kinase_signaling'="#B2DF8A"
,'hydrogen_peroxide_catabolic_process'="#FDBF6F"))

pheatmap(GeneSimilarity_Gruppo3,cluster_rows=T,cluster_cols=T,border_color=NA, annotation_col=annotation_gruppo3,annotation_colors = colors_gruppo3[1])
dev.off()


go

ha_column = rowAnnotation(annotation_gruppo3, 
col = colors_gruppo3)



Arguments:

terms: A data frame with columns for 'category', 'ID', 'term',
      adjusted p-value ('adj_pval') and 'genes'

genes: A data frame with columns for 'ID', 'logFC'

Gruppo1_GO=read.delim("GeneOntology_Gruppo1_DEGs.txt",head=T,stringsAsFactors=F,sep="\t",row.names=NULL)
Gruppo1_GO$ID=gsub("\\~.*","",Gruppo1_GO$Term)
Gruppo1_GO$Term2=gsub(".*~","",Gruppo1_GO$Term)
David_Gruppo1=Gruppo1_GO[,c(1,14,15,12,6)]
colnames(David_Gruppo1)=c("category","ID","Term","adj_pval","genes")

Genes_DEG_gruppo1=toptags_xypat_Kallisto_DEG_gruppo1[,c(1,2)]
colnames(Genes_DEG_gruppo1)=c("ID","logFC")
circ <- circle_dat(David_Gruppo1,Genes_DEG_gruppo1)


David_Gruppo1=Gruppo1_GO[,c(1,14,15,12,6)]
colnames(David_Gruppo3)=c("category","ID","Term","adj_pval","genes")
Genes_DEG_gruppo3=toptags_xypat_Kallisto_DEG_gruppo3[,c(1,2)]
colnames(Genes_DEG_gruppo3)=c("ID","logFC")



GO_Gruppo3_Gruppo2_Aggregati




#*************************************************************Vecchia analisi***********************************************************
load("Analisi_Kallisto_NUOVA.RData")



nrow(kallisto_cpm_DEGs)
[1] 609

pdf("FoldChanges_Kallisto_DEGs_15PAzienti.pdf")
breakList=seq(-5, 5, by = 0.1) 
pheatmap( toptags_xypat_Kallisto_cpm_DEG [,38:52], cluster_rows=T,cluster_cols=T,color=colorRampPalette( c("green", "black", "red"), space="rgb")(length(breakList)))
dev.off()


colnames(kallisto_cpm_DEGs)
[1] "genes"  "ALFE1"  "ALFE2"  "BESU1"  "BESU2"  "CALU1"  "CALU2"  "DEIV2"
[9] "DEIV2b" "DEST1"  "DEST2"  "DR41"   "DR42"   "DR51"   "DR52"   "FOCA1"
[17] "FOCA2"  "GAGRA1" "GAGRA2" "LUAN1"  "LUAN2"  "MABI1"  "MABI2"  "MOGE1"
[25] "MOGE2"  "PIAG1"  "PIAG2"  "PRELU1" "PRELU2" "logFC"  "logCPM" "LR"
[33] "PValue" "FDR"    "Qvalue"

kallisto_cpm_DEGs$ALFE_FC=kallisto_cpm_DEGs$ALFE2-kallisto_cpm_DEGs$ALFE1
kallisto_cpm_DEGs$BESU_FC=kallisto_cpm_DEGs$BESU2-kallisto_cpm_DEGs$BESU1
kallisto_cpm_DEGs$CALU_FC=kallisto_cpm_DEGs$CALU2-kallisto_cpm_DEGs$CALU1
kallisto_cpm_DEGs$DEIV_FC=kallisto_cpm_DEGs$DEIV2b-kallisto_cpm_DEGs$DEIV2
kallisto_cpm_DEGs$DEST_FC=kallisto_cpm_DEGs$DEST2-kallisto_cpm_DEGs$DEST1
kallisto_cpm_DEGs$DR4_FC=kallisto_cpm_DEGs$DR42-kallisto_cpm_DEGs$DR41
kallisto_cpm_DEGs$DR5_FC=kallisto_cpm_DEGs$DR52-kallisto_cpm_DEGs$DR51
kallisto_cpm_DEGs$FOCA_FC=kallisto_cpm_DEGs$FOCA2-kallisto_cpm_DEGs$FOCA1
kallisto_cpm_DEGs$GAGRA_FC=kallisto_cpm_DEGs$GAGRA2-kallisto_cpm_DEGs$GAGRA1
kallisto_cpm_DEGs$LUAN_FC=kallisto_cpm_DEGs$LUAN2-kallisto_cpm_DEGs$LUAN1
kallisto_cpm_DEGs$MABI_FC=kallisto_cpm_DEGs$MABI2-kallisto_cpm_DEGs$MABI1
kallisto_cpm_DEGs$MOGE_FC=kallisto_cpm_DEGs$MOGE2-kallisto_cpm_DEGs$MOGE1
kallisto_cpm_DEGs$PIAG_FC=kallisto_cpm_DEGs$PIAG2-kallisto_cpm_DEGs$PIAG1
kallisto_cpm_DEGs$PRELU_FC=kallisto_cpm_DEGs$PRELU2-kallisto_cpm_DEGs$PRELU1

colnames(kallisto_cpm_DEGs)
[1] "genes"    "ALFE1"    "ALFE2"    "BESU1"    "BESU2"    "CALU1"
[7] "CALU2"    "DEIV2"    "DEIV2b"   "DEST1"    "DEST2"    "DR41"
[13] "DR42"     "DR51"     "DR52"     "FOCA1"    "FOCA2"    "GAGRA1"
[19] "GAGRA2"   "LUAN1"    "LUAN2"    "MABI1"    "MABI2"    "MOGE1"
[25] "MOGE2"    "PIAG1"    "PIAG2"    "PRELU1"   "PRELU2"   "logFC"
[31] "logCPM"   "LR"       "PValue"   "FDR"      "Qvalue"   "ALFE_FC"
[37] "BESU_FC"  "CALU_FC"  "DEIV_FC"  "DEST_FC"  "DR4_FC"   "DR5_FC"
[43] "FOGA_FC"  "GAGRA_FC" "LUAN_FC"  "MABI_FC"  "MOGE_FC"  "PIAG_FC"
[49] "PRELU_FC"


ALFE_RNK=kallisto_cpm_DEGs[(kallisto_cpm_DEGs$ALFE2>=0.58&kallisto_cpm_DEGs$ALFE1>=0.58) & (kallisto_cpm_DEGs$ALFE_FC<=-0.58 | kallisto_cpm_DEGs$ALFE_FC>=0.58),c(1,2,3,36)]
BESU_RNK=kallisto_cpm_DEGs[(kallisto_cpm_DEGs$BESU2>=0.58&kallisto_cpm_DEGs$BESU1>=0.58) & (kallisto_cpm_DEGs$BESU_FC<=-0.58 | kallisto_cpm_DEGs$BESU_FC>=0.58),c(1,4,5,37)]
CALU_RNK=kallisto_cpm_DEGs[(kallisto_cpm_DEGs$CALU2>=0.58&kallisto_cpm_DEGs$CALU1>=0.58) & (kallisto_cpm_DEGs$CALU_FC<=-0.58 | kallisto_cpm_DEGs$CALU_FC>=0.58),c(1,6,7,38)]
DEIV_RNK=kallisto_cpm_DEGs[(kallisto_cpm_DEGs$DEIV2b>=0.58&kallisto_cpm_DEGs$DEIV2>=0.58) & (kallisto_cpm_DEGs$DEIV_FC<=-0.58 | kallisto_cpm_DEGs$DEIV_FC>=0.58),c(1,8,9,39)]
DEST_RNK=kallisto_cpm_DEGs[(kallisto_cpm_DEGs$DEST2>=0.58&kallisto_cpm_DEGs$DEST1>=0.58) & (kallisto_cpm_DEGs$DEST_FC<=-0.58 | kallisto_cpm_DEGs$DEST_FC>=0.58),c(1,10,11,40)]
DR4_RNK=kallisto_cpm_DEGs[(kallisto_cpm_DEGs$DR42>=0.58&kallisto_cpm_DEGs$DR41>=0.58) & (kallisto_cpm_DEGs$DR4_FC<=-0.58 | kallisto_cpm_DEGs$DR4_FC>=0.58),c(1,12,13,41)]
DR5_RNK=kallisto_cpm_DEGs[(kallisto_cpm_DEGs$DR52>=0.58&kallisto_cpm_DEGs$DR51>=0.58) & (kallisto_cpm_DEGs$DR5_FC<=-0.58 | kallisto_cpm_DEGs$DR5_FC>=0.58),c(1,14,15,42)]
FOCA_RNK=kallisto_cpm_DEGs[(kallisto_cpm_DEGs$FOCA2>=0.58&kallisto_cpm_DEGs$FOCA1>=0.58) & (kallisto_cpm_DEGs$FOCA_FC<=-0.58 | kallisto_cpm_DEGs$FOCA_FC>=0.58),c(1,16,17,43)]
GAGRA_RNK=kallisto_cpm_DEGs[(kallisto_cpm_DEGs$GAGRA2>=0.58&kallisto_cpm_DEGs$GAGRA1>=0.58) & (kallisto_cpm_DEGs$GAGRA_FC<=-0.58 | kallisto_cpm_DEGs$GAGRA_FC>=0.58),c(1,18,19,44)]
LUAN_RNK=kallisto_cpm_DEGs[(kallisto_cpm_DEGs$LUAN2>=0.58&kallisto_cpm_DEGs$LUAN1>=0.58) & (kallisto_cpm_DEGs$LUAN_FC<=-0.58 | kallisto_cpm_DEGs$LUAN_FC>=0.58),c(1,20,21,45)]
MABI_RNK=kallisto_cpm_DEGs[(kallisto_cpm_DEGs$MABI2>=0.58&kallisto_cpm_DEGs$MABI1>=0.58) & (kallisto_cpm_DEGs$MABI_FC<=-0.58 | kallisto_cpm_DEGs$MABI_FC>=0.58),c(1,22,23,46)]
MOGE_RNK=kallisto_cpm_DEGs[(kallisto_cpm_DEGs$MOGE2>=0.58&kallisto_cpm_DEGs$MOGE1>=0.58) & (kallisto_cpm_DEGs$MOGE_FC<=-0.58 | kallisto_cpm_DEGs$MOGE_FC>=0.58),c(1,24,25,47)]
PIAG_RNK=kallisto_cpm_DEGs[(kallisto_cpm_DEGs$PIAG2>=0.58&kallisto_cpm_DEGs$PIAG1>=0.58) & (kallisto_cpm_DEGs$PIAG_FC<=-0.58 | kallisto_cpm_DEGs$PIAG_FC>=0.58),c(1,26,27,48)]
PRELU_RNK=kallisto_cpm_DEGs[(kallisto_cpm_DEGs$PRELU2>=0.58&kallisto_cpm_DEGs$PRELU1>=0.58) & (kallisto_cpm_DEGs$PRELU_FC<=-0.58 | kallisto_cpm_DEGs$PRELU_FC>=0.58),c(1,28,29,49)]

head(ALFE_RNK[-order(ALFE_RNK$ALFE_FC),])

mds=plotMDS(countdata_filtDGE_keep,dim.plot=c(1,2))
mds2=data.frame(x=mds$x,y=mds$y)
mds2=cbind(mds2,as.data.frame(condition))
samples=c("UPN_8","UPN_8","UPN_9","UPN_9","UPN_1","UPN_1","UPN_12","UPN_12","UPN_2","UPN_2","UPN_14","UPN_14","UPN_15","UPN_15","UPN_3","UPN_3","UPN_4","UPN_4",   "UPN_13","UPN_13","UPN_5","UPN_5","UPN_6","UPN_6","UPN_11","UPN_11","UPN_7","UPN_7")

mds2=cbind(mds2,as.data.frame(samples))

pdf("PCA_Pazienti_NoDR1_Codici.pdf")
ggplot(mds2, aes(x, y,size=4)) + geom_point(aes(colour = condition)) + geom_text(label=mds2$samples,vjust=0,nudge_y=0.05,size=3)+ scale_y_reverse()+scale_x_continuous(limits = c(-3.5, 3.5))+ guides(colour = guide_legend(override.aes = list(size=4)))
dev.off()

library(biomaRt)
mart=useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="hsapiens_gene_ensembl")
TcellStim_Symbols <- getBM(filters= "entrezgene", attributes= c("hgnc_symbol","entrezgene"),values=TcellStim,mart=mart)
write.table(as.data.frame(G_list_gruppo3$entrezgene),"/home/fsantaniello/Background_Geni_Gruppo3.txt",sep="\t",col.names=F,row.names=F,quote=F)


TcellStim=c(3113, 3127, 3115, 3117, 3119, 100133941, 9020, 3123, 3122)

3113, 3127, 3115, 3117, 3119, 100133941, 9020, 3123, 3122

GO:003129





#Rifaccio analisi utilizzando RsubRead e un altro GTF alternativo

newcounts2=featureCounts(c('/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/Sample_ALFE1_Allfiles_sorted.bam',
'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/Sample_ALFE2_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/BESU1_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/BESU2_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/CALU1_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/CALU2_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/Sample_DEIV2_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/Sample_DEIV2b_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/Sample_DR1_1_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/Sample_DR1_2_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/DR41_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/DR42_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/DR51_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/DR52_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/Sample_DEST1_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/Sample_DEST2_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/FOCA1_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/FOCA2_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/GAGRA1_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/GAGRA2_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/Sample_LUAN1_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/Sample_LUAN2_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/MABI1_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/MABI2_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/MOGE1_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/MOGE2_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/PIAG1_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/PIAG2_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/PRELU1_Allfiles_sorted.bam'
,'/lustre2/scratch/bioinfotree/common/bioinfotree/prj/Vago_161_RNASeqVariantCalling_Hisat2/dataset/20170503/PRELU2_Allfiles_sorted.bam'),
annot.ext='/lustre1/genomes/hg38/annotation/gencode.v26.chr_patch_hapl_scaff.annotation.gtf',
isGTFAnnotationFile = TRUE,GTF.attrType='gene_name',chrAliases='Aliases_GTF.txt'),useMetaFeatures=T)


newcounts=readRDS('Conte_GabriLucia_NewGTF.rds')
colnames(newcounts2$counts)=c("ALFE1","ALFE2","BESU1","BESU2","CALU1","CALU2","DEIV2","DEIV2b","DR1_1" ,"DR1_2" ,"DR41",
"DR42","DR51","DR52","DEST1","DEST2","FOCA1","FOCA2","GAGRA1","GAGRA2","LUAN1",
"LUAN2","MABI1","MABI2","MOGE1","MOGE2","PIAG1","PIAG2","PRELU1","PRELU2")
condition_new=rep(c("Diagnosis","Relapse"),15)
countdata_new=newcounts2$counts
countdata_filt_new <- countdata_new[rowSums(countdata_new[, 1:30]) > 0, ]
countdata_filtDGE_new=DGEList(countdata_filt_new[,1:30],genes=rownames(countdata_filt_new),group=condition_new)
keep_new <- rowSums(cpm(countdata_filtDGE_new)>1) >= 8 #PRENDO SOLO I GENI CON CPM>1 IN ALMENO 8 CAMPIONI
countdata_filtDGE_keep_new<- countdata_filtDGE_new[keep_new, , keep.lib.sizes=FALSE]
countdata_filtDGE_keep_new=calcNormFactors(countdata_filtDGE_keep_new,method="TMM")
mds_new=plotMDS(countdata_filtDGE_keep_new,dim.plot=c(1,2))
mds2_new=data.frame(x=mds_new$x,y=mds_new$y)
samples=rownames(mds2_new)
mds2_new=cbind(mds2_new,as.data.frame(samples),as.data.frame(condition_new))
patient=rep(1:15,each=2)
mds2_new=cbind(mds2_new,as.data.frame(patient))
colnames(mds2_new)=c("x","y","samples","condition","patient")
design_new = model.matrix (~condition + x + y + as.factor(patient),data=mds2_new)
countdata_filtDGE_keep_new <- estimateDisp(countdata_filtDGE_keep_new,design_new)
fit_new <- glmFit(countdata_filtDGE_keep_new, design_new)
lrt_new <- glmLRT(fit_new,coef=2)
#Correggo per qvalue e seleziono solo i significativi per qvalue<=0.25
library(qvalue)
toptags_xypat_Kallisto_new=as.data.frame(topTags(lrt_new,sort.by="PValue",n="all"))
q_new = qvalue(toptags_xypat_Kallisto_new$PValue)
toptags_xypat_Kallisto_new$Qvalue = q_new$qvalues
toptags_xypat_Kallisto_DEG_new=toptags_xypat_Kallisto_new[toptags_xypat_Kallisto_new$Qvalue<=0.25,]
nrow(toptags_xypat_Kallisto_DEG_new)

All_Genes_NewLuciaGabri <- getBM(filters= "entrezgene", attributes= c("hgnc_symbol","entrezgene"),values=toptags_xypat_Kallisto,mart=mart)
write.table(as.data.frame(All_Genes_NewLuciaGabri$entrezgene),"/home/fsantaniello/Background_Geni_NewGabriLucia.txt",sep="\t",col.names=F,row.names=F,quote=F)

Deg_Genes_NewLuciaGabri <- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","entrezgene"),values=toptags_xypat_Kallisto_DEG$genes,mart=mart)
write.table(as.data.frame(Deg_Genes_NewLuciaGabri$entrezgene),"/home/fsantaniello/Significantly_DEGs_Geni_NewGabriLucia.txt",sep="\t",col.names=F,row.names=F,quote=F)


lsEOG<-list()

for (j in 2006:2012){
  z <- j
  sEOG <- paste("EOG", z, sep="")
  print(sEOG)
  dEOG <- get(paste("EOG", z, sep=""))
  print(dEOG)}
  lsEOG[[sEOG]] <-dEOG
}
   


cat /lustre2/scratch/gbucci/Vago/161/Exome/2017/BAM/CNS/PlayGround/*.segment-gainloss.bed | bedSort stdin stdout | awk -v OFS='\t' '{print $4,$4}' -  | sort | uniq > /lustre2/scratch/gbucci/Vago/161/Exome/2017/BAM/CNS/PlayGround/AllGenes_CNVs.txt

library(plyr)
All_Genes=read.table('/lustre2/scratch/gbucci/Vago/161/Exome/2017/BAM/CNS/PlayGround/AllGenes_CNVs.txt',head=T,stringsAsFactors=F,col.names=c('gene','gene_2'),colClasses(NA,'NULL'))
filenames=list.files(path='/lustre2/scratch/gbucci/Vago/161/Exome/2017/BAM/CNS/PlayGround/', ,pattern='segment-gainloss.bed',full.names=TRUE)
datalist = lapply(filenames, function(x){read.table(file=x,header=F,colClasses=c("NULL","NULL","NULL",NA,NA),col.names=c('chr','start','end','gene','logRatio'),stringsAsFactors=F)})
dataList=list(All_Genes,datalist)
names_file=gsub('/lustre2/scratch/gbucci/Vago/161/Exome/2017/BAM/CNS/PlayGround//','',filenames)
CNVs_Genes=join_all(prova,by='gene',match='first')
colnames(CNVs_Genes)=c('genes',names_file)

write.table(CNVs_Genes,'/lustre2/scratch/gbucci/Vago/161/Exome/2017/BAM/CNS/PlayGround/All_Genes_AllSamples_CNVsLogRatios.txt',sep='\t',col.names=T,row.names=F,quote=F)



geniImmuni_DaEleonora_X_Lucia=c('ABCF3',
'ABI1',
'ABL1',
'ABL2',
'ACIN1',
'ACTB',
'ACTG1',
'ACTN1',
'ACTR1A',
'ACTR1B',
'ACTR2',
'ACTR3',
'ADA',
'ADAM10',
'ADAM15',
'ADAM17',
'ADAM8',
'ADAR',
'ADD1',
'ADSS',
'AGBL5',
'AGPAT5',
'AHCY',
'AIF1',
'AIMP1',
'AKAP17A',
'AKAP8',
'AKIRIN2',
'AKT1',
'ALCAM',
'AMPD3',
'ANKHD1',
'ANKRD17',
'ANO6',
'ANXA1',
'ANXA2',
'AP1B1',
'AP1G1',
'AP1M1',
'AP1S1',
'AP1S2',
'AP2A1',
'AP2A2',
'AP2B1',
'AP2M1',
'AP2S1',
'AP3B1',
'AP3D1',
'APBB1IP',
'APOBEC3C',
'APOBEC3G',
'ARF1',
'ARHGEF2',
'ARHGEF7',
'ARID4A',
'ARIH2',
'ARMC6',
'ARPC1A',
'ARPC1B',
'ARPC2',
'ARPC3',
'ARPC4',
'ARPC5',
'ASH2L',
'ATG12',
'ATG5',
'ATG7',
'ATM',
'ATP11C',
'ATP1B1',
'ATP1B3',
'ATP6V0A2',
'ATPIF1',
'AZI2',
'B2M',
'B4GALT1',
'BAG6',
'BAIAP2',
'BAK1',
'BAX',
'BCAP31',
'BCL10',
'BCL3',
'BECN1',
'BIRC2',
'BIRC3',
'BMI1',
'BNIP3L',
'BRAF',
'BRK1',
'BSG',
'BST2',
'BTK',
'BTN3A2',
'C19orf66',
'C1QBP',
'C1RL',
'CALCOCO2',
'CALR',
'CAMK2G',
'CANX',
'CAPZA1',
'CAPZA2',
'CASP10',
'CASP3',
'CASP4',
'CASP8',
'CASP9',
'CBFA2T3',
'CBFB',
'CCL5',
'CCND3',
'CD109',
'CD151',
'CD164',
'CD244',
'CD300A',
'CD300LF',
'CD44',
'CD46',
'CD47',
'CD55',
'CD58',
'CD74',
'CD79B',
'CD84',
'CDC42',
'CDK13',
'CDK6',
'CDKN1C',
'CEBPA',
'CEBPB',
'CEBPG',
'CFP',
'CHD2',
'CHID1',
'CHUK',
'CIAPIN1',
'CIB1',
'CITED2',
'CKLF',
'CLTA',
'CLTC',
'CMTM7',
'CNPY3',
'COL24A1',
'COL4A3BP',
'CORO1A',
'CRCP',
'CREBBP',
'CRIP1',
'CRK',
'CRKL',
'CRTC3',
'CSF3R',
'CSK',
'CTC1',
'CTNNB1',
'CTNNBL1',
'CTSB',
'CTSC',
'CTSD',
'CTSS',
'CUL1',
'CUL4A',
'CXCL2',
'CXCR4',
'CYBA',
'CYFIP1',
'CYFIP2',
'CYLD',
'DAPK1',
'DAPK3',
'DBNL',
'DCLRE1C',
'DCTN1',
'DCTN2',
'DCTN3',
'DCTN4',
'DCTN5',
'DCTN6',
'DDIT4',
'DDOST',
'DDX3X',
'DDX41',
'DDX58',
'DENND1B',
'DHTKD1',
'DLG1',
'DNAJA3',
'DNAJC3',
'DNASE2',
'DNM2',
'DOCK2',
'DOCK8',
'DOK2',
'DPP8',
'DYNC1H1',
'DYNC1I2',
'DYNC1LI2',
'DYNC2LI1',
'DYNLL1',
'DYNLL2',
'EBP',
'ECSIT',
'EEF2',
'EGR1',
'EIF2AK1',
'EIF2AK2',
'EIF2AK4',
'ELF4',
'ELMO1',
'ELMO2',
'EP300',
'EPRS',
'ERAP1',
'ERAP2',
'ERCC1',
'ETV6',
'EXOSC3',
'EXOSC4',
'EXOSC6',
'EXOSC9',
'F11R',
'FADD',
'FAM111A',
'FARP2',
'FAS',
'FAU',
'FBXO9',
'FBXW11',
'FES',
'FKBP1A',
'FLCN',
'FLI1',
'FLT3',
'FLVCR1',
'FNIP1',
'FOXP1',
'FTH1',
'FYN',
'G6PD',
'GAB2',
'GAB3',
'GALNT2',
'GAPDH',
'GATA2',
'GBF1',
'GBP5',
'GCH1',
'GCNT1',
'GFI1',
'GLG1',
'GLO1',
'GLRX5',
'GNL1',
'GOLPH3',
'GPI',
'GPR183',
'GRB2',
'GSDMD',
'GTPBP1',
'HAVCR2',
'HCLS1',
'HDAC5',
'HHEX',
'HIF1A',
'HIPK1',
'HIST1H2BK',
'HIST2H2BE',
'HLA-A',
'HLA-B',
'HLA-C',
'HLA-DMA',
'HLA-DPA1',
'HLA-DPB1',
'HLA-DRA',
'HLA-DRB1',
'HLA-DRB5',
'HLA-E',
'HLA-F',
'HLA-H',
'HMGB1',
'HMGB2',
'HMGB3',
'HNRNPK',
'HPRT1',
'HRAS',
'HSH2D',
'HSP90AA1',
'HSP90AB1',
'HSP90B1',
'HSPD1',
'ICAM1',
'ICAM2',
'ICAM3',
'ICOSLG',
'ID2',
'IFI16',
'IFI30',
'IFI35',
'IFI6',
'IFIH1',
'IFITM2',
'IFNAR1',
'IFNAR2',
'IFNGR1',
'IFNGR2',
'IGBP1',
'IK',
'IKBKAP',
'IKBKB',
'IKBKG',
'IKZF1',
'IL10RB',
'IL18',
'IL1RAP',
'IL27RA',
'IL2RG',
'IL4R',
'IL6R',
'ILF2',
'ILF3',
'IMPDH1',
'IMPDH2',
'INPP5D',
'INPPL1',
'IP6K2',
'IPO7',
'IRAK1',
'IRAK3',
'IRAK4',
'IRF1',
'IRF2',
'IRF3',
'IRF5',
'IRF7',
'IRF9',
'ISG15',
'ITCH',
'ITFG2',
'ITGA4',
'ITGA5',
'ITGA6',
'ITGAL',
'ITGB1',
'ITGB2',
'ITPKB',
'JAGN1',
'JAK1',
'JAK2',
'JAK3',
'JARID2',
'JMJD6',
'JUN',
'JUNB',
'KAT6A',
'KAT8',
'KCNAB2',
'KCNN4',
'KDM1A',
'KDM6B',
'KIF13B',
'KIF22',
'KIF2A',
'KIF3B',
'KIFAP3',
'KIT',
'KLC1',
'KLF2',
'KLF6',
'KLHL6',
'KRAS',
'LAIR1',
'LAT2',
'LCP1',
'LCP2',
'LFNG',
'LGALS1',
'LGALS9',
'LIG1',
'LIG3',
'LIMK1',
'LNPEP',
'LRMP',
'LRRC8A',
'LRRK1',
'LSM14A',
'LST1',
'LTB',
'LTB4R',
'LTBR',
'LY75',
'LYL1',
'LYN',
'LYST',
'MAD1L1',
'MAEA',
'MALT1',
'MAP2K1',
'MAP2K2',
'MAP3K1',
'MAP3K14',
'MAP3K5',
'MAP3K7',
'MAP3K8',
'MAP4K2',
'MAPK1',
'MAPK14',
'MAPK3',
'MAPKAP1',
'MAPKAPK2',
'MAPKAPK3',
'MARCH8',
'MAVS',
'MB21D1',
'MBP',
'MCM3AP',
'MED1',
'MEF2C',
'MFNG',
'MIF',
'MILR1',
'MKNK2',
'MLH1',
'MLST8',
'MR1',
'MSH2',
'MSH3',
'MSH6',
'MSN',
'MTHFD1',
'MTOR',
'MX1',
'MYD88',
'MYH9',
'MYO1G',
'NAIP',
'NARFL',
'NBEAL2',
'NBN',
'NCF2',
'NCF4',
'NCK1',
'NCK2',
'NCKAP1L',
'NCOA6',
'NCSTN',
'NDRG1',
'NEDD4',
'NFE2L1',
'NFIL3',
'NFKB1',
'NFKBIA',
'NHEJ1',
'NKAP',
'NLRC3',
'NLRC5',
'NLRP3',
'NOTCH1',
'NOTCH2',
'NRAS',
'NUB1',
'NUDCD1',
'NUP85',
'ORAI1',
'OSBPL1A',
'OTUB1',
'PAFAH1B1',
'PAK1',
'PAK2',
'PAPD4',
'PARP1',
'PATZ1',
'PCBP2',
'PCID2',
'PDE4D',
'PDIA3',
'PDPK1',
'PECAM1',
'PIBF1',
'PIK3AP1',
'PIK3C3',
'PIK3CA',
'PIK3CB',
'PIK3CD',
'PIK3R1',
'PIK3R2',
'PIK3R4',
'PIP4K2A',
'PIP5K1C',
'PKN1',
'PKNOX1',
'PLA2G6',
'PLCG2',
'PLEK',
'PLSCR1',
'PMAIP1',
'PML',
'PNMA1',
'PNP',
'POLB',
'POLL',
'POLM',
'POLR3A',
'POLR3C',
'POLR3D',
'POLR3E',
'POLR3F',
'POLR3K',
'PPIA',
'PPIL2',
'PPP1R14B',
'PPP2R3C',
'PPP3CA',
'PPP3CB',
'PRDX1',
'PRDX3',
'PRELID1',
'PREX1',
'PRKACA',
'PRKACB',
'PRKCD',
'PRKCQ',
'PRKD2',
'PRKDC',
'PRKRA',
'PRKX',
'PRR5',
'PRRC2C',
'PSEN1',
'PSMA1',
'PSMA2',
'PSMA3',
'PSMA4',
'PSMA5',
'PSMA6',
'PSMA7',
'PSMB1',
'PSMB10',
'PSMB2',
'PSMB3',
'PSMB4',
'PSMB5',
'PSMB6',
'PSMB7',
'PSMB8',
'PSMB9',
'PSMC1',
'PSMC2',
'PSMC3',
'PSMC4',
'PSMC5',
'PSMC6',
'PSMD1',
'PSMD10',
'PSMD11',
'PSMD12',
'PSMD13',
'PSMD14',
'PSMD2',
'PSMD3',
'PSMD4',
'PSMD6',
'PSMD7',
'PSMD8',
'PSMD9',
'PSME1',
'PSME2',
'PSME3',
'PSME4',
'PSMF1',
'PTEN',
'PTGER4',
'PTK2B',
'PTPN11',
'PTPN2',
'PTPN22',
'PTPN6',
'PTPRC',
'PYCARD',
'RAB10',
'RAB32',
'RAB35',
'RAB4A',
'RAB5B',
'RAB6A',
'RAB7A',
'RAB8B',
'RAC1',
'RACGAP1',
'RAF1',
'RASGRP4',
'RB1',
'RBCK1',
'RBM15',
'RBPJ',
'RC3H1',
'RC3H2',
'REL',
'RELA',
'RELT',
'REST',
'RFTN1',
'RFX1',
'RGS1',
'RHOH',
'RICTOR',
'RILP',
'RIPK1',
'RIPK2',
'RNF125',
'RNF168',
'RNF19B',
'RNF31',
'RNF8',
'ROCK1',
'RPL13A',
'RPL22',
'RPL39',
'RPS14',
'RPS17',
'RPS19',
'RPS24',
'RPS27A',
'RPS6',
'RPS6KA3',
'RUNX1',
'RUNX3',
'S100A13',
'S1PR4',
'SAR1B',
'SART3',
'SATB1',
'SBDS',
'SBNO2',
'SEC13',
'SEC14L1',
'SEC23A',
'SEC24A',
'SEC24B',
'SEC24C',
'SEC24D',
'SEC31A',
'SEC61A1',
'SELPLG',
'SEMA4A',
'SEMA4D',
'SERINC3',
'SGPL1',
'SH2B3',
'SHC1',
'SHMT2',
'SIN3A',
'SIRPA',
'SIRT1',
'SIRT2',
'SIVA1',
'SKAP2',
'SKP1',
'SLA2',
'SLC16A3',
'SLC25A38',
'SLC39A3',
'SLC3A2',
'SLC7A5',
'SLC7A6',
'SLC7A6OS',
'SLFN11',
'SMAD3',
'SMAD5',
'SNAP23',
'SNCA',
'SNRK',
'SNX27',
'SOD1',
'SOS1',
'SOX4',
'SP100',
'SP2',
'SP3',
'SPG21',
'SPI1',
'SPN',
'SPNS2',
'SPPL3',
'SQSTM1',
'SRF',
'SRPK1',
'SRPK2',
'SSBP3',
'STAT1',
'STAT2',
'STAT5B',
'STAT6',
'STK11',
'STK3',
'STK4',
'STOML2',
'STXBP2',
'STXBP3',
'SWAP70',
'SYK',
'SYNCRIP',
'TAB1',
'TAB2',
'TAB3',
'TANK',
'TAP1',
'TAP2',
'TAPBP',
'TAPBPL',
'TAZ',
'TBK1',
'TBKBP1',
'TCF12',
'TCF3',
'TEC',
'TET2',
'TFE3',
'TFRC',
'TGFB1',
'TGFBR1',
'TGFBR2',
'THOC5',
'TIPARP',
'TMEM173',
'TMOD3',
'TNFAIP1',
'TNFAIP3',
'TNFAIP8L2',
'TNFRSF10B',
'TNFRSF10D',
'TNFRSF14',
'TNFRSF1A',
'TNFRSF1B',
'TNFSF13',
'TNIP1',
'TNIP2',
'TNK2',
'TOLLIP',
'TRAF3',
'TRAF3IP2',
'TRAF6',
'TRIM11',
'TRIM13',
'TRIM14',
'TRIM21',
'TRIM22',
'TRIM25',
'TRIM26',
'TRIM27',
'TRIM28',
'TRIM38',
'TRIM4',
'TRIM5',
'TRIM8',
'TROVE2',
'TSC1',
'TTC7A',
'TUBB',
'TUBB4B',
'TUSC2',
'TWSG1',
'TXLNA',
'TYK2',
'TYROBP',
'UBA52',
'UBB',
'UBC',
'UBE2D2',
'UBE2D3',
'UBE2N',
'UBE2V1',
'UNC13D',
'UNC93B1',
'UNG',
'VAMP2',
'VAMP7',
'VAMP8',
'VAV1',
'VAV3',
'VEGFA',
'VPS33A',
'VPS33B',
'WAS',
'WASF2',
'WIPF1',
'WIPF2',
'XBP1',
'XRCC5',
'YES1',
'YTHDF2',
'ZBTB1',
'ZBTB24',
'ZC3H12A',
'ZC3H8',
'ZC3HAV1',
'ZFP36L2',
'ZNF160',
'ZNF175',
'ZNF3',
'ZNF385A')





res2allgenes <- rcorr(as.matrix(t(Toptags_Immune[,38:52])))
cor.matrix=res2allgenes$r
colnames(cor.matrix)=Toptags_Immune$genes
rownames(cor.matrix)=Toptags_Immune$genes

cor.matrix_rows=apply(cor.matrix, 1, percent_rank)
cor.matrix_cols=apply(cor.matrix, 2, percent_rank)
cor.matrix_sqrt=sqrt(cor.matrix_rows*cor.matrix_cols)
library(reshape2)
my_cor_df <- melt(cor.matrix_sqrt)
my_cor_df= my_cor_df %>% filter(Var1 != Var2)
my_adj_list <- my_cor_df %>% filter(value > 0.995)
names(my_adj_list) <- c('from', 'to', 'weight') 
g=graph.data.frame(my_adj_list, directed=F)
communities=cluster_louvain(g, weights = NULL)
saveRDS(g,'Grafo_Geni_Immuni_15Pazienti.RDS')




kal=Toptags_Immune[Toptags_Immune$Community==1,];x=table(kal$QvaluePositive)[[2]];m=table(kal$QvaluePositive)[[2]]+table(Toptags_Immune$QvaluePositive)[[2]];n=table(kal$QvaluePositive)[[1]]+table(Toptags_Immune$QvaluePositive)[[1]];k=table(kal$QvaluePositive)[[2]]+table(kal$QvaluePositive)[[1]];print(paste('Enrichment of Significantly Deregulated Genes in Community ',1,' is: ',sep=''));print(phyper(x,m,n,k,lower.tail = F))
kal=Toptags_Immune[Toptags_Immune$Community==2,];x=table(kal$QvaluePositive)[[2]];m=table(kal$QvaluePositive)[[2]]+table(Toptags_Immune$QvaluePositive)[[2]];n=table(kal$QvaluePositive)[[1]]+table(Toptags_Immune$QvaluePositive)[[1]];k=table(kal$QvaluePositive)[[2]]+table(kal$QvaluePositive)[[1]];print(paste('Enrichment of Significantly Deregulated Genes in Community ',3,' is: ',sep=''));print(phyper(x,m,n,k,lower.tail = F))
kal=Toptags_Immune[Toptags_Immune$Community==3,];x=table(kal$QvaluePositive)[[2]];m=table(kal$QvaluePositive)[[2]]+table(Toptags_Immune$QvaluePositive)[[2]];n=table(kal$QvaluePositive)[[1]]+table(Toptags_Immune$QvaluePositive)[[1]];k=table(kal$QvaluePositive)[[2]]+table(kal$QvaluePositive)[[1]];print(paste('Enrichment of Significantly Deregulated Genes in Community ',3,' is: ',sep=''));print(phyper(x,m,n,k,lower.tail = F))
kal=Toptags_Immune[Toptags_Immune$Community==4,];x=table(kal$QvaluePositive)[[2]];m=table(kal$QvaluePositive)[[2]]+table(Toptags_Immune$QvaluePositive)[[2]];n=table(kal$QvaluePositive)[[1]]+table(Toptags_Immune$QvaluePositive)[[1]];k=table(kal$QvaluePositive)[[2]]+table(kal$QvaluePositive)[[1]];print(paste('Enrichment of Significantly Deregulated Genes in Community ',4,' is: ',sep=''));print(phyper(x,m,n,k,lower.tail = F))
kal=Toptags_Immune[Toptags_Immune$Community==5,];x=table(kal$QvaluePositive)[[2]];m=table(kal$QvaluePositive)[[2]]+table(Toptags_Immune$QvaluePositive)[[2]];n=table(kal$QvaluePositive)[[1]]+table(Toptags_Immune$QvaluePositive)[[1]];k=table(kal$QvaluePositive)[[2]]+table(kal$QvaluePositive)[[1]];print(paste('Enrichment of Significantly Deregulated Genes in Community ',5,' is: ',sep=''));print(phyper(x,m,n,k,lower.tail = F))
kal=Toptags_Immune[Toptags_Immune$Community==6,];x=table(kal$QvaluePositive)[[2]];m=table(kal$QvaluePositive)[[2]]+table(Toptags_Immune$QvaluePositive)[[2]];n=table(kal$QvaluePositive)[[1]]+table(Toptags_Immune$QvaluePositive)[[1]];k=table(kal$QvaluePositive)[[2]]+table(kal$QvaluePositive)[[1]];print(paste('Enrichment of Significantly Deregulated Genes in Community ',6,' is: ',sep=''));print(phyper(x,m,n,k,lower.tail = F))
kal=Toptags_Immune[Toptags_Immune$Community==7,];x=table(kal$QvaluePositive)[[2]];m=table(kal$QvaluePositive)[[2]]+table(Toptags_Immune$QvaluePositive)[[2]];n=table(kal$QvaluePositive)[[1]]+table(Toptags_Immune$QvaluePositive)[[1]];k=table(kal$QvaluePositive)[[2]]+table(kal$QvaluePositive)[[1]];print(paste('Enrichment of Significantly Deregulated Genes in Community ',7,' is: ',sep=''));print(phyper(x,m,n,k,lower.tail = F))
kal=Toptags_Immune[Toptags_Immune$Community==8,];x=table(kal$QvaluePositive)[[2]];m=table(kal$QvaluePositive)[[2]]+table(Toptags_Immune$QvaluePositive)[[2]];n=table(kal$QvaluePositive)[[1]]+table(Toptags_Immune$QvaluePositive)[[1]];k=table(kal$QvaluePositive)[[2]]+table(kal$QvaluePositive)[[1]];print(paste('Enrichment of Significantly Deregulated Genes in Community ',8,' is: ',sep=''));print(phyper(x,m,n,k,lower.tail = F))
kal=Toptags_Immune[Toptags_Immune$Community==9,];x=table(kal$QvaluePositive)[[2]];m=table(kal$QvaluePositive)[[2]]+table(Toptags_Immune$QvaluePositive)[[2]];n=table(kal$QvaluePositive)[[1]]+table(Toptags_Immune$QvaluePositive)[[1]];k=table(kal$QvaluePositive)[[2]]+table(kal$QvaluePositive)[[1]];print(paste('Enrichment of Significantly Deregulated Genes in Community ',9,' is: ',sep=''));print(phyper(x,m,n,k,lower.tail = F))
kal=Toptags_Immune[Toptags_Immune$Community==10,];x=table(kal$QvaluePositive)[[2]];m=table(kal$QvaluePositive)[[2]]+table(Toptags_Immune$QvaluePositive)[[2]];n=table(kal$QvaluePositive)[[1]]+table(Toptags_Immune$QvaluePositive)[[1]];k=table(kal$QvaluePositive)[[2]]+table(kal$QvaluePositive)[[1]];print(paste('Enrichment of Significantly Deregulated Genes in Community ',10,' is: ',sep=''));print(phyper(x,m,n,k,lower.tail = F))
kal=Toptags_Immune[Toptags_Immune$Community==11,];x=table(kal$QvaluePositive)[[2]];m=table(kal$QvaluePositive)[[2]]+table(Toptags_Immune$QvaluePositive)[[2]];n=table(kal$QvaluePositive)[[1]]+table(Toptags_Immune$QvaluePositive)[[1]];k=table(kal$QvaluePositive)[[2]]+table(kal$QvaluePositive)[[1]];print(paste('Enrichment of Significantly Deregulated Genes in Community ',11,' is: ',sep=''));print(phyper(x,m,n,k,lower.tail = F))
kal=Toptags_Immune[Toptags_Immune$Community==12,];x=table(kal$QvaluePositive)[[2]];m=table(kal$QvaluePositive)[[2]]+table(Toptags_Immune$QvaluePositive)[[2]];n=table(kal$QvaluePositive)[[1]]+table(Toptags_Immune$QvaluePositive)[[1]];k=table(kal$QvaluePositive)[[2]]+table(kal$QvaluePositive)[[1]];print(paste('Enrichment of Significantly Deregulated Genes in Community ',12,' is: ',sep=''));print(phyper(x,m,n,k,lower.tail = F))
kal=Toptags_Immune[Toptags_Immune$Community==13,];x=table(kal$QvaluePositive)[[2]];m=table(kal$QvaluePositive)[[2]]+table(Toptags_Immune$QvaluePositive)[[2]];n=table(kal$QvaluePositive)[[1]]+table(Toptags_Immune$QvaluePositive)[[1]];k=table(kal$QvaluePositive)[[2]]+table(kal$QvaluePositive)[[1]];print(paste('Enrichment of Significantly Deregulated Genes in Community ',13,' is: ',sep=''));print(phyper(x,m,n,k,lower.tail = F))
kal=Toptags_Immune[Toptags_Immune$Community==14,];x=table(kal$QvaluePositive)[[2]];m=table(kal$QvaluePositive)[[2]]+table(Toptags_Immune$QvaluePositive)[[2]];n=table(kal$QvaluePositive)[[1]]+table(Toptags_Immune$QvaluePositive)[[1]];k=table(kal$QvaluePositive)[[2]]+table(kal$QvaluePositive)[[1]];print(paste('Enrichment of Significantly Deregulated Genes in Community ',14,' is: ',sep=''));print(phyper(x,m,n,k,lower.tail = F))
kal=Toptags_Immune[Toptags_Immune$Community==15,];x=table(kal$QvaluePositive)[[2]];m=table(kal$QvaluePositive)[[2]]+table(Toptags_Immune$QvaluePositive)[[2]];n=table(kal$QvaluePositive)[[1]]+table(Toptags_Immune$QvaluePositive)[[1]];k=table(kal$QvaluePositive)[[2]]+table(kal$QvaluePositive)[[1]];print(paste('Enrichment of Significantly Deregulated Genes in Community ',15,' is: ',sep=''));print(phyper(x,m,n,k,lower.tail = F))
kal=Toptags_Immune[Toptags_Immune$Community==16,];x=table(kal$QvaluePositive)[[2]];m=table(kal$QvaluePositive)[[2]]+table(Toptags_Immune$QvaluePositive)[[2]];n=table(kal$QvaluePositive)[[1]]+table(Toptags_Immune$QvaluePositive)[[1]];k=table(kal$QvaluePositive)[[2]]+table(kal$QvaluePositive)[[1]];print(paste('Enrichment of Significantly Deregulated Genes in Community ',16,' is: ',sep=''));print(phyper(x,m,n,k,lower.tail = F))



library(RColorBrewer)


qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

color=sample(col_vector)

color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
par(mar=c(5.1,5.1,5.1,5.1))
#ets the bottom, left, top and right margins respectively of the plot region in number of lines of text.
color=sample(color)
plot(lines(loess.smooth(1:length(alfe_strong_filtered$IC50_Binder), y=alfe_strong_filtered$IC50_Binder)),xlim=c(0,150),ylim=c(0,50),col=color[1],xlab='Number of Strong Binders',ylab='IC50 (nM)')
lines(loess.smooth(1:length(alfe_strong_filtered$IC50_Binder), y=alfe_strong_filtered$IC50_Binder,degree = 2),xlim=c(0,150),ylim=c(0,50),col=color[1],xlab='Number of Strong Binders',ylab='IC50 (nM)')
par(new=TRUE)
lines(loess.smooth(1:length(besu_strong_filtered$IC50_Binder), y=besu_strong_filtered$IC50_Binder,degree = 2),xlim=c(0,150),ylim=c(0,50),col=color[2],xlab='Number of Strong Binders',ylab='IC50 (nM)')
par(new=TRUE)
lines(loess.smooth(1:length(calu_strong_filtered$IC50_Binder), y=calu_strong_filtered$IC50_Binder,degree = 2),xlim=c(0,150),ylim=c(0,50),col=color[3],xlab='Number of Strong Binders',ylab='IC50 (nM)')
par(new=TRUE)
lines(loess.smooth(1:length(deiv_strong_filtered$IC50_Binder), y=deiv_strong_filtered$IC50_Binder,degree = 2),xlim=c(0,150),ylim=c(0,50),col=color[4],xlab='Number of Strong Binders',ylab='IC50 (nM)')
par(new=TRUE)
lines(loess.smooth(1:length(dest_strong_filtered$IC50_Binder), y=dest_strong_filtered$IC50_Binder,degree = 2),xlim=c(0,150),ylim=c(0,50),col=color[5],xlab='Number of Strong Binders',ylab='IC50 (nM)')
par(new=TRUE)
lines(loess.smooth(1:length(dr1_strong_filtered$IC50_Binder), y=dr1_strong_filtered$IC50_Binder,degree = 2),xlim=c(0,150),ylim=c(0,50),col=color[6],xlab='Number of Strong Binders',ylab='IC50 (nM)')
par(new=TRUE)
lines(loess.smooth(1:length(dr4_strong_filtered$IC50_Binder), y=dr4_strong_filtered$IC50_Binder,degree = 2),xlim=c(0,150),ylim=c(0,50),col=color[7],xlab='Number of Strong Binders',ylab='IC50 (nM)')
par(new=TRUE)
lines(loess.smooth(1:length(dr5_strong_filtered$IC50_Binder), y=dr5_strong_filtered$IC50_Binder,degree = 2),xlim=c(0,150),ylim=c(0,50),col=color[8],xlab='Number of Strong Binders',ylab='IC50 (nM)')
par(new=TRUE)
lines(loess.smooth(1:length(foca_strong_filtered$IC50_Binder), y=foca_strong_filtered$IC50_Binder,degree = 2),xlim=c(0,150),ylim=c(0,50),col=color[9],xlab='Number of Strong Binders',ylab='IC50 (nM)')
par(new=TRUE)
lines(loess.smooth(1:length(gagra_strong_filtered$IC50_Binder), y=gagra_strong_filtered$IC50_Binder,degree = 2),xlim=c(0,150),ylim=c(0,50),col=color[10],xlab='Number of Strong Binders',ylab='IC50 (nM)')
par(new=TRUE)
lines(loess.smooth(1:length(luan_strong_filtered$IC50_Binder), y=luan_strong_filtered$IC50_Binder,degree = 2),xlim=c(0,150),ylim=c(0,50),col=color[11],xlab='Number of Strong Binders',ylab='IC50 (nM)')
par(new=TRUE)
lines(loess.smooth(1:length(mabi_strong_filtered$IC50_Binder), y=mabi_strong_filtered$IC50_Binder,degree = 2),xlim=c(0,150),ylim=c(0,50),col=color[12],xlab='Number of Strong Binders',ylab='IC50 (nM)')
par(new=TRUE)
lines(loess.smooth(1:length(moge_strong_filtered$IC50_Binder), y=moge_strong_filtered$IC50_Binder,degree = 2),xlim=c(0,150),ylim=c(0,50),col=color[13],xlab='Number of Strong Binders',ylab='IC50 (nM)')
par(new=TRUE)
lines(loess.smooth(1:length(piag_strong_filtered$IC50_Binder), y=piag_strong_filtered$IC50_Binder,degree = 2),xlim=c(0,150),ylim=c(0,50),col=color[14],xlab='Number of Strong Binders',ylab='IC50 (nM)')
par(new=TRUE)
lines(loess.smooth(1:length(prelu_strong_filtered$IC50_Binder), y=prelu_strong_filtered$IC50_Binder,degree = 2),xlim=c(0,150),ylim=c(0,50),col=color[15],xlab='Number of Strong Binders',ylab='IC50 (nM)')
legend("topright", inset=c(-0.2,0), c('alfe','besu','calu','deiv','dest','dr1','dr4','dr5','foca','gagra','luan','mabi','moge','piag','prelu'),col=color[1:15], title="Group")

