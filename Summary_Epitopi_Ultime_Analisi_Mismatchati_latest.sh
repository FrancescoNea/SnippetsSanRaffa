PRELU_StrongBinders_mhAGs_MHCII.sh
#PBS -l select=1:app=java:ncpus=4:mem=12gb
#PBS -P 161
#PBS -q workq
#PBS -m ae

python Summary_Epitopi_New_MHCII.py PRELU1_PRELU2_Somatic_MinorAGs_new_MHCII_PeptidesPrediction_MHCII.txt RESULTS_PRELU1_PRELU2_Somatic_MinorAGs_new_MHCII_PeptidesPrediction_MinorAntigens_MHCII.xls RESULTS_PRELU1_PRELU2_Somatic_MinorAGs_new_MHCII_strong_binders_MHC_II_2 PRELU_StrongBinders_Annotati.txt




#Da tenere a mente sempre quello che sta sotto

awk '{OFS="\t"; print $11":"$12"-"$12}'   /home/fsantaniello/Summaries_NeoEpitopi/MHCI/Summary_DEIV_RelOnly_NoHeaders.txt | sort | uniq  > Summary_DEIV_Annotato_NeoEpitopi_RelOnly_MHCI_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2_Allfiles_sorted.bam >> DEIV_NeoAGs_I_RelOnly_Binders_Regions_Diagnosis.txt; done < Summary_DEIV_Annotato_NeoEpitopi_RelOnly_MHCI_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2b_Allfiles_sorted.bam >> DEIV_NeoAGs_I_RelOnly_Binders_Regions_Relapse.txt; done < Summary_DEIV_Annotato_NeoEpitopi_RelOnly_MHCI_Regions.txt
 


deiv_NeoAG_mhci_RelOnly=read.table('/home/fsantaniello/Summaries_NeoEpitopi/MHCI/Summary_DEIV_RelOnly_NoHeaders.txt')
 colnames(deiv_NeoAG_mhci_RelOnly)=c("Gene_ID","Peptide","HLA","HLA","HLA","HLA","HLA","HLA","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA","IC50_HLA","IC,50_HLA","IC50_HLA","IC50_HLA","IC50_HLA")

deiv_NeoAG_mhci_RelOnly=merge(deiv_NeoAG_mhci_RelOnly,DEIV_RNK,by="genes",all.x=T)
relapse=read.table('DEIV_NeoAGs_I_RelOnly_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DEIV_NeoAGs_I_RelOnly_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
deiv_NeoAG_mhci_RelOnly=merge(deiv_NeoAG_mhci_RelOnly,relapse,by=c("chr","pos"),all.x=T)
deiv_NeoAG_mhci_RelOnly=merge(deiv_NeoAG_mhci_RelOnly,diagnosis,by=c("chr","pos"),all.x=T)
deiv_NeoAG_mhci_RelOnly$Rel_Mut=str_count(str_to_upper(deiv_NeoAG_mhci_RelOnly$rel_bases), as.character(deiv_NeoAG_mhci_RelOnly$mut_nt))
deiv_NeoAG_mhci_RelOnly$Rel_Wt=str_count(str_to_upper(deiv_NeoAG_mhci_RelOnly$rel_bases), as.character(deiv_NeoAG_mhci_RelOnly$wt_nt))
deiv_NeoAG_mhci_RelOnly$Dx_Mut=str_count(str_to_upper(deiv_NeoAG_mhci_RelOnly$dx_bases), as.character(deiv_NeoAG_mhci_RelOnly$mut_nt))
deiv_NeoAG_mhci_RelOnly$Dx_Wt=str_count(str_to_upper(deiv_NeoAG_mhci_RelOnly$dx_bases), as.character(deiv_NeoAG_mhci_RelOnly$wt_nt))
 write.table(unique(deiv_NeoAG_mhci_RelOnly[,c(1:26,35:38)]),"Summary_DEIV_Annotato_NeoAGs_MHCI_RelOnly_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)

mhag=read.table('Summary_DEIV_Annotato_NeoAGs_MHCI_RelOnly_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$DEIV2_kallisto[is.na(mhag$DEIV2b_kallisto)] <- 0
mhag$DEIV2_kallisto[is.na(mhag$DEIV2_kallisto)] <- 0

  nrow(mhag[apply(mhag[,17:22], MARGIN = 1, function(x) any(x<=500 & x> 50)), ])
[1] 10
  nrow(mhag[apply(mhag[,17:22], MARGIN = 1, function(x) any(x<=50 & x> 0)), ])
[1] 0

mhag=mhag[(mhag$DEIV2b_kallisto>=0.58 & mhag$DEIV2_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
8

> table(mhag$Rel_Mut>=3)

FALSE
   1
> table(mhag$Rel_Mut>mhag$Dx_Mut)

TRUE
   1

nrow(mhag[apply(mhag[,17:22], MARGIN = 1, function(x) any(x<=500 & x> 50)), ])
[1] 8
 nrow(mhag[apply(mhag[,17:22], MARGIN = 1, function(x) any(x<=50 & x> 0)), ])
[1] 0

#A livello di trascritto
 table(mhag$DEIV2b_kallisto > mhag$DEIV2_kallisto)

FALSE
    1
> table(mhag$DEIV2b_kallisto < mhag$DEIV2_kallisto)

TRUE
   1

wilcox.test(mhag$DEIV2b_kallisto, mhag$DEIV2_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$DEIV2_kallisto and mhag$DEIV1_kallisto
W = 17, p-value = 0.9497
alternative hypothesis: true location shift is greater than 0

Warning message:
In wilcox.test.default(mhag$DEIV2_kallisto, mhag$DEIV1_kallisto,  :
  cannot compute exact p-value with ties

table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
    4     4


wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 7, p-value = 0.9966
alternative hypothesis: true location shift is greater than 0

Warning message:
In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
  cannot compute exact p-value with ties


library(edgeR)
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


#Ricalcolo trascritti e numero reads alleli per NeoEpitopi ClasseI

#Calcolo reads su alleli per classe I


 for i in *Summary*txt; do echo $i; prefix=$(echo $i | sed 's/.txt/_NoHeaders.txt/g'); echo $prefix;  sed -e '/Strong/,+1d' $i | sed -e '/Weak/,+1d' - | grep -v WT_ > $prefix ; done

#ALFE RelDiag

 cat StrongBinders_Annotati_ALFE_RelDiag WeakBinders_Annotati_ALFE_RelDiag | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ > Summary_ALFE_Annotato_RelDiag_NoHeaders_NoMismatches.txt


ALFE1 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE1_Allfiles_sorted.bam
ALFE2 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE2_Allfiles_sorted.bam

awk '{OFS="\t"; print $8":"$9"-"$9}'   /home/fsantaniello/Summary_ALFE_Annotato_RelDiag_NoHeaders_NoMismatches.txt | sort | uniq  > Summary_ALFE_Annotato_RelDiag_NoHeaders_NoMismatches_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE1_Allfiles_sorted.bam >> ALFE_NeoAGs_I_RelDiag_Binders_Regions_Diagnosis_NoMismatches.txt; done < Summary_ALFE_Annotato_RelDiag_NoHeaders_NoMismatches_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE2_Allfiles_sorted.bam >> ALFE_NeoAGs_I_RelDiag_Binders_Regions_Relapse_NoMismatches.txt; done < Summary_ALFE_Annotato_RelDiag_NoHeaders_NoMismatches_Regions.txt


library(stringr)

alfe_NeoAG_mhci_RelDiag=read.table('Summary_ALFE_Annotato_RelDiag_NoHeaders_NoMismatches.txt')
 colnames(alfe_NeoAG_mhci_RelDiag)=c("Gene_ID","Peptide","HLA","HLA","HLA","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA","IC50_HLA","IC,50_HLA")

 alfe_NeoAG_mhci_RelDiag=merge(alfe_NeoAG_mhci_RelDiag,ALFE_RNK,by="genes",all.x=T)
relapse=read.table('ALFE_NeoAGs_I_RelDiag_Binders_Regions_Relapse_NoMismatches.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('ALFE_NeoAGs_I_RelDiag_Binders_Regions_Diagnosis_NoMismatches.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
alfe_NeoAG_mhci_RelDiag=merge(alfe_NeoAG_mhci_RelDiag,relapse,by=c("chr","pos"),all.x=T)
alfe_NeoAG_mhci_RelDiag=merge(alfe_NeoAG_mhci_RelDiag,diagnosis,by=c("chr","pos"),all.x=T)
alfe_NeoAG_mhci_RelDiag$Rel_Mut=str_count(str_to_upper(alfe_NeoAG_mhci_RelDiag$rel_bases), as.character(alfe_NeoAG_mhci_RelDiag$mut_nt))
alfe_NeoAG_mhci_RelDiag$Rel_Wt=str_count(str_to_upper(alfe_NeoAG_mhci_RelDiag$rel_bases), as.character(alfe_NeoAG_mhci_RelDiag$wt_nt))
alfe_NeoAG_mhci_RelDiag$Dx_Mut=str_count(str_to_upper(alfe_NeoAG_mhci_RelDiag$dx_bases), as.character(alfe_NeoAG_mhci_RelDiag$mut_nt))
alfe_NeoAG_mhci_RelDiag$Dx_Wt=str_count(str_to_upper(alfe_NeoAG_mhci_RelDiag$dx_bases), as.character(alfe_NeoAG_mhci_RelDiag$wt_nt))
 write.table(unique(alfe_NeoAG_mhci_RelDiag[,c(1:18,27:30)]),"../../Documents/Progetto_LuciaGabri/Summary_ALFE_Annotato_NeoAGs_MHCI_RelDiag_NoMismatches_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)


#ALFE RelOnly

cat StrongBinders_Annotati_ALFE_RelOnly WeakBinders_Annotati_ALFE_RelOnly | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ > Summary_ALFE_Annotato_RelOnly_NoHeaders_NoMismatches.txt


ALFE1 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE1_Allfiles_sorted.bam
ALFE2 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE2_Allfiles_sorted.bam

awk '{OFS="\t"; print $8":"$9"-"$9}'   /home/fsantaniello/Summary_ALFE_Annotato_RelOnly_NoHeaders_NoMismatches.txt | sort | uniq  > Summary_ALFE_Annotato_RelOnly_NoHeaders_NoMismatches_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE1_Allfiles_sorted.bam >> ALFE_NeoAGs_I_RelOnly_Binders_Regions_Diagnosis_NoMismatches.txt; done < Summary_ALFE_Annotato_RelOnly_NoHeaders_NoMismatches_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE2_Allfiles_sorted.bam >> ALFE_NeoAGs_I_RelOnly_Binders_Regions_Relapse_NoMismatches.txt; done < Summary_ALFE_Annotato_RelOnly_NoHeaders_NoMismatches_Regions.txt


library(stringr)

alfe_NeoAG_mhci_RelOnly=read.table('Summary_ALFE_Annotato_RelOnly_NoHeaders_NoMismatches.txt')
colnames(alfe_NeoAG_mhci_RelOnly)=c("Gene_ID","Peptide","HLA","HLA","HLA","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA","IC50_HLA","IC,50_HLA")

 alfe_NeoAG_mhci_RelOnly=merge(alfe_NeoAG_mhci_RelOnly,ALFE_RNK,by="genes",all.x=T)
relapse=read.table('ALFE_NeoAGs_I_RelOnly_Binders_Regions_Relapse_NoMismatches.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('ALFE_NeoAGs_I_RelOnly_Binders_Regions_Diagnosis_NoMismatches.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
alfe_NeoAG_mhci_RelOnly=merge(alfe_NeoAG_mhci_RelOnly,relapse,by=c("chr","pos"),all.x=T)
alfe_NeoAG_mhci_RelOnly=merge(alfe_NeoAG_mhci_RelOnly,diagnosis,by=c("chr","pos"),all.x=T)
alfe_NeoAG_mhci_RelOnly$Rel_Mut=str_count(str_to_upper(alfe_NeoAG_mhci_RelOnly$rel_bases), as.character(alfe_NeoAG_mhci_RelOnly$mut_nt))
alfe_NeoAG_mhci_RelOnly$Rel_Wt=str_count(str_to_upper(alfe_NeoAG_mhci_RelOnly$rel_bases), as.character(alfe_NeoAG_mhci_RelOnly$wt_nt))
alfe_NeoAG_mhci_RelOnly$Dx_Mut=str_count(str_to_upper(alfe_NeoAG_mhci_RelOnly$dx_bases), as.character(alfe_NeoAG_mhci_RelOnly$mut_nt))
alfe_NeoAG_mhci_RelOnly$Dx_Wt=str_count(str_to_upper(alfe_NeoAG_mhci_RelOnly$dx_bases), as.character(alfe_NeoAG_mhci_RelOnly$wt_nt))
 write.table(unique(alfe_NeoAG_mhci_RelOnly[,c(1:16,27:30)]),"../../Documents/Progetto_LuciaGabri/Summary_ALFE_Annotato_NeoAGs_MHCI_RelOnly_AllInfos_NoMismatches.txt",sep="\t",col.names=T,row.names=F,quote=F)
#ALFE NON HA MUTAZIONI REL ONLY

#ALFE DiagOnly

cat StrongBinders_Annotati_ALFE_DiagOnly WeakBinders_Annotati_ALFE_DiagOnly | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ > Summary_ALFE_Annotato_DiagOnly_NoHeaders_NoMismatches.txt


ALFE1 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE1_Allfiles_sorted.bam
ALFE2 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE2_Allfiles_sorted.bam

awk '{OFS="\t"; print $8":"$9"-"$9}'   ../../Documents/Progetto_LuciaGabri/Summary_ALFE_Annotato_DiagOnly_NoHeaders_NoMismatches.txt | sort | uniq  > Summary_ALFE_Annotato_DiagOnly_NoHeaders_NoMismatches_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE1_Allfiles_sorted.bam >> ALFE_NeoAGs_I_DiagOnly_Binders_Regions_Diagnosis_NoMismatches.txt; done < Summary_ALFE_Annotato_DiagOnly_NoHeaders_NoMismatches_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE2_Allfiles_sorted.bam >> ALFE_NeoAGs_I_DiagOnly_Binders_Regions_Relapse_NoMismatches.txt; done < Summary_ALFE_Annotato_DiagOnly_NoHeaders_NoMismatches_Regions.txt


library(stringr)

alfe_NeoAG_mhci_DiagOnly=read.table('../../Documents/Progetto_LuciaGabri/Summary_ALFE_Annotato_DiagOnly_NoHeaders_NoMismatches.txt')
 colnames(alfe_NeoAG_mhci_DiagOnly)=c("Gene_ID","Peptide","HLA","HLA","HLA","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA","IC50_HLA","IC,50_HLA")

 alfe_NeoAG_mhci_DiagOnly=merge(alfe_NeoAG_mhci_DiagOnly,ALFE_RNK,by="genes",all.x=T)
relapse=read.table('ALFE_NeoAGs_I_DiagOnly_Binders_Regions_Relapse_NoMismatches.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('ALFE_NeoAGs_I_DiagOnly_Binders_Regions_Diagnosis_NoMismatches.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
alfe_NeoAG_mhci_DiagOnly=merge(alfe_NeoAG_mhci_DiagOnly,relapse,by=c("chr","pos"),all.x=T)
alfe_NeoAG_mhci_DiagOnly=merge(alfe_NeoAG_mhci_DiagOnly,diagnosis,by=c("chr","pos"),all.x=T)
alfe_NeoAG_mhci_DiagOnly$Rel_Mut=str_count(str_to_upper(alfe_NeoAG_mhci_DiagOnly$rel_bases), as.character(alfe_NeoAG_mhci_DiagOnly$mut_nt))
alfe_NeoAG_mhci_DiagOnly$Rel_Wt=str_count(str_to_upper(alfe_NeoAG_mhci_DiagOnly$rel_bases), as.character(alfe_NeoAG_mhci_DiagOnly$wt_nt))
alfe_NeoAG_mhci_DiagOnly$Dx_Mut=str_count(str_to_upper(alfe_NeoAG_mhci_DiagOnly$dx_bases), as.character(alfe_NeoAG_mhci_DiagOnly$mut_nt))
alfe_NeoAG_mhci_DiagOnly$Dx_Wt=str_count(str_to_upper(alfe_NeoAG_mhci_DiagOnly$dx_bases), as.character(alfe_NeoAG_mhci_DiagOnly$wt_nt))
 write.table(unique(alfe_NeoAG_mhci_DiagOnly[,c(1:16,27:30)]),"Summary_ALFE_Annotato_NeoAGs_MHCI_DiagOnly_AllInfos_NoMismatches.txt",sep="\t",col.names=T,row.names=F,quote=F)
#ALFE NON HA MUTAZIONI DIAG ONLY




#BESU
#RelDiag
cat StrongBinders_Annotati_BESU_RelDiag WeakBinders_Annotati_BESU_RelDiag | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ > Summary_BESU_Annotato_RelDiag_NoHeaders_NoMismatches.txt

BESU1 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/BESU1_Allfiles_sorted.bam
BESU2 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/BESU2_Allfiles_sorted.bam


awk '{OFS="\t"; print $11":"$12"-"$12}'   Summary_BESU_Annotato_RelDiag_NoHeaders_NoMismatches.txt | sort | uniq  > Summary_BESU_Annotato_RelDiag_NoHeaders_NoMismatches_Regions.txt
while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/BESU1_Allfiles_sorted.bam >> BESU_NeoAGs_I_RelDiag_Binders_Regions_Diagnosis.txt; done < Summary_BESU_Annotato_RelDiag_NoHeaders_NoMismatches_Regions.txt
while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/BESU2_Allfiles_sorted.bam >> BESU_NeoAGs_I_RelDiag_Binders_Regions_Relapse.txt; done < Summary_BESU_Annotato_RelDiag_NoHeaders_NoMismatches_Regions.txt
 

besu_NeoAG_mhci_RelDiag=read.table('Summary_BESU_Annotato_RelDiag_NoHeaders_NoMismatches.txt')
colnames(besu_NeoAG_mhci_RelDiag)=c("Gene_ID","Peptide","HLA","HLA","HLA","HLA","HLA","HLA","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA","IC50_HLA","IC,50_HLA","IC50_HLA","IC50_HLA","IC50_HLA")

besu_NeoAG_mhci_RelDiag=merge(besu_NeoAG_mhci_RelDiag,BESU_RNK,by="genes",all.x=T)
relapse=read.table('BESU_NeoAGs_I_RelDiag_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('BESU_NeoAGs_I_RelDiag_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
besu_NeoAG_mhci_RelDiag=merge(besu_NeoAG_mhci_RelDiag,relapse,by=c("chr","pos"),all.x=T)
besu_NeoAG_mhci_RelDiag=merge(besu_NeoAG_mhci_RelDiag,diagnosis,by=c("chr","pos"),all.x=T)
besu_NeoAG_mhci_RelDiag$Rel_Mut=str_count(str_to_upper(besu_NeoAG_mhci_RelDiag$rel_bases), as.character(besu_NeoAG_mhci_RelDiag$mut_nt))
besu_NeoAG_mhci_RelDiag$Rel_Wt=str_count(str_to_upper(besu_NeoAG_mhci_RelDiag$rel_bases), as.character(besu_NeoAG_mhci_RelDiag$wt_nt))
besu_NeoAG_mhci_RelDiag$Dx_Mut=str_count(str_to_upper(besu_NeoAG_mhci_RelDiag$dx_bases), as.character(besu_NeoAG_mhci_RelDiag$mut_nt))
besu_NeoAG_mhci_RelDiag$Dx_Wt=str_count(str_to_upper(besu_NeoAG_mhci_RelDiag$dx_bases), as.character(besu_NeoAG_mhci_RelDiag$wt_nt))
 write.table(unique(besu_NeoAG_mhci_RelDiag[,c(1:24,33:36)]),"../../Documents/Progetto_LuciaGabri/Summary_BESU_Annotato_NeoAGs_MHCI_RelDiag_AllInfos_NoMismatches.txt",sep="\t",col.names=T,row.names=F,quote=F)


#RelOnly


#BESU non ha NeoAGs RelOnly



#BESU non ha NeoAGs DiagOnly 





#CALU
#RelDiag
cat StrongBinders_Annotati_CALU_RelDiag WeakBinders_Annotati_CALU_RelDiag | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ > Summary_CALU_Annotato_RelDiag_NoHeaders_NoMismatches.txt

#CALU
CALU1 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/CALU1_Allfiles_sorted.bam
CALU2 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/CALU2_Allfiles_sorted.bam



awk '{OFS="\t"; print $8":"$9"-"$9}'   Summary_CALU_Annotato_RelDiag_NoHeaders_NoMismatches.txt | sort | uniq  > Summary_CALU_Annotato_RelDiag_NoHeaders_NoMismatches_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/CALU1_Allfiles_sorted.bam >> CALU_NeoAGs_I_RelDiag_Binders_Regions_Diagnosis_NoMismatches.txt; done < Summary_CALU_Annotato_RelDiag_NoHeaders_NoMismatches_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/CALU2_Allfiles_sorted.bam >> CALU_NeoAGs_I_RelDiag_Binders_Regions_Relapse_NoMismatches.txt; done < Summary_CALU_Annotato_RelDiag_NoHeaders_NoMismatches_Regions.txt
 


calu_NeoAG_mhci_RelDiag=read.table('Summary_CALU_Annotato_RelDiag_NoHeaders_NoMismatches.txt')
colnames(calu_NeoAG_mhci_RelDiag)=c("Gene_ID","Peptide","HLA","HLA","HLA","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA","IC50_HLA","IC50_HLA")
calu_NeoAG_mhci_RelDiag=merge(calu_NeoAG_mhci_RelDiag,CALU_RNK,by="genes",all.x=T)
relapse=read.table('CALU_NeoAGs_I_RelDiag_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('CALU_NeoAGs_I_RelDiag_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
calu_NeoAG_mhci_RelDiag=merge(calu_NeoAG_mhci_RelDiag,relapse,by=c("chr","pos"),all.x=T)
calu_NeoAG_mhci_RelDiag=merge(calu_NeoAG_mhci_RelDiag,diagnosis,by=c("chr","pos"),all.x=T)
calu_NeoAG_mhci_RelDiag$Rel_Mut=str_count(str_to_upper(calu_NeoAG_mhci_RelDiag$rel_bases), as.character(calu_NeoAG_mhci_RelDiag$mut_nt))
calu_NeoAG_mhci_RelDiag$Rel_Wt=str_count(str_to_upper(calu_NeoAG_mhci_RelDiag$rel_bases), as.character(calu_NeoAG_mhci_RelDiag$wt_nt))
calu_NeoAG_mhci_RelDiag$Dx_Mut=str_count(str_to_upper(calu_NeoAG_mhci_RelDiag$dx_bases), as.character(calu_NeoAG_mhci_RelDiag$mut_nt))
calu_NeoAG_mhci_RelDiag$Dx_Wt=str_count(str_to_upper(calu_NeoAG_mhci_RelDiag$dx_bases), as.character(calu_NeoAG_mhci_RelDiag$wt_nt))
 write.table(unique(calu_NeoAG_mhci_RelDiag[,c(1:18,27:30)]),"../../Documents/Progetto_LuciaGabri/Summary_CALU_Annotato_NeoAGs_MHCI_RelDiag_NoMismatches_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)


#CALU non ha NeoAGs RelOnly



#CALU non ha NeoAGs DiagOnly 




#DEIV
#DEIV non ha NeoAGs RelDiag

cat StrongBinders_Annotati_DEIV_DiagOnly WeakBinders_Annotati_DEIV_DiagOnly | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ > Summary_DEIV_Annotato_DiagOnly_NoHeaders_NoMismatches.txt


#DiagOnly
DEIV2 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2_Allfiles_sorted.bam
DEIV2b = /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2b_Allfiles_sorted.bam


awk '{OFS="\t"; print $11":"$12"-"$12}'   Summary_DEIV_Annotato_DiagOnly_NoHeaders_NoMismatches.txt | sort | uniq  > Summary_DEIV_Annotato_NeoEpitopi_DiagOnly_MHCI_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2_Allfiles_sorted.bam >> DEIV_NeoAGs_I_DiagOnly_Binders_Regions_Diagnosis.txt; done < Summary_DEIV_Annotato_NeoEpitopi_DiagOnly_MHCI_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2b_Allfiles_sorted.bam >> DEIV_NeoAGs_I_DiagOnly_Binders_Regions_Relapse.txt; done < Summary_DEIV_Annotato_NeoEpitopi_DiagOnly_MHCI_Regions.txt
 


deiv_NeoAG_mhci_DiagOnly=read.table('Summary_DEIV_Annotato_DiagOnly_NoHeaders_NoMismatches.txt')
colnames(deiv_NeoAG_mhci_DiagOnly)=c("Gene_ID","Peptide","HLA","HLA","HLA","HLA","HLA","HLA","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA","IC50_HLA","IC,50_HLA","IC50_HLA","IC50_HLA","IC50_HLA")

deiv_NeoAG_mhci_DiagOnly=merge(deiv_NeoAG_mhci_DiagOnly,DEIV_RNK,by="genes",all.x=T)
relapse=read.table('DEIV_NeoAGs_I_DiagOnly_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DEIV_NeoAGs_I_DiagOnly_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
deiv_NeoAG_mhci_DiagOnly=merge(deiv_NeoAG_mhci_DiagOnly,relapse,by=c("chr","pos"),all.x=T)
deiv_NeoAG_mhci_DiagOnly=merge(deiv_NeoAG_mhci_DiagOnly,diagnosis,by=c("chr","pos"),all.x=T)
deiv_NeoAG_mhci_DiagOnly$Rel_Mut=str_count(str_to_upper(deiv_NeoAG_mhci_DiagOnly$rel_bases), as.character(deiv_NeoAG_mhci_DiagOnly$mut_nt))
deiv_NeoAG_mhci_DiagOnly$Rel_Wt=str_count(str_to_upper(deiv_NeoAG_mhci_DiagOnly$rel_bases), as.character(deiv_NeoAG_mhci_DiagOnly$wt_nt))
deiv_NeoAG_mhci_DiagOnly$Dx_Mut=str_count(str_to_upper(deiv_NeoAG_mhci_DiagOnly$dx_bases), as.character(deiv_NeoAG_mhci_DiagOnly$mut_nt))
deiv_NeoAG_mhci_DiagOnly$Dx_Wt=str_count(str_to_upper(deiv_NeoAG_mhci_DiagOnly$dx_bases), as.character(deiv_NeoAG_mhci_DiagOnly$wt_nt))
 write.table(unique(deiv_NeoAG_mhci_DiagOnly[,c(1:24,33:36)]),"../../Documents/Progetto_LuciaGabri/Summary_DEIV_Annotato_NeoAGs_MHCI_DiagOnly_NoMismatches_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)



#DEIV RelOnly
cat StrongBinders_Annotati_DEIV_RelOnly WeakBinders_Annotati_DEIV_RelOnly | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ > Summary_DEIV_Annotato_RelOnly_NoHeaders_NoMismatches.txt

DEIV2 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2_Allfiles_sorted.bam
DEIV2b = /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2b_Allfiles_sorted.bam


awk '{OFS="\t"; print $11":"$12"-"$12}'   Summary_DEIV_Annotato_RelOnly_NoHeaders_NoMismatches.txt | sort | uniq  > Summary_DEIV_Annotato_NeoEpitopi_RelOnly_MHCI_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2_Allfiles_sorted.bam >> DEIV_NeoAGs_I_RelOnly_Binders_Regions_Diagnosis.txt; done < Summary_DEIV_Annotato_NeoEpitopi_RelOnly_MHCI_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2b_Allfiles_sorted.bam >> DEIV_NeoAGs_I_RelOnly_Binders_Regions_Relapse.txt; done < Summary_DEIV_Annotato_NeoEpitopi_RelOnly_MHCI_Regions.txt
 


deiv_NeoAG_mhci_RelOnly=read.table('Summary_DEIV_Annotato_RelOnly_NoHeaders_NoMismatches.txt')
colnames(deiv_NeoAG_mhci_RelOnly)=c("Gene_ID","Peptide","HLA","HLA","HLA","HLA","HLA","HLA","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA","IC50_HLA","IC,50_HLA","IC50_HLA","IC50_HLA","IC50_HLA")

deiv_NeoAG_mhci_RelOnly=merge(deiv_NeoAG_mhci_RelOnly,DEIV_RNK,by="genes",all.x=T)
relapse=read.table('DEIV_NeoAGs_I_RelOnly_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DEIV_NeoAGs_I_RelOnly_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
deiv_NeoAG_mhci_RelOnly=merge(deiv_NeoAG_mhci_RelOnly,relapse,by=c("chr","pos"),all.x=T)
deiv_NeoAG_mhci_RelOnly=merge(deiv_NeoAG_mhci_RelOnly,diagnosis,by=c("chr","pos"),all.x=T)
deiv_NeoAG_mhci_RelOnly$Rel_Mut=str_count(str_to_upper(deiv_NeoAG_mhci_RelOnly$rel_bases), as.character(deiv_NeoAG_mhci_RelOnly$mut_nt))
deiv_NeoAG_mhci_RelOnly$Rel_Wt=str_count(str_to_upper(deiv_NeoAG_mhci_RelOnly$rel_bases), as.character(deiv_NeoAG_mhci_RelOnly$wt_nt))
deiv_NeoAG_mhci_RelOnly$Dx_Mut=str_count(str_to_upper(deiv_NeoAG_mhci_RelOnly$dx_bases), as.character(deiv_NeoAG_mhci_RelOnly$mut_nt))
deiv_NeoAG_mhci_RelOnly$Dx_Wt=str_count(str_to_upper(deiv_NeoAG_mhci_RelOnly$dx_bases), as.character(deiv_NeoAG_mhci_RelOnly$wt_nt))
 write.table(unique(deiv_NeoAG_mhci_RelOnly[,c(1:24,33:36)]),"../../Documents/Progetto_LuciaGabri/Summary_DEIV_Annotato_NeoAGs_MHCI_RelOnly_NoMismatches_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)




#DEST
#Dest non ha neoAGs DiagOnly
#Dest RelDiag
cat StrongBinders_Annotati_DEST_RelDiag WeakBinders_Annotati_DEST_RelDiag | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ > Summary_DEST_Annotato_RelDiag_NoHeaders_NoMismatches.txt



awk '{OFS="\t"; print $8":"$9"-"$9}'   Summary_DEST_Annotato_RelDiag_NoHeaders_NoMismatches.txt | sort | uniq  > Summary_DEST_Annotato_NeoEpitopi_RelDiag_MHCI_Regions.txt

DEST1 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST1_Allfiles_sorted.bam
DEST2 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST2_Allfiles_sorted.bam


 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST1_Allfiles_sorted.bam >> DEST_NeoAGs_I_RelDiag_Binders_Regions_Diagnosis.txt; done < Summary_DEST_Annotato_NeoEpitopi_RelDiag_MHCI_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST2_Allfiles_sorted.bam >> DEST_NeoAGs_I_RelDiag_Binders_Regions_Relapse.txt; done < Summary_DEST_Annotato_NeoEpitopi_RelDiag_MHCI_Regions.txt
 


dest_NeoAG_mhci_RelDiag=read.table('Summary_DEST_Annotato_RelDiag_NoHeaders_NoMismatches.txt')
 colnames(dest_NeoAG_mhci_RelDiag)=c("Gene_ID","Peptide","HLA","HLA","HLA","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA","IC50_HLA","IC,50_HLA")

dest_NeoAG_mhci_RelDiag=merge(dest_NeoAG_mhci_RelDiag,DEST_RNK,by="genes",all.x=T)
relapse=read.table('DEST_NeoAGs_I_RelDiag_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DEST_NeoAGs_I_RelDiag_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
dest_NeoAG_mhci_RelDiag=merge(dest_NeoAG_mhci_RelDiag,relapse,by=c("chr","pos"),all.x=T)
dest_NeoAG_mhci_RelDiag=merge(dest_NeoAG_mhci_RelDiag,diagnosis,by=c("chr","pos"),all.x=T)
dest_NeoAG_mhci_RelDiag$Rel_Mut=str_count(str_to_upper(dest_NeoAG_mhci_RelDiag$rel_bases), as.character(dest_NeoAG_mhci_RelDiag$mut_nt))
dest_NeoAG_mhci_RelDiag$Rel_Wt=str_count(str_to_upper(dest_NeoAG_mhci_RelDiag$rel_bases), as.character(dest_NeoAG_mhci_RelDiag$wt_nt))
dest_NeoAG_mhci_RelDiag$Dx_Mut=str_count(str_to_upper(dest_NeoAG_mhci_RelDiag$dx_bases), as.character(dest_NeoAG_mhci_RelDiag$mut_nt))
dest_NeoAG_mhci_RelDiag$Dx_Wt=str_count(str_to_upper(dest_NeoAG_mhci_RelDiag$dx_bases), as.character(dest_NeoAG_mhci_RelDiag$wt_nt))
 write.table(unique(dest_NeoAG_mhci_RelDiag[,c(1:18,27:30)]),"../../Documents/Progetto_LuciaGabri/Summary_DEST_Annotato_NeoAGs_MHCI_RelDiag_NoMismatches_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)



#RelOnly
#Dest RelOnly
cat StrongBinders_Annotati_DEST_RelOnly WeakBinders_Annotati_DEST_RelOnly | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ > Summary_DEST_Annotato_RelOnly_NoHeaders_NoMismatches.txt



awk '{OFS="\t"; print $8":"$9"-"$9}'   Summary_DEST_Annotato_RelOnly_NoHeaders_NoMismatches.txt | sort | uniq  > Summary_DEST_Annotato_NeoEpitopi_RelOnly_MHCI_Regions.txt

DEST1 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST1_Allfiles_sorted.bam
DEST2 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST2_Allfiles_sorted.bam


 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST1_Allfiles_sorted.bam >> DEST_NeoAGs_I_RelOnly_Binders_Regions_Diagnosis.txt; done < Summary_DEST_Annotato_NeoEpitopi_RelOnly_MHCI_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST2_Allfiles_sorted.bam >> DEST_NeoAGs_I_RelOnly_Binders_Regions_Relapse.txt; done < Summary_DEST_Annotato_NeoEpitopi_RelOnly_MHCI_Regions.txt
 


dest_NeoAG_mhci_RelOnly=read.table('Summary_DEST_Annotato_RelOnly_NoHeaders_NoMismatches.txt')
colnames(dest_NeoAG_mhci_RelOnly)=c("Gene_ID","Peptide","HLA","HLA","HLA","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA","IC50_HLA","IC,50_HLA")

dest_NeoAG_mhci_RelOnly=merge(dest_NeoAG_mhci_RelOnly,DEST_RNK,by="genes",all.x=T)
relapse=read.table('DEST_NeoAGs_I_RelOnly_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DEST_NeoAGs_I_RelOnly_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
dest_NeoAG_mhci_RelOnly=merge(dest_NeoAG_mhci_RelOnly,relapse,by=c("chr","pos"),all.x=T)
dest_NeoAG_mhci_RelOnly=merge(dest_NeoAG_mhci_RelOnly,diagnosis,by=c("chr","pos"),all.x=T)
dest_NeoAG_mhci_RelOnly$Rel_Mut=str_count(str_to_upper(dest_NeoAG_mhci_RelOnly$rel_bases), as.character(dest_NeoAG_mhci_RelOnly$mut_nt))
dest_NeoAG_mhci_RelOnly$Rel_Wt=str_count(str_to_upper(dest_NeoAG_mhci_RelOnly$rel_bases), as.character(dest_NeoAG_mhci_RelOnly$wt_nt))
dest_NeoAG_mhci_RelOnly$Dx_Mut=str_count(str_to_upper(dest_NeoAG_mhci_RelOnly$dx_bases), as.character(dest_NeoAG_mhci_RelOnly$mut_nt))
dest_NeoAG_mhci_RelOnly$Dx_Wt=str_count(str_to_upper(dest_NeoAG_mhci_RelOnly$dx_bases), as.character(dest_NeoAG_mhci_RelOnly$wt_nt))
 write.table(unique(dest_NeoAG_mhci_RelOnly[,c(1:18,27:30)]),"../../Documents/Progetto_LuciaGabri/Summary_DEST_Annotato_NeoAGs_MHCI_RelOnly_NoMismatches_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)


#
#mhag=read.table('Summary_DEST_Annotato_NeoAGs_MHCI_RelOnly_AllInfos.txt',head=T)
#mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
#mhag$logFC[is.na(mhag$logFC)] <- 0
#mhag$DEST2_kallisto[is.na(mhag$DEST2_kallisto)] <- 0
#mhag$DEST1_kallisto[is.na(mhag$DEST1_kallisto)] <- 0
#
#   row
#[1] 6
#
#  nrow(mhag[apply(mhag[,17:22], MARGIN = 1, function(x) any(x<=50 & x> 0)), ])
#[1] 1
#
#  
#mhag=mhag[(mhag$DEST2_kallisto>=0.58 & mhag$DEST1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
#nrow(mhag)
#4
#
#
#
#table(mhag$Rel_Mut>=3)
#
#TRUE
#   8
#table(mhag$Dx_Mut>=3)
#
#   nrow(mhag[apply(mhag[,17:22], MARGIN = 1, function(x) any(x<=500 & x> 50)), ])
#[1] 3
#
#  nrow(mhag[apply(mhag[,17:22], MARGIN = 1, function(x) any(x<=50 & x> 0)), ])
#[1] 1
#  
#  
##A livello di trascritto
# table(mhag$DEST2_kallisto > mhag$DEST1_kallisto)
#
#FALSE
#    4
#
# wilcox.test(mhag$DEST2_kallisto, mhag$DEST1_kallisto,alternative="greater")
#
#        Wilcoxon rank sum test with continuity correction
#
#data:  mhag$DEST2_kallisto and mhag$DEST1_kallisto
#W = 0, p-value = 0.9975
#alternative hypothesis: true location shift is greater than 0
#
#Warning message:
#In wilcox.test.default(mhag$DEST2_kallisto, mhag$DEST1_kallisto,  :
#  cannot compute exact p-value with ties
#
#table(mhag$Rel_Mut>mhag$Dx_Mut)
#
#FALSE  TRUE
#    4     0
#
#
#wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")
#
#        Wilcoxon rank sum test with continuity correction
#
#data:  mhag$Rel_Mut and mhag$Dx_Mut
#W = 7, p-value = 0.9966
#alternative hypothesis: true location shift is greater than 0
#
#Warning message:
#In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
#  cannot compute exact p-value with ties


#DR1 RelOnly

sh Launcher_snp2epi_Paper_Annotation.sh DR1_Somatic_Rerun.tsv DR11_DR12_MHCI_NoMismatchHla HLA-A01:01,HLA-A03:01,HLA-B35:01,HLA-B15:01,HLA-C04:01,HLA-C03:04 DR1_Expression_Kallisto_NeoEpitopes.txt

       #Summary Binders Diagnosis Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A01:01\tHLA-A03:01\tHLA-B35:01\tHLA-B15:01\tHLA-C04:01\tHLA-C03:04" )
        more +2 ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_weak_binders_DiagOnly
#         ../../Documents/Progetto_LuciaGabri/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt ../../Documents/Progetto_LuciaGabri/weak_binders_DiagOnly ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_weak_binders_annotated_DiagOnly
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A01:01\tHLA-A03:01\tHLA-B35:01\tHLA-B15:01\tHLA-C04:01\tHLA-C03:04" )
        more +2 ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls |   awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' -  | (echo -e $strong_binders; cat -) > ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_strong_binders_DiagOnly
 #        ../../Documents/Progetto_LuciaGabri/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt ../../Documents/Progetto_LuciaGabri/strong_binders_DiagOnly ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls ../../Documents/Progetto_LuciaGabri/strong_binders_annotated_DiagOnly

        #cat ../../Documents/Progetto_LuciaGabri/strong_binders_annotated_DiagOnly ../../Documents/Progetto_LuciaGabri/weak_binders_annotated_DiagOnly > ./Summary_Binders_RESULTS_DR11_DR12_MHCI_NoMismatchHla__PeptidesPrediction_DiagnosisOnly.txt



        #Summary Binders Relapse Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A01:01\tHLA-A03:01\tHLA-B35:01\tHLA-B15:01\tHLA-C04:01\tHLA-C03:04" )
        more +2 ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_weak_binders_RelOnly
        # ../../Documents/Progetto_LuciaGabri/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt ../../Documents/Progetto_LuciaGabri/weak_binders_RelOnly ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls ../../Documents/Progetto_LuciaGabri/weak_binders_annotated_RelOnly

        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A01:01\tHLA-A03:01\tHLA-B35:01\tHLA-B15:01\tHLA-C04:01\tHLA-C03:04" )
        more +2 ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' -  | (echo -e $strong_binders; cat -) > ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_strong_binders_RelOnly
         #../../Documents/Progetto_LuciaGabri/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt ../../Documents/Progetto_LuciaGabri/strong_binders_RelOnly ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls ../../Documents/Progetto_LuciaGabri/strong_binders_annotated_RelOnly
        
        #cat ../../Documents/Progetto_LuciaGabri/strong_binders_annotated_RelOnly ../../Documents/Progetto_LuciaGabri/weak_binders_annotated_RelOnly > ./Summary_Binders_RESULTS_DR11_DR12_MHCI_NoMismatchHla__PeptidesPrediction_RelapseOnly.txt
        
        
        #Summary Binders Common Diagnosis Relapse
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A01:01\tHLA-A03:01\tHLA-B35:01\tHLA-B15:01\tHLA-C04:01\tHLA-C03:04" )
        more +2 ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_weak_binders_RelDiag
  #       ../../Documents/Progetto_LuciaGabri/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt ../../Documents/Progetto_LuciaGabri/weak_binders_RelDiag ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls ../../Documents/Progetto_LuciaGabri/weak_binders_annotated_RelDiag
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A01:01\tHLA-A03:01\tHLA-B35:01\tHLA-B15:01\tHLA-C04:01\tHLA-C03:04" )
        more +2 ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' -  | (echo -e $strong_binders; cat -) > ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_strong_binders_RelDiag
   #      ../../Documents/Progetto_LuciaGabri/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt ../../Documents/Progetto_LuciaGabri/strong_binders_RelDiag ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls ../../Documents/Progetto_LuciaGabri/strong_binders_annotated_RelDiag
       
        #cat ../../Documents/Progetto_LuciaGabri/strong_binders_annotated_RelDiag ../../Documents/Progetto_LuciaGabri/weak_binders_annotated_RelDiag > ./Summary_Binders_RESULTS_DR11_DR12_MHCI_NoMismatchHla__PeptidesPrediction_CommonRelapseDiagnosisOnly.txt


#DR1
#DR1 non ha NeoAGs DiagOnly
#DR1 RelDiag
cat StrongBinders_Annotati_DR1_RelDiag WeakBinders_Annotati_DR1_RelDiag | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ > Summary_DR1_Annotato_RelDiag_NoHeaders_NoMismatches.txt



awk '{OFS="\t"; print $11":"$12"-"$12}'   Summary_DR1_Annotato_RelDiag_NoHeaders_NoMismatches.txt | sort | uniq  > Summary_DR1_Annotato_NeoEpitopi_RelDiag_MHCI_Regions.txt




DR1_1 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DR1_1_Allfiles_sorted.bam
DR1_2 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DR1_2_Allfiles_sorted.bam


 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DR1_1_Allfiles_sorted.bam >> DR1_NeoAGs_I_RelDiag_Binders_Regions_Relapse.txt; done < Summary_DR1_Annotato_NeoEpitopi_RelDiag_MHCI_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DR1_2_Allfiles_sorted.bam >> DR1_NeoAGs_I_RelDiag_Binders_Regions_Relapse.txt; done < Summary_DR1_Annotato_NeoEpitopi_RelDiag_MHCI_Regions.txt
 


dr1_NeoAG_mhci_RelDiag=read.table('Summary_DR1_Annotato_RelDiag_NoHeaders_NoMismatches.txt')
colnames(dr1_NeoAG_mhci_RelDiag)=c("Gene_ID","Peptide","HLA","HLA","HLA","HLA","HLA","HLA","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA","IC50_HLA","IC,50_HLA","IC50_HLA","IC50_HLA","IC,50_HLA")

dr1_NeoAG_mhci_RelDiag=merge(dr1_NeoAG_mhci_RelDiag,DR1_RNK,by="genes",all.x=T)
relapse=read.table('DR1_NeoAGs_I_RelDiag_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DR1_NeoAGs_I_RelDiag_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
dr1_NeoAG_mhci_RelDiag=merge(dr1_NeoAG_mhci_RelDiag,relapse,by=c("chr","pos"),all.x=T)
dr1_NeoAG_mhci_RelDiag=merge(dr1_NeoAG_mhci_RelDiag,diagnosis,by=c("chr","pos"),all.x=T)
dr1_NeoAG_mhci_RelDiag$Rel_Mut=str_count(str_to_upper(dr1_NeoAG_mhci_RelDiag$rel_bases), as.character(dr1_NeoAG_mhci_RelDiag$mut_nt))
dr1_NeoAG_mhci_RelDiag$Rel_Wt=str_count(str_to_upper(dr1_NeoAG_mhci_RelDiag$rel_bases), as.character(dr1_NeoAG_mhci_RelDiag$wt_nt))
dr1_NeoAG_mhci_RelDiag$Dx_Mut=str_count(str_to_upper(dr1_NeoAG_mhci_RelDiag$dx_bases), as.character(dr1_NeoAG_mhci_RelDiag$mut_nt))
dr1_NeoAG_mhci_RelDiag$Dx_Wt=str_count(str_to_upper(dr1_NeoAG_mhci_RelDiag$dx_bases), as.character(dr1_NeoAG_mhci_RelDiag$wt_nt))
write.table(unique(dr1_NeoAG_mhci_RelDiag[,c(1:24,33:36)]),"../../Documents/Progetto_LuciaGabri/Summary_DR1_Annotato_NeoAGs_MHCI_RelDiag_NoMismatches_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)


#DR1 RelOnly
cat StrongBinders_Annotati_DR1_RelOnly WeakBinders_Annotati_DR1_RelOnly | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ > Summary_DR1_Annotato_RelOnly_NoHeaders_NoMismatches.txt



awk '{OFS="\t"; print $11":"$12"-"$12}'   Summary_DR1_Annotato_RelOnly_NoHeaders_NoMismatches.txt | sort | uniq  > Summary_DR1_Annotato_NeoEpitopi_RelOnly_MHCI_Regions.txt




DR1_1 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DR1_1_Allfiles_sorted.bam
DR1_2 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DR1_2_Allfiles_sorted.bam


 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DR1_1_Allfiles_sorted.bam >> DR1_NeoAGs_I_RelOnly_Binders_Regions_Relapse.txt; done < Summary_DR1_Annotato_NeoEpitopi_RelOnly_MHCI_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DR1_2_Allfiles_sorted.bam >> DR1_NeoAGs_I_RelOnly_Binders_Regions_Relapse.txt; done < Summary_DR1_Annotato_NeoEpitopi_RelOnly_MHCI_Regions.txt
 


dr1_NeoAG_mhci_RelOnly=read.table('Summary_DR1_Annotato_RelOnly_NoHeaders_NoMismatches.txt')
colnames(dr1_NeoAG_mhci_RelOnly)=c("Gene_ID","Peptide","HLA","HLA","HLA","HLA","HLA","HLA","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA","IC50_HLA","IC,50_HLA","IC50_HLA","IC50_HLA","IC,50_HLA")

dr1_NeoAG_mhci_RelOnly=merge(dr1_NeoAG_mhci_RelOnly,DR1_RNK,by="genes",all.x=T)
relapse=read.table('DR1_NeoAGs_I_RelOnly_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DR1_NeoAGs_I_RelOnly_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
dr1_NeoAG_mhci_RelOnly=merge(dr1_NeoAG_mhci_RelOnly,relapse,by=c("chr","pos"),all.x=T)
dr1_NeoAG_mhci_RelOnly=merge(dr1_NeoAG_mhci_RelOnly,diagnosis,by=c("chr","pos"),all.x=T)
dr1_NeoAG_mhci_RelOnly$Rel_Mut=str_count(str_to_upper(dr1_NeoAG_mhci_RelOnly$rel_bases), as.character(dr1_NeoAG_mhci_RelOnly$mut_nt))
dr1_NeoAG_mhci_RelOnly$Rel_Wt=str_count(str_to_upper(dr1_NeoAG_mhci_RelOnly$rel_bases), as.character(dr1_NeoAG_mhci_RelOnly$wt_nt))
dr1_NeoAG_mhci_RelOnly$Dx_Mut=str_count(str_to_upper(dr1_NeoAG_mhci_RelOnly$dx_bases), as.character(dr1_NeoAG_mhci_RelOnly$mut_nt))
dr1_NeoAG_mhci_RelOnly$Dx_Wt=str_count(str_to_upper(dr1_NeoAG_mhci_RelOnly$dx_bases), as.character(dr1_NeoAG_mhci_RelOnly$wt_nt))
write.table(unique(dr1_NeoAG_mhci_RelOnly[,c(1:24,33:36)]),"../../Documents/Progetto_LuciaGabri/Summary_DR1_Annotato_NeoAGs_MHCI_RelOnly_NoMismatches_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)

##DR1
#mhag=read.table('Summary_DR1_Annotato_NeoAGs_MHCI_RelDiag_AllInfos.txt',head=T)
#mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
#mhag$logFC[is.na(mhag$logFC)] <- 0
#mhag$Expr_Relapse[is.na(mhag$Expr_Relapse)] <- 0
#mhag$Expr_Diag[is.na(mhag$Expr_Diag)] <- 0
##Weak
# nrow(mhag[apply(mhag[,16:20], MARGIN = 1, function(x) any(x<=500 & x> 50)), ])
#[1] 7
##Strong 
# nrow(mhag[apply(mhag[,16:20], MARGIN = 1, function(x) any(x<=50 & x> 0)), ])
#
#mhag=mhag[(mhag$ALFE2_kallisto>=0.58 & mhag$ALFE1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
#nrow(mhag)
#8
#
#table(mhag$Rel_Mut>=3)
#
#FALSE  TRUE
#    1     7
#
#table(mhag$Rel_Mut>mhag$Dx_Mut)
#
#FALSE
#    8
#
#
##Weak
# nrow(mhag[apply(mhag[,16:20], MARGIN = 1, function(x) any(x<=500 & x> 50)), ])
#[1] 5
##Strong 
# nrow(mhag[apply(mhag[,16:20], MARGIN = 1, function(x) any(x<=50 & x> 0)), ])
# nrow(mhag[apply(mhag[,16:20], MARGIN = 1, function(x) any(x<=50 & x> 0)), ])
#[1] 3
#
##A livello di trascritto
#table(mhag$ALFE2_kallisto > mhag$ALFE1_kallisto)
#
#FALSE  TRUE
#    5     3
#
#table(mhag$ALFE2_kallisto < mhag$ALFE1_kallisto)
#
#FALSE  TRUE
#    3     5
#
#
#wilcox.test(mhag$ALFE2_kallisto, mhag$ALFE1_kallisto,alternative="greater")
#
#        Wilcoxon rank sum test with continuity correction
#
#data:  mhag$ALFE2_kallisto and mhag$ALFE1_kallisto
#W = 27, p-value = 0.7194
#alternative hypothesis: true location shift is greater than 0
#
#Warning message:
#In wilcox.test.default(mhag$ALFE2_kallisto, mhag$ALFE1_kallisto,  :
#  cannot compute exact p-value with ties
#
#
##A livello di trascritto
# table(mhag$Rel_Mut>mhag$Dx_Mut)
#
#FALSE
#    8
#
#table(mhag$Rel_Mut<mhag$Dx_Mut)
#
#TRUE
#   8
#
#wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")
#
#        Wilcoxon rank sum test with continuity correction
#
#data:  mhag$Rel_Mut and mhag$Dx_Mut
#W = 7, p-value = 0.9966
#alternative hypothesis: true location shift is greater than 0
#
#Warning message:
#In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
#  cannot compute exact p-value with ties




#DR4



sh Launcher_snp2epi_Paper_Annotation.sh DR41_DR42_Somatic.Anno.tsv DR41_DR42_MHCI_NoMismatchHla HLA-A29:02,HLA-A31:01,HLA-B07:02,HLA-C07:02,HLA-C03:04 DR4_Expression_Kallisto_NeoEpitopes.txt

       #Summary Binders Diagnosis Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A29:02\tHLA-A31:01\tHLA-B07:02\tHLA-C07:02\tHLA-C03:04" )
        more +2 ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_weak_binders_DiagOnly
#         ../../Documents/Progetto_LuciaGabri/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt ../../Documents/Progetto_LuciaGabri/weak_binders_DiagOnly ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_weak_binders_annotated_DiagOnly
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A29:02\tHLA-A31:01\tHLA-B07:02\tHLA-C07:02\tHLA-C03:04" )
        more +2 ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls |   awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_strong_binders_DiagOnly
 #        ../../Documents/Progetto_LuciaGabri/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt ../../Documents/Progetto_LuciaGabri/strong_binders_DiagOnly ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls ../../Documents/Progetto_LuciaGabri/strong_binders_annotated_DiagOnly

        #cat ../../Documents/Progetto_LuciaGabri/strong_binders_annotated_DiagOnly ../../Documents/Progetto_LuciaGabri/weak_binders_annotated_DiagOnly > ./Summary_Binders_RESULTS_DR41_DR42_MHCI_NoMismatchHla__PeptidesPrediction_DiagnosisOnly.txt



        #Summary Binders Relapse Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A29:02\tHLA-A31:01\tHLA-B07:02\tHLA-C07:02\tHLA-C03:04" )
        more +2 ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1}' -| (echo -e $weak_binders; cat -) > ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_weak_binders_RelOnly
        # ../../Documents/Progetto_LuciaGabri/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt ../../Documents/Progetto_LuciaGabri/weak_binders_RelOnly ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls ../../Documents/Progetto_LuciaGabri/weak_binders_annotated_RelOnly

        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A29:02\tHLA-A31:01\tHLA-B07:02\tHLA-C07:02\tHLA-C03:04" )
        more +2 ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_strong_binders_RelOnly
         #../../Documents/Progetto_LuciaGabri/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt ../../Documents/Progetto_LuciaGabri/strong_binders_RelOnly ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls ../../Documents/Progetto_LuciaGabri/strong_binders_annotated_RelOnly
        
        #cat ../../Documents/Progetto_LuciaGabri/strong_binders_annotated_RelOnly ../../Documents/Progetto_LuciaGabri/weak_binders_annotated_RelOnly > ./Summary_Binders_RESULTS_DR41_DR42_MHCI_NoMismatchHla__PeptidesPrediction_RelapseOnly.txt
        
        
        #Summary Binders Common Diagnosis Relapse
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A29:02\tHLA-A31:01\tHLA-B07:02\tHLA-C07:02\tHLA-C03:04" )
        more +2 ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_weak_binders_RelDiag
  #       ../../Documents/Progetto_LuciaGabri/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt ../../Documents/Progetto_LuciaGabri/weak_binders_RelDiag ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls ../../Documents/Progetto_LuciaGabri/weak_binders_annotated_RelDiag
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A29:02\tHLA-A31:01\tHLA-B07:02\tHLA-C07:02\tHLA-C03:04" )
        more +2 ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_strong_binders_RelDiag
   #      ../../Documents/Progetto_LuciaGabri/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt ../../Documents/Progetto_LuciaGabri/strong_binders_RelDiag ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls ../../Documents/Progetto_LuciaGabri/strong_binders_annotated_RelDiag
       
        #cat ../../Documents/Progetto_LuciaGabri/strong_binders_annotated_RelDiag ../../Documents/Progetto_LuciaGabri/weak_binders_annotated_RelDiag > ./Summary_Binders_RESULTS_DR41_DR42_MHCI_NoMismatchHla__PeptidesPrediction_CommonRelapseDiagnosisOnly.txt
        


#DR4 RelDiag
cat StrongBinders_Annotati_DR4_RelDiag WeakBinders_Annotati_DR4_RelDiag | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ > Summary_DR4_Annotato_RelDiag_NoHeaders_NoMismatches.txt


DR41 =  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR41_Allfiles_sorted.bam
DR42 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR42_Allfiles_sorted.bam

 
awk '{OFS="\t"; print $10":"$11"-"$11}'   Summary_DR4_Annotato_RelDiag_NoHeaders_NoMismatches.txt | sort | uniq  > Summary_DR4_Annotato_NeoEpitopi_RelDiag_MHCI_Regions.txt

while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR41_Allfiles_sorted.bam >> DR4_NeoAGs_I_RelDiag_Binders_Regions_Diagnosis.txt; done < Summary_DR4_Annotato_NeoEpitopi_RelDiag_MHCI_Regions.txt
while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR42_Allfiles_sorted.bam >> DR4_NeoAGs_I_RelDiag_Binders_Regions_Relapse.txt; done < Summary_DR4_Annotato_NeoEpitopi_RelDiag_MHCI_Regions.txt


dr4_NeoAG_mhci_RelDiag=read.table('Summary_DR4_Annotato_RelDiag_NoHeaders_NoMismatches.txt')
colnames(dr4_NeoAG_mhci_RelDiag)=c("Gene_ID","Peptide","HLA","HLA","HLA","HLA","HLA","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA","IC50_HLA","IC,50_HLA","IC50_HLA","IC50_HLA")

dr4_NeoAG_mhci_RelDiag=merge(dr4_NeoAG_mhci_RelDiag,DR4_RNK,by="genes",all.x=T)
relapse=read.table('DR4_NeoAGs_I_RelDiag_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DR4_NeoAGs_I_RelDiag_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
dr4_NeoAG_mhci_RelDiag=merge(dr4_NeoAG_mhci_RelDiag,relapse,by=c("chr","pos"),all.x=T)
dr4_NeoAG_mhci_RelDiag=merge(dr4_NeoAG_mhci_RelDiag,diagnosis,by=c("chr","pos"),all.x=T)
dr4_NeoAG_mhci_RelDiag$Rel_Mut=str_count(str_to_upper(dr4_NeoAG_mhci_RelDiag$rel_bases), as.character(dr4_NeoAG_mhci_RelDiag$mut_nt))
dr4_NeoAG_mhci_RelDiag$Rel_Wt=str_count(str_to_upper(dr4_NeoAG_mhci_RelDiag$rel_bases), as.character(dr4_NeoAG_mhci_RelDiag$wt_nt))
dr4_NeoAG_mhci_RelDiag$Dx_Mut=str_count(str_to_upper(dr4_NeoAG_mhci_RelDiag$dx_bases), as.character(dr4_NeoAG_mhci_RelDiag$mut_nt))
dr4_NeoAG_mhci_RelDiag$Dx_Wt=str_count(str_to_upper(dr4_NeoAG_mhci_RelDiag$dx_bases), as.character(dr4_NeoAG_mhci_RelDiag$wt_nt))
 write.table(unique(dr4_NeoAG_mhci_RelDiag[,c(1:22,31:34)]),"../../Documents/Progetto_LuciaGabri/Summary_DR4_Annotato_NeoAGs_MHCI_RelDiag_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)



#DR4 RelOnly
cat StrongBinders_Annotati_DR4_RelOnly WeakBinders_Annotati_DR4_RelOnly | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ > Summary_DR4_Annotato_RelOnly_NoHeaders_NoMismatches.txt


DR41 =  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR41_Allfiles_sorted.bam
DR42 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR42_Allfiles_sorted.bam

 
awk '{OFS="\t"; print $10":"$11"-"$11}'   Summary_DR4_Annotato_RelOnly_NoHeaders_NoMismatches.txt | sort | uniq  > Summary_DR4_Annotato_NeoEpitopi_RelOnly_MHCI_Regions.txt

while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR41_Allfiles_sorted.bam >> DR4_NeoAGs_I_RelOnly_Binders_Regions_Diagnosis.txt; done < Summary_DR4_Annotato_NeoEpitopi_RelOnly_MHCI_Regions.txt
while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR42_Allfiles_sorted.bam >> DR4_NeoAGs_I_RelOnly_Binders_Regions_Relapse.txt; done < Summary_DR4_Annotato_NeoEpitopi_RelOnly_MHCI_Regions.txt


dr4_NeoAG_mhci_RelOnly=read.table('Summary_DR4_Annotato_RelOnly_NoHeaders_NoMismatches.txt')
colnames(dr4_NeoAG_mhci_RelOnly)=c("Gene_ID","Peptide","HLA","HLA","HLA","HLA","HLA","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA","IC50_HLA","IC,50_HLA","IC50_HLA","IC50_HLA")

dr4_NeoAG_mhci_RelOnly=merge(dr4_NeoAG_mhci_RelOnly,DR4_RNK,by="genes",all.x=T)
relapse=read.table('DR4_NeoAGs_I_RelOnly_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DR4_NeoAGs_I_RelOnly_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
dr4_NeoAG_mhci_RelOnly=merge(dr4_NeoAG_mhci_RelOnly,relapse,by=c("chr","pos"),all.x=T)
dr4_NeoAG_mhci_RelOnly=merge(dr4_NeoAG_mhci_RelOnly,diagnosis,by=c("chr","pos"),all.x=T)
dr4_NeoAG_mhci_RelOnly$Rel_Mut=str_count(str_to_upper(dr4_NeoAG_mhci_RelOnly$rel_bases), as.character(dr4_NeoAG_mhci_RelOnly$mut_nt))
dr4_NeoAG_mhci_RelOnly$Rel_Wt=str_count(str_to_upper(dr4_NeoAG_mhci_RelOnly$rel_bases), as.character(dr4_NeoAG_mhci_RelOnly$wt_nt))
dr4_NeoAG_mhci_RelOnly$Dx_Mut=str_count(str_to_upper(dr4_NeoAG_mhci_RelOnly$dx_bases), as.character(dr4_NeoAG_mhci_RelOnly$mut_nt))
dr4_NeoAG_mhci_RelOnly$Dx_Wt=str_count(str_to_upper(dr4_NeoAG_mhci_RelOnly$dx_bases), as.character(dr4_NeoAG_mhci_RelOnly$wt_nt))
write.table(unique(dr4_NeoAG_mhci_RelOnly[,c(1:22,31:34)]),"../../Documents/Progetto_LuciaGabri/Summary_DR4_Annotato_NeoAGs_MHCI_RelOnly_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)


#
#mhag=read.table('Summary_DR4_Annotato_NeoAGs_MHCI_RelDiag_AllInfos.txt',head=T)
#mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
#mhag$logFC[is.na(mhag$logFC)] <- 0
#mhag$DR42_kallisto[is.na(mhag$DR42_kallisto)] <- 0
#mhag$DR41_kallisto[is.na(mhag$DR41_kallisto)] <- 0
#
#   nrow(mhag[apply(mhag[,16:20], MARGIN = 1, function(x) any(x<=500 & x> 50)), ])
#[1] 6
#
#  nrow(mhag[apply(mhag[,16:20], MARGIN = 1, function(x) any(x<=50 & x> 0)), ])
#[1] 1
#  
#  
#  
#mhag=mhag[(mhag$DR42_kallisto>=0.58 & mhag$DR41_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
#nrow(mhag)
#4
#
#table(mhag$Rel_Mut>=3)
#
#
#   nrow(mhag[apply(mhag[,16:20], MARGIN = 1, function(x) any(x<=500 & x> 50)), ])
#[1] 4
#
#  nrow(mhag[apply(mhag[,16:20], MARGIN = 1, function(x) any(x<=50 & x> 0)), ])
#[1] 0
#  
##A livello di trascritto
# table(mhag$DR42_kallisto > mhag$DR41_kallisto)
#
#FALSE  TRUE
#    3     1
#>  wilcox.test(mhag$DR42_kallisto, mhag$DR41_kallisto,alternative="greater")
#
#        Wilcoxon rank sum test with continuity correction
#
#data:  mhag$DR42_kallisto and mhag$DR41_kallisto
#W = 6, p-value = 0.7674
#alternative hypothesis: true location shift is greater than 0
#
#Warning message:
#In wilcox.test.default(mhag$DR42_kallisto, mhag$DR41_kallisto, alternative = "greater") :
#  cannot compute exact p-value with ties
#
#table(mhag$Rel_Mut>mhag$Dx_Mut)
#
#FALSE  TRUE
#    2     2
#
#wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")
#
#        Wilcoxon rank sum test with continuity correction
#
#data:  mhag$Rel_Mut and mhag$Dx_Mut
#W = 7, p-value = 0.6724
#alternative hypothesis: true location shift is greater than 0
#
#Warning message:
#In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
#  cannot compute exact p-value with ties
#



#DR5

#DR5


#DR5 RelDiag
cat StrongBinders_Annotati_DR5_RelDiag WeakBinders_Annotati_DR5_RelDiag | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ > Summary_DR5_Annotato_RelDiag_MCHI_NoHeaders_NoMismatches.txt



DR51 =  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR51_Allfiles_sorted.bam
DR52 =  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR52_Allfiles_sorted.bam
 
awk '{OFS="\t"; print $11":"$12"-"$12}'   Summary_DR5_Annotato_RelDiag_MCHI_NoHeaders_NoMismatches.txt | sort | uniq  > Summary_DR5_Annotato_RelDiag_MCHI_NoHeaders_NoMismatches_Regions.txt

while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR51_Allfiles_sorted.bam >> DR5_MHCI__Binders_Regions_Diagnosis.txt; done < Summary_DR5_Annotato_RelDiag_MCHI_NoHeaders_NoMismatches_Regions.txt
while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR52_Allfiles_sorted.bam >> DR5_MHCI__Binders_Regions_Relapse.txt; done < Summary_DR5_Annotato_RelDiag_MCHI_NoHeaders_NoMismatches_Regions.txt
 

dr5_NeoAG_mhci_RelDiag=read.table('Summary_DR5_Annotato_RelDiag_MCHI_NoHeaders_NoMismatches.txt')
colnames(dr5_NeoAG_mhci_RelDiag)=c("Gene_ID","Peptide","HLA","HLA","HLA","HLA","HLA","HLA","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA","IC50_HLA","IC,50_HLA","IC50_HLA","IC50_HLA","IC50_HLA")

dr5_NeoAG_mhci_RelDiag=merge(dr5_NeoAG_mhci_RelDiag,DR5_RNK,by="genes",all.x=T)
relapse=read.table('DR5_MHCI__Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DR5_MHCI__Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
dr5_NeoAG_mhci_RelDiag=merge(dr5_NeoAG_mhci_RelDiag,relapse,by=c("chr","pos"),all.x=T)
dr5_NeoAG_mhci_RelDiag=merge(dr5_NeoAG_mhci_RelDiag,diagnosis,by=c("chr","pos"),all.x=T)
dr5_NeoAG_mhci_RelDiag$Rel_Mut=str_count(str_to_upper(dr5_NeoAG_mhci_RelDiag$rel_bases), as.character(dr5_NeoAG_mhci_RelDiag$mut_nt))
dr5_NeoAG_mhci_RelDiag$Rel_Wt=str_count(str_to_upper(dr5_NeoAG_mhci_RelDiag$rel_bases), as.character(dr5_NeoAG_mhci_RelDiag$wt_nt))
dr5_NeoAG_mhci_RelDiag$Dx_Mut=str_count(str_to_upper(dr5_NeoAG_mhci_RelDiag$dx_bases), as.character(dr5_NeoAG_mhci_RelDiag$mut_nt))
dr5_NeoAG_mhci_RelDiag$Dx_Wt=str_count(str_to_upper(dr5_NeoAG_mhci_RelDiag$dx_bases), as.character(dr5_NeoAG_mhci_RelDiag$wt_nt))
write.table(unique(dr5_NeoAG_mhci_RelDiag[,c(1:24,33:36)]),"../../Documents/Progetto_LuciaGabri/Summary_DR5_Annotato_NeoAGs_MHCI_RelDiag_NoMismatches_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)


#FOCA
#RelDiag
cat StrongBinders_Annotati_FOCA_RelDiag WeakBinders_Annotati_FOCA_RelDiag | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ > Summary_FOCA_Annotato_RelDiag_NoHeaders_NoMismatches.txt



FOCA1 =  /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/FOCA1_Allfiles_sorted.bam
FOCA2 =  /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/FOCA2_Allfiles_sorted.bam

awk '{OFS="\t"; print $11":"$12"-"$12}' Summary_FOCA_Annotato_RelDiag_NoHeaders_NoMismatches.txt | sort | uniq  > Summary_FOCA_Annotato_RelDiag_NoHeaders_NoMismatches_Regions.txt
while read region;  do samtools mpileup -r $region    /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/FOCA1_Allfiles_sorted.bam >> FOCA_NeoAGs_I_RelDiag_Binders_Regions_Diagnosis.txt; done < Summary_FOCA_Annotato_RelDiag_NoHeaders_NoMismatches_Regions.txt
while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/FOCA2_Allfiles_sorted.bam >> FOCA_NeoAGs_I_RelDiag_Binders_Regions_Relapse.txt; done < Summary_FOCA_Annotato_RelDiag_NoHeaders_NoMismatches_Regions.txt



foca_NeoAG_mhci_RelDiag=read.table('Summary_FOCA_Annotato_MHCI_RelDiag_NoHeaders_NoMismatches.txt')
colnames(foca_NeoAG_mhci_RelDiag)=c("Gene_ID","Peptide","HLA","HLA","HLA","HLA","HLA","HLA","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA","IC50_HLA","IC,50_HLA","IC50_HLA","IC50_HLA","IC50_HLA")

foca_NeoAG_mhci_RelDiag=merge(foca_NeoAG_mhci_RelDiag,FOCA_RNK,by="genes",all.x=T)
relapse=read.table('FOCA_NeoAGs_I_RelDiag_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('FOCA_NeoAGs_I_RelDiag_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
foca_NeoAG_mhci_RelDiag=merge(foca_NeoAG_mhci_RelDiag,relapse,by=c("chr","pos"),all.x=T)
foca_NeoAG_mhci_RelDiag=merge(foca_NeoAG_mhci_RelDiag,diagnosis,by=c("chr","pos"),all.x=T)
foca_NeoAG_mhci_RelDiag$Rel_Mut=str_count(str_to_upper(foca_NeoAG_mhci_RelDiag$rel_bases), as.character(foca_NeoAG_mhci_RelDiag$mut_nt))
foca_NeoAG_mhci_RelDiag$Rel_Wt=str_count(str_to_upper(foca_NeoAG_mhci_RelDiag$rel_bases), as.character(foca_NeoAG_mhci_RelDiag$wt_nt))
foca_NeoAG_mhci_RelDiag$Dx_Mut=str_count(str_to_upper(foca_NeoAG_mhci_RelDiag$dx_bases), as.character(foca_NeoAG_mhci_RelDiag$mut_nt))
foca_NeoAG_mhci_RelDiag$Dx_Wt=str_count(str_to_upper(foca_NeoAG_mhci_RelDiag$dx_bases), as.character(foca_NeoAG_mhci_RelDiag$wt_nt))
 write.table(unique(foca_NeoAG_mhci_RelDiag[,c(1:24,33:36)]),"../../Documents/Progetto_LuciaGabri/Summary_FOCA_Annotato_NeoAGs_MHCI_RelDiag_NoMismatches_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)

#mhag=read.table('Summary_FOCA_Annotato_NeoAGs_MHCI_RelDiag_AllInfos.txt',head=T)
#mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
#mhag$logFC[is.na(mhag$logFC)] <- 0
#mhag$FOCA2_kallisto[is.na(mhag$FOCA2_kallisto)] <- 0
#mhag$FOCA1_kallisto[is.na(mhag$FOCA1_kallisto)] <- 0
#
#   nrow(mhag[apply(mhag[,17:22], MARGIN = 1, function(x) any(x<=500 & x> 50)), ])
#[1] 4
#
#  nrow(mhag[apply(mhag[,17:22], MARGIN = 1, function(x) any(x<=50 & x> 0)), ])
#[1] 0
#  
#mhag=mhag[(mhag$FOCA2_kallisto>=0.58 & mhag$FOCA1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
#nrow(mhag)
#3
#
#table(mhag$Rel_Mut>=3)
#
#   nrow(mhag[apply(mhag[,17:22], MARGIN = 1, function(x) any(x<=500 & x> 50)), ])
#[1] 4
#
#  nrow(mhag[apply(mhag[,17:22], MARGIN = 1, function(x) any(x<=50 & x> 0)), ])
#[1] 0
#  
#  
##A livello di trascritto
# table(mhag$FOCA2_kallisto > mhag$FOCA1_kallisto)
#
#TRUE  
#    3  
#wilcox.test(mhag$FOCA2_kallisto, mhag$FOCA1_kallisto,alternative="greater")
#
#        Wilcoxon rank sum test with continuity correction
#
#data:  mhag$FOCA2_kallisto and mhag$FOCA1_kallisto
#W = 9, p-value = 0.0361
#alternative hypothesis: true location shift is greater than 0
#
#Warning message:
#In wilcox.test.default(mhag$FOCA2_kallisto, mhag$FOCA1_kallisto,  :
#  cannot compute exact p-value with ties
#
#
#table(mhag$Rel_Mut>mhag$Dx_Mut)
#
#FALSE
#   3
#
#
#wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")
#
#        Wilcoxon rank sum test with continuity correction
#
#data:  mhag$Rel_Mut and mhag$Dx_Mut
#W = 4, p-value = 0.09697
#alternative hypothesis: true location shift is greater than 0
#
#Warning message:
#In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
#  cannot compute exact p-value with ties
#
#


#GAGRA
#Reldiag
cat StrongBinders_Annotati_GAGRA_RelDiag WeakBinders_Annotati_GAGRA_RelDiag | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ > Summary_GAGRA_Annotato_RelDiag_MHCI_NoHeaders_NoMismatches.txt


GAGRA1 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/GAGRA1_Allfiles_sorted.bam
GAGRA2 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/GAGRA2_Allfiles_sorted.bam

awk '{OFS="\t"; print $11":"$12"-"$12}' Summary_GAGRA_Annotato_RelDiag_MHCI_NoHeaders_NoMismatches.txt | sort | uniq  > Summary_GAGRA_Annotato_RelDiag_MHCI_NoHeaders_NoMismatches_Regions.txt
while read region;  do samtools mpileup -r $region    /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/GAGRA1_Allfiles_sorted.bam >> GAGRA_NeoAGs_I_RelDiag_Binders_Regions_Diagnosis.txt; done < Summary_GAGRA_Annotato_RelDiag_MHCI_NoHeaders_NoMismatches_Regions.txt
while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/GAGRA2_Allfiles_sorted.bam >> GAGRA_NeoAGs_I_RelDiag_Binders_Regions_Relapse.txt; done < Summary_GAGRA_Annotato_RelDiag_MHCI_NoHeaders_NoMismatches_Regions.txt



gagra_NeoAG_mhci_RelDiag=read.table('Summary_GAGRA_Annotato_RelDiag_MHCI_NoHeaders_NoMismatches.txt')
 colnames(gagra_NeoAG_mhci_RelDiag)=c("Gene_ID","Peptide","HLA","HLA","HLA","HLA","HLA","HLA","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA","IC50_HLA","IC,50_HLA","IC50_HLA","IC50_HLA","IC50_HLA")

gagra_NeoAG_mhci_RelDiag=merge(gagra_NeoAG_mhci_RelDiag,GAGRA_RNK,by="genes",all.x=T)
relapse=read.table('GAGRA_NeoAGs_I_RelDiag_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('GAGRA_NeoAGs_I_RelDiag_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
gagra_NeoAG_mhci_RelDiag=merge(gagra_NeoAG_mhci_RelDiag,relapse,by=c("chr","pos"),all.x=T)
gagra_NeoAG_mhci_RelDiag=merge(gagra_NeoAG_mhci_RelDiag,diagnosis,by=c("chr","pos"),all.x=T)
gagra_NeoAG_mhci_RelDiag$Rel_Mut=str_count(str_to_upper(gagra_NeoAG_mhci_RelDiag$rel_bases), as.character(gagra_NeoAG_mhci_RelDiag$mut_nt))
gagra_NeoAG_mhci_RelDiag$Rel_Wt=str_count(str_to_upper(gagra_NeoAG_mhci_RelDiag$rel_bases), as.character(gagra_NeoAG_mhci_RelDiag$wt_nt))
gagra_NeoAG_mhci_RelDiag$Dx_Mut=str_count(str_to_upper(gagra_NeoAG_mhci_RelDiag$dx_bases), as.character(gagra_NeoAG_mhci_RelDiag$mut_nt))
gagra_NeoAG_mhci_RelDiag$Dx_Wt=str_count(str_to_upper(gagra_NeoAG_mhci_RelDiag$dx_bases), as.character(gagra_NeoAG_mhci_RelDiag$wt_nt))
write.table(unique(gagra_NeoAG_mhci_RelDiag[,c(1:24,33:36)]),"../../Documents/Progetto_LuciaGabri/Summary_GAGRA_Annotato_NeoAGs_MHCI_RelDiag_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)

#GAGRA non ha NeoAGs RelOnly
#GAGRA non ha NeaoAGs DiagOnly



#mhag=read.table('Summary_GAGRA_Annotato_NeoAGs_MHCI_RelDiag_AllInfos.txt',head=T)
#mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
#mhag$logFC[is.na(mhag$logFC)] <- 0
#mhag$GAGRA2_kallisto[is.na(mhag$GAGRA2_kallisto)] <- 0
#mhag$GAGRA1_kallisto[is.na(mhag$GAGRA1_kallisto)] <- 0
#
#
#   nrow(mhag[apply(mhag[,17:22], MARGIN = 1, function(x) any(x<=500 & x> 50)), ])
#[1] 4
#
#  nrow(mhag[apply(mhag[,17:22], MARGIN = 1, function(x) any(x<=50 & x> 0)), ])
#[1] 0
#  
#  
#  
#mhag=mhag[(mhag$GAGRA2_kallisto>=0.58 & mhag$GAGRA1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
#nrow(mhag)
#4
#
# table(mhag$Dx_Mut>=3)
#
#
#
##A livello di trascritto
# table(mhag$GAGRA2_kallisto > mhag$GAGRA1_kallisto)
#
#FALSE  
#    4
#
#wilcox.test(mhag$GAGRA2_kallisto, mhag$GAGRA1_kallisto,alternative="greater")
#
#        Wilcoxon rank sum test with continuity correction
#
#data:  mhag$GAGRA2_kallisto and mhag$GAGRA1_kallisto
#W = 5, p-value = 0.8468
#alternative hypothesis: true location shift is greater than 0
#
#Warning message:
#In wilcox.test.default(mhag$GAGRA2_kallisto, mhag$GAGRA1_kallisto,  :
#  cannot compute exact p-value with ties
#
#
#
# table(mhag$Rel_Mut>mhag$Dx_Mut)
#
#FALSE  TRUE
#    1     3
# table(mhag$Rel_Mut<mhag$Dx_Mut)
#
#FALSE  TRUE
#    3     1
#
#
#wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")
#
#        Wilcoxon rank sum test with continuity correction
#
#data:  mhag$Rel_Mut and mhag$Dx_Mut
#W = 10, p-value = 0.3306
#alternative hypothesis: true location shift is greater than 0
#
#Warning message:
#In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
#  cannot compute exact p-value with ties


#LUAN
#LUAN RelDiag
cat StrongBinders_Annotati_LUAN_RelDiag WeakBinders_Annotati_LUAN_RelDiag | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ > Summary_LUAN_Annotato_RelDiag_MHCI_NoHeaders_NoMismatches.txt

LUAN1 = /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN1_Allfiles_sorted.bam
LUAN2 = /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN2_Allfiles_sorted.bam

awk '{OFS="\t"; print $10":"$11"-"$11}' Summary_LUAN_Annotato_RelDiag_MHCI_NoHeaders_NoMismatches.txt | sort | uniq  >  Summary_LUAN_Annotato_RelDiag_MHCI_NoHeaders_NoMismatches_Regions.txt
 while read region;  do samtools mpileup -r $region    /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN1_Allfiles_sorted.bam >> LUAN_NeoAGs_I_RelDiag_Binders_Regions_Diagnosis.txt; done < Summary_LUAN_Annotato_RelDiag_MHCI_NoHeaders_NoMismatches_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN2_Allfiles_sorted.bam >> LUAN_NeoAGs_I_RelDiag_Binders_Regions_Relapse.txt; done < Summary_LUAN_Annotato_RelDiag_MHCI_NoHeaders_NoMismatches_Regions.txt


luan_NeoAG_mhci_RelDiag=read.table('Summary_LUAN_Annotato_RelDiag_MHCI_NoHeaders_NoMismatches.txt')
colnames(luan_NeoAG_mhci_RelDiag)=c("Gene_ID","Peptide","HLA","HLA","HLA","HLA","HLA","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA","IC50_HLA","IC,50_HLA","IC50_HLA","IC50_HLA")

luan_NeoAG_mhci_RelDiag=merge(luan_NeoAG_mhci_RelDiag,LUAN_RNK,by="genes",all.x=T)
relapse=read.table('LUAN_NeoAGs_I_RelDiag_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('LUAN_NeoAGs_I_RelDiag_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
luan_NeoAG_mhci_RelDiag=merge(luan_NeoAG_mhci_RelDiag,relapse,by=c("chr","pos"),all.x=T)
luan_NeoAG_mhci_RelDiag=merge(luan_NeoAG_mhci_RelDiag,diagnosis,by=c("chr","pos"),all.x=T)
luan_NeoAG_mhci_RelDiag$Rel_Mut=str_count(str_to_upper(luan_NeoAG_mhci_RelDiag$rel_bases), as.character(luan_NeoAG_mhci_RelDiag$mut_nt))
luan_NeoAG_mhci_RelDiag$Rel_Wt=str_count(str_to_upper(luan_NeoAG_mhci_RelDiag$rel_bases), as.character(luan_NeoAG_mhci_RelDiag$wt_nt))
luan_NeoAG_mhci_RelDiag$Dx_Mut=str_count(str_to_upper(luan_NeoAG_mhci_RelDiag$dx_bases), as.character(luan_NeoAG_mhci_RelDiag$mut_nt))
luan_NeoAG_mhci_RelDiag$Dx_Wt=str_count(str_to_upper(luan_NeoAG_mhci_RelDiag$dx_bases), as.character(luan_NeoAG_mhci_RelDiag$wt_nt))
write.table(unique(luan_NeoAG_mhci_RelDiag[,c(1:22,31:34)]),"../../Documents/Progetto_LuciaGabri/Summary_LUAN_Annotato_NeoAGs_MHCI_RelDiag_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)

#LUAN RelOnly
cat StrongBinders_Annotati_LUAN_RelOnly WeakBinders_Annotati_LUAN_RelOnly | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ > Summary_LUAN_Annotato_RelOnly_MHCI_NoHeaders_NoMismatches.txt

LUAN1 = /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN1_Allfiles_sorted.bam
LUAN2 = /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN2_Allfiles_sorted.bam

awk '{OFS="\t"; print $10":"$11"-"$11}' Summary_LUAN_Annotato_RelOnly_MHCI_NoHeaders_NoMismatches.txt | sort | uniq  >  Summary_LUAN_Annotato_RelOnly_MHCI_NoHeaders_NoMismatches_Regions.txt
 while read region;  do samtools mpileup -r $region    /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN1_Allfiles_sorted.bam >> LUAN_NeoAGs_I_RelOnly_Binders_Regions_Diagnosis.txt; done < Summary_LUAN_Annotato_RelOnly_MHCI_NoHeaders_NoMismatches_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN2_Allfiles_sorted.bam >> LUAN_NeoAGs_I_RelOnly_Binders_Regions_Relapse.txt; done < Summary_LUAN_Annotato_RelOnly_MHCI_NoHeaders_NoMismatches_Regions.txt


luan_NeoAG_mhci_RelOnly=read.table('Summary_LUAN_Annotato_RelOnly_MHCI_NoHeaders_NoMismatches.txt')
colnames(luan_NeoAG_mhci_RelOnly)=c("Gene_ID","Peptide","HLA","HLA","HLA","HLA","HLA","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA","IC50_HLA","IC,50_HLA","IC50_HLA","IC50_HLA")

luan_NeoAG_mhci_RelOnly=merge(luan_NeoAG_mhci_RelOnly,LUAN_RNK,by="genes",all.x=T)
relapse=read.table('LUAN_NeoAGs_I_RelOnly_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('LUAN_NeoAGs_I_RelOnly_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
luan_NeoAG_mhci_RelOnly=merge(luan_NeoAG_mhci_RelOnly,relapse,by=c("chr","pos"),all.x=T)
luan_NeoAG_mhci_RelOnly=merge(luan_NeoAG_mhci_RelOnly,diagnosis,by=c("chr","pos"),all.x=T)
luan_NeoAG_mhci_RelOnly$Rel_Mut=str_count(str_to_upper(luan_NeoAG_mhci_RelOnly$rel_bases), as.character(luan_NeoAG_mhci_RelOnly$mut_nt))
luan_NeoAG_mhci_RelOnly$Rel_Wt=str_count(str_to_upper(luan_NeoAG_mhci_RelOnly$rel_bases), as.character(luan_NeoAG_mhci_RelOnly$wt_nt))
luan_NeoAG_mhci_RelOnly$Dx_Mut=str_count(str_to_upper(luan_NeoAG_mhci_RelOnly$dx_bases), as.character(luan_NeoAG_mhci_RelOnly$mut_nt))
luan_NeoAG_mhci_RelOnly$Dx_Wt=str_count(str_to_upper(luan_NeoAG_mhci_RelOnly$dx_bases), as.character(luan_NeoAG_mhci_RelOnly$wt_nt))
write.table(unique(luan_NeoAG_mhci_RelOnly[,c(1:22,31:34)]),"../../Documents/Progetto_LuciaGabri/Summary_LUAN_Annotato_NeoAGs_MHCI_RelOnly_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)


#LUAN non ha NeoAGs DiagOnly

#mhag=read.table('Summary_LUAN_Annotato_NeoAGs_MHCI_RelOnly_AllInfos.txt',head=T)
#mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
#mhag$logFC[is.na(mhag$logFC)] <- 0
#mhag$LUAN2_kallisto[is.na(mhag$LUAN2_kallisto)] <- 0
#mhag$LUAN1_kallisto[is.na(mhag$LUAN1_kallisto)] <- 0
#
# nrow(mhag[apply(mhag[,16:20], MARGIN = 1, function(x) any(x<=500 & x> 50)), ])
#[1] 7
# nrow(mhag[apply(mhag[,16:20], MARGIN = 1, function(x) any(x<=50 & x> 0)), ])
#[1] 1
#
#  
#  
#mhag=mhag[(mhag$LUAN2_kallisto>=0.58 & mhag$LUAN1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
#nrow(mhag)
#4
#
# table(mhag$Dx_Mut>=3)
#
#
#
# nrow(mhag[apply(mhag[,16:20], MARGIN = 1, function(x) any(x<=500 & x> 50)), ])
#[1] 7
# nrow(mhag[apply(mhag[,16:20], MARGIN = 1, function(x) any(x<=50 & x> 0)), ])
#[1] 1
# 
##A livello di trascritto
# table(mhag$LUAN2_kallisto > mhag$LUAN1_kallisto)
#
#FALSE  TRUE
#    2     1
# table(mhag$LUAN2_kallisto > mhag$LUAN1_kallisto)
#
#FALSE  TRUE
#    2     1
#
#wilcox.test(mhag$LUAN2_kallisto, mhag$LUAN1_kallisto,alternative="greater")
#
#        Wilcoxon rank sum test with continuity correction
#
#data:  mhag$LUAN2_kallisto and mhag$LUAN1_kallisto
#W = 3, p-value = 0.8157
#alternative hypothesis: true location shift is greater than 0
#
#Warning message:
#In wilcox.test.default(mhag$LUAN2_kallisto, mhag$LUAN1_kallisto,  :
#  cannot compute exact p-value with ties
#
#
#
#table(mhag$Rel_Mut>mhag$Dx_Mut)
#
#TRUE	
#   3
#
#wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")
#
#        Wilcoxon rank sum test with continuity correction
#
#data:  mhag$Rel_Mut and mhag$Dx_Mut
#W = 9, p-value = 0.02967
#alternative hypothesis: true location shift is greater than 0
#
#Warning message:
#In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
#  cannot compute exact p-value with ties
#
#MABI
#MABI RelDiag

cat StrongBinders_Annotati_MABI_RelDiag WeakBinders_Annotati_MABI_RelDiag | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ > Summary_MABI_Annotato_RelDiag_MHCI_NoHeaders_NoMismatches.txt


MABI1 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/MABI1_Allfiles_sorted.bam
MABI2 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/MABI2_Allfiles_sorted.bam

awk '{OFS="\t"; print $8":"$9"-"$9}' Summary_MABI_Annotato_RelDiag_MHCI_NoHeaders_NoMismatches.txt |  sort | uniq  > Summary_MABI_Annotato_RelDiag_MHCI_NoHeaders_NoMismatches_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/MABI1_Allfiles_sorted.bam >> MABI_NeoAGs_I_RelDiag_Binders_Regions_Diagnosis.txt; done < Summary_MABI_Annotato_RelDiag_MHCI_NoHeaders_NoMismatches_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/MABI2_Allfiles_sorted.bam >> MABI_NeoAGs_I_RelDiag_Binders_Regions_Relapse.txt; done < Summary_MABI_Annotato_RelDiag_MHCI_NoHeaders_NoMismatches_Regions.txt


mabi_NeoAG_mhci_RelDiag=read.table('Summary_MABI_Annotato_RelDiag_MHCI_NoHeaders_NoMismatches.txt')
colnames(mabi_NeoAG_mhci_RelDiag)=c("Gene_ID","Peptide","HLA","HLA","HLA","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA","IC50_HLA","IC,50_HLA")

mabi_NeoAG_mhci_RelDiag=merge(mabi_NeoAG_mhci_RelDiag,MABI_RNK,by="genes",all.x=T)
relapse=read.table('MABI_NeoAGs_I_RelDiag_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('MABI_NeoAGs_I_RelDiag_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
mabi_NeoAG_mhci_RelDiag=merge(mabi_NeoAG_mhci_RelDiag,relapse,by=c("chr","pos"),all.x=T)
mabi_NeoAG_mhci_RelDiag=merge(mabi_NeoAG_mhci_RelDiag,diagnosis,by=c("chr","pos"),all.x=T)
mabi_NeoAG_mhci_RelDiag$Rel_Mut=str_count(str_to_upper(mabi_NeoAG_mhci_RelDiag$rel_bases), as.character(mabi_NeoAG_mhci_RelDiag$mut_nt))
mabi_NeoAG_mhci_RelDiag$Rel_Wt=str_count(str_to_upper(mabi_NeoAG_mhci_RelDiag$rel_bases), as.character(mabi_NeoAG_mhci_RelDiag$wt_nt))
mabi_NeoAG_mhci_RelDiag$Dx_Mut=str_count(str_to_upper(mabi_NeoAG_mhci_RelDiag$dx_bases), as.character(mabi_NeoAG_mhci_RelDiag$mut_nt))
mabi_NeoAG_mhci_RelDiag$Dx_Wt=str_count(str_to_upper(mabi_NeoAG_mhci_RelDiag$dx_bases), as.character(mabi_NeoAG_mhci_RelDiag$wt_nt))
write.table(unique(mabi_NeoAG_mhci_RelDiag[,c(1:18,27:30)]),"../../Documents/Progetto_LuciaGabri/Summary_MABI_Annotato_NeoAGs_MHCI_RelDiag_NoMismatches_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)

#
#mhag=read.table('Summary_MABI_Annotato_NeoAGs_MHCI_RelDiag_AllInfos.txt',head=T)
#mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
#mhag$logFC[is.na(mhag$logFC)] <- 0
#mhag$MABI2_kallisto[is.na(mhag$MABI2_kallisto)] <- 0
#mhag$MABI1_kallisto[is.na(mhag$MABI1_kallisto)] <- 0
#
# nrow(mhag[apply(mhag[,16:20], MARGIN = 1, function(x) any(x<=500 & x> 50)), ])
#[1] 7
# nrow(mhag[apply(mhag[,16:20], MARGIN = 1, function(x) any(x<=50 & x> 0)), ])
#[1] 1
# 
#mhag=mhag[(mhag$MABI2_kallisto>=0.58 & mhag$MABI1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
#nrow(mhag)
#12
#
# table(mhag$Dx_Mut>=3)
#
#
#
##A livello di trascritto
# table(mhag$MABI2_kallisto > mhag$MABI1_kallisto)
#
#
#FALSE  TRUE
#    6     6
#	
#wilcox.test(mhag$MABI2_kallisto, mhag$MABI1_kallisto,alternative="greater")
#
#        Wilcoxon rank sum test with continuity correction
#
#data:  mhag$MABI2_kallisto and mhag$MABI1_kallisto
#W = 64, p-value = 0.6918
#alternative hypothesis: true location shift is greater than 0
#
#Warning message:
#In wilcox.test.default(mhag$MABI2_kallisto, mhag$MABI1_kallisto,  :
#  cannot compute exact p-value with ties
#
#
#table(mhag$Rel_Mut>mhag$Dx_Mut)
#
#FALSE  TRUE
#    7     5
#
#
#wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")
#
#        Wilcoxon rank sum test with continuity correction
#
#data:  mhag$Rel_Mut and mhag$Dx_Mut
#W = 81, p-value = 0.2991
#alternative hypothesis: true location shift is greater than 0
#
#Warning message:
#In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
#  cannot compute exact p-value with ties
#

#MOGE
#MOGE RelDiag

cat StrongBinders_Annotati_MOGE_RelDiag WeakBinders_Annotati_MOGE_RelDiag | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ > Summary_MOGE_Annotato_RelDiag_MHCI_NoHeaders_NoMismatches.txt


MOGE1 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE1_Allfiles_sorted.bam
MOGE2 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE2_Allfiles_sorted.bam

awk '{OFS="\t"; print $10":"$11"-"$11}'  Summary_MOGE_Annotato_RelDiag_MHCI_NoHeaders_NoMismatches.txt  |  sort | uniq  > Summary_MOGE_Annotato_RelDiag_MHCI_NoHeaders_NoMismatches_Regions.txt
 while read region;  do samtools mpileup -r $region   /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE1_Allfiles_sorted.bam >> MOGE_NeoAGs_I_RelDiag_Binders_Regions_Diagnosis.txt; done < Summary_MOGE_Annotato_RelDiag_MHCI_NoHeaders_NoMismatches_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE2_Allfiles_sorted.bam >> MOGE_NeoAGs_I_RelDiag_Binders_Regions_Relapse.txt; done < Summary_MOGE_Annotato_RelDiag_MHCI_NoHeaders_NoMismatches_Regions.txt

moge_NeoAG_mhci_RelDiag=read.table('Summary_MOGE_Annotato_RelDiag_MHCI_NoHeaders_NoMismatches.txt')
colnames(moge_NeoAG_mhci_RelDiag)=c("Gene_ID","Peptide","HLA","HLA","HLA","HLA","HLA","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA","IC50_HLA","IC,50_HLA","IC50_HLA","IC50_HLA")

moge_NeoAG_mhci_RelDiag=merge(moge_NeoAG_mhci_RelDiag,MOGE_RNK,by="genes",all.x=T)
moge_NeoAG_mhci_RelDiag$mut_nt="T"


relapse=read.table('MOGE_NeoAGs_I_RelDiag_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('MOGE_NeoAGs_I_RelDiag_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
moge_NeoAG_mhci_RelDiag=merge(moge_NeoAG_mhci_RelDiag,relapse,by=c("chr","pos"),all.x=T)
moge_NeoAG_mhci_RelDiag=merge(moge_NeoAG_mhci_RelDiag,diagnosis,by=c("chr","pos"),all.x=T)
moge_NeoAG_mhci_RelDiag$Rel_Mut=str_count(str_to_upper(moge_NeoAG_mhci_RelDiag$rel_bases), as.character(moge_NeoAG_mhci_RelDiag$mut_nt))
moge_NeoAG_mhci_RelDiag$Rel_Wt=str_count(str_to_upper(moge_NeoAG_mhci_RelDiag$rel_bases), as.character(moge_NeoAG_mhci_RelDiag$wt_nt))
moge_NeoAG_mhci_RelDiag$Dx_Mut=str_count(str_to_upper(moge_NeoAG_mhci_RelDiag$dx_bases), as.character(moge_NeoAG_mhci_RelDiag$mut_nt))
moge_NeoAG_mhci_RelDiag$Dx_Wt=str_count(str_to_upper(moge_NeoAG_mhci_RelDiag$dx_bases), as.character(moge_NeoAG_mhci_RelDiag$wt_nt))
 write.table(unique(moge_NeoAG_mhci_RelDiag[,c(1:22,31:34)]),"../../Documents/Progetto_LuciaGabri/Summary_MOGE_Annotato_NeoAGs_MHCI_RelDiag_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)


#MOGE RelOnly

cat StrongBinders_Annotati_MOGE_RelOnly WeakBinders_Annotati_MOGE_RelOnly | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ > Summary_MOGE_Annotato_RelOnly_MHCI_NoHeaders_NoMismatches.txt


MOGE1 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE1_Allfiles_sorted.bam
MOGE2 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE2_Allfiles_sorted.bam

awk '{OFS="\t"; print $10":"$11"-"$11}'  Summary_MOGE_Annotato_RelOnly_MHCI_NoHeaders_NoMismatches.txt  |  sort | uniq  > Summary_MOGE_Annotato_RelOnly_MHCI_NoHeaders_NoMismatches_Regions.txt
 while read region;  do samtools mpileup -r $region   /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE1_Allfiles_sorted.bam >> MOGE_NeoAGs_I_RelOnly_Binders_Regions_Diagnosis.txt; done < Summary_MOGE_Annotato_RelOnly_MHCI_NoHeaders_NoMismatches_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE2_Allfiles_sorted.bam >> MOGE_NeoAGs_I_RelOnly_Binders_Regions_Relapse.txt; done < Summary_MOGE_Annotato_RelOnly_MHCI_NoHeaders_NoMismatches_Regions.txt

moge_NeoAG_mhci_RelOnly=read.table('Summary_MOGE_Annotato_RelOnly_MHCI_NoHeaders_NoMismatches.txt')
colnames(moge_NeoAG_mhci_RelOnly)=c("Gene_ID","Peptide","HLA","HLA","HLA","HLA","HLA","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA","IC50_HLA","IC,50_HLA","IC50_HLA","IC50_HLA")

moge_NeoAG_mhci_RelOnly=merge(moge_NeoAG_mhci_RelOnly,MOGE_RNK,by="genes",all.x=T)
moge_NeoAG_mhci_RelOnly$mut_nt="T"


relapse=read.table('MOGE_NeoAGs_I_RelOnly_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('MOGE_NeoAGs_I_RelOnly_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
moge_NeoAG_mhci_RelOnly=merge(moge_NeoAG_mhci_RelOnly,relapse,by=c("chr","pos"),all.x=T)
moge_NeoAG_mhci_RelOnly=merge(moge_NeoAG_mhci_RelOnly,diagnosis,by=c("chr","pos"),all.x=T)
moge_NeoAG_mhci_RelOnly$Rel_Mut=str_count(str_to_upper(moge_NeoAG_mhci_RelOnly$rel_bases), as.character(moge_NeoAG_mhci_RelOnly$mut_nt))
moge_NeoAG_mhci_RelOnly$Rel_Wt=str_count(str_to_upper(moge_NeoAG_mhci_RelOnly$rel_bases), as.character(moge_NeoAG_mhci_RelOnly$wt_nt))
moge_NeoAG_mhci_RelOnly$Dx_Mut=str_count(str_to_upper(moge_NeoAG_mhci_RelOnly$dx_bases), as.character(moge_NeoAG_mhci_RelOnly$mut_nt))
moge_NeoAG_mhci_RelOnly$Dx_Wt=str_count(str_to_upper(moge_NeoAG_mhci_RelOnly$dx_bases), as.character(moge_NeoAG_mhci_RelOnly$wt_nt))
 write.table(unique(moge_NeoAG_mhci_RelOnly[,c(1:22,31:34)]),"../../Documents/Progetto_LuciaGabri/Summary_MOGE_Annotato_NeoAGs_MHCI_RelOnly_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)


#
#mhag=read.table('Summary_MOGE_Annotato_NeoAGs_MHCI_RelDiag_AllInfos.txt',head=T)
#mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
#mhag$logFC[is.na(mhag$logFC)] <- 0
#mhag$MOGE2_kallisto[is.na(mhag$MOGE2_kallisto)] <- 0
#mhag$MOGE1_kallisto[is.na(mhag$MOGE1_kallisto)] <- 0
#
# nrow(mhag[apply(mhag[,16:20], MARGIN = 1, function(x) any(x<=500 & x> 50)), ])
#[1] 7
# nrow(mhag[apply(mhag[,16:20], MARGIN = 1, function(x) any(x<=50 & x> 0)), ])
#[1] 1
# 
# 
#mhag=mhag[(mhag$MOGE2_kallisto>=0.58 & mhag$MOGE1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
#nrow(mhag)
#1
#
# table(mhag$Dx_Mut>=3)
#
#
#
##A livello di trascritto
# table(mhag$MOGE2_kallisto > mhag$MOGE1_kallisto)
#
#
#FALSE  
#    1  
#	
#wilcox.test(mhag$MOGE2_kallisto, mhag$MOGE1_kallisto,alternative="greater")
#
#        Wilcoxon rank sum test
#
#data:  mhag$MOGE2_kallisto and mhag$MOGE1_kallisto
#W = 0, p-value = 1
#alternative hypothesis: true location shift is greater than 0
#
# table(mhag$Rel_Mut>mhag$Dx_Mut)
#
#FALSE
#    1
#
#
#wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")
#
#        Wilcoxon rank sum test
#
#data:  mhag$Rel_Mut and mhag$Dx_Mut
#W = 0, p-value = 1
#alternative hypothesis: true location shift is greater than 0
#
#



#PIAG
#Piag non ha NeoAGs espressi



#PRELU
#PRELU non ha NeoAGs espressi 



#!/home/fsantaniello/.ctgb-sistemistica/bin/python
import sys
Mutationi=open(str(sys.argv[1])).readlines()
Excel=open(str(sys.argv[2])).readlines()
Summary=open(str(sys.argv[3])).readlines()
out=open(str(sys.argv[4]),'w')
header=Summary[0:2]

#Mutazioni=open('BESU1_BESU2_Somatic_MinorAGs_new_MHCII_PeptidesPrediction_MHCII.txt').readlines()
Mutazioni_mut= [i for i in Mutazioni if str(i.split('\t')[6]) !="synonymous_variant" ]

#Summary=open('RESULTS_BESU1_BESU2_Somatic_MinorAGs_new_MHCII_weak_binders_MHC_II_2').readlines()
Summary_mut= [i for i in Summary if "WT_" not in str(i.split('\t')[0])]
#Excel=open('RESULTS_BESU1_BESU2_Somatic_MinorAGs_new_MHCII_PeptidesPrediction_MinorAntigens_MHCII.xls').readlines()
Excel_mut= [i for i in Excel if "WT_"  not in str(i.split('\t')[2])]

#out=open('BESU_WeakBinders_Annotati.txt','w')
header=Summary[0:2]

Excel_mut= [i for i in Excel if "WT_"  not in str(i.split('\t')[2])]

if str(header).count("DRB1")==2 :
	print(''.join(header).replace(' DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'\t'.join(''.join(header).replace(' DRB1','\tDRB1').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')


binders=[]
for j in Mutazioni_mut:
	for i in Summary_mut[2:]:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0].rstrip()):
			binders.append(str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(i.split('\t')[0])+'\t'+str(i.split('\t')[1])+'\t'+str(i.split('\t')[2])+'\t'+str(i.split('\t')[3]).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip())
for i in binders:
	for j in Excel_mut[2:]:
		if str(i.split('\t')[3])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))

out.close()

elif str(header).count("DRB1")==1:

	print(''.join(header).replace(' DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'\t'.join(''.join(header).replace(' DRB1','\tDRB1').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')


binders=[]
for j in Mutazioni_mut:
	for i in Summary_mut[2:]:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0].rstrip()):
			binders.append(str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(i.split('\t')[0])+'\t'+str(i.split('\t')[1])+'\t'+str(i.split('\t')[2])+'\t'+str(i.split('\t')[3]).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip())
for i in binders:
	for j in Excel_mut[2:]:
		if str(i.split('\t')[3])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\n'))

out.close()


#Antigeni Minori MHCI
 for i in  Summary_*_Annotato_MinorAntigens.txt; do echo $i; prefix=$(echo $i | sed 's/.txt/_MHCI.txt/g'); echo $prefix;  sed -e '/Strong/,+1d' $i | sed -e '/Weak/,+1d' - | grep -v WT_ > $prefix ; done


#Antigeni Minori MHCII
 cat ALFE_StrongBinders_Annotati.txt ALFE_WeakBinders_Annotati.txt > Summary_ALFE_Annotato_MinorAntigens_MHCII.txt
 cat BESU_StrongBinders_Annotati.txt BESU_WeakBinders_Annotati.txt > Summary_BESU_Annotato_MinorAntigens_MHCII.txt
 cat CALU_StrongBinders_Annotati.txt CALU_WeakBinders_Annotati.txt > Summary_CALU_Annotato_MinorAntigens_MHCII.txt
 cat DEST_StrongBinders_Annotati.txt DEST_WeakBinders_Annotati.txt > Summary_DEST_Annotato_MinorAntigens_MHCII.txt
 cat DR4_StrongBinders_Annotati.txt DR4_WeakBinders_Annotati.txt > Summary_DR4_Annotato_MinorAntigens_MHCII.txt
 cat DR5_StrongBinders_Annotati.txt DR5_WeakBinders_Annotati.txt > Summary_DR5_Annotato_MinorAntigens_MHCII.txt
 cat FOCA_StrongBinders_Annotati.txt FOCA_WeakBinders_Annotati.txt > Summary_FOCA_Annotato_MinorAntigens_MHCII.txt
 cat GAGRA_StrongBinders_Annotati.txt GAGRA_WeakBinders_Annotati.txt > Summary_GAGRA_Annotato_MinorAntigens_MHCII.txt
 cat LUAN_StrongBinders_Annotati.txt LUAN_WeakBinders_Annotati.txt > Summary_LUAN_Annotato_MinorAntigens_MHCII.txt
 cat MABI_StrongBinders_Annotati.txt MABI_WeakBinders_Annotati.txt > Summary_MABI_Annotato_MinorAntigens_MHCII.txt
 cat MOGE_StrongBinders_Annotati.txt MOGE_WeakBinders_Annotati.txt > Summary_MOGE_Annotato_MinorAntigens_MHCII.txt
 cat PIAG_StrongBinders_Annotati.txt PIAG_WeakBinders_Annotati.txt > Summary_PIAG_Annotato_MinorAntigens_MHCII.txt
 cat PRELU_StrongBinders_Annotati.txt PRELU_WeakBinders_Annotati.txt > Summary_PRELU_Annotato_MinorAntigens_MHCII.txt
cat DEIV_StrongBinders_Annotati.txt DEIV_WeakBinders_Annotati.txt > Summary_DEIV_Annotato_MinorAntigens_MHCII.txt



#Calcolo reads su alleli per classe I


ALFE1 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE1_Allfiles_sorted.bam
ALFE2 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE2_Allfiles_sorted.bam

awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_ALFE_Annotato_MinorAntigens_MHCI_reformatted.txt | sort | uniq  > Summary_ALFE_Annotato_MinorAntigens_MHCI_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE1_Allfiles_sorted.bam >> ALFE_MHCI_Binders_Regions_Diagnosis.txt; done < Summary_ALFE_Annotato_MinorAntigens_MHCI_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE2_Allfiles_sorted.bam >> ALFE_MHCI_Binders_Regions_Relapse.txt; done < Summary_ALFE_Annotato_MinorAntigens_MHCI_Regions.txt
 

alfe_mhci=read.table('Summary_ALFE_Annotato_MinorAntigens_MHCI_reformatted.txt')
 colnames(alfe_mhci)=c("chr","pos","id","peptide","hla","hla","hla","hla","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50","ic50","ic50","ic50","ic50")
  alfe_mhci=merge(alfe_mhci,ALFE_RNK,by="genes",all.x=T)
relapse=read.table('ALFE_MHCI_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('ALFE_MHCI_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
alfe_mhci=merge(alfe_mhci,relapse,by=c("chr","pos"),all.x=T)
alfe_mhci=merge(alfe_mhci,diagnosis,by=c("chr","pos"),all.x=T)
alfe_mhci$Rel_Mut=str_count(str_to_upper(alfe_mhci$rel_bases), as.character(alfe_mhci$mut_nt))
alfe_mhci$Rel_Wt=str_count(str_to_upper(alfe_mhci$rel_bases), as.character(alfe_mhci$wt_nt))
alfe_mhci$Dx_Mut=str_count(str_to_upper(alfe_mhci$dx_bases), as.character(alfe_mhci$mut_nt))
alfe_mhci$Dx_Wt=str_count(str_to_upper(alfe_mhci$dx_bases), as.character(alfe_mhci$wt_nt))
 write.table(unique(alfe_mhci[,c(1:24,33:36)]),"Summary_ALFE_Annotato_MinorAntigens_MHCI_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)



BESU1 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/BESU1_Allfiles_sorted.bam
BESU2 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/BESU2_Allfiles_sorted.bam

awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_BESU_Annotato_MinorAntigens_MHCI_reformatted.txt  | sort | uniq  > Summary_BESU_Annotato_MinorAntigens_MHCI__Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/BESU1_Allfiles_sorted.bam >> BESU_MHCI__Binders_Regions_Diagnosis.txt; done < Summary_BESU_Annotato_MinorAntigens_MHCI__Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/BESU2_Allfiles_sorted.bam >> BESU_MHCI__Binders_Regions_Relapse.txt; done < Summary_BESU_Annotato_MinorAntigens_MHCI__Regions.txt


besu_mhci=read.table('Summary_BESU_Annotato_MinorAntigens_MHCI_reformatted.txt')
 colnames(besu_mhci)=c("chr","pos","id","peptide","hla","hla","hla","hla","hla","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50","ic50","ic50","ic50","ic50","ic50")
  besu_mhci=merge(besu_mhci,BESU_RNK,by="genes",all.x=T)
relapse=read.table('BESU_MHCI__Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('BESU_MHCI__Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
besu_mhci=merge(besu_mhci,relapse,by=c("chr","pos"),all.x=T)
besu_mhci=merge(besu_mhci,diagnosis,by=c("chr","pos"),all.x=T)
besu_mhci$Rel_Mut=str_count(str_to_upper(besu_mhci$rel_bases), as.character(besu_mhci$mut_nt))
besu_mhci$Rel_Wt=str_count(str_to_upper(besu_mhci$rel_bases), as.character(besu_mhci$wt_nt))
besu_mhci$Dx_Mut=str_count(str_to_upper(besu_mhci$dx_bases), as.character(besu_mhci$mut_nt))
besu_mhci$Dx_Wt=str_count(str_to_upper(besu_mhci$dx_bases), as.character(besu_mhci$wt_nt))
 write.table(unique(besu_mhci[,c(1:26,35:38)]),"Summary_BESU_Annotato_MinorAntigens_MHCI_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)

 
 

 
CALU1 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/CALU1_Allfiles_sorted.bam
CALU2 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/CALU2_Allfiles_sorted.bam


awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_CALU_Annotato_MinorAntigens_MHCI_reformatted.txt | sort | uniq  > Summary_CALU_Annotato_MinorAntigens_MHCI__Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/CALU1_Allfiles_sorted.bam >> CALU_MHCI__Binders_Regions_Diagnosis.txt; done < Summary_CALU_Annotato_MinorAntigens_MHCI__Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/CALU2_Allfiles_sorted.bam >> CALU_MHCI__Binders_Regions_Relapse.txt; done < Summary_CALU_Annotato_MinorAntigens_MHCI__Regions.txt
 
 calu_mhci=read.table('Summary_CALU_Annotato_MinorAntigens_MHCI_reformatted.txt')
 colnames(calu_mhci)=c("chr","pos","id","peptide","hla","hla","hla","hla","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50","ic50","ic50","ic50","ic50")
calu_mhci=merge(calu_mhci,CALU_RNK,by="genes",all.x=T)

relapse=read.table('CALU_MHCI__Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('CALU_MHCI__Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
calu_mhci=merge(calu_mhci,relapse,by=c("chr","pos"),all.x=T)
calu_mhci=merge(calu_mhci,diagnosis,by=c("chr","pos"),all.x=T)
calu_mhci$Rel_Mut=str_count(str_to_upper(calu_mhci$rel_bases), as.character(calu_mhci$mut_nt))
calu_mhci$Rel_Wt=str_count(str_to_upper(calu_mhci$rel_bases), as.character(calu_mhci$wt_nt))
calu_mhci$Dx_Mut=str_count(str_to_upper(calu_mhci$dx_bases), as.character(calu_mhci$mut_nt))
calu_mhci$Dx_Wt=str_count(str_to_upper(calu_mhci$dx_bases), as.character(calu_mhci$wt_nt))
 write.table(unique(calu_mhci[,c(1:24,33:36)]),"Summary_CALU_Annotato_MinorAntigens_MHCI_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)

#CALU
mhag=read.table('Summary_CALU_Annotato_MinorAntigens_MHCI_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$CALU2_kallisto[is.na(mhag$CALU2_kallisto)] <- 0
mhag$CALU1_kallisto[is.na(mhag$CALU1_kallisto)] <- 0

mhag=mhag[(mhag$CALU2_kallisto>=0.58 & mhag$CALU1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
1134
#A livello di trascritto
table(mhag$CALU2_kallisto > mhag$CALU1_kallisto)

FALSE  TRUE
  535   589
table(mhag$CALU2_kallisto < mhag$CALU1_kallisto)

FALSE  TRUE
  589   535


 wilcox.test(mhag$CALU2_kallisto, mhag$CALU1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$CALU2_kallisto and mhag$CALU1_kallisto
W = 641847, p-value = 0.2546
alternative hypothesis: true location shift is greater than 0



#A livello di trascritto
table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
  844   280

table(mhag$Rel_Mut<mhag$Dx_Mut)

 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 631688, p-value = 0.5
alternative hypothesis: true location shift is greater than 0

library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,CALU2_kallisto,CALU1_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,CALU2_kallisto,CALU1_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,CALU2_kallisto,CALU1_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,CALU2_kallisto,CALU1_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)
	

#DEIV
 awk '{OFS="\t"; print $11,$12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$14,$13,$15,$16,$17,$18,$19,$20,$21,$22}' Summary_DEIV_Annotato_MinorAntigens_MHCI.txt | sed '/^\s*$/d' > Summary_DEIV_Annotato_MinorAntigens_MHCI_reformatted.txt

DEIV2 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2_Allfiles_sorted.bam
DEIV2b = /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2b_Allfiles_sorted.bam


awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_DEIV_Annotato_MinorAntigens_MHCI_reformatted.txt | sort | uniq  > Summary_DEIV_Annotato_MinorAntigens_MHCI__Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2_Allfiles_sorted.bam >> DEIV_MHCI__Binders_Regions_Diagnosis.txt; done < Summary_DEIV_Annotato_MinorAntigens_MHCI__Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2b_Allfiles_sorted.bam >> DEIV_MHCI__Binders_Regions_Relapse.txt; done < Summary_DEIV_Annotato_MinorAntigens_MHCI__Regions.txt
 
 deiv_mhci=read.table('Summary_DEIV_Annotato_MinorAntigens_MHCI_reformatted.txt')
 colnames(deiv_mhci)=c("chr","pos","id","peptide","hla","hla","hla","hla","hla","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50","ic50","ic50","ic50","ic50","ic50")
deiv_mhci=merge(deiv_mhci,DEIV_RNK,by="genes",all.x=T)
relapse=read.table('DEIV_MHCI__Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DEIV_MHCI__Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
deiv_mhci=merge(deiv_mhci,relapse,by=c("chr","pos"),all.x=T)
deiv_mhci=merge(deiv_mhci,diagnosis,by=c("chr","pos"),all.x=T)
deiv_mhci$Rel_Mut=str_count(str_to_upper(deiv_mhci$rel_bases), as.character(deiv_mhci$mut_nt))
deiv_mhci$Rel_Wt=str_count(str_to_upper(deiv_mhci$rel_bases), as.character(deiv_mhci$wt_nt))
deiv_mhci$Dx_Mut=str_count(str_to_upper(deiv_mhci$dx_bases), as.character(deiv_mhci$mut_nt))
deiv_mhci$Dx_Wt=str_count(str_to_upper(deiv_mhci$dx_bases), as.character(deiv_mhci$wt_nt))
 write.table(unique(deiv_mhci[,c(1:26,35:38)]),"Summary_DEIV_Annotato_MinorAntigens_MHCI_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)

 
 
DEST1 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST1_Allfiles_sorted.bam
DEST2 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST2_Allfiles_sorted.bam

awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_DEST_Annotato_MinorAntigens_MHCI_reformatted.txt | sort | uniq  > Summary_DEST_Annotato_MinorAntigens_MHCI__Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST1_Allfiles_sorted.bam >> DEST_MHCI__Binders_Regions_Diagnosis.txt; done < Summary_DEST_Annotato_MinorAntigens_MHCI__Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST2_Allfiles_sorted.bam >> DEST_MHCI__Binders_Regions_Relapse.txt; done < Summary_DEST_Annotato_MinorAntigens_MHCI__Regions.txt
 
  dest_mhci=read.table('Summary_DEST_Annotato_MinorAntigens_MHCI_reformatted.txt')
 colnames(dest_mhci)=c("chr","pos","id","peptide","hla","hla","hla","hla","hla","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50","ic50","ic50","ic50","ic50","ic50")
dest_mhci=merge(dest_mhci,DEST_RNK,by="genes",all.x=T)
relapse=read.table('DEST_MHCI__Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DEST_MHCI__Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
dest_mhci=merge(dest_mhci,relapse,by=c("chr","pos"),all.x=T)
dest_mhci=merge(dest_mhci,diagnosis,by=c("chr","pos"),all.x=T)
dest_mhci$Rel_Mut=str_count(str_to_upper(dest_mhci$rel_bases), as.character(dest_mhci$mut_nt))
dest_mhci$Rel_Wt=str_count(str_to_upper(dest_mhci$rel_bases), as.character(dest_mhci$wt_nt))
dest_mhci$Dx_Mut=str_count(str_to_upper(dest_mhci$dx_bases), as.character(dest_mhci$mut_nt))
dest_mhci$Dx_Wt=str_count(str_to_upper(dest_mhci$dx_bases), as.character(dest_mhci$wt_nt))
 write.table(unique(dest_mhci[,c(1:26,35:38)]),"Summary_DEST_Annotato_MinorAntigens_MHCI_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)



 
 
# DR1_1 = /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DR1_1_Allfiles_sorted.bam
# DR1_2 = /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DR1_2_Allfiles_sorted.bam
# 
# 
# awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_ALFE_Annotato_MinorAntigens_MHCI_.txt | sort | uniq  > Summary_ALFE_Annotato_MinorAntigens_MHCI__Regions.txt
#  while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE1_Allfiles_sorted.bam >> ALFE_MHCI__Binders_Regions_Diagnosis.txt; done < Summary_ALFE_Annotato_MinorAntigens_MHCI__Regions.txt
#  while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE2_Allfiles_sorted.bam >> ALFE_MHCI__Binders_Regions_Relapse.txt; done < Summary_ALFE_Annotato_MinorAntigens_MHCI__Regions.txt
#  
#  
DR41 =  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR41_Allfiles_sorted.bam
DR42 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR42_Allfiles_sorted.bam

 
 awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_DR4_Annotato_MinorAntigens_MHCI_reformatted.txt | sort | uniq  > Summary_DR4_Annotato_MinorAntigens_MHCI__Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR41_Allfiles_sorted.bam>> DR4_MHCI__Binders_Regions_Diagnosis.txt; done < Summary_DR4_Annotato_MinorAntigens_MHCI__Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR42_Allfiles_sorted.bam >> DR4_MHCI__Binders_Regions_Relapse.txt; done < Summary_DR4_Annotato_MinorAntigens_MHCI__Regions.txt
 
  dr4_mhci=read.table('Summary_DR4_Annotato_MinorAntigens_MHCI_reformatted.txt')
 colnames(dr4_mhci)=c("chr","pos","id","peptide","hla","hla","hla","hla","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50","ic50","ic50","ic50","ic50")
dr4_mhci=merge(dr4_mhci,DR4_RNK,by="genes",all.x=T)
relapse=read.table('DR4_MHCI__Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DR4_MHCI__Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
dr4_mhci=merge(dr4_mhci,relapse,by=c("chr","pos"),all.x=T)
dr4_mhci=merge(dr4_mhci,diagnosis,by=c("chr","pos"),all.x=T)
dr4_mhci$Rel_Mut=str_count(str_to_upper(dr4_mhci$rel_bases), as.character(dr4_mhci$mut_nt))
dr4_mhci$Rel_Wt=str_count(str_to_upper(dr4_mhci$rel_bases), as.character(dr4_mhci$wt_nt))
dr4_mhci$Dx_Mut=str_count(str_to_upper(dr4_mhci$dx_bases), as.character(dr4_mhci$mut_nt))
dr4_mhci$Dx_Wt=str_count(str_to_upper(dr4_mhci$dx_bases), as.character(dr4_mhci$wt_nt))
 write.table(unique(dr4_mhci[,c(1:24,33:36)]),"Summary_DR4_Annotato_MinorAntigens_MHCI_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)






DR51 =  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR51_Allfiles_sorted.bam
DR52 =  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR52_Allfiles_sorted.bam
 
 awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_DR5_Annotato_MinorAntigens_MHCI_reformatted.txt | sort | uniq  > Summary_DR5_Annotato_MinorAntigens_MHCI__Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR51_Allfiles_sorted.bam >> DR5_MHCI__Binders_Regions_Diagnosis.txt; done < Summary_DR5_Annotato_MinorAntigens_MHCI__Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR52_Allfiles_sorted.bam >> DR5_MHCI__Binders_Regions_Relapse.txt; done < Summary_DR5_Annotato_MinorAntigens_MHCI__Regions.txt
 
 
  dr5_mhci=read.table('Summary_DR5_Annotato_MinorAntigens_MHCI_reformatted.txt')
 colnames(dr5_mhci)=c("chr","pos","id","peptide","hla","hla","hla","hla","hla","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50","ic50","ic50","ic50","ic50","ic50")
dr5_mhci=merge(dr5_mhci,DR5_RNK,by="genes",all.x=T)
relapse=read.table('DR5_MHCI__Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DR5_MHCI__Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
dr5_mhci=merge(dr5_mhci,relapse,by=c("chr","pos"),all.x=T)
dr5_mhci=merge(dr5_mhci,diagnosis,by=c("chr","pos"),all.x=T)
dr5_mhci$Rel_Mut=str_count(str_to_upper(dr5_mhci$rel_bases), as.character(dr5_mhci$mut_nt))
dr5_mhci$Rel_Wt=str_count(str_to_upper(dr5_mhci$rel_bases), as.character(dr5_mhci$wt_nt))
dr5_mhci$Dx_Mut=str_count(str_to_upper(dr5_mhci$dx_bases), as.character(dr5_mhci$mut_nt))
dr5_mhci$Dx_Wt=str_count(str_to_upper(dr5_mhci$dx_bases), as.character(dr5_mhci$wt_nt))
 write.table(unique(dr5_mhci[,c(1:26,35:38)]),"Summary_DR5_Annotato_MinorAntigens_MHCI_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)

 


 
FOCA1 =  /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/FOCA1_Allfiles_sorted.bam
FOCA2 =  /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/FOCA2_Allfiles_sorted.bam

 
 awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_FOCA_Annotato_MinorAntigens_MHCI_reformatted.txt | sort | uniq  > Summary_FOCA_Annotato_MinorAntigens_MHCI__Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/FOCA1_Allfiles_sorted.bam >> FOCA_MHCI__Binders_Regions_Diagnosis.txt; done < Summary_FOCA_Annotato_MinorAntigens_MHCI__Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/FOCA2_Allfiles_sorted.bam >> FOCA_MHCI__Binders_Regions_Relapse.txt; done < Summary_FOCA_Annotato_MinorAntigens_MHCI__Regions.txt

  foca_mhci=read.table('Summary_FOCA_Annotato_MinorAntigens_MHCI_reformatted.txt')
 colnames(foca_mhci)=c("chr","pos","id","peptide","hla","hla","hla","hla","hla","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50","ic50","ic50","ic50","ic50","ic50")
foca_mhci=merge(foca_mhci,FOCA_RNK,by="genes",all.x=T)
relapse=read.table('FOCA_MHCI__Binders_Regions_Relapse.txt',head=F,sep="\t")
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
#100152307

diagnosis=read.table('FOCA_MHCI__Binders_Regions_Diagnosis.txt',head=F,sep="\t")
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
foca_mhci=merge(foca_mhci,relapse,by=c("chr","pos"),all.x=T)
foca_mhci=merge(foca_mhci,diagnosis,by=c("chr","pos"),all.x=T)
foca_mhci$Rel_Mut=str_count(str_to_upper(foca_mhci$rel_bases), as.character(foca_mhci$mut_nt))
foca_mhci$Rel_Wt=str_count(str_to_upper(foca_mhci$rel_bases), as.character(foca_mhci$wt_nt))
foca_mhci$Dx_Mut=str_count(str_to_upper(foca_mhci$dx_bases), as.character(foca_mhci$mut_nt))
foca_mhci$Dx_Wt=str_count(str_to_upper(foca_mhci$dx_bases), as.character(foca_mhci$wt_nt))
 write.table(unique(foca_mhci[,c(1:26,35:38)]),"Summary_FOCA_Annotato_MinorAntigens_MHCI_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)



 
GAGRA1 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/GAGRA1_Allfiles_sorted.bam
GAGRA2 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/GAGRA2_Allfiles_sorted.bam

 
 awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_GAGRA_Annotato_MinorAntigens_MHCI_reformatted.txt | sort | uniq  > Summary_GAGRA_Annotato_MinorAntigens_MHCI__Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/GAGRA1_Allfiles_sorted.bam >> GAGRA_MHCI__Binders_Regions_Diagnosis.txt; done < Summary_GAGRA_Annotato_MinorAntigens_MHCI__Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/GAGRA2_Allfiles_sorted.bam >> GAGRA_MHCI__Binders_Regions_Relapse.txt; done < Summary_GAGRA_Annotato_MinorAntigens_MHCI__Regions.txt


  gagra_mhci=read.table('Summary_GAGRA_Annotato_MinorAntigens_MHCI_reformatted.txt')
 colnames(gagra_mhci)=c("chr","pos","id","peptide","hla","hla","hla","hla","hla","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50","ic50","ic50","ic50","ic50","ic50")
gagra_mhci=merge(gagra_mhci,GAGRA_RNK,by="genes",all.x=T)
relapse=read.table('GAGRA_MHCI__Binders_Regions_Relapse.txt',sep="\t",head=F)
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('GAGRA_MHCI__Binders_Regions_Diagnosis.txt',sep="\t",head=F)
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
gagra_mhci=merge(gagra_mhci,relapse,by=c("chr","pos"),all.x=T)
gagra_mhci=merge(gagra_mhci,diagnosis,by=c("chr","pos"),all.x=T)
gagra_mhci$Rel_Mut=str_count(str_to_upper(gagra_mhci$rel_bases), as.character(gagra_mhci$mut_nt))
gagra_mhci$Rel_Wt=str_count(str_to_upper(gagra_mhci$rel_bases), as.character(gagra_mhci$wt_nt))
gagra_mhci$Dx_Mut=str_count(str_to_upper(gagra_mhci$dx_bases), as.character(gagra_mhci$mut_nt))
gagra_mhci$Dx_Wt=str_count(str_to_upper(gagra_mhci$dx_bases), as.character(gagra_mhci$wt_nt))
 write.table(unique(gagra_mhci[,c(1:26,35:38)]),"Summary_GAGRA_Annotato_MinorAntigens_MHCI_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)


 
LUAN1 = /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN1_Allfiles_sorted.bam
LUAN2 = /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN2_Allfiles_sorted.bam

 
 awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_LUAN_Annotato_MinorAntigens_MHCI_reformatted.txt | sort | uniq  > Summary_LUAN_Annotato_MinorAntigens_MHCI__Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN1_Allfiles_sorted.bam >> LUAN_MHCI__Binders_Regions_Diagnosis.txt; done < Summary_LUAN_Annotato_MinorAntigens_MHCI__Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN2_Allfiles_sorted.bam >> LUAN_MHCI__Binders_Regions_Relapse.txt; done < Summary_LUAN_Annotato_MinorAntigens_MHCI__Regions.txt
 
  luan_mhci=read.table('Summary_LUAN_Annotato_MinorAntigens_MHCI_reformatted.txt')
 colnames(luan_mhci)=c("chr","pos","id","peptide","hla","hla","hla","hla","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50","ic50","ic50","ic50","ic50")
luan_mhci=merge(luan_mhci,LUAN_RNK,by="genes",all.x=T)
relapse=read.table('LUAN_MHCI__Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('LUAN_MHCI__Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
luan_mhci=merge(luan_mhci,relapse,by=c("chr","pos"),all.x=T)
luan_mhci=merge(luan_mhci,diagnosis,by=c("chr","pos"),all.x=T)
luan_mhci$Rel_Mut=str_count(str_to_upper(luan_mhci$rel_bases), as.character(luan_mhci$mut_nt))
luan_mhci$Rel_Wt=str_count(str_to_upper(luan_mhci$rel_bases), as.character(luan_mhci$wt_nt))
luan_mhci$Dx_Mut=str_count(str_to_upper(luan_mhci$dx_bases), as.character(luan_mhci$mut_nt))
luan_mhci$Dx_Wt=str_count(str_to_upper(luan_mhci$dx_bases), as.character(luan_mhci$wt_nt))
 write.table(unique(luan_mhci[,c(1:24,33:36)]),"Summary_LUAN_Annotato_MinorAntigens_MHCI_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)





MABI1 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/MABI1_Allfiles_sorted.bam
MABI2 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/MABI2_Allfiles_sorted.bam

 
 awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_MABI_Annotato_MinorAntigens_MHCI_reformatted.txt | sort | uniq  > Summary_MABI_Annotato_MinorAntigens_MHCI__Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/MABI1_Allfiles_sorted.bam >> MABI_MHCI__Binders_Regions_Diagnosis.txt; done < Summary_MABI_Annotato_MinorAntigens_MHCI__Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/MABI2_Allfiles_sorted.bam >> MABI_MHCI__Binders_Regions_Relapse.txt; done < Summary_MABI_Annotato_MinorAntigens_MHCI__Regions.txt
 
  mabi_mhci=read.table('Summary_MABI_Annotato_MinorAntigens_MHCI_reformatted.txt')
 colnames(mabi_mhci)=c("chr","pos","id","peptide","hla","hla","hla","hla","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50","ic50","ic50","ic50","ic50")
mabi_mhci=merge(mabi_mhci,MABI_RNK,by="genes",all.x=T)
relapse=read.table('MABI_MHCI__Binders_Regions_Relapse.txt',sep="\t")
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('MABI_MHCI__Binders_Regions_Diagnosis.txt',sep="\t")
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
mabi_mhci=merge(mabi_mhci,relapse,by=c("chr","pos"),all.x=T)
mabi_mhci=merge(mabi_mhci,diagnosis,by=c("chr","pos"),all.x=T)
mabi_mhci$Rel_Mut=str_count(str_to_upper(mabi_mhci$rel_bases), as.character(mabi_mhci$mut_nt))
mabi_mhci$Rel_Wt=str_count(str_to_upper(mabi_mhci$rel_bases), as.character(mabi_mhci$wt_nt))
mabi_mhci$Dx_Mut=str_count(str_to_upper(mabi_mhci$dx_bases), as.character(mabi_mhci$mut_nt))
mabi_mhci$Dx_Wt=str_count(str_to_upper(mabi_mhci$dx_bases), as.character(mabi_mhci$wt_nt))
 write.table(unique(mabi_mhci[,c(1:24,33:36)]),"Summary_MABI_Annotato_MinorAntigens_MHCI_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)
 
 
 
 
 
MOGE1 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE1_Allfiles_sorted.bam
MOGE2 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE2_Allfiles_sorted.bam

 
 awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_MOGE_Annotato_MinorAntigens_MHCI_reformatted.txt | sort | uniq  > Summary_MOGE_Annotato_MinorAntigens_MHCI__Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE1_Allfiles_sorted.bam >> MOGE_MHCI__Binders_Regions_Diagnosis.txt; done < Summary_MOGE_Annotato_MinorAntigens_MHCI__Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE2_Allfiles_sorted.bam >> MOGE_MHCI__Binders_Regions_Relapse.txt; done < Summary_MOGE_Annotato_MinorAntigens_MHCI__Regions.txt
 
 
  moge_mhci=read.table('Summary_MOGE_Annotato_MinorAntigens_MHCI_reformatted.txt')
 colnames(moge_mhci)=c("chr","pos","id","peptide","hla","hla","hla","hla","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50","ic50","ic50","ic50","ic50")
moge_mhci=merge(moge_mhci,MOGE_RNK,by="genes",all.x=T)
relapse=read.table('MOGE_MHCI__Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('MOGE_MHCI__Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
moge_mhci=merge(moge_mhci,relapse,by=c("chr","pos"),all.x=T)
moge_mhci=merge(moge_mhci,diagnosis,by=c("chr","pos"),all.x=T)
moge_mhci$Rel_Mut=str_count(str_to_upper(moge_mhci$rel_bases), as.character(moge_mhci$mut_nt))
moge_mhci$Rel_Wt=str_count(str_to_upper(moge_mhci$rel_bases), as.character(moge_mhci$wt_nt))
moge_mhci$Dx_Mut=str_count(str_to_upper(moge_mhci$dx_bases), as.character(moge_mhci$mut_nt))
moge_mhci$Dx_Wt=str_count(str_to_upper(moge_mhci$dx_bases), as.character(moge_mhci$wt_nt))
 write.table(unique(moge_mhci[,c(1:24,33:36)]),"Summary_MOGE_Annotato_MinorAntigens_MHCI_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)

 
 

 
PIAG1 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/PIAG1_Allfiles_sorted.bam
PIAG2 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/PIAG2_Allfiles_sorted.bam

 
 awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_PIAG_Annotato_MinorAntigens_MHCI_reformatted.txt | sort | uniq  > Summary_PIAG_Annotato_MinorAntigens_MHCI__Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/PIAG1_Allfiles_sorted.bam >> PIAG_MHCI__Binders_Regions_Diagnosis.txt; done < Summary_PIAG_Annotato_MinorAntigens_MHCI__Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/PIAG2_Allfiles_sorted.bam >> PIAG_MHCI__Binders_Regions_Relapse.txt; done < Summary_PIAG_Annotato_MinorAntigens_MHCI__Regions.txt

  piag_mhci=read.table('Summary_PIAG_Annotato_MinorAntigens_MHCI_reformatted.txt')
 colnames(piag_mhci)=c("chr","pos","id","peptide","hla","hla","hla","hla","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50","ic50","ic50","ic50","ic50")
piag_mhci=merge(piag_mhci,PIAG_RNK,by="genes",all.x=T)
relapse=read.table('PIAG_MHCI__Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('PIAG_MHCI__Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
piag_mhci=merge(piag_mhci,relapse,by=c("chr","pos"),all.x=T)
piag_mhci=merge(piag_mhci,diagnosis,by=c("chr","pos"),all.x=T)
piag_mhci$Rel_Mut=str_count(str_to_upper(piag_mhci$rel_bases), as.character(piag_mhci$mut_nt))
piag_mhci$Rel_Wt=str_count(str_to_upper(piag_mhci$rel_bases), as.character(piag_mhci$wt_nt))
piag_mhci$Dx_Mut=str_count(str_to_upper(piag_mhci$dx_bases), as.character(piag_mhci$mut_nt))
piag_mhci$Dx_Wt=str_count(str_to_upper(piag_mhci$dx_bases), as.character(piag_mhci$wt_nt))
 write.table(unique(piag_mhci[,c(1:24,33:36)]),"Summary_PIAG_Annotato_MinorAntigens_MHCI_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)

 




 
PRELU1 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/PRELU1_Allfiles_sorted.bam
PRELU2 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/PRELU2_Allfiles_sorted.bam

 
 awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_PRELU_Annotato_MinorAntigens_MHCI_reformatted.txt | sort | uniq  > Summary_PRELU_Annotato_MinorAntigens_MHCI__Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/PRELU1_Allfiles_sorted.bam >> PRELU_MHCI__Binders_Regions_Diagnosis.txt; done < Summary_PRELU_Annotato_MinorAntigens_MHCI__Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/PRELU2_Allfiles_sorted.bam >> PRELU_MHCI__Binders_Regions_Relapse.txt; done < Summary_PRELU_Annotato_MinorAntigens_MHCI__Regions.txt
 


prelu_mhci=read.table('Summary_PRELU_Annotato_MinorAntigens_MHCI_reformatted.txt')
colnames(prelu_mhci)=c("chr","pos","id","peptide","hla","hla","hla","hla","hla","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50","ic50","ic50","ic50","ic50","ic50")
prelu_mhci=merge(prelu_mhci,PRELU_RNK,by="genes",all.x=T)
relapse=read.table('PRELU_MHCI__Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('PRELU_MHCI__Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
prelu_mhci=merge(prelu_mhci,relapse,by=c("chr","pos"),all.x=T)
prelu_mhci=merge(prelu_mhci,diagnosis,by=c("chr","pos"),all.x=T)
prelu_mhci$Rel_Mut=str_count(str_to_upper(prelu_mhci$rel_bases), as.character(prelu_mhci$mut_nt))
prelu_mhci$Rel_Wt=str_count(str_to_upper(prelu_mhci$rel_bases), as.character(prelu_mhci$wt_nt))
prelu_mhci$Dx_Mut=str_count(str_to_upper(prelu_mhci$dx_bases), as.character(prelu_mhci$mut_nt))
prelu_mhci$Dx_Wt=str_count(str_to_upper(prelu_mhci$dx_bases), as.character(prelu_mhci$wt_nt))
 write.table(unique(prelu_mhci[,c(1:26,35:38)]),"Summary_PRELU_Annotato_MinorAntigens_MHCI_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)


#mhAGs Classe I

#ALFE
mhag=read.table('Summary_ALFE_Annotato_MinorAntigens_MHCI_AllInfos.txt',head=T)

mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$ALFE2_kallisto[is.na(mhag$ALFE2_kallisto)] <- 0
mhag$ALFE1_kallisto[is.na(mhag$ALFE1_kallisto)] <- 0
mhag$Rel_Mut[is.na(mhag$Rel_Mut)] <- 0
mhag$Dx_Mut[is.na(mhag$Dx_Mut)] <- 0
mhag=mhag[mhag$ALFE2_kallisto>=0.58  & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
392



library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,ALFE2_kallisto,ALFE1_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,ALFE2_kallisto,ALFE1_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,ALFE2_kallisto,ALFE1_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,ALFE2_kallisto,ALFE1_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)

nrow(result)
	
table(result$ALFE2_kallisto<result$ALFE1_kallisto)

table(result$ALFE2_kallisto>result$ALFE1_kallisto)

table(result$Rel_Mut<result$Dx_Mut)

table(result$Rel_Mut>result$Dx_Mut)
	
#A livello di trascritto
 table(mhag$ALFE2_kallisto > mhag$ALFE1_kallisto)

FALSE  TRUE
  211   143
  
 table(mhag$ALFE2_kallisto < mhag$ALFE1_kallisto)


 wilcox.test(mhag$ALFE2_kallisto, mhag$ALFE1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$ALFE2_kallisto and mhag$ALFE1_kallisto
W = 61288, p-value = 0.6928
alternative hypothesis: true location shift is greater than 0



#A livello di trascritto
table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
  269    77

table(mhag$Rel_Mut<mhag$Dx_Mut)
FALSE  TRUE
   85   261

 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 43840.5, p-value = 1
alternative hypothesis: true location shift is greater than 0



 
#BESU
mhag=read.table('Summary_BESU_Annotato_MinorAntigens_MHCI_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$BESU2_kallisto[is.na(mhag$BESU2_kallisto)] <- 0
mhag$BESU1_kallisto[is.na(mhag$BESU1_kallisto)] <- 0
mhag$Rel_Mut[is.na(mhag$Rel_Mut)] <- 0
mhag$Dx_Mut[is.na(mhag$Dx_Mut)] <- 0

mhag=mhag[mhag$BESU2_kallisto>=0.58 & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
774


library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,BESU2_kallisto,BESU1_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,BESU2_kallisto,BESU1_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,BESU2_kallisto,BESU1_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,BESU2_kallisto,BESU1_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)
	
nrow(result)
table(result$BESU2_kallisto<result$BESU1_kallisto)

table(result$BESU2_kallisto>result$BESU1_kallisto)

table(result$Rel_Mut<result$Dx_Mut)

table(result$Rel_Mut>result$Dx_Mut)



#A livello di trascritto
 table(mhag$BESU2_kallisto > mhag$BESU1_kallisto)

FALSE  TRUE
  294   417

 table(mhag$BESU2_kallisto < mhag$BESU1_kallisto)

table(mhag$BESU2_kallisto < mhag$BESU1_kallisto)

FALSE  TRUE
  417   294



 wilcox.test(mhag$BESU2_kallisto, mhag$BESU1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$BESU2_kallisto and mhag$BESU1_kallisto
W = 261696, p-value = 0.1242
alternative hypothesis: true location shift is greater than 0




#A livello di trascritto
 table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
  290   413

 table(mhag$Rel_Mut<mhag$Dx_Mut)

FALSE  TRUE
  445   258


 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 272560, p-value = 0.001589
alternative hypothesis: true location shift is greater than 0




#CALU
mhag=read.table('Summary_CALU_Annotato_MinorAntigens_MHCI_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$CALU2_kallisto[is.na(mhag$CALU2_kallisto)] <- 0
mhag$CALU1_kallisto[is.na(mhag$CALU1_kallisto)] <- 0
mhag$Rel_Mut[is.na(mhag$Rel_Mut)] <- 0
mhag$Dx_Mut[is.na(mhag$Dx_Mut)] <- 0
mhag=mhag[mhag$CALU2_kallisto>=0.58 & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
1134
#A livello di trascritto
table(mhag$CALU2_kallisto > mhag$CALU1_kallisto)

FALSE  TRUE
  535   589
table(mhag$CALU2_kallisto < mhag$CALU1_kallisto)

FALSE  TRUE
  589   535


 wilcox.test(mhag$CALU2_kallisto, mhag$CALU1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$CALU2_kallisto and mhag$CALU1_kallisto
W = 641847, p-value = 0.2546
alternative hypothesis: true location shift is greater than 0



#A livello di trascritto
table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
  844   280

table(mhag$Rel_Mut<mhag$Dx_Mut)

 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 631688, p-value = 0.5
alternative hypothesis: true location shift is greater than 0



library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,CALU2_kallisto,CALU1_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,CALU2_kallisto,CALU1_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,CALU2_kallisto,CALU1_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,CALU2_kallisto,CALU1_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)

nrow(result)
table(result$CALU2_kallisto<result$CALU1_kallisto)

table(result$CALU2_kallisto>result$CALU1_kallisto)

table(result$Rel_Mut<result$Dx_Mut)


table(result$Rel_Mut>result$Dx_Mut)


	

#DEIV


#DEIV
mhag=read.table('Summary_DEIV_Annotato_MinorAntigens_MHCI_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$DEIV2_kallisto[is.na(mhag$DEIV2_kallisto)] <- 0
mhag$DEIV2b_kallisto[is.na(mhag$DEIV2b_kallisto)] <- 0
mhag$Rel_Mut[is.na(mhag$Rel_Mut)] <- 0
mhag$Dx_Mut[is.na(mhag$Dx_Mut)] <- 0

mhag=mhag[mhag$DEIV2_kallisto>=0.58  & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
661
#A livello di trascritto
table(mhag$DEIV2b_kallisto > mhag$DEIV2_kallisto)

FALSE  TRUE
  258   358

 table(mhag$DEIV2b_kallisto < mhag$DEIV2_kallisto)

FALSE  TRUE
  358   258



 wilcox.test(mhag$DEIV2b_kallisto, mhag$DEIV2_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$DEIV2b_kallisto and mhag$DEIV2_kallisto
W = 200422, p-value = 0.0434
alternative hypothesis: true location shift is greater than 0



#A livello di trascritto
table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
  202   411


table(mhag$Rel_Mut<mhag$Dx_Mut)

FALSE  TRUE
  438   175

wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 220947.5, p-value = 1.177e-07
alternative hypothesis: true location shift is greater than 0




library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,DEIV2b_kallisto,DEIV2_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,DEIV2b_kallisto,DEIV2_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,DEIV2b_kallisto,DEIV2_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,DEIV2b_kallisto,DEIV2_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)

nrow(result)
	
table(result$DEIV2b_kallisto<result$DEIV2_kallisto)

table(result$DEIV2b_kallisto>result$DEIV2_kallisto)

table(result$Rel_Mut<result$Dx_Mut)


table(result$Rel_Mut>result$Dx_Mut)

	
#DEST
mhag=read.table('Summary_DEST_Annotato_MinorAntigens_MHCI_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$DEST2_kallisto[is.na(mhag$DEST2_kallisto)] <- 0
mhag$DEST1_kallisto[is.na(mhag$DEST1_kallisto)] <- 0
mhag$Rel_Mut[is.na(mhag$Rel_Mut)] <- 0
mhag$Dx_Mut[is.na(mhag$Dx_Mut)] <- 0

mhag=mhag[mhag$DEST2_kallisto>=0.58  & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
496
#A livello di trascritto
table(mhag$DEST2_kallisto > mhag$DEST1_kallisto)

FALSE  TRUE
  219   234
table(mhag$DEST2_kallisto < mhag$DEST1_kallisto)

FALSE  TRUE
  234   219


wilcox.test(mhag$DEST2_kallisto, mhag$DEST1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$DEST2_kallisto and mhag$DEST1_kallisto
W = 101763, p-value = 0.5846
alternative hypothesis: true location shift is greater than 0



#A livello di trascritto
table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
  181   268
table(mhag$Rel_Mut<mhag$Dx_Mut)

FALSE  TRUE
  283   166



 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 108636.5, p-value = 0.03803
alternative hypothesis: true location shift is greater than 0


library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,DEST2_kallisto,DEST1_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,DEST2_kallisto,DEST1_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,DEST2_kallisto,DEST1_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,DEST2_kallisto,DEST1_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)
	
nrow(result)

table(result$DEST2_kallisto<result$DEST1_kallisto)

table(result$DEST2_kallisto>result$DEST1_kallisto)

table(result$Rel_Mut<result$Dx_Mut)


table(result$Rel_Mut>result$Dx_Mut)




#DR4
mhag=read.table('Summary_DR4_Annotato_MinorAntigens_MHCI_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$DR42_kallisto[is.na(mhag$DR42_kallisto)] <- 0
mhag$DR41_kallisto[is.na(mhag$DR41_kallisto)] <- 0
mhag$Rel_Mut[is.na(mhag$Rel_Mut)] <- 0
mhag$Dx_Mut[is.na(mhag$Dx_Mut)] <- 0

mhag=mhag[mhag$DR42_kallisto>=0.58 & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
465
#A livello di trascritto
 table(mhag$DR42_kallisto > mhag$DR41_kallisto)

FALSE  TRUE
  170   211

 table(mhag$DR42_kallisto < mhag$DR41_kallisto)

FALSE  TRUE
  211   170


wilcox.test(mhag$DR42_kallisto, mhag$DR41_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$DR42_kallisto and mhag$DR41_kallisto
W = 72523, p-value = 0.5076
alternative hypothesis: true location shift is greater than 0


#A livello di trascritto
   table(mhag$Rel_Mut>mhag$Dx_Mut)

table(mhag$Rel_Mut<mhag$Dx_Mut)

FALSE  TRUE
  176   189



 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 66839, p-value = 0.822
alternative hypothesis: true location shift is greater than 0



library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,DR42_kallisto,DR41_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,DR42_kallisto,DR41_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,DR42_kallisto,DR41_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,DR42_kallisto,DR41_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)
nrow(result)
	table(result$DR42_kallisto<result$DR41_kallisto)

table(result$DR42_kallisto>result$DR41_kallisto)

table(result$Rel_Mut<result$Dx_Mut)


table(result$Rel_Mut>result$Dx_Mut)





#DR5
mhag=read.table('Summary_DR5_Annotato_MinorAntigens_MHCI_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$DR52_kallisto[is.na(mhag$DR52_kallisto)] <- 0
mhag$DR51_kallisto[is.na(mhag$DR51_kallisto)] <- 0
mhag$Rel_Mut[is.na(mhag$Rel_Mut)] <- 0
mhag$Dx_Mut[is.na(mhag$Dx_Mut)] <- 0

mhag=mhag[mhag$DR52_kallisto>=0.58  & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
305
#A livello di trascritto
 table(mhag$DR52_kallisto > mhag$DR51_kallisto)

FALSE  TRUE
  134   119

table(mhag$DR52_kallisto < mhag$DR51_kallisto)

FALSE  TRUE
  119   134



 wilcox.test(mhag$DR52_kallisto, mhag$DR51_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$DR52_kallisto and mhag$DR51_kallisto
W = 30378, p-value = 0.8388
alternative hypothesis: true location shift is greater than 0


#A livello di trascritto
 table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
   75   160


 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 34677, p-value = 0.0007908
alternative hypothesis: true location shift is greater than 0



library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,DR52_kallisto,DR51_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,DR52_kallisto,DR51_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,DR52_kallisto,DR51_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,DR52_kallisto,DR51_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)
	
nrow(result)
table(result$DR52_kallisto<result$DR51_kallisto)

table(result$DR52_kallisto>result$DR51_kallisto)

table(result$Rel_Mut<result$Dx_Mut)


table(result$Rel_Mut>result$Dx_Mut)


  
#FOCA
mhag=read.table('Summary_FOCA_Annotato_MinorAntigens_MHCI_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$FOCA2_kallisto[is.na(mhag$FOCA2_kallisto)] <- 0
mhag$FOCA1_kallisto[is.na(mhag$FOCA1_kallisto)] <- 0
mhag$Rel_Mut[is.na(mhag$Rel_Mut)] <- 0
mhag$Dx_Mut[is.na(mhag$Dx_Mut)] <- 0

mhag=mhag[mhag$FOCA2_kallisto>=0.58  & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
550
#A livello di trascritto
table(mhag$FOCA2_kallisto > mhag$FOCA1_kallisto)

FALSE  TRUE
    7     4

wilcox.test(mhag$FOCA2_kallisto, mhag$FOCA1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$FOCA2_kallisto and mhag$FOCA1_kallisto
W = 56, p-value = 0.6287
alternative hypothesis: true location shift is greater than 0

Warning message:
In wilcox.test.default(mhag$FOCA2_kallisto, mhag$FOCA1_kallisto,  :
  cannot compute exact p-value with ties


#A livello di trascritto
table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
    8     2
table(mhag$Rel_Mut<mhag$Dx_Mut)

FALSE  TRUE
    3     7


 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 46.5, p-value = 0.7376
alternative hypothesis: true location shift is greater than 0

Warning message:
In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
  cannot compute exact p-value with ties

library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,FOCA2_kallisto,FOCA1_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,FOCA2_kallisto,FOCA1_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,FOCA2_kallisto,FOCA1_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,FOCA2_kallisto,FOCA1_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)
	
	
table(result$FOCA2_kallisto<result$FOCA1_kallisto)

table(result$FOCA2_kallisto>result$FOCA1_kallisto)

table(result$Rel_Mut<result$Dx_Mut)


table(result$Rel_Mut>result$Dx_Mut)



#GAGRA
mhag=read.table('Summary_GAGRA_Annotato_MinorAntigens_MHCI_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$GAGRA2_kallisto[is.na(mhag$GAGRA2_kallisto)] <- 0
mhag$GAGRA1_kallisto[is.na(mhag$GAGRA1_kallisto)] <- 0
mhag$Rel_Mut[is.na(mhag$Rel_Mut)] <- 0
mhag$Dx_Mut[is.na(mhag$Dx_Mut)] <- 0

mhag=mhag[mhag$GAGRA2_kallisto>=0.58  & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
777



library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,GAGRA2_kallisto,GAGRA1_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,GAGRA2_kallisto,GAGRA1_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,GAGRA2_kallisto,GAGRA1_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,GAGRA2_kallisto,GAGRA1_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)
	
nrow(result)

table(result$GAGRA2_kallisto<result$GAGRA1_kallisto)

table(result$GAGRA2_kallisto>result$GAGRA1_kallisto)

table(result$Rel_Mut<result$Dx_Mut)


table(result$Rel_Mut>result$Dx_Mut)



#A livello di trascritto
 table(mhag$GAGRA2_kallisto > mhag$GAGRA1_kallisto)

FALSE  TRUE
    4    12


 wilcox.test(mhag$GAGRA2_kallisto, mhag$GAGRA1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$GAGRA2_kallisto and mhag$GAGRA1_kallisto
W = 138, p-value = 0.36
alternative hypothesis: true location shift is greater than 0

Warning message:
In wilcox.test.default(mhag$GAGRA2_kallisto, mhag$GAGRA1_kallisto,  :
  cannot compute exact p-value with ties

table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
    6    10



 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 143.5, p-value = 0.2854
alternative hypothesis: true location shift is greater than 0

Warning message:
In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
  cannot compute exact p-value with ties



#LUAN
mhag=read.table('Summary_LUAN_Annotato_MinorAntigens_MHCI_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$LUAN2_kallisto[is.na(mhag$LUAN2_kallisto)] <- 0
mhag$LUAN1_kallisto[is.na(mhag$LUAN1_kallisto)] <- 0
mhag$Rel_Mut[is.na(mhag$Rel_Mut)] <- 0
mhag$Dx_Mut[is.na(mhag$Dx_Mut)] <- 0

mhag=mhag[mhag$LUAN2_kallisto>=0.58 & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
290


library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,LUAN2_kallisto,LUAN1_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,LUAN2_kallisto,LUAN1_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,LUAN2_kallisto,LUAN1_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,LUAN2_kallisto,LUAN1_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)
	
	table(result$LUAN2_kallisto<result$LUAN1_kallisto)

table(result$LUAN2_kallisto>result$LUAN1_kallisto)

table(result$Rel_Mut<result$Dx_Mut)


table(result$Rel_Mut>result$Dx_Mut)


#A livello di trascritto
 table(mhag$LUAN2_kallisto > mhag$LUAN1_kallisto)

FALSE  TRUE
  143   111

 wilcox.test(mhag$LUAN2_kallisto, mhag$LUAN1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$LUAN2_kallisto and mhag$LUAN1_kallisto
W = 32591, p-value = 0.4203
alternative hypothesis: true location shift is greater than 0

 table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
  127   122

table(mhag$Rel_Mut<mhag$Dx_Mut)

FALSE  TRUE
  133   116



wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 32818.5, p-value = 0.2316
alternative hypothesis: true location shift is greater than 0



#MABI
mhag=read.table('Summary_MABI_Annotato_MinorAntigens_MHCI_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$MABI2_kallisto[is.na(mhag$MABI2_kallisto)] <- 0
mhag$MABI1_kallisto[is.na(mhag$MABI1_kallisto)] <- 0
mhag$Rel_Mut[is.na(mhag$Rel_Mut)] <- 0
mhag$Dx_Mut[is.na(mhag$Dx_Mut)] <- 0

mhag=mhag[mhag$MABI2_kallisto>=0.58  & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
518




library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,MABI2_kallisto,MABI1_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,MABI2_kallisto,MABI1_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,MABI2_kallisto,MABI1_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,MABI2_kallisto,MABI1_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)

nrow(result)
	table(result$MABI2_kallisto<result$MABI1_kallisto)

table(result$MABI2_kallisto>result$MABI1_kallisto)

table(result$Rel_Mut<result$Dx_Mut)


table(result$Rel_Mut>result$Dx_Mut)



#A livello di trascritto
 table(mhag$MABI2_kallisto > mhag$MABI1_kallisto)

FALSE  TRUE
  272   198

  
 wilcox.test(mhag$MABI2_kallisto, mhag$MABI1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$MABI2_kallisto and mhag$MABI1_kallisto
W = 106240, p-value = 0.8441
alternative hypothesis: true location shift is greater than 0


 table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
  236   203


 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 103291.5, p-value = 0.4873
alternative hypothesis: true location shift is greater than 0

  
#MOGE
mhag=read.table('Summary_MOGE_Annotato_MinorAntigens_MHCI_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$MOGE2_kallisto[is.na(mhag$MOGE2_kallisto)] <- 0
mhag$MOGE1_kallisto[is.na(mhag$MOGE1_kallisto)] <- 0
mhag$Rel_Mut[is.na(mhag$Rel_Mut)] <- 0
mhag$Dx_Mut[is.na(mhag$Dx_Mut)] <- 0

mhag=mhag[mhag$MOGE2_kallisto>=0.58 & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
241





library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,MOGE2_kallisto,MOGE1_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,MOGE2_kallisto,MOGE1_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,MOGE2_kallisto,MOGE1_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,MOGE2_kallisto,MOGE1_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)

nrow(result)
	table(result$MOGE2_kallisto<result$MOGE1_kallisto)

table(result$MOGE2_kallisto>result$MOGE1_kallisto)

table(result$Rel_Mut<result$Dx_Mut)


table(result$Rel_Mut>result$Dx_Mut)





#A livello di trascritto
table(mhag$MOGE2_kallisto > mhag$MOGE1_kallisto)

FALSE  TRUE
  122   107


 wilcox.test(mhag$MOGE2_kallisto, mhag$MOGE1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$MOGE2_kallisto and mhag$MOGE1_kallisto
W = 26269, p-value = 0.4865
alternative hypothesis: true location shift is greater than 0


table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
  105   116

table(mhag$Rel_Mut<mhag$Dx_Mut)

FALSE  TRUE
  126    95

 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 26529, p-value = 0.1887
alternative hypothesis: true location shift is greater than 0




#PIAG
mhag=read.table('Summary_PIAG_Annotato_MinorAntigens_MHCI_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$PIAG2_kallisto[is.na(mhag$PIAG2_kallisto)] <- 0
mhag$PIAG1_kallisto[is.na(mhag$PIAG1_kallisto)] <- 0
mhag$Rel_Mut[is.na(mhag$Rel_Mut)] <- 0
mhag$Dx_Mut[is.na(mhag$Dx_Mut)] <- 0
mhag=mhag[mhag$PIAG2_kallisto>=0.58 & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
283



library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,PIAG2_kallisto,PIAG1_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,PIAG2_kallisto,PIAG1_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,PIAG2_kallisto,PIAG1_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,logFC,Qvalue,PIAG2_kallisto,PIAG1_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)
nrow(result)
	table(result$PIAG2_kallisto<result$PIAG1_kallisto)

table(result$PIAG2_kallisto>result$PIAG1_kallisto)

table(result$Rel_Mut<result$Dx_Mut)


table(result$Rel_Mut>result$Dx_Mut)





#A livello di trascritto
table(mhag$PIAG2_kallisto > mhag$PIAG1_kallisto)

FALSE  TRUE
  107    90

 wilcox.test(mhag$PIAG2_kallisto, mhag$PIAG1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$PIAG2_kallisto and mhag$PIAG1_kallisto
W = 18649, p-value = 0.7482
alternative hypothesis: true location shift is greater than 0


table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
    4     9


 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 1327, p-value = 0.4141
alternative hypothesis: true location shift is greater than 0


#PRELU
mhag=read.table('Summary_PRELU_Annotato_MinorAntigens_MHCI_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$PRELU2_kallisto[is.na(mhag$PRELU2_kallisto)] <- 0
mhag$PRELU1_kallisto[is.na(mhag$PRELU1_kallisto)] <- 0
mhag$Rel_Mut[is.na(mhag$Rel_Mut)] <- 0
mhag$Dx_Mut[is.na(mhag$Dx_Mut)] <- 0
mhag=mhag[mhag$PRELU2_kallisto>=0.58  & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
323



library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,PRELU2_kallisto,PRELU1_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,PRELU2_kallisto,PRELU1_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,PRELU2_kallisto,PRELU1_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,hla.2,hla.3,hla.4,hla.5,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,ic50.2,ic50.3,ic50.4,ic50.5,logFC,Qvalue,PRELU2_kallisto,PRELU1_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)
	
nrow(result)
table(result$PRELU2_kallisto<result$PRELU1_kallisto)

table(result$PRELU2_kallisto>result$PRELU1_kallisto)

table(result$Rel_Mut<result$Dx_Mut)


table(result$Rel_Mut>result$Dx_Mut)




#A livello di trascritto
table(mhag$PRELU2_kallisto > mhag$PRELU1_kallisto)

FALSE  TRUE
   84   159


 wilcox.test(mhag$PRELU2_kallisto, mhag$PRELU1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$PRELU2_kallisto and mhag$PRELU1_kallisto
W = 32081, p-value = 0.04936
alternative hypothesis: true location shift is greater than 0



table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
    8     5


 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 2133, p-value = 0.01668
alternative hypothesis: true location shift is greater than 0



#Calcolo reads su alleli per classe II


ALFE1 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE1_Allfiles_sorted.bam
ALFE2 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE2_Allfiles_sorted.bam

awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_ALFE_Annotato_MinorAntigens_MHCII.txt | sort | uniq  > Summary_ALFE_Annotato_MinorAntigens_MHCII_Regions.txt
while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE1_Allfiles_sorted.bam >> ALFE_MHCII_Binders_Regions_Diagnosis.txt; done < Summary_ALFE_Annotato_MinorAntigens_MHCII_Regions.txt
while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE2_Allfiles_sorted.bam >> ALFE_MHCII_Binders_Regions_Relapse.txt; done < Summary_ALFE_Annotato_MinorAntigens_MHCII_Regions.txt

alfe_mhcii=read.table('Summary_ALFE_Annotato_MinorAntigens_MHCII.txt')
colnames(alfe_mhcii)=c("chr","pos","id","peptide","hla","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50","ic50")
alfe_mhcii=merge(alfe_mhcii,ALFE_RNK,by="genes",all.x=T)
relapse=read.table('ALFE_MHCII_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('ALFE_MHCII_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
alfe_mhcii=merge(alfe_mhcii,relapse,by=c("chr","pos"),all.x=T)
alfe_mhcii=merge(alfe_mhcii,diagnosis,by=c("chr","pos"),all.x=T)
alfe_mhcii$Rel_Mut=str_count(str_to_upper(alfe_mhcii$rel_bases), as.character(alfe_mhcii$mut_nt))
alfe_mhcii$Rel_Wt=str_count(str_to_upper(alfe_mhcii$rel_bases), as.character(alfe_mhcii$wt_nt))
alfe_mhcii$Dx_Mut=str_count(str_to_upper(alfe_mhcii$dx_bases), as.character(alfe_mhcii$mut_nt))
alfe_mhcii$Dx_Wt=str_count(str_to_upper(alfe_mhcii$dx_bases), as.character(alfe_mhcii$wt_nt))
 write.table(unique(alfe_mhcii[,c(1:18,27:30)]),"Summary_ALFE_Annotato_MinorAntigens_MHCII_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)



BESU1 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/BESU1_Allfiles_sorted.bam
BESU2 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/BESU2_Allfiles_sorted.bam

awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_BESU_Annotato_MinorAntigens_MHCII.txt | sort | uniq  > Summary_BESU_Annotato_MinorAntigens_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/BESU1_Allfiles_sorted.bam >> BESU_MHCII_Binders_Regions_Diagnosis.txt; done < Summary_BESU_Annotato_MinorAntigens_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/BESU2_Allfiles_sorted.bam >> BESU_MHCII_Binders_Regions_Relapse.txt; done < Summary_BESU_Annotato_MinorAntigens_MHCII_Regions.txt
 

besu_mhcii=read.table('Summary_BESU_Annotato_MinorAntigens_MHCII.txt')
 colnames(besu_mhcii)=c("chr","pos","id","peptide","hla","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50","ic50")
  besu_mhcii=merge(besu_mhcii,BESU_RNK,by="genes",all.x=T)
relapse=read.table('BESU_MHCII_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('BESU_MHCII_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
besu_mhcii=merge(besu_mhcii,relapse,by=c("chr","pos"),all.x=T)
besu_mhcii=merge(besu_mhcii,diagnosis,by=c("chr","pos"),all.x=T)
besu_mhcii$Rel_Mut=str_count(str_to_upper(besu_mhcii$rel_bases), as.character(besu_mhcii$mut_nt))
besu_mhcii$Rel_Wt=str_count(str_to_upper(besu_mhcii$rel_bases), as.character(besu_mhcii$wt_nt))
besu_mhcii$Dx_Mut=str_count(str_to_upper(besu_mhcii$dx_bases), as.character(besu_mhcii$mut_nt))
besu_mhcii$Dx_Wt=str_count(str_to_upper(besu_mhcii$dx_bases), as.character(besu_mhcii$wt_nt))
 write.table(unique(besu_mhcii[,c(1:18,27:30)]),"Summary_BESU_Annotato_MinorAntigens_MHCII_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)




CALU1 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/CALU1_Allfiles_sorted.bam
CALU2 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/CALU2_Allfiles_sorted.bam


awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_CALU_Annotato_MinorAntigens_MHCII.txt | sort | uniq  > Summary_CALU_Annotato_MinorAntigens_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/CALU1_Allfiles_sorted.bam >> CALU_MHCII_Binders_Regions_Diagnosis.txt; done < Summary_CALU_Annotato_MinorAntigens_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/CALU2_Allfiles_sorted.bam >> CALU_MHCII_Binders_Regions_Relapse.txt; done < Summary_CALU_Annotato_MinorAntigens_MHCII_Regions.txt



calu_mhcii=read.table('Summary_CALU_Annotato_MinorAntigens_MHCII.txt')
 colnames(calu_mhcii)=c("chr","pos","id","peptide","hla","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50","ic50")
  calu_mhcii=merge(calu_mhcii,CALU_RNK,by="genes",all.x=T)
relapse=read.table('CALU_MHCII_Binders_Regions_Relapse.txt',sep="\t")
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('CALU_MHCII_Binders_Regions_Diagnosis.txt',sep="\t")
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
calu_mhcii=merge(calu_mhcii,relapse,by=c("chr","pos"),all.x=T)
calu_mhcii=merge(calu_mhcii,diagnosis,by=c("chr","pos"),all.x=T)
calu_mhcii$Rel_Mut=str_count(str_to_upper(calu_mhcii$rel_bases), as.character(calu_mhcii$mut_nt))
calu_mhcii$Rel_Wt=str_count(str_to_upper(calu_mhcii$rel_bases), as.character(calu_mhcii$wt_nt))
calu_mhcii$Dx_Mut=str_count(str_to_upper(calu_mhcii$dx_bases), as.character(calu_mhcii$mut_nt))
calu_mhcii$Dx_Wt=str_count(str_to_upper(calu_mhcii$dx_bases), as.character(calu_mhcii$wt_nt))
 write.table(unique(calu_mhcii[,c(1:18,27:30)]),"Summary_CALU_Annotato_MinorAntigens_MHCII_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)




#DEIV
DEIV2 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2_Allfiles_sorted.bam
DEIV2b = /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2b_Allfiles_sorted.bam



awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_DEIV_Annotato_MinorAntigens_MHCII.txt | sort | uniq  > Summary_DEIV_Annotato_MinorAntigens_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2_Allfiles_sorted.bam >> DEIV_MHCII_Binders_Regions_Diagnosis.txt; done < Summary_DEIV_Annotato_MinorAntigens_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2b_Allfiles_sorted.bam >> DEIV_MHCII_Binders_Regions_Relapse.txt; done < Summary_DEIV_Annotato_MinorAntigens_MHCII_Regions.txt



deiv_mhcii=read.table('Summary_DEIV_Annotato_MinorAntigens_MHCII.txt')
 colnames(deiv_mhcii)=c("chr","pos","id","peptide","hla","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50","ic50")
  deiv_mhcii=merge(deiv_mhcii,DEIV_RNK,by="genes",all.x=T)
relapse=read.table('DEIV_MHCII_Binders_Regions_Relapse.txt',sep="\t")
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DEIV_MHCII_Binders_Regions_Diagnosis.txt',sep="\t")
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
deiv_mhcii=merge(deiv_mhcii,relapse,by=c("chr","pos"),all.x=T)
deiv_mhcii=merge(deiv_mhcii,diagnosis,by=c("chr","pos"),all.x=T)
deiv_mhcii$Rel_Mut=str_count(str_to_upper(deiv_mhcii$rel_bases), as.character(deiv_mhcii$mut_nt))
deiv_mhcii$Rel_Wt=str_count(str_to_upper(deiv_mhcii$rel_bases), as.character(deiv_mhcii$wt_nt))
deiv_mhcii$Dx_Mut=str_count(str_to_upper(deiv_mhcii$dx_bases), as.character(deiv_mhcii$mut_nt))
deiv_mhcii$Dx_Wt=str_count(str_to_upper(deiv_mhcii$dx_bases), as.character(deiv_mhcii$wt_nt))
 write.table(unique(deiv_mhcii[,c(1:18,27:30)]),"Summary_DEIV_Annotato_MinorAntigens_MHCII_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)



#DEST

 
DEST1 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST1_Allfiles_sorted.bam
DEST2 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST2_Allfiles_sorted.bam

awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_DEST_Annotato_MinorAntigens_MHCII.txt | sort | uniq  > Summary_DEST_Annotato_MinorAntigens_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST1_Allfiles_sorted.bam >> DEST_MHCII_Binders_Regions_Diagnosis.txt; done < Summary_DEST_Annotato_MinorAntigens_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST2_Allfiles_sorted.bam >> DEST_MHCII_Binders_Regions_Relapse.txt; done < Summary_DEST_Annotato_MinorAntigens_MHCII_Regions.txt




dest_mhcii=read.table('Summary_DEST_Annotato_MinorAntigens_MHCII.txt')
 colnames(dest_mhcii)=c("chr","pos","id","peptide","hla","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50","ic50")
  dest_mhcii=merge(dest_mhcii,DEST_RNK,by="genes",all.x=T)
relapse=read.table('DEST_MHCII_Binders_Regions_Relapse.txt',sep="\t")
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DEST_MHCII_Binders_Regions_Diagnosis.txt',sep="\t")
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
dest_mhcii=merge(dest_mhcii,relapse,by=c("chr","pos"),all.x=T)
dest_mhcii=merge(dest_mhcii,diagnosis,by=c("chr","pos"),all.x=T)
dest_mhcii$Rel_Mut=str_count(str_to_upper(dest_mhcii$rel_bases), as.character(dest_mhcii$mut_nt))
dest_mhcii$Rel_Wt=str_count(str_to_upper(dest_mhcii$rel_bases), as.character(dest_mhcii$wt_nt))
dest_mhcii$Dx_Mut=str_count(str_to_upper(dest_mhcii$dx_bases), as.character(dest_mhcii$mut_nt))
dest_mhcii$Dx_Wt=str_count(str_to_upper(dest_mhcii$dx_bases), as.character(dest_mhcii$wt_nt))
 write.table(unique(dest_mhcii[,c(1:18,27:30)]),"Summary_DEST_Annotato_MinorAntigens_MHCII_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)



 
 
# DR1_1 = /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DR1_1_Allfiles_sorted.bam
# DR1_2 = /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DR1_2_Allfiles_sorted.bam
# 
# 
# awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_ALFE_Annotato_MinorAntigens_MHCII.txt | sort | uniq  > Summary_ALFE_Annotato_MinorAntigens_MHCII_Regions.txt
#  while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE1_Allfiles_sorted.bam >> ALFE_MHCII_Binders_Regions_Diagnosis.txt; done < Summary_ALFE_Annotato_MinorAntigens_MHCII_Regions.txt
#  while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE2_Allfiles_sorted.bam >> ALFE_MHCII_Binders_Regions_Relapse.txt; done < Summary_ALFE_Annotato_MinorAntigens_MHCII_Regions.txt
#  
#  
DR41 =  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR41_Allfiles_sorted.bam
DR42 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR42_Allfiles_sorted.bam

 
 awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_DR4_Annotato_MinorAntigens_MHCII.txt | sort | uniq  > Summary_DR4_Annotato_MinorAntigens_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR41_Allfiles_sorted.bam>> DR4_MHCII_Binders_Regions_Diagnosis.txt; done < Summary_DR4_Annotato_MinorAntigens_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR42_Allfiles_sorted.bam >> DR4_MHCII_Binders_Regions_Relapse.txt; done < Summary_DR4_Annotato_MinorAntigens_MHCII_Regions.txt




dr4_mhcii=read.table('Summary_DR4_Annotato_MinorAntigens_MHCII.txt')
 colnames(dr4_mhcii)=c("chr","pos","id","peptide","hla","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50","ic50")
  dr4_mhcii=merge(dr4_mhcii,DR4_RNK,by="genes",all.x=T)
relapse=read.table('DR4_MHCII_Binders_Regions_Relapse.txt',sep="\t")
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DR4_MHCII_Binders_Regions_Diagnosis.txt',sep="\t")
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
dr4_mhcii=merge(dr4_mhcii,relapse,by=c("chr","pos"),all.x=T)
dr4_mhcii=merge(dr4_mhcii,diagnosis,by=c("chr","pos"),all.x=T)
dr4_mhcii$Rel_Mut=str_count(str_to_upper(dr4_mhcii$rel_bases), as.character(dr4_mhcii$mut_nt))
dr4_mhcii$Rel_Wt=str_count(str_to_upper(dr4_mhcii$rel_bases), as.character(dr4_mhcii$wt_nt))
dr4_mhcii$Dx_Mut=str_count(str_to_upper(dr4_mhcii$dx_bases), as.character(dr4_mhcii$mut_nt))
dr4_mhcii$Dx_Wt=str_count(str_to_upper(dr4_mhcii$dx_bases), as.character(dr4_mhcii$wt_nt))
 write.table(unique(dr4_mhcii[,c(1:18,27:30)]),"Summary_DR4_Annotato_MinorAntigens_MHCII_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)




 
DR51 =  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR51_Allfiles_sorted.bam
DR52 =  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR52_Allfiles_sorted.bam
 
 awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_DR5_Annotato_MinorAntigens_MHCII.txt | sort | uniq  > Summary_DR5_Annotato_MinorAntigens_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR51_Allfiles_sorted.bam >> DR5_MHCII_Binders_Regions_Diagnosis.txt; done < Summary_DR5_Annotato_MinorAntigens_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR52_Allfiles_sorted.bam >> DR5_MHCII_Binders_Regions_Relapse.txt; done < Summary_DR5_Annotato_MinorAntigens_MHCII_Regions.txt
 
dr5_mhcii=read.table('Summary_DR5_Annotato_MinorAntigens_MHCII.txt')
 colnames(dr5_mhcii)=c("chr","pos","id","peptide","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50")
  dr5_mhcii=merge(dr5_mhcii,DR5_RNK,by="genes",all.x=T)
relapse=read.table('DR5_MHCII_Binders_Regions_Relapse.txt',sep="\t")
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DR5_MHCII_Binders_Regions_Diagnosis.txt',sep="\t")
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
dr5_mhcii=merge(dr5_mhcii,relapse,by=c("chr","pos"),all.x=T)
dr5_mhcii=merge(dr5_mhcii,diagnosis,by=c("chr","pos"),all.x=T)
dr5_mhcii$Rel_Mut=str_count(str_to_upper(dr5_mhcii$rel_bases), as.character(dr5_mhcii$mut_nt))
dr5_mhcii$Rel_Wt=str_count(str_to_upper(dr5_mhcii$rel_bases), as.character(dr5_mhcii$wt_nt))
dr5_mhcii$Dx_Mut=str_count(str_to_upper(dr5_mhcii$dx_bases), as.character(dr5_mhcii$mut_nt))
dr5_mhcii$Dx_Wt=str_count(str_to_upper(dr5_mhcii$dx_bases), as.character(dr5_mhcii$wt_nt))
 write.table(unique(dr5_mhcii[,c(1:16,25:28)]),"Summary_DR5_Annotato_MinorAntigens_MHCII_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)






FOCA1 =  /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/FOCA1_Allfiles_sorted.bam
FOCA2 =  /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/FOCA2_Allfiles_sorted.bam

 
 awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_FOCA_Annotato_MinorAntigens_MHCII.txt | sort | uniq  > Summary_FOCA_Annotato_MinorAntigens_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/FOCA1_Allfiles_sorted.bam >> FOCA_MHCII_Binders_Regions_Diagnosis.txt; done < Summary_FOCA_Annotato_MinorAntigens_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/FOCA2_Allfiles_sorted.bam >> FOCA_MHCII_Binders_Regions_Relapse.txt; done < Summary_FOCA_Annotato_MinorAntigens_MHCII_Regions.txt
 
foca_mhcii=read.table('Summary_FOCA_Annotato_MinorAntigens_MHCII.txt')
 colnames(foca_mhcii)=c("chr","pos","id","peptide","hla","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50","ic50")
  foca_mhcii=merge(foca_mhcii,FOCA_RNK,by="genes",all.x=T)
relapse=read.table('FOCA_MHCII_Binders_Regions_Relapse.txt',sep="\t")
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('FOCA_MHCII_Binders_Regions_Diagnosis.txt',sep="\t")
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
foca_mhcii=merge(foca_mhcii,relapse,by=c("chr","pos"),all.x=T)
foca_mhcii=merge(foca_mhcii,diagnosis,by=c("chr","pos"),all.x=T)
foca_mhcii$Rel_Mut=str_count(str_to_upper(foca_mhcii$rel_bases), as.character(foca_mhcii$mut_nt))
foca_mhcii$Rel_Wt=str_count(str_to_upper(foca_mhcii$rel_bases), as.character(foca_mhcii$wt_nt))
foca_mhcii$Dx_Mut=str_count(str_to_upper(foca_mhcii$dx_bases), as.character(foca_mhcii$mut_nt))
foca_mhcii$Dx_Wt=str_count(str_to_upper(foca_mhcii$dx_bases), as.character(foca_mhcii$wt_nt))
 write.table(unique(foca_mhcii[,c(1:18,27:30)]),"Summary_FOCA_Annotato_MinorAntigens_MHCII_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)

 
 


GAGRA1 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/GAGRA1_Allfiles_sorted.bam
GAGRA2 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/GAGRA2_Allfiles_sorted.bam

 
 awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_GAGRA_Annotato_MinorAntigens_MHCII.txt | sort | uniq  > Summary_GAGRA_Annotato_MinorAntigens_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/GAGRA1_Allfiles_sorted.bam >> GAGRA_MHCII_Binders_Regions_Diagnosis.txt; done < Summary_GAGRA_Annotato_MinorAntigens_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/GAGRA2_Allfiles_sorted.bam >> GAGRA_MHCII_Binders_Regions_Relapse.txt; done < Summary_GAGRA_Annotato_MinorAntigens_MHCII_Regions.txt



gagra_mhcii=read.table('Summary_GAGRA_Annotato_MinorAntigens_MHCII.txt')
 colnames(gagra_mhcii)=c("chr","pos","id","peptide","hla","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50","ic50")
  gagra_mhcii=merge(gagra_mhcii,GAGRA_RNK,by="genes",all.x=T)
relapse=read.table('GAGRA_MHCII_Binders_Regions_Relapse.txt',sep="\t")
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('GAGRA_MHCII_Binders_Regions_Diagnosis.txt',sep="\t")
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
gagra_mhcii=merge(gagra_mhcii,relapse,by=c("chr","pos"),all.x=T)
gagra_mhcii=merge(gagra_mhcii,diagnosis,by=c("chr","pos"),all.x=T)
gagra_mhcii$Rel_Mut=str_count(str_to_upper(gagra_mhcii$rel_bases), as.character(gagra_mhcii$mut_nt))
gagra_mhcii$Rel_Wt=str_count(str_to_upper(gagra_mhcii$rel_bases), as.character(gagra_mhcii$wt_nt))
gagra_mhcii$Dx_Mut=str_count(str_to_upper(gagra_mhcii$dx_bases), as.character(gagra_mhcii$mut_nt))
gagra_mhcii$Dx_Wt=str_count(str_to_upper(gagra_mhcii$dx_bases), as.character(gagra_mhcii$wt_nt))
 write.table(unique(gagra_mhcii[,c(1:18,27:30)]),"Summary_GAGRA_Annotato_MinorAntigens_MHCII_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)

 
 
 
 
LUAN1 = /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN1_Allfiles_sorted.bam
LUAN2 = /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN2_Allfiles_sorted.bam

 
 awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_LUAN_Annotato_MinorAntigens_MHCII.txt | sort | uniq  > Summary_LUAN_Annotato_MinorAntigens_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN1_Allfiles_sorted.bam >> LUAN_MHCII_Binders_Regions_Diagnosis.txt; done < Summary_LUAN_Annotato_MinorAntigens_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN2_Allfiles_sorted.bam >> LUAN_MHCII_Binders_Regions_Relapse.txt; done < Summary_LUAN_Annotato_MinorAntigens_MHCII_Regions.txt
 


luan_mhcii=read.table('Summary_LUAN_Annotato_MinorAntigens_MHCII.txt')
 colnames(luan_mhcii)=c("chr","pos","id","peptide","hla","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50","ic50")
  luan_mhcii=merge(luan_mhcii,LUAN_RNK,by="genes",all.x=T)
relapse=read.table('LUAN_MHCII_Binders_Regions_Relapse.txt',sep="\t")
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('LUAN_MHCII_Binders_Regions_Diagnosis.txt',sep="\t")
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
luan_mhcii=merge(luan_mhcii,relapse,by=c("chr","pos"),all.x=T)
luan_mhcii=merge(luan_mhcii,diagnosis,by=c("chr","pos"),all.x=T)
luan_mhcii$Rel_Mut=str_count(str_to_upper(luan_mhcii$rel_bases), as.character(luan_mhcii$mut_nt))
luan_mhcii$Rel_Wt=str_count(str_to_upper(luan_mhcii$rel_bases), as.character(luan_mhcii$wt_nt))
luan_mhcii$Dx_Mut=str_count(str_to_upper(luan_mhcii$dx_bases), as.character(luan_mhcii$mut_nt))
luan_mhcii$Dx_Wt=str_count(str_to_upper(luan_mhcii$dx_bases), as.character(luan_mhcii$wt_nt))
 write.table(unique(luan_mhcii[,c(1:18,27:30)]),"Summary_LUAN_Annotato_MinorAntigens_MHCII_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)

 
 
 
MABI1 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/MABI1_Allfiles_sorted.bam
MABI2 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/MABI2_Allfiles_sorted.bam

 
 awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_MABI_Annotato_MinorAntigens_MHCII.txt | sort | uniq  > Summary_MABI_Annotato_MinorAntigens_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/MABI1_Allfiles_sorted.bam >> MABI_MHCII_Binders_Regions_Diagnosis.txt; done < Summary_MABI_Annotato_MinorAntigens_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/MABI2_Allfiles_sorted.bam >> MABI_MHCII_Binders_Regions_Relapse.txt; done < Summary_MABI_Annotato_MinorAntigens_MHCII_Regions.txt
 

mabi_mhcii=read.table('Summary_MABI_Annotato_MinorAntigens_MHCII.txt')
 colnames(mabi_mhcii)=c("chr","pos","id","peptide","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50")
  mabi_mhcii=merge(mabi_mhcii,MABI_RNK,by="genes",all.x=T)
relapse=read.table('MABI_MHCII_Binders_Regions_Relapse.txt',sep="\t")
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('MABI_MHCII_Binders_Regions_Diagnosis.txt',sep="\t")
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
mabi_mhcii=merge(mabi_mhcii,relapse,by=c("chr","pos"),all.x=T)
mabi_mhcii=merge(mabi_mhcii,diagnosis,by=c("chr","pos"),all.x=T)
mabi_mhcii$Rel_Mut=str_count(str_to_upper(mabi_mhcii$rel_bases), as.character(mabi_mhcii$mut_nt))
mabi_mhcii$Rel_Wt=str_count(str_to_upper(mabi_mhcii$rel_bases), as.character(mabi_mhcii$wt_nt))
mabi_mhcii$Dx_Mut=str_count(str_to_upper(mabi_mhcii$dx_bases), as.character(mabi_mhcii$mut_nt))
mabi_mhcii$Dx_Wt=str_count(str_to_upper(mabi_mhcii$dx_bases), as.character(mabi_mhcii$wt_nt))
write.table(unique(mabi_mhcii[,c(1:16,25:28)]),"Summary_MABI_Annotato_MinorAntigens_MHCII_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)



MOGE1 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE1_Allfiles_sorted.bam
MOGE2 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE2_Allfiles_sorted.bam

 
 awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_MOGE_Annotato_MinorAntigens_MHCII.txt | sort | uniq  > Summary_MOGE_Annotato_MinorAntigens_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE1_Allfiles_sorted.bam >> MOGE_MHCII_Binders_Regions_Diagnosis.txt; done < Summary_MOGE_Annotato_MinorAntigens_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE2_Allfiles_sorted.bam >> MOGE_MHCII_Binders_Regions_Relapse.txt; done < Summary_MOGE_Annotato_MinorAntigens_MHCII_Regions.txt
 

moge_mhcii=read.table('Summary_MOGE_Annotato_MinorAntigens_MHCII.txt')
 colnames(moge_mhcii)=c("chr","pos","id","peptide","hla","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50","ic50")
  moge_mhcii=merge(moge_mhcii,MOGE_RNK,by="genes",all.x=T)
relapse=read.table('MOGE_MHCII_Binders_Regions_Relapse.txt',sep="\t")
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('MOGE_MHCII_Binders_Regions_Diagnosis.txt',sep="\t")
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
moge_mhcii=merge(moge_mhcii,relapse,by=c("chr","pos"),all.x=T)
moge_mhcii=merge(moge_mhcii,diagnosis,by=c("chr","pos"),all.x=T)
moge_mhcii$Rel_Mut=str_count(str_to_upper(moge_mhcii$rel_bases), as.character(moge_mhcii$mut_nt))
moge_mhcii$Rel_Wt=str_count(str_to_upper(moge_mhcii$rel_bases), as.character(moge_mhcii$wt_nt))
moge_mhcii$Dx_Mut=str_count(str_to_upper(moge_mhcii$dx_bases), as.character(moge_mhcii$mut_nt))
moge_mhcii$Dx_Wt=str_count(str_to_upper(moge_mhcii$dx_bases), as.character(moge_mhcii$wt_nt))
write.table(unique(moge_mhcii[,c(1:18,27:30)]),"Summary_MOGE_Annotato_MinorAntigens_MHCII_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)


 
PIAG1 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/PIAG1_Allfiles_sorted.bam
PIAG2 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/PIAG2_Allfiles_sorted.bam



 
 awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_PIAG_Annotato_MinorAntigens_MHCII.txt | sort | uniq  > Summary_PIAG_Annotato_MinorAntigens_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/PIAG1_Allfiles_sorted.bam >> PIAG_MHCII_Binders_Regions_Diagnosis.txt; done < Summary_PIAG_Annotato_MinorAntigens_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/PIAG2_Allfiles_sorted.bam >> PIAG_MHCII_Binders_Regions_Relapse.txt; done < Summary_PIAG_Annotato_MinorAntigens_MHCII_Regions.txt


piag_mhcii=read.table('Summary_PIAG_Annotato_MinorAntigens_MHCII.txt')
 colnames(piag_mhcii)=c("chr","pos","id","peptide","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50")
  piag_mhcii=merge(piag_mhcii,PIAG_RNK,by="genes",all.x=T)
relapse=read.table('PIAG_MHCII_Binders_Regions_Relapse.txt',sep="\t")
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('PIAG_MHCII_Binders_Regions_Diagnosis.txt',sep="\t")
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
piag_mhcii=merge(piag_mhcii,relapse,by=c("chr","pos"),all.x=T)
piag_mhcii=merge(piag_mhcii,diagnosis,by=c("chr","pos"),all.x=T)
piag_mhcii$Rel_Mut=str_count(str_to_upper(piag_mhcii$rel_bases), as.character(piag_mhcii$mut_nt))
piag_mhcii$Rel_Wt=str_count(str_to_upper(piag_mhcii$rel_bases), as.character(piag_mhcii$wt_nt))
piag_mhcii$Dx_Mut=str_count(str_to_upper(piag_mhcii$dx_bases), as.character(piag_mhcii$mut_nt))
piag_mhcii$Dx_Wt=str_count(str_to_upper(piag_mhcii$dx_bases), as.character(piag_mhcii$wt_nt))
write.table(unique(piag_mhcii[,c(1:16,25:28)]),"Summary_PIAG_Annotato_MinorAntigens_MHCII_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)




 
PRELU1 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/PRELU1_Allfiles_sorted.bam
PRELU2 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/PRELU2_Allfiles_sorted.bam

 
 awk '{OFS="\t"; print $1":"$2"-"$2}'   Summary_PRELU_Annotato_MinorAntigens_MHCII.txt | sort | uniq  > Summary_PRELU_Annotato_MinorAntigens_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/PRELU1_Allfiles_sorted.bam >> PRELU_MHCII_Binders_Regions_Diagnosis.txt; done < Summary_PRELU_Annotato_MinorAntigens_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/PRELU2_Allfiles_sorted.bam >> PRELU_MHCII_Binders_Regions_Relapse.txt; done < Summary_PRELU_Annotato_MinorAntigens_MHCII_Regions.txt
 

prelu_mhcii=read.table('Summary_PRELU_Annotato_MinorAntigens_MHCII.txt')
 colnames(prelu_mhcii)=c("chr","pos","id","peptide","hla","hla","effect","genes","mut_nt","wt_nt","relapse_tx","diagnosis_tx","ic50","ic50")
  prelu_mhcii=merge(prelu_mhcii,PRELU_RNK,by="genes",all.x=T)
relapse=read.table('PRELU_MHCII_Binders_Regions_Relapse.txt',sep="\t")
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('PRELU_MHCII_Binders_Regions_Diagnosis.txt',sep="\t")
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
prelu_mhcii=merge(prelu_mhcii,relapse,by=c("chr","pos"),all.x=T)
prelu_mhcii=merge(prelu_mhcii,diagnosis,by=c("chr","pos"),all.x=T)
prelu_mhcii$Rel_Mut=str_count(str_to_upper(prelu_mhcii$rel_bases), as.character(prelu_mhcii$mut_nt))
prelu_mhcii$Rel_Wt=str_count(str_to_upper(prelu_mhcii$rel_bases), as.character(prelu_mhcii$wt_nt))
prelu_mhcii$Dx_Mut=str_count(str_to_upper(prelu_mhcii$dx_bases), as.character(prelu_mhcii$mut_nt))
prelu_mhcii$Dx_Wt=str_count(str_to_upper(prelu_mhcii$dx_bases), as.character(prelu_mhcii$wt_nt))
write.table(unique(prelu_mhcii[,c(1:18,27:30)]),"Summary_PRELU_Annotato_MinorAntigens_MHCII_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)



#mhAGs Classe II

#ALFE
mhag=read.table('Summary_ALFE_Annotato_MinorAntigens_MHCII_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$ALFE2_kallisto[is.na(mhag$ALFE2_kallisto)] <- 0
mhag$ALFE1_kallisto[is.na(mhag$ALFE1_kallisto)] <- 0

mhag=mhag[(mhag$ALFE2_kallisto>=0.58 & mhag$ALFE1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
6462


library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,ALFE2_kallisto,ALFE1_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,ALFE2_kallisto,ALFE1_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,ALFE2_kallisto,ALFE1_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,ALFE2_kallisto,ALFE1_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)

nrow(result)
	
table(result$ALFE2_kallisto<result$ALFE1_kallisto)

table(result$ALFE2_kallisto>result$ALFE1_kallisto)

table(result$Rel_Mut<result$Dx_Mut)


table(result$Rel_Mut>result$Dx_Mut)







#A lvello di trascritto
 table(mhag$ALFE2_kallisto > mhag$ALFE1_kallisto)

FALSE  TRUE
 3380  2285


 wilcox.test(mhag$ALFE2_kallisto, mhag$ALFE1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$ALFE2_kallisto and mhag$ALFE1_kallisto
W = 15209356, p-value = 1
alternative hypothesis: true location shift is greater than 0



#A livello di trascritto
 table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
 4315  1181

table(mhag$Rel_Mut<mhag$Dx_Mut)

FALSE  TRUE
 1331  4165


 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 10810557, p-value = 1
alternative hypothesis: true location shift is greater than 0



 
#BESU
mhag=read.table('Summary_BESU_Annotato_MinorAntigens_MHCII_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$BESU2_kallisto[is.na(mhag$BESU2_kallisto)] <- 0
mhag$BESU1_kallisto[is.na(mhag$BESU1_kallisto)] <- 0
mhag$Rel_Mut[is.na(mhag$Rel_Mut)] <- 0
mhag$Dx_Mut[is.na(mhag$Dx_Mut)] <- 0

mhag=mhag[(mhag$BESU2_kallisto>=0.58 & mhag$BESU1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
28986


library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,BESU2_kallisto,BESU1_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,BESU2_kallisto,BESU1_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,BESU2_kallisto,BESU1_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,BESU2_kallisto,BESU1_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)

nrow(result)
	
table(result$BESU2_kallisto<result$BESU1_kallisto)

table(result$BESU2_kallisto>result$BESU1_kallisto)

table(result$Rel_Mut<result$Dx_Mut)


table(result$Rel_Mut>result$Dx_Mut)





#A livello di trascritto
 table(mhag$BESU2_kallisto > mhag$BESU1_kallisto)

FALSE  TRUE
10640 15331



 wilcox.test(mhag$BESU2_kallisto, mhag$BESU1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$BESU2_kallisto and mhag$BESU1_kallisto
W = 344963810, p-value = 3.141e-06
alternative hypothesis: true location shift is greater than 0




#A livello di trascritto
 table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
10634 14897



 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 354436363, p-value < 2.2e-16
alternative hypothesis: true location shift is greater than 0


#CALU
mhag=read.table('Summary_CALU_Annotato_MinorAntigens_MHCII_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$CALU2_kallisto[is.na(mhag$CALU2_kallisto)] <- 0
mhag$CALU1_kallisto[is.na(mhag$CALU1_kallisto)] <- 0
mhag$Rel_Mut[is.na(mhag$Rel_Mut)] <- 0
mhag$Dx_Mut[is.na(mhag$Dx_Mut)] <- 0

mhag=mhag[(mhag$CALU2_kallisto>=0.58 & mhag$CALU1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
3864



library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,CALU2_kallisto,CALU1_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,CALU2_kallisto,CALU1_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,CALU2_kallisto,CALU1_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,CALU2_kallisto,CALU1_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)

nrow(result)
	
table(result$CALU2_kallisto<result$CALU1_kallisto)

table(result$CALU2_kallisto>result$CALU1_kallisto)

table(result$Rel_Mut<result$Dx_Mut)


table(result$Rel_Mut>result$Dx_Mut)




#A livello di trascritto
 table(mhag$CALU2_kallisto > mhag$CALU1_kallisto)

FALSE  TRUE
 1633  1698




 wilcox.test(mhag$CALU2_kallisto, mhag$CALU1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$CALU2_kallisto and mhag$CALU1_kallisto
W = 5596524, p-value = 0.2673
alternative hypothesis: true location shift is greater than 0




#A livello di trascritto
 table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
 2722    24

table(mhag$Rel_Mut<mhag$Dx_Mut)

FALSE  TRUE
 2724    22


 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 4547278, p-value = 0.6499
alternative hypothesis: true location shift is greater than 0




#DEIV
mhag=read.table('Summary_DEIV_Annotato_MinorAntigens_MHCII_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$DEIV2b_kallisto[is.na(mhag$DEIV2b_kallisto)] <- 0
mhag$DEIV2_kallisto[is.na(mhag$DEIV2_kallisto)] <- 0
mhag$Rel_Mut[is.na(mhag$Rel_Mut)] <- 0
mhag$Dx_Mut[is.na(mhag$Dx_Mut)] <- 0

mhag=mhag[(mhag$DEIV2b_kallisto>=0.58 & mhag$DEIV2_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
30677


library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,DEIV2b_kallisto,DEIV2_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,DEIV2b_kallisto,DEIV2_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,DEIV2b_kallisto,DEIV2_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,DEIV2b_kallisto,DEIV2_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)

nrow(result)
	
table(result$DEIV2b_kallisto<result$DEIV2_kallisto)

table(result$DEIV2b_kallisto>result$DEIV2_kallisto)

table(result$Rel_Mut<result$Dx_Mut)


table(result$Rel_Mut>result$Dx_Mut)





#A livello di trascritto
table(mhag$DEIV2b_kallisto > mhag$DEIV2_kallisto)

FALSE  TRUE
 6387  9564

table(mhag$DEIV2b_kallisto < mhag$DEIV2_kallisto)

FALSE  TRUE
 9564  6387




wilcox.test(mhag$DEIV2b_kallisto, mhag$DEIV2_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$DEIV2b_kallisto and mhag$DEIV2_kallisto
W = 137509365, p-value < 2.2e-16
alternative hypothesis: true location shift is greater than 0



#A livello di trascritto
table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
 1528  4960


table(mhag$Rel_Mut<mhag$Dx_Mut)

FALSE  TRUE
 5269  1219



wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 79189336, p-value < 2.2e-16
alternative hypothesis: true location shift is greater than 0

#DEST
mhag=read.table('Summary_DEST_Annotato_MinorAntigens_MHCII_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$DEST2_kallisto[is.na(mhag$DEST2_kallisto)] <- 0
mhag$DEST1_kallisto[is.na(mhag$DEST1_kallisto)] <- 0
mhag$Rel_Mut[is.na(mhag$Rel_Mut)] <- 0
mhag$Dx_Mut[is.na(mhag$Dx_Mut)] <- 0

mhag=mhag[(mhag$DEST2_kallisto>=0.58 & mhag$DEST1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
5675


library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,DEST2_kallisto,DEST1_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,DEST2_kallisto,DEST1_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,DEST2_kallisto,DEST1_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,DEST2_kallisto,DEST1_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)

nrow(result)
	
table(result$DEST2_kallisto<result$DEST1_kallisto)

table(result$DEST2_kallisto>result$DEST1_kallisto)

table(result$Rel_Mut<result$Dx_Mut)


table(result$Rel_Mut>result$Dx_Mut)







#A livello di trascritto
 table(mhag$DEST2_kallisto > mhag$DEST1_kallisto)

FALSE  TRUE
 1555  1808


wilcox.test(mhag$DEST2_kallisto, mhag$DEST1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$DEST2_kallisto and mhag$DEST1_kallisto
W = 5756134, p-value = 0.1018
alternative hypothesis: true location shift is greater than 0


#A livello di trascritto
table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
  271   492


 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 2011216, p-value = 0.1825
alternative hypothesis: true location shift is greater than 0

#DR4
mhag=read.table('Summary_DR4_Annotato_MinorAntigens_MHCII_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$DR42_kallisto[is.na(mhag$DR42_kallisto)] <- 0
mhag$DR41_kallisto[is.na(mhag$DR41_kallisto)] <- 0
mhag$Rel_Mut[is.na(mhag$Rel_Mut)] <- 0
mhag$Dx_Mut[is.na(mhag$Dx_Mut)] <- 0
mhag=mhag[(mhag$DR42_kallisto>=0.58 & mhag$DR41_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
5681




library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,DR42_kallisto,DR41_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,DR42_kallisto,DR41_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,DR42_kallisto,DR41_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,DR42_kallisto,DR41_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)

nrow(result)
	
table(result$DR42_kallisto<result$DR41_kallisto)

table(result$DR42_kallisto>result$DR41_kallisto)

table(result$Rel_Mut<result$Dx_Mut)


table(result$Rel_Mut>result$Dx_Mut)





#A livello di trascritto
table(mhag$DR42_kallisto > mhag$DR41_kallisto)

FALSE  TRUE
 2110  2361



wilcox.test(mhag$DR42_kallisto, mhag$DR41_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$DR42_kallisto and mhag$DR41_kallisto
W = 9872534, p-value = 0.842
alternative hypothesis: true location shift is greater than 0



#A livello di trascritto
table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
 2429  1815

wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 8837292, p-value = 1
alternative hypothesis: true location shift is greater than 0


#DR5
mhag=read.table('Summary_DR5_Annotato_MinorAntigens_MHCII_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$DR52_kallisto[is.na(mhag$DR52_kallisto)] <- 0
mhag$DR51_kallisto[is.na(mhag$DR51_kallisto)] <- 0
mhag$Rel_Mut[is.na(mhag$Rel_Mut)] <- 0
mhag$Dx_Mut[is.na(mhag$Dx_Mut)] <- 0
mhag=mhag[(mhag$DR52_kallisto>=0.58 & mhag$DR51_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
4200



library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,logFC,Qvalue,DR52_kallisto,DR51_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,logFC,Qvalue,DR52_kallisto,DR51_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,logFC,Qvalue,DR52_kallisto,DR51_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,logFC,Qvalue,DR52_kallisto,DR51_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)

nrow(result)
	
table(result$DR52_kallisto<result$DR51_kallisto)

table(result$DR52_kallisto>result$DR51_kallisto)

table(result$Rel_Mut<result$Dx_Mut)


table(result$Rel_Mut>result$Dx_Mut)




#A livello di trascritto
table(mhag$DR52_kallisto > mhag$DR51_kallisto)

FALSE  TRUE
 1696  1623


wilcox.test(mhag$DR52_kallisto, mhag$DR51_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$DR52_kallisto and mhag$DR51_kallisto
W = 5158598, p-value = 1
alternative hypothesis: true location shift is greater than 0




#A livello di trascritto
 table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
 1064  1645


 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 5146724, p-value < 2.2e-16
alternative hypothesis: true location shift is greater than 0



#FOCA
mhag=read.table('Summary_FOCA_Annotato_MinorAntigens_MHCII_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$FOCA2_kallisto[is.na(mhag$FOCA2_kallisto)] <- 0
mhag$FOCA1_kallisto[is.na(mhag$FOCA1_kallisto)] <- 0
mhag$Rel_Mut[is.na(mhag$Rel_Mut)] <- 0
mhag$Dx_Mut[is.na(mhag$Dx_Mut)] <- 0
mhag=mhag[(mhag$FOCA2_kallisto>=0.58 & mhag$FOCA1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
13563


library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,FOCA2_kallisto,FOCA1_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,FOCA2_kallisto,FOCA1_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,FOCA2_kallisto,FOCA1_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,FOCA2_kallisto,FOCA1_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)

nrow(result)
	
table(result$FOCA2_kallisto<result$FOCA1_kallisto)

table(result$FOCA2_kallisto>result$FOCA1_kallisto)

table(result$Rel_Mut<result$Dx_Mut)


table(result$Rel_Mut>result$Dx_Mut)




#A livello di trascritto
table(mhag$FOCA2_kallisto > mhag$FOCA1_kallisto)

FALSE  TRUE
 4406  5077


 wilcox.test(mhag$FOCA2_kallisto, mhag$FOCA1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$FOCA2_kallisto and mhag$FOCA1_kallisto
W = 46236940, p-value = 0.0003658
alternative hypothesis: true location shift is greater than 0

#A livello di trascritto
 table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
 2033  2125



 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 19511078, p-value = 0.9591
alternative hypothesis: true location shift is greater than 0


#GAGRA
mhag=read.table('Summary_GAGRA_Annotato_MinorAntigens_MHCII_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$GAGRA2_kallisto[is.na(mhag$GAGRA2_kallisto)] <- 0
mhag$GAGRA1_kallisto[is.na(mhag$GAGRA1_kallisto)] <- 0
mhag$Rel_Mut[is.na(mhag$Rel_Mut)] <- 0
mhag$Dx_Mut[is.na(mhag$Dx_Mut)] <- 0
mhag=mhag[(mhag$GAGRA2_kallisto>=0.58 & mhag$GAGRA1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
22628



library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,GAGRA2_kallisto,GAGRA1_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,GAGRA2_kallisto,GAGRA1_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,GAGRA2_kallisto,GAGRA1_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,GAGRA2_kallisto,GAGRA1_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)

nrow(result)
	
table(result$GAGRA2_kallisto<result$GAGRA1_kallisto)

table(result$GAGRA2_kallisto>result$GAGRA1_kallisto)

table(result$Rel_Mut<result$Dx_Mut)


table(result$Rel_Mut>result$Dx_Mut)






#A livello di trascritto
table(mhag$GAGRA2_kallisto > mhag$GAGRA1_kallisto)

FALSE  TRUE
 5985 10004


wilcox.test(mhag$GAGRA2_kallisto, mhag$GAGRA1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$GAGRA2_kallisto and mhag$GAGRA1_kallisto
W = 140314548, p-value < 2.2e-16
alternative hypothesis: true location shift is greater than 0


#A livello di trascritto
table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
 1170  1139

table(mhag$Rel_Mut<mhag$Dx_Mut)

FALSE  TRUE
 1334   975


 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 25174962, p-value = 0.4445
alternative hypothesis: true location shift is greater than 0


#LUAN
mhag=read.table('Summary_LUAN_Annotato_MinorAntigens_MHCII_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$LUAN2_kallisto[is.na(mhag$LUAN2_kallisto)] <- 0
mhag$LUAN1_kallisto[is.na(mhag$LUAN1_kallisto)] <- 0
mhag$Rel_Mut[is.na(mhag$Rel_Mut)] <- 0
mhag$Dx_Mut[is.na(mhag$Dx_Mut)] <- 0
mhag=mhag[(mhag$LUAN2_kallisto>=0.58 & mhag$LUAN1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
4242



library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,LUAN2_kallisto,LUAN1_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,LUAN2_kallisto,LUAN1_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,LUAN2_kallisto,LUAN1_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,LUAN2_kallisto,LUAN1_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)

nrow(result)
	
table(result$LUAN2_kallisto<result$LUAN1_kallisto)

table(result$LUAN2_kallisto>result$LUAN1_kallisto)

table(result$Rel_Mut<result$Dx_Mut)


table(result$Rel_Mut>result$Dx_Mut)





#A livello di trascritto
 table(mhag$LUAN2_kallisto > mhag$LUAN1_kallisto)

FALSE  TRUE
 2297  1177



 wilcox.test(mhag$LUAN2_kallisto, mhag$LUAN1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$LUAN2_kallisto and mhag$LUAN1_kallisto
W = 5894897, p-value = 0.9524
alternative hypothesis: true location shift is greater than 0



#A livello di trascritto
 table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
 1132   712

 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 3238198, p-value = 0.2545
alternative hypothesis: true location shift is greater than 0


#MABI
mhag=read.table('Summary_MABI_Annotato_MinorAntigens_MHCII_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$MABI2_kallisto[is.na(mhag$MABI2_kallisto)] <- 0
mhag$MABI1_kallisto[is.na(mhag$MABI1_kallisto)] <- 0
mhag$Rel_Mut[is.na(mhag$Rel_Mut)] <- 0
mhag$Dx_Mut[is.na(mhag$Dx_Mut)] <- 0
mhag=mhag[(mhag$MABI2_kallisto>=0.58 & mhag$MABI1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
2801	




library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,logFC,Qvalue,MABI2_kallisto,MABI1_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,logFC,Qvalue,MABI2_kallisto,MABI1_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,logFC,Qvalue,MABI2_kallisto,MABI1_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,logFC,Qvalue,MABI2_kallisto,MABI1_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)

nrow(result)
	
table(result$MABI2_kallisto<result$MABI1_kallisto)

table(result$MABI2_kallisto>result$MABI1_kallisto)

table(result$Rel_Mut<result$Dx_Mut)


table(result$Rel_Mut>result$Dx_Mut)





#A livello di trascritto
 table(mhag$MABI2_kallisto > mhag$MABI1_kallisto)

FALSE  TRUE
 1737  1064


 wilcox.test(mhag$MABI2_kallisto, mhag$MABI1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$MABI2_kallisto and mhag$MABI1_kallisto
W = 3513375, p-value = 1
alternative hypothesis: true location shift is greater than 0



#A livello di trascritto
 table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
 1513  1153

 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 3556325, p-value = 0.9989
alternative hypothesis: true location shift is greater than 0

#MOGE
mhag=read.table('Summary_MOGE_Annotato_MinorAntigens_MHCII_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$MOGE2_kallisto[is.na(mhag$MOGE2_kallisto)] <- 0
mhag$MOGE1_kallisto[is.na(mhag$MOGE1_kallisto)] <- 0
mhag$Rel_Mut[is.na(mhag$Rel_Mut)] <- 0
mhag$Dx_Mut[is.na(mhag$Dx_Mut)] <- 0
mhag=mhag[(mhag$MOGE2_kallisto>=0.58 & mhag$MOGE1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
8838


library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,MOGE2_kallisto,MOGE1_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,MOGE2_kallisto,MOGE1_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,MOGE2_kallisto,MOGE1_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,MOGE2_kallisto,MOGE1_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)

nrow(result)
	
table(result$MOGE2_kallisto<result$MOGE1_kallisto)

table(result$MOGE2_kallisto>result$MOGE1_kallisto)

table(result$Rel_Mut<result$Dx_Mut)


table(result$Rel_Mut>result$Dx_Mut)





#A livello di trascritto
table(mhag$MOGE2_kallisto > mhag$MOGE1_kallisto)

FALSE  TRUE
 2586  2401


wilcox.test(mhag$MOGE2_kallisto, mhag$MOGE1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$MOGE2_kallisto and mhag$MOGE1_kallisto
W = 12578836, p-value = 0.1587
alternative hypothesis: true location shift is greater than 0



#A livello di trascritto
 table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
 1650  1698


wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 8980556, p-value = 6.103e-09
alternative hypothesis: true location shift is greater than 0

#PIAG
mhag=read.table('Summary_PIAG_Annotato_MinorAntigens_MHCII_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$PIAG2_kallisto[is.na(mhag$PIAG2_kallisto)] <- 0
mhag$PIAG1_kallisto[is.na(mhag$PIAG1_kallisto)] <- 0
mhag$Rel_Mut[is.na(mhag$Rel_Mut)] <- 0
mhag$Dx_Mut[is.na(mhag$Dx_Mut)] <- 0
mhag=mhag[(mhag$PIAG2_kallisto>=0.58 & mhag$PIAG1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
12791


library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,logFC,Qvalue,PIAG2_kallisto,PIAG1_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,logFC,Qvalue,PIAG2_kallisto,PIAG1_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,logFC,Qvalue,PIAG2_kallisto,PIAG1_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,logFC,Qvalue,PIAG2_kallisto,PIAG1_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)

nrow(result)
	
table(result$PIAG2_kallisto<result$PIAG1_kallisto)

table(result$PIAG2_kallisto>result$PIAG1_kallisto)

table(result$Rel_Mut<result$Dx_Mut)


table(result$Rel_Mut>result$Dx_Mut)



#A livello di trascritto
 table(mhag$PIAG2_kallisto > mhag$PIAG1_kallisto)

FALSE  TRUE
 5278  3802



wilcox.test(mhag$PIAG2_kallisto, mhag$PIAG1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$PIAG2_kallisto and mhag$PIAG1_kallisto
W = 38607283, p-value = 1
alternative hypothesis: true location shift is greater than 0


#A livello di trascritto
table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
 4610  1795



 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 24979378, p-value = 1
alternative hypothesis: true location shift is greater than 0

#PRELU
mhag=read.table('Summary_PRELU_Annotato_MinorAntigens_MHCII_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$PRELU2_kallisto[is.na(mhag$PRELU2_kallisto)] <- 0
mhag$PRELU1_kallisto[is.na(mhag$PRELU1_kallisto)] <- 0

mhag=mhag[(mhag$PRELU2_kallisto>=0.58 & mhag$PRELU1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
18143




library(dplyr)
result <- mhag %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,PRELU2_kallisto,PRELU1_kallisto,Dx_Mut) %>%
    filter(Rel_Mut == max(Rel_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,PRELU2_kallisto,PRELU1_kallisto,Rel_Mut,Dx_Mut)
result=as.data.frame(result)

result <- result %>% 
    group_by(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,PRELU2_kallisto,PRELU1_kallisto,Rel_Mut) %>%
    filter(Dx_Mut == max(Dx_Mut)) %>%
    arrange(chr,pos,genes,id,peptide,hla,hla.1,effect,mut_nt,wt_nt,relapse_tx,diagnosis_tx,ic50,ic50.1,logFC,Qvalue,PRELU2_kallisto,PRELU1_kallisto,Rel_Mut,Dx_Mut)

result=as.data.frame(result)

nrow(result)
	
table(result$PRELU2_kallisto<result$PRELU1_kallisto)

table(result$PRELU2_kallisto>result$PRELU1_kallisto)

table(result$Rel_Mut<result$Dx_Mut)


table(result$Rel_Mut>result$Dx_Mut)




#A livello di trascritto
	table(mhag$PRELU2_kallisto > mhag$PRELU1_kallisto)

FALSE  TRUE
 6024  8910


wilcox.test(mhag$PRELU2_kallisto, mhag$PRELU1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$PRELU2_kallisto and mhag$PRELU1_kallisto
W = 115022845, p-value = 1.227e-06
alternative hypothesis: true location shift is greater than 0


#A livello di trascritto
 table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
 4941  7569



wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 99869526, p-value < 2.2e-16
alternative hypothesis: true location shift is greater than 0





#Calcolo reads su alleli per classe II


 for i in *Summary*_formatted.txt; do echo $i; prefix=$(echo $i | sed 's/.txt/_NoHeaders.txt/g'); echo $prefix;  sed -e '/Strong/d' $i | sed -e '/Weak/d' - | grep -v WT_ > $prefix ; done

#ALFE RelDiag
ALFE1 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE1_Allfiles_sorted.bam
ALFE2 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE2_Allfiles_sorted.bam

awk '{OFS="\t"; print $7":"$8"-"$8}'   /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_ALFE_RelDiag_MHCII_formatted_NoHeaders.txt | sort | uniq  > Summary_ALFE_NeoEpitopi_RelDiag_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE1_Allfiles_sorted.bam >> ALFE_NeoAGs_MHCII_RelDiag_Binders_Regions_Diagnosis.txt; done < Summary_ALFE_NeoEpitopi_RelDiag_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE2_Allfiles_sorted.bam >> ALFE_NeoAGs_MHCII_RelDiag_Binders_Regions_Relapse.txt; done < Summary_ALFE_NeoEpitopi_RelDiag_MHCII_Regions.txt
 

library(stringr)

alfe_NeoAG_mhcii_RelDiag=read.table('/home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_ALFE_RelDiag_MHCII_formatted_NoHeaders.txt')
 colnames(alfe_NeoAG_mhcii_RelDiag)=c("Gene","Peptide","DRB1","DRB1","Mutation","genes","chr","pos","wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_DRB1","IC50_DRB1")
alfe_NeoAG_mhcii_RelDiag=merge(alfe_NeoAG_mhcii_RelDiag,ALFE_RNK,by="genes",all.x=T)
relapse=read.table('ALFE_NeoAGs_MHCII_RelDiag_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('ALFE_NeoAGs_MHCII_RelDiag_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
alfe_NeoAG_mhcii_RelDiag=merge(alfe_NeoAG_mhcii_RelDiag,relapse,by=c("chr","pos"),all.x=T)
alfe_NeoAG_mhcii_RelDiag=merge(alfe_NeoAG_mhcii_RelDiag,diagnosis,by=c("chr","pos"),all.x=T)
alfe_NeoAG_mhcii_RelDiag$Rel_Mut=str_count(str_to_upper(alfe_NeoAG_mhcii_RelDiag$rel_bases), as.character(alfe_NeoAG_mhcii_RelDiag$mut_nt))
alfe_NeoAG_mhcii_RelDiag$Rel_Wt=str_count(str_to_upper(alfe_NeoAG_mhcii_RelDiag$rel_bases), as.character(alfe_NeoAG_mhcii_RelDiag$wt_nt))
alfe_NeoAG_mhcii_RelDiag$Dx_Mut=str_count(str_to_upper(alfe_NeoAG_mhcii_RelDiag$dx_bases), as.character(alfe_NeoAG_mhcii_RelDiag$mut_nt))
alfe_NeoAG_mhcii_RelDiag$Dx_Wt=str_count(str_to_upper(alfe_NeoAG_mhcii_RelDiag$dx_bases), as.character(alfe_NeoAG_mhcii_RelDiag$wt_nt))
 write.table(unique(alfe_NeoAG_mhcii_RelDiag[,c(1:18,26:30)]),"Summary_ALFE_Annotato_NeoAGs_MHCII_RelDiag_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)

#ALFE
mhag=read.table('Summary_ALFE_Annotato_NeoAGs_MHCII_RelDiag_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$ALFE2_kallisto[is.na(mhag$ALFE2_kallisto)] <- 0
mhag$ALFE1_kallisto[is.na(mhag$ALFE1_kallisto)] <- 0


mhag=mhag[(mhag$ALFE2_kallisto>=0.58 & mhag$ALFE1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
112

table(mhag$Dx_Mut>=3)
table(mhag$Rel_Mut>=3)



#A livello di trascritto
> table(mhag$ALFE2_kallisto > mhag$ALFE1_kallisto)

FALSE  TRUE
   90    22
> table(mhag$ALFE2_kallisto < mhag$ALFE1_kallisto)

FALSE  TRUE
   22    90
> wilcox.test(mhag$ALFE2_kallisto, mhag$ALFE1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$ALFE2_kallisto and mhag$ALFE1_kallisto
W = 3589, p-value = 1
alternative hypothesis: true location shift is greater than 0

#A livello di trascritto
>  table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE
  112
> table(mhag$Rel_Mut<mhag$Dx_Mut)

TRUE
 112


wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 7, p-value = 0.9966
alternative hypothesis: true location shift is greater than 0

Warning message:
In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
  cannot compute exact p-value with ties




#BESU
#RealDiag

/home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_BESU_RelDiag_MHCII_formatted_NoHeaders.txt

BESU1 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/BESU1_Allfiles_sorted.bam
BESU2 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/BESU2_Allfiles_sorted.bam

awk '{OFS="\t"; print $7":"$8"-"$8}'   /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_BESU_RelDiag_MHCII_formatted_NoHeaders.txt | sort | uniq  > Summary_BESU_NeoEpitopi_RelDiag_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/BESU1_Allfiles_sorted.bam >> BESU_NeoAGs_MHCII_RelDiag_Binders_Regions_Diagnosis.txt; done < Summary_BESU_NeoEpitopi_RelDiag_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/BESU2_Allfiles_sorted.bam >> BESU_NeoAGs_MHCII_RelDiag_Binders_Regions_Relapse.txt; done < Summary_BESU_NeoEpitopi_RelDiag_MHCII_Regions.txt
 

besu_NeoAG_mhcii_RelDiag=read.table('/home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_BESU_RelDiag_MHCII_formatted_NoHeaders.txt')
 colnames(besu_NeoAG_mhcii_RelDiag)=c("Gene","Peptide","DRB1","DRB1","Mutation","genes","chr","pos","wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_DRB1","IC50_DRB1")
besu_NeoAG_mhcii_RelDiag=merge(besu_NeoAG_mhcii_RelDiag,BESU_RNK,by="genes",all.x=T)
relapse=read.table('BESU_NeoAGs_MHCII_RelDiag_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('BESU_NeoAGs_MHCII_RelDiag_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
besu_NeoAG_mhcii_RelDiag=merge(besu_NeoAG_mhcii_RelDiag,relapse,by=c("chr","pos"),all.x=T)
besu_NeoAG_mhcii_RelDiag=merge(besu_NeoAG_mhcii_RelDiag,diagnosis,by=c("chr","pos"),all.x=T)
besu_NeoAG_mhcii_RelDiag$Rel_Mut=str_count(str_to_upper(besu_NeoAG_mhcii_RelDiag$rel_bases), as.character(besu_NeoAG_mhcii_RelDiag$mut_nt))
besu_NeoAG_mhcii_RelDiag$Rel_Wt=str_count(str_to_upper(besu_NeoAG_mhcii_RelDiag$rel_bases), as.character(besu_NeoAG_mhcii_RelDiag$wt_nt))
besu_NeoAG_mhcii_RelDiag$Dx_Mut=str_count(str_to_upper(besu_NeoAG_mhcii_RelDiag$dx_bases), as.character(besu_NeoAG_mhcii_RelDiag$mut_nt))
besu_NeoAG_mhcii_RelDiag$Dx_Wt=str_count(str_to_upper(besu_NeoAG_mhcii_RelDiag$dx_bases), as.character(besu_NeoAG_mhcii_RelDiag$wt_nt))
 write.table(unique(besu_NeoAG_mhcii_RelDiag[,c(1:18,26:30)]),"Summary_BESU_Annotato_NeoAGs_MHCII_RelDiag_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)


mhag=read.table('Summary_BESU_Annotato_NeoAGs_MHCII_RelDiag_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$BESU2_kallisto[is.na(mhag$BESU2_kallisto)] <- 0
mhag$BESU1_kallisto[is.na(mhag$BESU1_kallisto)] <- 0


mhag=mhag[(mhag$BESU2_kallisto>=0.58 & mhag$BESU1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
112

table(mhag$Dx_Mut>=3)
table(mhag$Rel_Mut>=3)


#A livello di trascritto
  table(mhag$BESU2_kallisto > mhag$BESU1_kallisto)

FALSE  TRUE
   60    49
  table(mhag$BESU2_kallisto < mhag$BESU1_kallisto)

FALSE  TRUE
   49    60

wilcox.test(mhag$BESU2_kallisto, mhag$BESU1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$BESU2_kallisto and mhag$BESU1_kallisto
W = 3589, p-value = 1
alternative hypothesis: true location shift is greater than 0

#A livello di trascritto
>  table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
   60    49

table(mhag$Rel_Mut<mhag$Dx_Mut)

TRUE
 112


wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 7, p-value = 0.9966
alternative hypothesis: true location shift is greater than 0

Warning message:
In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
  cannot compute exact p-value with ties




#DiagOnly
BESU1 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/BESU1_Allfiles_sorted.bam
BESU2 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/BESU2_Allfiles_sorted.bam

awk '{OFS="\t"; print $7":"$8"-"$8}'   /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_BESU_DiagOnly_MHCII_formatted_NoHeaders.txt | sort | uniq  > Summary_BESU_NeoEpitopi_DiagOnly_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/BESU1_Allfiles_sorted.bam >> BESU_NeoAGs_MHCII_DiagOnly_Binders_Regions_Diagnosis.txt; done < Summary_BESU_NeoEpitopi_DiagOnly_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/BESU2_Allfiles_sorted.bam >> BESU_NeoAGs_MHCII_DiagOnly_Binders_Regions_Relapse.txt; done < Summary_BESU_NeoEpitopi_DiagOnly_MHCII_Regions.txt
 

besu_NeoAG_mhcii_DiagOnly=read.table('/home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_BESU_DiagOnly_MHCII_formatted_NoHeaders.txt')
 colnames(besu_NeoAG_mhcii_DiagOnly)=c("Gene","Peptide","DRB1","DRB1","Mutation","genes","chr","pos","wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_DRB1","IC50_DRB1")

#Besu DiagOnly non ha NeoAGs espressi per Classe II




#CALU

#Calu RelDiag
CALU1 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/CALU1_Allfiles_sorted.bam
CALU2 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/CALU2_Allfiles_sorted.bam


awk '{OFS="\t"; print $7":"$8"-"$8}'   /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_CALU_RelDiag_MHCII_formatted_NoHeaders.txt | sort | uniq  > Summary_CALU_NeoEpitopi_RelDiag_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/CALU1_Allfiles_sorted.bam >> CALU_NeoAGs_MHCII_RelDiag_Binders_Regions_Diagnosis.txt; done < Summary_CALU_NeoEpitopi_RelDiag_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/CALU2_Allfiles_sorted.bam >> CALU_NeoAGs_MHCII_RelDiag_Binders_Regions_Relapse.txt; done < Summary_CALU_NeoEpitopi_RelDiag_MHCII_Regions.txt
 

calu_NeoAG_mhcii_RelDiag=read.table('/home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_CALU_RelDiag_MHCII_formatted_NoHeaders.txt')
 colnames(calu_NeoAG_mhcii_RelDiag)=c("Gene","Peptide","DRB1","DRB1","Mutation","genes","chr","pos","wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_DRB1","IC50_DRB1")
calu_NeoAG_mhcii_RelDiag=merge(calu_NeoAG_mhcii_RelDiag,CALU_RNK,by="genes",all.x=T)
relapse=read.table('CALU_NeoAGs_MHCII_RelDiag_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('CALU_NeoAGs_MHCII_RelDiag_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
calu_NeoAG_mhcii_RelDiag=merge(calu_NeoAG_mhcii_RelDiag,relapse,by=c("chr","pos"),all.x=T)
calu_NeoAG_mhcii_RelDiag=merge(calu_NeoAG_mhcii_RelDiag,diagnosis,by=c("chr","pos"),all.x=T)
calu_NeoAG_mhcii_RelDiag$Rel_Mut=str_count(str_to_upper(calu_NeoAG_mhcii_RelDiag$rel_bases), as.character(calu_NeoAG_mhcii_RelDiag$mut_nt))
calu_NeoAG_mhcii_RelDiag$Rel_Wt=str_count(str_to_upper(calu_NeoAG_mhcii_RelDiag$rel_bases), as.character(calu_NeoAG_mhcii_RelDiag$wt_nt))
calu_NeoAG_mhcii_RelDiag$Dx_Mut=str_count(str_to_upper(calu_NeoAG_mhcii_RelDiag$dx_bases), as.character(calu_NeoAG_mhcii_RelDiag$mut_nt))
calu_NeoAG_mhcii_RelDiag$Dx_Wt=str_count(str_to_upper(calu_NeoAG_mhcii_RelDiag$dx_bases), as.character(calu_NeoAG_mhcii_RelDiag$wt_nt))
 write.table(unique(calu_NeoAG_mhcii_RelDiag[,c(1:18,26:30)]),"Summary_CALU_Annotato_NeoAGs_MHCII_RelDiag_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)


mhag=read.table('Summary_CALU_Annotato_NeoAGs_MHCII_RelDiag_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$CALU2_kallisto[is.na(mhag$CALU2_kallisto)] <- 0
mhag$CALU1_kallisto[is.na(mhag$CALU1_kallisto)] <- 0


mhag=mhag[(mhag$CALU2_kallisto>=0.58 & mhag$CALU1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
112

table(mhag$Dx_Mut>=3)
table(mhag$Rel_Mut>=3)

#A livello di trascritto
table(mhag$CALU2_kallisto > mhag$CALU1_kallisto)

FALSE  TRUE
    5     2
table(mhag$CALU2_kallisto < mhag$CALU1_kallisto)

FALSE  TRUE
    2     5


wilcox.test(mhag$CALU2_kallisto, mhag$CALU1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$CALU2_kallisto and mhag$CALU1_kallisto
W = 3589, p-value = 1
alternative hypothesis: true location shift is greater than 0

#A livello di trascritto
 table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
   60    49

table(mhag$Rel_Mut<mhag$Dx_Mut)

TRUE
 112


wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 7, p-value = 0.9966
alternative hypothesis: true location shift is greater than 0

Warning message:
In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
  cannot compute exact p-value with ties



#DEIV

 sed -e '/Strong/d' Summary_DEIV_RelOnly_MHCII.txt  | sed -e '/Weak/d' - | awk '{OFS="\t"; print $1"_"$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' | sed 's/\"//g' | grep -v WT_  > Summary_DEIV_RelOnly_MHCII_formatted_NoHeaders.txt
sed -e '/Strong/d' Summary_DEIV_DiagOnly_MHCII.txt  | sed -e '/Weak/d' - | awk '{OFS="\t"; print $1"_"$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' | sed 's/\"//g' | grep -v WT_  > Summary_DEIV_DiagOnly_MHCII_formatted_NoHeaders.txt




#RelOnly
DEIV2 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2_Allfiles_sorted.bam
DEIV2b = /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2b_Allfiles_sorted.bam


awk '{OFS="\t"; print $7":"$8"-"$8}'   /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_DEIV_RelOnly_MHCII_formatted_NoHeaders.txt | sort | uniq  > Summary_DEIV_Annotato_NeoEpitopi_RelOnly_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2_Allfiles_sorted.bam >> DEIV_NeoAGs_MHCII_RelOnly_Binders_Regions_Diagnosis.txt; done < Summary_DEIV_Annotato_NeoEpitopi_RelOnly_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2b_Allfiles_sorted.bam >> DEIV_NeoAGs_MHCII_RelOnly_Binders_Regions_Relapse.txt; done < Summary_DEIV_Annotato_NeoEpitopi_RelOnly_MHCII_Regions.txt
 

deiv_NeoAG_mhcii_RelOnly=read.table('/home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_DEIV_RelOnly_MHCII_formatted_NoHeaders.txt')
 colnames(deiv_NeoAG_mhcii_RelOnly)=c("Gene","Peptide","DRB1","DRB1","Mutation","genes","chr","pos","wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_DRB1","IC50_DRB1")
deiv_NeoAG_mhcii_RelOnly=merge(deiv_NeoAG_mhcii_RelOnly,DEIV_RNK,by="genes",all.x=T)
relapse=read.table('DEIV_NeoAGs_MHCII_RelOnly_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DEIV_NeoAGs_MHCII_RelOnly_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
deiv_NeoAG_mhcii_RelOnly=merge(deiv_NeoAG_mhcii_RelOnly,relapse,by=c("chr","pos"),all.x=T)
deiv_NeoAG_mhcii_RelOnly=merge(deiv_NeoAG_mhcii_RelOnly,diagnosis,by=c("chr","pos"),all.x=T)
deiv_NeoAG_mhcii_RelOnly$Rel_Mut=str_count(str_to_upper(deiv_NeoAG_mhcii_RelOnly$rel_bases), as.character(deiv_NeoAG_mhcii_RelOnly$mut_nt))
deiv_NeoAG_mhcii_RelOnly$Rel_Wt=str_count(str_to_upper(deiv_NeoAG_mhcii_RelOnly$rel_bases), as.character(deiv_NeoAG_mhcii_RelOnly$wt_nt))
deiv_NeoAG_mhcii_RelOnly$Dx_Mut=str_count(str_to_upper(deiv_NeoAG_mhcii_RelOnly$dx_bases), as.character(deiv_NeoAG_mhcii_RelOnly$mut_nt))
deiv_NeoAG_mhcii_RelOnly$Dx_Wt=str_count(str_to_upper(deiv_NeoAG_mhcii_RelOnly$dx_bases), as.character(deiv_NeoAG_mhcii_RelOnly$wt_nt))
 write.table(unique(deiv_NeoAG_mhcii_RelOnly[,c(1:18,27:30)]),"Summary_DEIV_Annotato_NeoAGs_MHCII_RelOnly_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)


mhag=read.table('Summary_DEIV_Annotato_NeoAGs_MHCII_RelOnly_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$DEIV2_kallisto[is.na(mhag$DEIV2b_kallisto)] <- 0
mhag$DEIV1_kallisto[is.na(mhag$DEIV2_kallisto)] <- 0


mhag=mhag[(mhag$DEIV2b_kallisto>=0.58 & mhag$DEIV2_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
112

table(mhag$Dx_Mut>=3)
table(mhag$Rel_Mut>=3)

#A livello di trascritto
table(mhag$DEIV2b_kallisto > mhag$DEIV2_kallisto)

FALSE  TRUE
    5     2
table(mhag$DEIV2_kallisto < mhag$DEIV1_kallisto)

FALSE  TRUE
    2     5


wilcox.test(mhag$DEIV2b_kallisto, mhag$DEIV2_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$DEIV2_kallisto and mhag$DEIV1_kallisto
W = 3589, p-value = 1
alternative hypothesis: true location shift is greater than 0

#A livello di trascritto
 table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
   60    49

table(mhag$Rel_Mut<mhag$Dx_Mut)

TRUE
 112


wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 7, p-value = 0.9966
alternative hypothesis: true location shift is greater than 0

Warning message:
In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
  cannot compute exact p-value with ties



#DiagOnly
DEIV2 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2_Allfiles_sorted.bam
DEIV2b = /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2b_Allfiles_sorted.bam


awk '{OFS="\t"; print $7":"$8"-"$8}'   /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_DEIV_DiagOnly_MHCII_formatted_NoHeaders.txt | sort | uniq  > Summary_DEIV_Annotato_NeoEpitopi_DiagOnly_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2_Allfiles_sorted.bam >> DEIV_NeoAGs_MHCII_DiagOnly_Binders_Regions_Diagnosis.txt; done < Summary_DEIV_Annotato_NeoEpitopi_DiagOnly_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2b_Allfiles_sorted.bam >> DEIV_NeoAGs_MHCII_DiagOnly_Binders_Regions_Relapse.txt; done < Summary_DEIV_Annotato_NeoEpitopi_DiagOnly_MHCII_Regions.txt
 

deiv_NeoAG_mhcii_DiagOnly=read.table('/home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_DEIV_DiagOnly_MHCII_formatted_NoHeaders.txt')
 colnames(deiv_NeoAG_mhcii_DiagOnly)=c("Gene","Peptide","DRB1","DRB1","Mutation","genes","chr","pos","wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_DRB1","IC50_DRB1")
deiv_NeoAG_mhcii_DiagOnly=merge(deiv_NeoAG_mhcii_DiagOnly,DEIV_RNK,by="genes",all.x=T)
relapse=read.table('DEIV_NeoAGs_MHCII_DiagOnly_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DEIV_NeoAGs_MHCII_DiagOnly_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
deiv_NeoAG_mhcii_DiagOnly=merge(deiv_NeoAG_mhcii_DiagOnly,relapse,by=c("chr","pos"),all.x=T)
deiv_NeoAG_mhcii_DiagOnly=merge(deiv_NeoAG_mhcii_DiagOnly,diagnosis,by=c("chr","pos"),all.x=T)
deiv_NeoAG_mhcii_DiagOnly$Rel_Mut=str_count(str_to_upper(deiv_NeoAG_mhcii_DiagOnly$rel_bases), as.character(deiv_NeoAG_mhcii_DiagOnly$mut_nt))
deiv_NeoAG_mhcii_DiagOnly$Rel_Wt=str_count(str_to_upper(deiv_NeoAG_mhcii_DiagOnly$rel_bases), as.character(deiv_NeoAG_mhcii_DiagOnly$wt_nt))
deiv_NeoAG_mhcii_DiagOnly$Dx_Mut=str_count(str_to_upper(deiv_NeoAG_mhcii_DiagOnly$dx_bases), as.character(deiv_NeoAG_mhcii_DiagOnly$mut_nt))
deiv_NeoAG_mhcii_DiagOnly$Dx_Wt=str_count(str_to_upper(deiv_NeoAG_mhcii_DiagOnly$dx_bases), as.character(deiv_NeoAG_mhcii_DiagOnly$wt_nt))
 write.table(unique(deiv_NeoAG_mhcii_DiagOnly[,c(1:18,27:30)]),"Summary_DEIV_Annotato_NeoAGs_MHCII_DiagOnly_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)


mhag=read.table('Summary_DEIV_Annotato_NeoAGs_MHCII_DiagOnly_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$DEIV2_kallisto[is.na(mhag$DEIV2b_kallisto)] <- 0
mhag$DEIV1_kallisto[is.na(mhag$DEIV2_kallisto)] <- 0


mhag=mhag[(mhag$DEIV2b_kallisto>=0.58 & mhag$DEIV2_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
112


table(mhag$Dx_Mut>=3)
table(mhag$Rel_Mut>=3)

#A livello di trascritto
table(mhag$DEIV2b_kallisto > mhag$DEIV2_kallisto)

FALSE  TRUE
    5     2
table(mhag$DEIV2_kallisto < mhag$DEIV1_kallisto)

FALSE  TRUE
    2     5


wilcox.test(mhag$DEIV2b_kallisto, mhag$DEIV2_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$DEIV2_kallisto and mhag$DEIV1_kallisto
W = 3589, p-value = 1
alternative hypothesis: true location shift is greater than 0

#A livello di trascritto
 table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
   60    49

table(mhag$Rel_Mut<mhag$Dx_Mut)

TRUE
 112


wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 7, p-value = 0.9966
alternative hypothesis: true location shift is greater than 0

Warning message:
In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
  cannot compute exact p-value with ties



#DEST
DEST1 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST1_Allfiles_sorted.bam
DEST2 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST2_Allfiles_sorted.bam

awk '{OFS="\t"; print $7":"$8"-"$8}'   /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_DEST_RelDiag_MHCII_formatted_NoHeaders.txt | sort | uniq  > Summary_DEST_NeoEpitopi_RelDiag_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST1_Allfiles_sorted.bam >> DEST_NeoAGs_MHCII_RelDiag_Binders_Regions_Diagnosis.txt; done < Summary_DEST_NeoEpitopi_RelDiag_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST2_Allfiles_sorted.bam >> DEST_NeoAGs_MHCII_RelDiag_Binders_Regions_Relapse.txt; done < Summary_DEST_NeoEpitopi_RelDiag_MHCII_Regions.txt
 

dest_NeoAG_mhcii_RelDiag=read.table('/home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_DEST_RelDiag_MHCII_formatted_NoHeaders.txt')
 colnames(dest_NeoAG_mhcii_RelDiag)=c("Gene","Peptide","DRB1","DRB1","Mutation","genes","chr","pos","wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_DRB1","IC50_DRB1")
dest_NeoAG_mhcii_RelDiag=merge(dest_NeoAG_mhcii_RelDiag,DEST_RNK,by="genes",all.x=T)
relapse=read.table('DEST_NeoAGs_MHCII_RelDiag_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DEST_NeoAGs_MHCII_RelDiag_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
dest_NeoAG_mhcii_RelDiag=merge(dest_NeoAG_mhcii_RelDiag,relapse,by=c("chr","pos"),all.x=T)
dest_NeoAG_mhcii_RelDiag=merge(dest_NeoAG_mhcii_RelDiag,diagnosis,by=c("chr","pos"),all.x=T)
dest_NeoAG_mhcii_RelDiag$Rel_Mut=str_count(str_to_upper(dest_NeoAG_mhcii_RelDiag$rel_bases), as.character(dest_NeoAG_mhcii_RelDiag$mut_nt))
dest_NeoAG_mhcii_RelDiag$Rel_Wt=str_count(str_to_upper(dest_NeoAG_mhcii_RelDiag$rel_bases), as.character(dest_NeoAG_mhcii_RelDiag$wt_nt))
dest_NeoAG_mhcii_RelDiag$Dx_Mut=str_count(str_to_upper(dest_NeoAG_mhcii_RelDiag$dx_bases), as.character(dest_NeoAG_mhcii_RelDiag$mut_nt))
dest_NeoAG_mhcii_RelDiag$Dx_Wt=str_count(str_to_upper(dest_NeoAG_mhcii_RelDiag$dx_bases), as.character(dest_NeoAG_mhcii_RelDiag$wt_nt))
 write.table(unique(dest_NeoAG_mhcii_RelDiag[,c(1:18,26:30)]),"Summary_DEST_Annotato_NeoAGs_MHCII_RelDiag_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)

mhag=read.table('Summary_DEST_Annotato_NeoAGs_MHCII_RelDiag_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$DEST2_kallisto[is.na(mhag$DEST2_kallisto)] <- 0
mhag$DEST1_kallisto[is.na(mhag$DEST1_kallisto)] <- 0


mhag=mhag[(mhag$DEST2_kallisto>=0.58 & mhag$DEST1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
31
table(mhag$Dx_Mut>=3)
table(mhag$Rel_Mut>=3)



#A livello di trascritto
table(mhag$DEST2_kallisto > mhag$DEST1_kallisto)

FALSE  TRUE
    5     2
table(mhag$DEST2_kallisto < mhag$DEST1_kallisto)

FALSE  TRUE
    2     5


wilcox.test(mhag$DEST2_kallisto, mhag$DEST1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$DEST2_kallisto and mhag$DEST1_kallisto
W = 3589, p-value = 1
alternative hypothesis: true location shift is greater than 0

#A livello di trascritto
table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
   16    15


table(mhag$Rel_Mut<mhag$Dx_Mut)

FALSE  TRUE
   15    16



 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 465, p-value = 0.592
alternative hypothesis: true location shift is greater than 0

Warning message:
In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
  cannot compute exact p-value with ties



DEST1 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST1_Allfiles_sorted.bam
DEST2 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST2_Allfiles_sorted.bam



awk '{OFS="\t"; print $7":"$8"-"$8}'   /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_DEST_RelOnly_MHCII_formatted_NoHeaders.txt | sort | uniq  > Summary_DEST_NeoEpitopi_RelOnly_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST1_Allfiles_sorted.bam >> DEST_NeoAGs_MHCII_RelOnly_Binders_Regions_Diagnosis.txt; done < Summary_DEST_NeoEpitopi_RelOnly_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST2_Allfiles_sorted.bam >> DEST_NeoAGs_MHCII_RelOnly_Binders_Regions_Relapse.txt; done < Summary_DEST_NeoEpitopi_RelOnly_MHCII_Regions.txt
 

dest_NeoAG_mhcii_RelOnly=read.table('/home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_DEST_RelOnly_MHCII_formatted_NoHeaders.txt')
 colnames(dest_NeoAG_mhcii_RelOnly)=c("Gene","Peptide","DRB1","DRB1","Mutation","genes","chr","pos","wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_DRB1","IC50_DRB1")
dest_NeoAG_mhcii_RelOnly=merge(dest_NeoAG_mhcii_RelOnly,DEST_RNK,by="genes",all.x=T)
relapse=read.table('DEST_NeoAGs_MHCII_RelOnly_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DEST_NeoAGs_MHCII_RelOnly_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
dest_NeoAG_mhcii_RelOnly=merge(dest_NeoAG_mhcii_RelOnly,relapse,by=c("chr","pos"),all.x=T)
dest_NeoAG_mhcii_RelOnly=merge(dest_NeoAG_mhcii_RelOnly,diagnosis,by=c("chr","pos"),all.x=T)
dest_NeoAG_mhcii_RelOnly$Rel_Mut=str_count(str_to_upper(dest_NeoAG_mhcii_RelOnly$rel_bases), as.character(dest_NeoAG_mhcii_RelOnly$mut_nt))
dest_NeoAG_mhcii_RelOnly$Rel_Wt=str_count(str_to_upper(dest_NeoAG_mhcii_RelOnly$rel_bases), as.character(dest_NeoAG_mhcii_RelOnly$wt_nt))
dest_NeoAG_mhcii_RelOnly$Dx_Mut=str_count(str_to_upper(dest_NeoAG_mhcii_RelOnly$dx_bases), as.character(dest_NeoAG_mhcii_RelOnly$mut_nt))
dest_NeoAG_mhcii_RelOnly$Dx_Wt=str_count(str_to_upper(dest_NeoAG_mhcii_RelOnly$dx_bases), as.character(dest_NeoAG_mhcii_RelOnly$wt_nt))
 write.table(unique(dest_NeoAG_mhcii_RelOnly[,c(1:18,26:30)]),"Summary_DEST_Annotato_NeoAGs_MHCII_RelOnly_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)

mhag=read.table('Summary_DEST_Annotato_NeoAGs_MHCII_RelOnly_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$DEST2_kallisto[is.na(mhag$DEST2_kallisto)] <- 0
mhag$DEST1_kallisto[is.na(mhag$DEST1_kallisto)] <- 0


mhag=mhag[(mhag$DEST2_kallisto>=0.58 & mhag$DEST1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
10
table(mhag$Dx_Mut>=3)
table(mhag$Rel_Mut>=3)


#A livello di trascritto
table(mhag$DEST2_kallisto > mhag$DEST1_kallisto)

FALSE  TRUE
    5     2
table(mhag$DEST2_kallisto < mhag$DEST1_kallisto)

FALSE  TRUE
    2     5


wilcox.test(mhag$DEST2_kallisto, mhag$DEST1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$DEST2_kallisto and mhag$DEST1_kallisto
W = 3589, p-value = 1
alternative hypothesis: true location shift is greater than 0

#A livello di trascritto
table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
   16    15


table(mhag$Rel_Mut<mhag$Dx_Mut)

FALSE  TRUE
   15    16



 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 465, p-value = 0.592
alternative hypothesis: true location shift is greater than 0

Warning message:
In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
  cannot compute exact p-value with ties



#DR4
#DR4 RelDiag
DR41 =  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR41_Allfiles_sorted.bam
DR42 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR42_Allfiles_sorted.bam

 awk '{OFS="\t"; print $7":"$8"-"$8}'   /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_DR4_RelDiag_MHCII_formatted_NoHeaders.txt | sort | uniq  > Summary_DR4_NeoEpitopi_RelDiag_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR41_Allfiles_sorted.bam >> DR4_NeoAGs_MHCII_RelDiag_Binders_Regions_Diagnosis.txt; done < Summary_DR4_NeoEpitopi_RelDiag_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR42_Allfiles_sorted.bam >> DR4_NeoAGs_MHCII_RelDiag_Binders_Regions_Relapse.txt; done < Summary_DR4_NeoEpitopi_RelDiag_MHCII_Regions.txt
 

dr4_NeoAG_mhcii_RelDiag=read.table('/home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_DR4_RelDiag_MHCII_formatted_NoHeaders.txt')
 colnames(dr4_NeoAG_mhcii_RelDiag)=c("Gene","Peptide","DRB1","DRB1","Mutation","genes","chr","pos","wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_DRB1","IC50_DRB1")
dr4_NeoAG_mhcii_RelDiag=merge(dr4_NeoAG_mhcii_RelDiag,DR4_RNK,by="genes",all.x=T)
relapse=read.table('DR4_NeoAGs_MHCII_RelDiag_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DR4_NeoAGs_MHCII_RelDiag_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
dr4_NeoAG_mhcii_RelDiag=merge(dr4_NeoAG_mhcii_RelDiag,relapse,by=c("chr","pos"),all.x=T)
dr4_NeoAG_mhcii_RelDiag=merge(dr4_NeoAG_mhcii_RelDiag,diagnosis,by=c("chr","pos"),all.x=T)
dr4_NeoAG_mhcii_RelDiag$Rel_Mut=str_count(str_to_upper(dr4_NeoAG_mhcii_RelDiag$rel_bases), as.character(dr4_NeoAG_mhcii_RelDiag$mut_nt))
dr4_NeoAG_mhcii_RelDiag$Rel_Wt=str_count(str_to_upper(dr4_NeoAG_mhcii_RelDiag$rel_bases), as.character(dr4_NeoAG_mhcii_RelDiag$wt_nt))
dr4_NeoAG_mhcii_RelDiag$Dx_Mut=str_count(str_to_upper(dr4_NeoAG_mhcii_RelDiag$dx_bases), as.character(dr4_NeoAG_mhcii_RelDiag$mut_nt))
dr4_NeoAG_mhcii_RelDiag$Dx_Wt=str_count(str_to_upper(dr4_NeoAG_mhcii_RelDiag$dx_bases), as.character(dr4_NeoAG_mhcii_RelDiag$wt_nt))
 write.table(unique(dr4_NeoAG_mhcii_RelDiag[,c(1:18,26:30)]),"Summary_DR4_Annotato_NeoAGs_MHCII_RelDiag_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)

mhag=read.table('Summary_DR4_Annotato_NeoAGs_MHCII_RelDiag_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$DR42_kallisto[is.na(mhag$DR42_kallisto)] <- 0
mhag$DR41_kallisto[is.na(mhag$DR41_kallisto)] <- 0
33

mhag=mhag[(mhag$DR42_kallisto>=0.58 & mhag$DR41_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
10

table(mhag$Dx_Mut>=3)
table(mhag$Rel_Mut>=3)

#A livello di trascritto
table(mhag$DR42_kallisto > mhag$DR41_kallisto)

FALSE  TRUE
    5     2
table(mhag$DR42_kallisto < mhag$DR41_kallisto)

FALSE  TRUE
    2     5


wilcox.test(mhag$DR42_kallisto, mhag$DR41_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$DR42_kallisto and mhag$DR41_kallisto
W = 3589, p-value = 1
alternative hypothesis: true location shift is greater than 0

#A livello di trascritto
table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
   16    15


table(mhag$Rel_Mut<mhag$Dx_Mut)

FALSE  TRUE
   15    16



 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 465, p-value = 0.592
alternative hypothesis: true location shift is greater than 0

Warning message:
In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
  cannot compute exact p-value with ties



#DR4 REL ONLY
 awk '{OFS="\t"; print $7":"$8"-"$8}'   /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_DR4_RelOnly_MHCII_formatted_NoHeaders.txt | sort | uniq  > Summary_DR4_NeoEpitopi_RelOnly_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR41_Allfiles_sorted.bam >> DR4_NeoAGs_MHCII_RelOnly_Binders_Regions_Diagnosis.txt; done < Summary_DR4_NeoEpitopi_RelOnly_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR42_Allfiles_sorted.bam >> DR4_NeoAGs_MHCII_RelOnly_Binders_Regions_Relapse.txt; done < Summary_DR4_NeoEpitopi_RelOnly_MHCII_Regions.txt
 

dr4_NeoAG_mhcii_RelOnly=read.table('/home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_DR4_RelOnly_MHCII_formatted_NoHeaders.txt')
 colnames(dr4_NeoAG_mhcii_RelOnly)=c("Gene","Peptide","DRB1","DRB1","Mutation","genes","chr","pos","wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_DRB1","IC50_DRB1")
dr4_NeoAG_mhcii_RelOnly=merge(dr4_NeoAG_mhcii_RelOnly,DR4_RNK,by="genes",all.x=T)
relapse=read.table('DR4_NeoAGs_MHCII_RelOnly_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DR4_NeoAGs_MHCII_RelOnly_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
dr4_NeoAG_mhcii_RelOnly=merge(dr4_NeoAG_mhcii_RelOnly,relapse,by=c("chr","pos"),all.x=T)
dr4_NeoAG_mhcii_RelOnly=merge(dr4_NeoAG_mhcii_RelOnly,diagnosis,by=c("chr","pos"),all.x=T)
dr4_NeoAG_mhcii_RelOnly$Rel_Mut=str_count(str_to_upper(dr4_NeoAG_mhcii_RelOnly$rel_bases), as.character(dr4_NeoAG_mhcii_RelOnly$mut_nt))
dr4_NeoAG_mhcii_RelOnly$Rel_Wt=str_count(str_to_upper(dr4_NeoAG_mhcii_RelOnly$rel_bases), as.character(dr4_NeoAG_mhcii_RelOnly$wt_nt))
dr4_NeoAG_mhcii_RelOnly$Dx_Mut=str_count(str_to_upper(dr4_NeoAG_mhcii_RelOnly$dx_bases), as.character(dr4_NeoAG_mhcii_RelOnly$mut_nt))
dr4_NeoAG_mhcii_RelOnly$Dx_Wt=str_count(str_to_upper(dr4_NeoAG_mhcii_RelOnly$dx_bases), as.character(dr4_NeoAG_mhcii_RelOnly$wt_nt))
 write.table(unique(dr4_NeoAG_mhcii_RelOnly[,c(1:18,26:30)]),"Summary_DR4_Annotato_NeoAGs_MHCII_RelOnly_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)

mhag=read.table('Summary_DR4_Annotato_NeoAGs_MHCII_RelOnly_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$DR42_kallisto[is.na(mhag$DR42_kallisto)] <- 0
mhag$DR41_kallisto[is.na(mhag$DR41_kallisto)] <- 0
51

mhag=mhag[(mhag$DR42_kallisto>=0.58 & mhag$DR41_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
14
table(mhag$Dx_Mut>=3)
table(mhag$Rel_Mut>=3)


#A livello di trascritto
table(mhag$DR42_kallisto > mhag$DR41_kallisto)

FALSE  TRUE
    5     2
table(mhag$DR42_kallisto < mhag$DR41_kallisto)

FALSE  TRUE
    2     5


wilcox.test(mhag$DR42_kallisto, mhag$DR41_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$DR42_kallisto and mhag$DR41_kallisto
W = 3589, p-value = 1
alternative hypothesis: true location shift is greater than 0

#A livello di trascritto
table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
   16    15


table(mhag$Rel_Mut<mhag$Dx_Mut)

FALSE  TRUE
   15    16



 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 465, p-value = 0.592
alternative hypothesis: true location shift is greater than 0

Warning message:
In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
  cannot compute exact p-value with ties


#DR5
#DR5 RelDiag
DR51 =  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR51_Allfiles_sorted.bam
DR52 =  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR52_Allfiles_sorted.bam
 

 awk '{OFS="\t"; print $6":"$7"-"$7}'   /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_DR5_RelDiag_MHCII_formatted_NoHeaders.txt | sort | uniq  > Summary_DR5_NeoEpitopi_RelDiag_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR51_Allfiles_sorted.bam >> DR5_NeoAGs_MHCII_RelDiag_Binders_Regions_Diagnosis.txt; done < Summary_DR5_NeoEpitopi_RelDiag_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR52_Allfiles_sorted.bam >> DR5_NeoAGs_MHCII_RelDiag_Binders_Regions_Relapse.txt; done < Summary_DR5_NeoEpitopi_RelDiag_MHCII_Regions.txt
 

dr5_NeoAG_mhcii_RelDiag=read.table('/home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_DR5_RelDiag_MHCII_formatted_NoHeaders.txt')
 colnames(dr5_NeoAG_mhcii_RelDiag)=c("Gene","Peptide","DRB1","Mutation","genes","chr","pos","wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_DRB1")
dr5_NeoAG_mhcii_RelDiag=merge(dr5_NeoAG_mhcii_RelDiag,DR5_RNK,by="genes",all.x=T)
relapse=read.table('DR5_NeoAGs_MHCII_RelDiag_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DR5_NeoAGs_MHCII_RelDiag_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
dr5_NeoAG_mhcii_RelDiag=merge(dr5_NeoAG_mhcii_RelDiag,relapse,by=c("chr","pos"),all.x=T)
dr5_NeoAG_mhcii_RelDiag=merge(dr5_NeoAG_mhcii_RelDiag,diagnosis,by=c("chr","pos"),all.x=T)
dr5_NeoAG_mhcii_RelDiag$Rel_Mut=str_count(str_to_upper(dr5_NeoAG_mhcii_RelDiag$rel_bases), as.character(dr5_NeoAG_mhcii_RelDiag$mut_nt))
dr5_NeoAG_mhcii_RelDiag$Rel_Wt=str_count(str_to_upper(dr5_NeoAG_mhcii_RelDiag$rel_bases), as.character(dr5_NeoAG_mhcii_RelDiag$wt_nt))
dr5_NeoAG_mhcii_RelDiag$Dx_Mut=str_count(str_to_upper(dr5_NeoAG_mhcii_RelDiag$dx_bases), as.character(dr5_NeoAG_mhcii_RelDiag$mut_nt))
dr5_NeoAG_mhcii_RelDiag$Dx_Wt=str_count(str_to_upper(dr5_NeoAG_mhcii_RelDiag$dx_bases), as.character(dr5_NeoAG_mhcii_RelDiag$wt_nt))
 write.table(unique(dr5_NeoAG_mhcii_RelDiag[,c(1:16,25:28)]),"Summary_DR5_Annotato_NeoAGs_MHCII_RelDiag_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)


mhag=read.table('Summary_DR5_Annotato_NeoAGs_MHCII_RelDiag_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$DR52_kallisto[is.na(mhag$DR52_kallisto)] <- 0
mhag$DR51_kallisto[is.na(mhag$DR51_kallisto)] <- 0
63

mhag=mhag[(mhag$DR52_kallisto>=0.58 & mhag$DR51_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
12

table(mhag$Dx_Mut>=3)
table(mhag$Rel_Mut>=3)

#A livello di trascritto
table(mhag$DR52_kallisto > mhag$DR51_kallisto)

FALSE  TRUE
    5     2
table(mhag$DR52_kallisto < mhag$DR51_kallisto)

FALSE  TRUE
    2     5


wilcox.test(mhag$DR52_kallisto, mhag$DR51_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$DR52_kallisto and mhag$DR51_kallisto
W = 3589, p-value = 1
alternative hypothesis: true location shift is greater than 0

#A livello di trascritto
table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
   16    15


table(mhag$Rel_Mut<mhag$Dx_Mut)

FALSE  TRUE
   15    16



 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 465, p-value = 0.592
alternative hypothesis: true location shift is greater than 0

Warning message:
In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
  cannot compute exact p-value with ties



#FOCA
FOCA1 =  /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/FOCA1_Allfiles_sorted.bam
FOCA2 =  /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/FOCA2_Allfiles_sorted.bam

 awk '{OFS="\t"; print $7":"$8"-"$8}'   /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_FOCA_RelDiag_MHCII_formatted_NoHeaders.txt | sort | uniq  > Summary_FOCA_NeoEpitopi_RelDiag_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/FOCA1_Allfiles_sorted.bam >> FOCA_NeoAGs_MHCII_RelDiag_Binders_Regions_Diagnosis.txt; done < Summary_FOCA_NeoEpitopi_RelDiag_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/FOCA2_Allfiles_sorted.bam >> FOCA_NeoAGs_MHCII_RelDiag_Binders_Regions_Relapse.txt; done < Summary_FOCA_NeoEpitopi_RelDiag_MHCII_Regions.txt
 
foca_NeoAG_mhcii_RelDiag=read.table('/home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_FOCA_RelDiag_MHCII_formatted_NoHeaders.txt')
 colnames(foca_NeoAG_mhcii_RelDiag)=c("Gene","Peptide","DRB1","DRB1","Mutation","genes","chr","pos","wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_DRB1","IC50_DRB1")
foca_NeoAG_mhcii_RelDiag=merge(foca_NeoAG_mhcii_RelDiag,FOCA_RNK,by="genes",all.x=T)
relapse=read.table('FOCA_NeoAGs_MHCII_RelDiag_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('FOCA_NeoAGs_MHCII_RelDiag_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
foca_NeoAG_mhcii_RelDiag=merge(foca_NeoAG_mhcii_RelDiag,relapse,by=c("chr","pos"),all.x=T)
foca_NeoAG_mhcii_RelDiag=merge(foca_NeoAG_mhcii_RelDiag,diagnosis,by=c("chr","pos"),all.x=T)
foca_NeoAG_mhcii_RelDiag$Rel_Mut=str_count(str_to_upper(foca_NeoAG_mhcii_RelDiag$rel_bases), as.character(foca_NeoAG_mhcii_RelDiag$mut_nt))
foca_NeoAG_mhcii_RelDiag$Rel_Wt=str_count(str_to_upper(foca_NeoAG_mhcii_RelDiag$rel_bases), as.character(foca_NeoAG_mhcii_RelDiag$wt_nt))
foca_NeoAG_mhcii_RelDiag$Dx_Mut=str_count(str_to_upper(foca_NeoAG_mhcii_RelDiag$dx_bases), as.character(foca_NeoAG_mhcii_RelDiag$mut_nt))
foca_NeoAG_mhcii_RelDiag$Dx_Wt=str_count(str_to_upper(foca_NeoAG_mhcii_RelDiag$dx_bases), as.character(foca_NeoAG_mhcii_RelDiag$wt_nt))
 write.table(unique(foca_NeoAG_mhcii_RelDiag[,c(1:18,27:30)]),"Summary_FOCA_Annotato_NeoAGs_MHCII_RelDiag_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)


mhag=read.table('Summary_FOCA_Annotato_NeoAGs_MHCII_RelDiag_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$FOCA2_kallisto[is.na(mhag$FOCA2_kallisto)] <- 0
mhag$FOCA1_kallisto[is.na(mhag$FOCA1_kallisto)] <- 0
52

mhag=mhag[(mhag$FOCA2_kallisto>=0.58 & mhag$FOCA1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
12

table(mhag$Dx_Mut>=3)
table(mhag$Rel_Mut>=3)

#A livello di trascritto
table(mhag$FOCA2_kallisto > mhag$FOCA1_kallisto)

FALSE  TRUE
    5     2
table(mhag$FOCA2_kallisto < mhag$FOCA1_kallisto)

FALSE  TRUE
    2     5


wilcox.test(mhag$FOCA2_kallisto, mhag$FOCA1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$FOCA2_kallisto and mhag$FOCA1_kallisto
W = 3589, p-value = 1
alternative hypothesis: true location shift is greater than 0

#A livello di trascritto
table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
   16    15


table(mhag$Rel_Mut<mhag$Dx_Mut)

FALSE  TRUE
   15    16



 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 465, p-value = 0.592
alternative hypothesis: true location shift is greater than 0

Warning message:
In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
  cannot compute exact p-value with ties



#GAGRA
GAGRA1 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/GAGRA1_Allfiles_sorted.bam
GAGRA2 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/GAGRA2_Allfiles_sorted.bam

awk '{OFS="\t"; print $7":"$8"-"$8}'   /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_GAGRA_RelDiag_MHCII_formatted_NoHeaders.txt | sort | uniq  > Summary_GAGRA_NeoEpitopi_RelDiag_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/GAGRA1_Allfiles_sorted.bam >> GAGRA_NeoAGs_MHCII_RelDiag_Binders_Regions_Diagnosis.txt; done < Summary_GAGRA_NeoEpitopi_RelDiag_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/GAGRA2_Allfiles_sorted.bam >> GAGRA_NeoAGs_MHCII_RelDiag_Binders_Regions_Relapse.txt; done < Summary_GAGRA_NeoEpitopi_RelDiag_MHCII_Regions.txt
 
gagra_NeoAG_mhcii_RelDiag=read.table('/home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_GAGRA_RelDiag_MHCII_formatted_NoHeaders.txt')
 colnames(gagra_NeoAG_mhcii_RelDiag)=c("Gene","Peptide","DRB1","DRB1","Mutation","genes","chr","pos","wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_DRB1","IC50_DRB1")
gagra_NeoAG_mhcii_RelDiag=merge(gagra_NeoAG_mhcii_RelDiag,GAGRA_RNK,by="genes",all.x=T)
relapse=read.table('GAGRA_NeoAGs_MHCII_RelDiag_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('GAGRA_NeoAGs_MHCII_RelDiag_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
gagra_NeoAG_mhcii_RelDiag=merge(gagra_NeoAG_mhcii_RelDiag,relapse,by=c("chr","pos"),all.x=T)
gagra_NeoAG_mhcii_RelDiag=merge(gagra_NeoAG_mhcii_RelDiag,diagnosis,by=c("chr","pos"),all.x=T)
gagra_NeoAG_mhcii_RelDiag$Rel_Mut=str_count(str_to_upper(gagra_NeoAG_mhcii_RelDiag$rel_bases), as.character(gagra_NeoAG_mhcii_RelDiag$mut_nt))
gagra_NeoAG_mhcii_RelDiag$Rel_Wt=str_count(str_to_upper(gagra_NeoAG_mhcii_RelDiag$rel_bases), as.character(gagra_NeoAG_mhcii_RelDiag$wt_nt))
gagra_NeoAG_mhcii_RelDiag$Dx_Mut=str_count(str_to_upper(gagra_NeoAG_mhcii_RelDiag$dx_bases), as.character(gagra_NeoAG_mhcii_RelDiag$mut_nt))
gagra_NeoAG_mhcii_RelDiag$Dx_Wt=str_count(str_to_upper(gagra_NeoAG_mhcii_RelDiag$dx_bases), as.character(gagra_NeoAG_mhcii_RelDiag$wt_nt))
 write.table(unique(gagra_NeoAG_mhcii_RelDiag[,c(1:18,27:30)]),"Summary_GAGRA_Annotato_NeoAGs_MHCII_RelDiag_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)

mhag=read.table('Summary_GAGRA_Annotato_NeoAGs_MHCII_RelDiag_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$GAGRA2_kallisto[is.na(mhag$GAGRA2_kallisto)] <- 0
mhag$GAGRA1_kallisto[is.na(mhag$GAGRA1_kallisto)] <- 0
52

mhag=mhag[(mhag$GAGRA2_kallisto>=0.58 & mhag$GAGRA1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
12


table(mhag$Dx_Mut>=3)
table(mhag$Rel_Mut>=3)

#A livello di trascritto
table(mhag$GAGRA2_kallisto > mhag$GAGRA1_kallisto)

FALSE  TRUE
    5     2
table(mhag$GAGRA2_kallisto < mhag$GAGRA1_kallisto)

FALSE  TRUE
    2     5


wilcox.test(mhag$GAGRA2_kallisto, mhag$GAGRA1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$GAGRA2_kallisto and mhag$GAGRA1_kallisto
W = 3589, p-value = 1
alternative hypothesis: true location shift is greater than 0

#A livello di trascritto
table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
   16    15


table(mhag$Rel_Mut<mhag$Dx_Mut)

FALSE  TRUE
   15    16



 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 465, p-value = 0.592
alternative hypothesis: true location shift is greater than 0

Warning message:
In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
  cannot compute exact p-value with ties


#LUAN
#LUAN RelDiag
LUAN1 = /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN1_Allfiles_sorted.bam
LUAN2 = /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN2_Allfiles_sorted.bam

awk '{OFS="\t"; print $7":"$8"-"$8}'   /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_LUAN_RelDiag_MHCII_formatted_NoHeaders.txt | sort | uniq  > Summary_LUAN_NeoEpitopi_RelDiag_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN1_Allfiles_sorted.bam >> LUAN_NeoAGs_MHCII_RelDiag_Binders_Regions_Diagnosis.txt; done < Summary_LUAN_NeoEpitopi_RelDiag_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN2_Allfiles_sorted.bam >> LUAN_NeoAGs_MHCII_RelDiag_Binders_Regions_Relapse.txt; done < Summary_LUAN_NeoEpitopi_RelDiag_MHCII_Regions.txt
 
luan_NeoAG_mhcii_RelDiag=read.table('/home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_LUAN_RelDiag_MHCII_formatted_NoHeaders.txt')
 colnames(luan_NeoAG_mhcii_RelDiag)=c("Gene","Peptide","DRB1","DRB1","Mutation","genes","chr","pos","wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_DRB1","IC50_DRB1")
luan_NeoAG_mhcii_RelDiag=merge(luan_NeoAG_mhcii_RelDiag,LUAN_RNK,by="genes",all.x=T)
relapse=read.table('LUAN_NeoAGs_MHCII_RelDiag_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('LUAN_NeoAGs_MHCII_RelDiag_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
luan_NeoAG_mhcii_RelDiag=merge(luan_NeoAG_mhcii_RelDiag,relapse,by=c("chr","pos"),all.x=T)
luan_NeoAG_mhcii_RelDiag=merge(luan_NeoAG_mhcii_RelDiag,diagnosis,by=c("chr","pos"),all.x=T)
luan_NeoAG_mhcii_RelDiag$Rel_Mut=str_count(str_to_upper(luan_NeoAG_mhcii_RelDiag$rel_bases), as.character(luan_NeoAG_mhcii_RelDiag$mut_nt))
luan_NeoAG_mhcii_RelDiag$Rel_Wt=str_count(str_to_upper(luan_NeoAG_mhcii_RelDiag$rel_bases), as.character(luan_NeoAG_mhcii_RelDiag$wt_nt))
luan_NeoAG_mhcii_RelDiag$Dx_Mut=str_count(str_to_upper(luan_NeoAG_mhcii_RelDiag$dx_bases), as.character(luan_NeoAG_mhcii_RelDiag$mut_nt))
luan_NeoAG_mhcii_RelDiag$Dx_Wt=str_count(str_to_upper(luan_NeoAG_mhcii_RelDiag$dx_bases), as.character(luan_NeoAG_mhcii_RelDiag$wt_nt))
 write.table(unique(luan_NeoAG_mhcii_RelDiag[,c(1:18,27:30)]),"Summary_LUAN_Annotato_NeoAGs_MHCII_RelDiag_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)

mhag=read.table('Summary_LUAN_Annotato_NeoAGs_MHCII_RelDiag_AllInfos.txt',head=T)
mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
mhag$logFC[is.na(mhag$logFC)] <- 0
mhag$LUAN2_kallisto[is.na(mhag$LUAN2_kallisto)] <- 0
mhag$LUAN1_kallisto[is.na(mhag$LUAN1_kallisto)] <- 0
52

mhag=mhag[(mhag$LUAN2_kallisto>=0.58 & mhag$LUAN1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
12
table(mhag$Dx_Mut>=3)
table(mhag$Rel_Mut>=3)


#A livello di trascritto
table(mhag$LUAN2_kallisto > mhag$LUAN1_kallisto)

FALSE  TRUE
    5     2
table(mhag$LUAN2_kallisto < mhag$LUAN1_kallisto)

FALSE  TRUE
    2     5


wilcox.test(mhag$LUAN2_kallisto, mhag$LUAN1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$LUAN2_kallisto and mhag$LUAN1_kallisto
W = 3589, p-value = 1
alternative hypothesis: true location shift is greater than 0

#A livello di trascritto
table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
   16    15


table(mhag$Rel_Mut<mhag$Dx_Mut)

FALSE  TRUE
   15    16



 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 465, p-value = 0.592
alternative hypothesis: true location shift is greater than 0

Warning message:
In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
  cannot compute exact p-value with ties


#Luan RelOnly
LUAN1 = /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN1_Allfiles_sorted.bam
LUAN2 = /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN2_Allfiles_sorted.bam

awk '{OFS="\t"; print $7":"$8"-"$8}'   /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_LUAN_RelOnly_MHCII_formatted_NoHeaders.txt | sort | uniq  > Summary_LUAN_NeoEpitopi_RelOnly_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN1_Allfiles_sorted.bam >> LUAN_NeoAGs_MHCII_RelOnly_Binders_Regions_Diagnosis.txt; done < Summary_LUAN_NeoEpitopi_RelOnly_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN2_Allfiles_sorted.bam >> LUAN_NeoAGs_MHCII_RelOnly_Binders_Regions_Relapse.txt; done < Summary_LUAN_NeoEpitopi_RelOnly_MHCII_Regions.txt
 
luan_NeoAG_mhcii_RelOnly=read.table('/home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_LUAN_RelOnly_MHCII_formatted_NoHeaders.txt')
 colnames(luan_NeoAG_mhcii_RelOnly)=c("Gene","Peptide","DRB1","DRB1","Mutation","genes","chr","pos","wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_DRB1","IC50_DRB1")
luan_NeoAG_mhcii_RelOnly=merge(luan_NeoAG_mhcii_RelOnly,LUAN_RNK,by="genes",all.x=T)
relapse=read.table('LUAN_NeoAGs_MHCII_RelOnly_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('LUAN_NeoAGs_MHCII_RelOnly_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
luan_NeoAG_mhcii_RelOnly=merge(luan_NeoAG_mhcii_RelOnly,relapse,by=c("chr","pos"),all.x=T)
luan_NeoAG_mhcii_RelOnly=merge(luan_NeoAG_mhcii_RelOnly,diagnosis,by=c("chr","pos"),all.x=T)
luan_NeoAG_mhcii_RelOnly$Rel_Mut=str_count(str_to_upper(luan_NeoAG_mhcii_RelOnly$rel_bases), as.character(luan_NeoAG_mhcii_RelOnly$mut_nt))
luan_NeoAG_mhcii_RelOnly$Rel_Wt=str_count(str_to_upper(luan_NeoAG_mhcii_RelOnly$rel_bases), as.character(luan_NeoAG_mhcii_RelOnly$wt_nt))
luan_NeoAG_mhcii_RelOnly$Dx_Mut=str_count(str_to_upper(luan_NeoAG_mhcii_RelOnly$dx_bases), as.character(luan_NeoAG_mhcii_RelOnly$mut_nt))
luan_NeoAG_mhcii_RelOnly$Dx_Wt=str_count(str_to_upper(luan_NeoAG_mhcii_RelOnly$dx_bases), as.character(luan_NeoAG_mhcii_RelOnly$wt_nt))
 write.table(unique(luan_NeoAG_mhcii_RelOnly[,c(1:18,27:30)]),"Summary_LUAN_Annotato_NeoAGs_MHCII_RelOnly_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)

	mhag=read.table('Summary_LUAN_Annotato_NeoAGs_MHCII_RelOnly_AllInfos.txt',head=T)
	mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
	mhag$logFC[is.na(mhag$logFC)] <- 0
	mhag$LUAN2_kallisto[is.na(mhag$LUAN2_kallisto)] <- 0
	mhag$LUAN1_kallisto[is.na(mhag$LUAN1_kallisto)] <- 0
52

mhag=mhag[(mhag$LUAN2_kallisto>=0.58 & mhag$LUAN1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
12

table(mhag$Dx_Mut>=3)
table(mhag$Rel_Mut>=3)

#A livello di trascritto
table(mhag$LUAN2_kallisto > mhag$LUAN1_kallisto)

FALSE  TRUE
    5     2
table(mhag$LUAN2_kallisto < mhag$LUAN1_kallisto)

FALSE  TRUE
    2     5


wilcox.test(mhag$LUAN2_kallisto, mhag$LUAN1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$LUAN2_kallisto and mhag$LUAN1_kallisto
W = 3589, p-value = 1
alternative hypothesis: true location shift is greater than 0

#A livello di trascritto
table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
   16    15


table(mhag$Rel_Mut<mhag$Dx_Mut)

FALSE  TRUE
   15    16



 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 465, p-value = 0.592
alternative hypothesis: true location shift is greater than 0

Warning message:
In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
  cannot compute exact p-value with ties

#MABI
#MABI RelDiag
MABI1 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/MABI1_Allfiles_sorted.bam
MABI2 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/MABI2_Allfiles_sorted.bam

awk '{OFS="\t"; print $6":"$7"-"$7}'   /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_MABI_RelDiag_MHCII_formatted_NoHeaders.txt | sort | uniq  > Summary_MABI_NeoEpitopi_RelDiag_MHCII_Regions.txt

 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/MABI1_Allfiles_sorted.bam >> MABI_NeoAGs_MHCII_RelDiag_Binders_Regions_Diagnosis.txt; done < Summary_MABI_NeoEpitopi_RelDiag_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/MABI2_Allfiles_sorted.bam >> MABI_NeoAGs_MHCII_RelDiag_Binders_Regions_Relapse.txt; done < Summary_MABI_NeoEpitopi_RelDiag_MHCII_Regions.txt
 
mabi_NeoAG_mhcii_RelDiag=read.table('/home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_MABI_RelDiag_MHCII_formatted_NoHeaders.txt')
 colnames(mabi_NeoAG_mhcii_RelDiag)=c("Gene","Peptide","DRB1","Mutation","genes","chr","pos","wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_DRB1")
mabi_NeoAG_mhcii_RelDiag=merge(mabi_NeoAG_mhcii_RelDiag,MABI_RNK,by="genes",all.x=T)
relapse=read.table('MABI_NeoAGs_MHCII_RelDiag_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('MABI_NeoAGs_MHCII_RelDiag_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
mabi_NeoAG_mhcii_RelDiag=merge(mabi_NeoAG_mhcii_RelDiag,relapse,by=c("chr","pos"),all.x=T)
mabi_NeoAG_mhcii_RelDiag=merge(mabi_NeoAG_mhcii_RelDiag,diagnosis,by=c("chr","pos"),all.x=T)
mabi_NeoAG_mhcii_RelDiag$Rel_Mut=str_count(str_to_upper(mabi_NeoAG_mhcii_RelDiag$rel_bases), as.character(mabi_NeoAG_mhcii_RelDiag$mut_nt))
mabi_NeoAG_mhcii_RelDiag$Rel_Wt=str_count(str_to_upper(mabi_NeoAG_mhcii_RelDiag$rel_bases), as.character(mabi_NeoAG_mhcii_RelDiag$wt_nt))
mabi_NeoAG_mhcii_RelDiag$Dx_Mut=str_count(str_to_upper(mabi_NeoAG_mhcii_RelDiag$dx_bases), as.character(mabi_NeoAG_mhcii_RelDiag$mut_nt))
mabi_NeoAG_mhcii_RelDiag$Dx_Wt=str_count(str_to_upper(mabi_NeoAG_mhcii_RelDiag$dx_bases), as.character(mabi_NeoAG_mhcii_RelDiag$wt_nt))
 write.table(unique(mabi_NeoAG_mhcii_RelDiag[,c(1:16,25:28)]),"Summary_MABI_Annotato_NeoAGs_MHCII_RelDiag_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)

	mhag=read.table('Summary_MABI_Annotato_NeoAGs_MHCII_RelDiag_AllInfos.txt',head=T)
	mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
	mhag$logFC[is.na(mhag$logFC)] <- 0
	mhag$MABI2_kallisto[is.na(mhag$MABI2_kallisto)] <- 0
	mhag$MABI1_kallisto[is.na(mhag$MABI1_kallisto)] <- 0
52

mhag=mhag[(mhag$MABI2_kallisto>=0.58 & mhag$MABI1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
12

table(mhag$Dx_Mut>=3)
table(mhag$Rel_Mut>=3)

#A livello di trascritto
table(mhag$MABI2_kallisto > mhag$MABI1_kallisto)

FALSE  TRUE
    5     2
table(mhag$MABI2_kallisto < mhag$MABI1_kallisto)

FALSE  TRUE
    2     5


wilcox.test(mhag$MABI2_kallisto, mhag$MABI1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$MABI2_kallisto and mhag$MABI1_kallisto
W = 3589, p-value = 1
alternative hypothesis: true location shift is greater than 0

#A livello di trascritto
table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
   16    15


table(mhag$Rel_Mut<mhag$Dx_Mut)

FALSE  TRUE
   15    16



 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 465, p-value = 0.592
alternative hypothesis: true location shift is greater than 0

Warning message:
In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
  cannot compute exact p-value with ties
  

#MOGE
#MOGE RelDiag
MOGE1 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE1_Allfiles_sorted.bam
MOGE2 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE2_Allfiles_sorted.bam

awk '{OFS="\t"; print $7":"$8"-"$8}'   /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_MOGE_RelDiag_MHCII_formatted_NoHeaders.txt | sort | uniq  > Summary_MOGE_NeoEpitopi_RelDiag_MHCII_Regions.txt

 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE1_Allfiles_sorted.bam >> MOGE_NeoAGs_MHCII_RelDiag_Binders_Regions_Diagnosis.txt; done < Summary_MOGE_NeoEpitopi_RelDiag_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE2_Allfiles_sorted.bam >> MOGE_NeoAGs_MHCII_RelDiag_Binders_Regions_Relapse.txt; done < Summary_MOGE_NeoEpitopi_RelDiag_MHCII_Regions.txt
 
moge_NeoAG_mhcii_RelDiag=read.table('/home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_MOGE_RelDiag_MHCII_formatted_NoHeaders.txt')
 colnames(moge_NeoAG_mhcii_RelDiag)=c("Gene","Peptide","DRB1","DRB1","Mutation","genes","chr","pos","wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_DRB1","IC50_DRB1")
moge_NeoAG_mhcii_RelDiag=merge(moge_NeoAG_mhcii_RelDiag,MOGE_RNK,by="genes",all.x=T)
relapse=read.table('MOGE_NeoAGs_MHCII_RelDiag_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('MOGE_NeoAGs_MHCII_RelDiag_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
moge_NeoAG_mhcii_RelDiag=merge(moge_NeoAG_mhcii_RelDiag,relapse,by=c("chr","pos"),all.x=T)
moge_NeoAG_mhcii_RelDiag=merge(moge_NeoAG_mhcii_RelDiag,diagnosis,by=c("chr","pos"),all.x=T)
moge_NeoAG_mhcii_RelDiag$Rel_Mut=str_count(str_to_upper(moge_NeoAG_mhcii_RelDiag$rel_bases), as.character(moge_NeoAG_mhcii_RelDiag$mut_nt))
moge_NeoAG_mhcii_RelDiag$Rel_Wt=str_count(str_to_upper(moge_NeoAG_mhcii_RelDiag$rel_bases), as.character(moge_NeoAG_mhcii_RelDiag$wt_nt))
moge_NeoAG_mhcii_RelDiag$Dx_Mut=str_count(str_to_upper(moge_NeoAG_mhcii_RelDiag$dx_bases), as.character(moge_NeoAG_mhcii_RelDiag$mut_nt))
moge_NeoAG_mhcii_RelDiag$Dx_Wt=str_count(str_to_upper(moge_NeoAG_mhcii_RelDiag$dx_bases), as.character(moge_NeoAG_mhcii_RelDiag$wt_nt))
 write.table(unique(moge_NeoAG_mhcii_RelDiag[,c(1:18,27:30)]),"Summary_MOGE_Annotato_NeoAGs_MHCII_RelDiag_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)

	mhag=read.table('Summary_MOGE_Annotato_NeoAGs_MHCII_RelDiag_AllInfos.txt',head=T)
	mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
	mhag$logFC[is.na(mhag$logFC)] <- 0
	mhag$MOGE2_kallisto[is.na(mhag$MOGE2_kallisto)] <- 0
	mhag$MOGE1_kallisto[is.na(mhag$MOGE1_kallisto)] <- 0
52

mhag=mhag[(mhag$MOGE2_kallisto>=0.58 & mhag$MOGE1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
12


table(mhag$Dx_Mut>=3)
table(mhag$Rel_Mut>=3)

#A livello di trascritto
table(mhag$MOGE2_kallisto > mhag$MOGE1_kallisto)

FALSE  TRUE
    5     2
table(mhag$MOGE2_kallisto < mhag$MOGE1_kallisto)

FALSE  TRUE
    2     5


wilcox.test(mhag$MOGE2_kallisto, mhag$MOGE1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$MOGE2_kallisto and mhag$MOGE1_kallisto
W = 3589, p-value = 1
alternative hypothesis: true location shift is greater than 0

#A livello di trascritto
table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
   16    15


table(mhag$Rel_Mut<mhag$Dx_Mut)

FALSE  TRUE
   15    16



 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 465, p-value = 0.592
alternative hypothesis: true location shift is greater than 0

Warning message:
In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
  cannot compute exact p-value with ties
  

#MOGE RelOnly
MOGE1 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE1_Allfiles_sorted.bam
MOGE2 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE2_Allfiles_sorted.bam

 sed -e '/Strong/,+1d' /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_MOGE_RelOnly_MHCII.txt  | sed -e '/Weak/,+1d' - | grep -v WT_ > tt && mv tt /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_MOGE_RelOnly_MHCII_formatted_NoHeaders.txt
 
awk '{OFS="\t"; print $7":"$8"-"$8}'   /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_MOGE_RelOnly_MHCII_formatted_NoHeaders.txt | sort | uniq  > Summary_MOGE_NeoEpitopi_RelOnly_MHCII_Regions.txt

 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE1_Allfiles_sorted.bam >> MOGE_NeoAGs_MHCII_RelOnly_Binders_Regions_Diagnosis.txt; done < Summary_MOGE_NeoEpitopi_RelOnly_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE2_Allfiles_sorted.bam >> MOGE_NeoAGs_MHCII_RelOnly_Binders_Regions_Relapse.txt; done < Summary_MOGE_NeoEpitopi_RelOnly_MHCII_Regions.txt
 
moge_NeoAG_mhcii_RelOnly=read.table('/home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_MOGE_RelOnly_MHCII_formatted_NoHeaders.txt')
 colnames(moge_NeoAG_mhcii_RelOnly)=c("Gene","Peptide","DRB1","DRB1","Mutation","genes","chr","pos","wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_DRB1","IC50_DRB1")
moge_NeoAG_mhcii_RelOnly=merge(moge_NeoAG_mhcii_RelOnly,MOGE_RNK,by="genes",all.x=T)
relapse=read.table('MOGE_NeoAGs_MHCII_RelOnly_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('MOGE_NeoAGs_MHCII_RelOnly_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
moge_NeoAG_mhcii_RelOnly=merge(moge_NeoAG_mhcii_RelOnly,relapse,by=c("chr","pos"),all.x=T)
moge_NeoAG_mhcii_RelOnly=merge(moge_NeoAG_mhcii_RelOnly,diagnosis,by=c("chr","pos"),all.x=T)
moge_NeoAG_mhcii_RelOnly$mut_nt="T"
moge_NeoAG_mhcii_RelOnly$Rel_Mut=str_count(str_to_upper(moge_NeoAG_mhcii_RelOnly$rel_bases), as.character(moge_NeoAG_mhcii_RelOnly$mut_nt))
moge_NeoAG_mhcii_RelOnly$Rel_Wt=str_count(str_to_upper(moge_NeoAG_mhcii_RelOnly$rel_bases), as.character(moge_NeoAG_mhcii_RelOnly$wt_nt))
moge_NeoAG_mhcii_RelOnly$Dx_Mut=str_count(str_to_upper(moge_NeoAG_mhcii_RelOnly$dx_bases), as.character(moge_NeoAG_mhcii_RelOnly$mut_nt))
moge_NeoAG_mhcii_RelOnly$Dx_Wt=str_count(str_to_upper(moge_NeoAG_mhcii_RelOnly$dx_bases), as.character(moge_NeoAG_mhcii_RelOnly$wt_nt))
 write.table(unique(moge_NeoAG_mhcii_RelOnly[,c(1:18,27:30)]),"Summary_MOGE_Annotato_NeoAGs_MHCII_RelOnly_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)

	mhag=read.table('Summary_MOGE_Annotato_NeoAGs_MHCII_RelOnly_AllInfos.txt',head=T)
	mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
	mhag$logFC[is.na(mhag$logFC)] <- 0
	mhag$MOGE2_kallisto[is.na(mhag$MOGE2_kallisto)] <- 0
	mhag$MOGE1_kallisto[is.na(mhag$MOGE1_kallisto)] <- 0
52

mhag=mhag[(mhag$MOGE2_kallisto>=0.58 & mhag$MOGE1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
12

table(mhag$Dx_Mut>=3)
table(mhag$Rel_Mut>=3)

#A livello di trascritto
table(mhag$MOGE2_kallisto > mhag$MOGE1_kallisto)

FALSE  TRUE
    5     2
table(mhag$MOGE2_kallisto < mhag$MOGE1_kallisto)

FALSE  TRUE
    2     5


wilcox.test(mhag$MOGE2_kallisto, mhag$MOGE1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$MOGE2_kallisto and mhag$MOGE1_kallisto
W = 3589, p-value = 1
alternative hypothesis: true location shift is greater than 0

#A livello di trascritto
table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
   16    15


table(mhag$Rel_Mut<mhag$Dx_Mut)

FALSE  TRUE
   15    16



 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 465, p-value = 0.592
alternative hypothesis: true location shift is greater than 0

Warning message:
In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
  cannot compute exact p-value with ties
  

#PIAG
#PIAG RelOnly
PIAG1 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/PIAG1_Allfiles_sorted.bam
PIAG2 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/PIAG2_Allfiles_sorted.bam

 sed -e '/Strong/,+1d' /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_PIAG_RelOnly_MHCII.txt  | sed -e '/Weak/,+1d' - | grep -v WT_ > tt && mv tt /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_PIAG_RelOnly_MHCII_formatted_NoHeaders.txt

awk '{OFS="\t"; print $6":"$7"-"$7}'   /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_PIAG_RelOnly_MHCII_formatted_NoHeaders.txt | sort | uniq  > Summary_PIAG_NeoEpitopi_RelOnly_MHCII_Regions.txt

 while read region;  do samtools mpileup -r $region  //lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/PIAG1_Allfiles_sorted.bam >> PIAG_NeoAGs_MHCII_RelOnly_Binders_Regions_Diagnosis.txt; done < Summary_PIAG_NeoEpitopi_RelOnly_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/PIAG2_Allfiles_sorted.bam >> PIAG_NeoAGs_MHCII_RelOnly_Binders_Regions_Relapse.txt; done < Summary_PIAG_NeoEpitopi_RelOnly_MHCII_Regions.txt
 
piag_NeoAG_mhcii_RelOnly=read.table('/home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_PIAG_RelOnly_MHCII_formatted_NoHeaders.txt')
 colnames(piag_NeoAG_mhcii_RelOnly)=c("Gene","Peptide","DRB1","Mutation","genes","chr","pos","wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_DRB1")
piag_NeoAG_mhcii_RelOnly=merge(piag_NeoAG_mhcii_RelOnly,PIAG_RNK,by="genes",all.x=T)
relapse=read.table('PIAG_NeoAGs_MHCII_RelOnly_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('PIAG_NeoAGs_MHCII_RelOnly_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
piag_NeoAG_mhcii_RelOnly=merge(piag_NeoAG_mhcii_RelOnly,relapse,by=c("chr","pos"),all.x=T)
piag_NeoAG_mhcii_RelOnly=merge(piag_NeoAG_mhcii_RelOnly,diagnosis,by=c("chr","pos"),all.x=T)
piag_NeoAG_mhcii_RelOnly$Rel_Mut=str_count(str_to_upper(piag_NeoAG_mhcii_RelOnly$rel_bases), as.character(piag_NeoAG_mhcii_RelOnly$mut_nt))
piag_NeoAG_mhcii_RelOnly$Rel_Wt=str_count(str_to_upper(piag_NeoAG_mhcii_RelOnly$rel_bases), as.character(piag_NeoAG_mhcii_RelOnly$wt_nt))
piag_NeoAG_mhcii_RelOnly$Dx_Mut=str_count(str_to_upper(piag_NeoAG_mhcii_RelOnly$dx_bases), as.character(piag_NeoAG_mhcii_RelOnly$mut_nt))
piag_NeoAG_mhcii_RelOnly$Dx_Wt=str_count(str_to_upper(piag_NeoAG_mhcii_RelOnly$dx_bases), as.character(piag_NeoAG_mhcii_RelOnly$wt_nt))
 write.table(unique(piag_NeoAG_mhcii_RelOnly[,c(1:16,25:28)]),"Summary_PIAG_Annotato_NeoAGs_MHCII_RelOnly_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)


	mhag=read.table('Summary_PIAG_Annotato_NeoAGs_MHCII_RelOnly_AllInfos.txt',head=T)
	mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
	mhag$logFC[is.na(mhag$logFC)] <- 0
	mhag$PIAG2_kallisto[is.na(mhag$PIAG2_kallisto)] <- 0
	mhag$PIAG1_kallisto[is.na(mhag$PIAG1_kallisto)] <- 0
52

mhag=mhag[(mhag$PIAG2_kallisto>=0.58 & mhag$PIAG1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
12

table(mhag$Dx_Mut>=3)
table(mhag$Rel_Mut>=3)

#A livello di trascritto
table(mhag$PIAG2_kallisto > mhag$PIAG1_kallisto)

FALSE  TRUE
    5     2
table(mhag$PIAG2_kallisto < mhag$PIAG1_kallisto)

FALSE  TRUE
    2     5


wilcox.test(mhag$PIAG2_kallisto, mhag$PIAG1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$PIAG2_kallisto and mhag$PIAG1_kallisto
W = 3589, p-value = 1
alternative hypothesis: true location shift is greater than 0

#A livello di trascritto
table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
   16    15


table(mhag$Rel_Mut<mhag$Dx_Mut)

FALSE  TRUE
   15    16



 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 465, p-value = 0.592
alternative hypothesis: true location shift is greater than 0

Warning message:
In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
  cannot compute exact p-value with ties
  




#PIAG DiagOnly

 sed -e '/Strong/,+1d' /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_PIAG_DiagOnly_MHCII.txt  | sed -e '/Weak/,+1d' - | grep -v WT_ > tt && mv tt /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_PIAG_DiagOnly_MHCII_formatted_NoHeaders.txt


awk '{OFS="\t"; print $6":"$7"-"$7}'   /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_PIAG_DiagOnly_MHCII_formatted_NoHeaders.txt | sort | uniq  > Summary_PIAG_NeoEpitopi_DiagOnly_MHCII_Regions.txt

 while read region;  do samtools mpileup -r $region  //lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/PIAG1_Allfiles_sorted.bam >> PIAG_NeoAGs_MHCII_DiagOnly_Binders_Regions_Diagnosis.txt; done < Summary_PIAG_NeoEpitopi_DiagOnly_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/PIAG2_Allfiles_sorted.bam >> PIAG_NeoAGs_MHCII_DiagOnly_Binders_Regions_Relapse.txt; done < Summary_PIAG_NeoEpitopi_DiagOnly_MHCII_Regions.txt
 
piag_NeoAG_mhcii_DiagOnly=read.table('/home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_PIAG_DiagOnly_MHCII_formatted_NoHeaders.txt')
 colnames(piag_NeoAG_mhcii_DiagOnly)=c("Gene","Peptide","DRB1","Mutation","genes","chr","pos","wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_DRB1")
piag_NeoAG_mhcii_DiagOnly=merge(piag_NeoAG_mhcii_DiagOnly,PIAG_RNK,by="genes",all.x=T)

#PIAG DiagOnly non ha niente di espresso


	mhag=read.table('/home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_PIAG_DiagOnly_MHCII_formatted_NoHeaders.txt',head=T)
	mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
	mhag$logFC[is.na(mhag$logFC)] <- 0
	mhag$PIAG2_kallisto[is.na(mhag$PIAG2_kallisto)] <- 0
	mhag$PIAG1_kallisto[is.na(mhag$PIAG1_kallisto)] <- 0
52

mhag=mhag[(mhag$PIAG2_kallisto>=0.58 & mhag$PIAG1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
12

table(mhag$Dx_Mut>=3)
table(mhag$Rel_Mut>=3)

#A livello di trascritto
table(mhag$PIAG2_kallisto > mhag$PIAG1_kallisto)

FALSE  TRUE
    5     2
table(mhag$PIAG2_kallisto < mhag$PIAG1_kallisto)

FALSE  TRUE
    2     5


wilcox.test(mhag$PIAG2_kallisto, mhag$PIAG1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$PIAG2_kallisto and mhag$PIAG1_kallisto
W = 3589, p-value = 1
alternative hypothesis: true location shift is greater than 0

#A livello di trascritto
table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
   16    15


table(mhag$Rel_Mut<mhag$Dx_Mut)

FALSE  TRUE
   15    16



 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 465, p-value = 0.592
alternative hypothesis: true location shift is greater than 0

Warning message:
In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
  cannot compute exact p-value with ties
  



#PRELU
PRELU1 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/PRELU1_Allfiles_sorted.bam
PRELU2 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/PRELU2_Allfiles_sorted.bam

awk '{OFS="\t"; print $7":"$8"-"$8}'   /home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_PRELU_RelDiag_MHCII_formatted_NoHeaders.txt | sort | uniq  > Summary_PRELU_NeoEpitopi_RelDiag_MHCII_Regions.txt

 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/PRELU1_Allfiles_sorted.bam >> PRELU_NeoAGs_MHCII_RelDiag_Binders_Regions_Diagnosis.txt; done < Summary_PRELU_NeoEpitopi_RelDiag_MHCII_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/PRELU2_Allfiles_sorted.bam >> PRELU_NeoAGs_MHCII_RelDiag_Binders_Regions_Relapse.txt; done < Summary_PRELU_NeoEpitopi_RelDiag_MHCII_Regions.txt
 
prelu_NeoAG_mhcii_RelDiag=read.table('/home/fsantaniello/Summaries_NeoEpitopi/MHCII/Def/Summary_PRELU_RelDiag_MHCII_formatted_NoHeaders.txt')
 colnames(prelu_NeoAG_mhcii_RelDiag)=c("Gene","Peptide","DRB1","DRB1","Mutation","genes","chr","pos","wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_DRB1","IC50_DRB1")
prelu_NeoAG_mhcii_RelDiag=merge(prelu_NeoAG_mhcii_RelDiag,PRELU_RNK,by="genes",all.x=T)
relapse=read.table('PRELU_NeoAGs_MHCII_RelDiag_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('PRELU_NeoAGs_MHCII_RelDiag_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
prelu_NeoAG_mhcii_RelDiag=merge(prelu_NeoAG_mhcii_RelDiag,relapse,by=c("chr","pos"),all.x=T)
prelu_NeoAG_mhcii_RelDiag=merge(prelu_NeoAG_mhcii_RelDiag,diagnosis,by=c("chr","pos"),all.x=T)
prelu_NeoAG_mhcii_RelDiag$Rel_Mut=str_count(str_to_upper(prelu_NeoAG_mhcii_RelDiag$rel_bases), as.character(prelu_NeoAG_mhcii_RelDiag$mut_nt))
prelu_NeoAG_mhcii_RelDiag$Rel_Wt=str_count(str_to_upper(prelu_NeoAG_mhcii_RelDiag$rel_bases), as.character(prelu_NeoAG_mhcii_RelDiag$wt_nt))
prelu_NeoAG_mhcii_RelDiag$Dx_Mut=str_count(str_to_upper(prelu_NeoAG_mhcii_RelDiag$dx_bases), as.character(prelu_NeoAG_mhcii_RelDiag$mut_nt))
prelu_NeoAG_mhcii_RelDiag$Dx_Wt=str_count(str_to_upper(prelu_NeoAG_mhcii_RelDiag$dx_bases), as.character(prelu_NeoAG_mhcii_RelDiag$wt_nt))
 write.table(unique(prelu_NeoAG_mhcii_RelDiag[,c(1:18,27:30)]),"Summary_PRELU_Annotato_NeoAGs_MHCII_RelDiag_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)
 
 	mhag=read.table('Summary_PRELU_Annotato_NeoAGs_MHCII_RelDiag_AllInfos.txt',head=T)
	mhag$Qvalue[is.na(mhag$Qvalue)] <- 1
	mhag$logFC[is.na(mhag$logFC)] <- 0
	mhag$PRELU2_kallisto[is.na(mhag$PRELU2_kallisto)] <- 0
	mhag$PRELU1_kallisto[is.na(mhag$PRELU1_kallisto)] <- 0
52

mhag=mhag[(mhag$PRELU2_kallisto>=0.58 & mhag$PRELU1_kallisto>=0.58) & (mhag$Rel_Mut>=3 | mhag$Dx_Mut>=3 ),]
nrow(mhag)
12

table(mhag$Dx_Mut>=3)
table(mhag$Rel_Mut>=3)


#A livello di trascritto
table(mhag$PRELU2_kallisto > mhag$PRELU1_kallisto)

FALSE  TRUE
    5     2
table(mhag$PRELU2_kallisto < mhag$PRELU1_kallisto)

FALSE  TRUE
    2     5


wilcox.test(mhag$PRELU2_kallisto, mhag$PRELU1_kallisto,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$PRELU2_kallisto and mhag$PRELU1_kallisto
W = 3589, p-value = 1
alternative hypothesis: true location shift is greater than 0

#A livello di trascritto
table(mhag$Rel_Mut>mhag$Dx_Mut)

FALSE  TRUE
   16    15


table(mhag$Rel_Mut<mhag$Dx_Mut)

FALSE  TRUE
   15    16



 wilcox.test(mhag$Rel_Mut, mhag$Dx_Mut,alternative="greater")

        Wilcoxon rank sum test with continuity correction

data:  mhag$Rel_Mut and mhag$Dx_Mut
W = 465, p-value = 0.592
alternative hypothesis: true location shift is greater than 0

Warning message:
In wilcox.test.default(mhag$Rel_Mut, mhag$Dx_Mut, alternative = "greater") :
  cannot compute exact p-value with ties
  
 
 
 Deiv_MinorAntigens_MHCII.sh
 
#PBS -l select=1:app=java:ncpus=1:mem=40gb
#PBS -P 161
#PBS -q workq
#PBS -m ae

condactivate
cd /home/fsantaniello
sh Launcher_snp2epi_Paper_Annotation_MinorAGs_MHCII.sh DEIV_Germline_vs_Donor.tsv DEIV1_DEIV2_DEIV2b_MinorAGs_new_MHCII DRB1_0102,DRB1_1501 DEIV_Expression_NeoEpitopes.txt
condeactivate


DEIV_Germline_vs_Donor.tsv

DEIV_Germline_vs_Donor.tsv DEIV_MinorAGs_new HLA-A11:01,HLA-A24:02,HLA-B14:02,HLA-B35:03,HLA-C08:02,HLA-C12:03 DEIV_Expression_NeoEpitopes.txt

DEIV_WeakBinders_mhAGs_MHCII.sh

#PBS -l select=1:app=java:ncpus=4:mem=12gb
#PBS -P 161
#PBS -q workq
#PBS -m ae

python Summary_Epitopi_New_MHCII.py DEIV1_DEIV2_DEIV2b_MinorAGs_new_MHCII_PeptidesPrediction_MHCII.txt RESULTS_DEIV1_DEIV2_DEIV2b_MinorAGs_new_MHCII_PeptidesPrediction_MinorAntigens_MHCII.xls RESULTS_DEIV2_DEIV2b_Somatic_MinorAGs_new_MHCII_weak_binders_MHC_II_2 DEIV_WeakBinders_Annotati.txt



