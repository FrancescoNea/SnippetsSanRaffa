#PBS -l select=1:app=java:ncpus=4:mem=12gb
#PBS -P 161
#PBS -q workq
#PBS -m ae


#*************************************************************************************************************************************************
#*************************************************************************************************************************************************
#***************************************************************HLA classe I**********************************************************************
#*************************************************************************************************************************************************
#*************************************************************************************************************************************************





#ALFE
sh Launcher_snp2epi_Paper_Annotation.sh ALFE_Somatic_Rerun.tsv ALFE1_ALFE2_MHCI_NoMismatchHla HLA-A24:02,HLA-B35:01,HLA-C04:01 ALFE_Expression_Kallisto_NeoEpitopes.txt

        #Summary Binders Diagnosis Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A24:02\tHLA-B35:01\tHLA-C04:01" )
        more +2 ./RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0; if ($11>50&&$11<=500) print $3,$2,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_weak_binders_DiagOnly
#         /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/weak_binders_DiagOnly ./RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls ./RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_weak_binders_annotated_DiagOnly
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A24:02\tHLA-B35:01\tHLA-C04:01")
        more +2 ./RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls |   awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0; if ($11>0&&$11<=50) print $3,$2,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_strong_binders_DiagOnly
 #        /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/strong_binders_DiagOnly ./RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls /home/fsantaniello/strong_binders_annotated_DiagOnly

        #cat /home/fsantaniello/strong_binders_annotated_DiagOnly /home/fsantaniello/weak_binders_annotated_DiagOnly > ./Summary_Binders_RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla__PeptidesPrediction_DiagnosisOnly.txt



        #Summary Binders Relapse Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A24:02\tHLA-B35:01\tHLA-C04:01")
        more +2 ./RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0; if ($11>50&&$11<=500) print $3,$2,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_weak_binders_RelOnly
        # /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/weak_binders_RelOnly ./RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/weak_binders_annotated_RelOnly

        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A24:02\tHLA-B35:01\tHLA-C04:01" )
        more +2 ./RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0; if ($11>0&&$11<=50) print $3,$2,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_strong_binders_RelOnly
         #/home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/strong_binders_RelOnly ./RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/strong_binders_annotated_RelOnly
        
        #cat /home/fsantaniello/strong_binders_annotated_RelOnly /home/fsantaniello/weak_binders_annotated_RelOnly > ./Summary_Binders_RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla__PeptidesPrediction_RelapseOnly.txt
        
        
        #Summary Binders Common Diagnosis Relapse
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A24:02\tHLA-B35:01\tHLA-C04:01")
        more +2 ./RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0; if ($11>50&&$11<=500) print $3,$2,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_weak_binders_RelDiag
  #       /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/weak_binders_RelDiag ./RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/weak_binders_annotated_RelDiag
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A24:02\tHLA-B35:01\tHLA-C04:01")
        more +2 ./RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0; if ($11>0&&$11<=50) print $3,$2,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_strong_binders_RelDiag
   #      /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/strong_binders_RelDiag ./RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/strong_binders_annotated_RelDiag
       
        #cat /home/fsantaniello/strong_binders_annotated_RelDiag /home/fsantaniello/weak_binders_annotated_RelDiag > ./Summary_Binders_RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla__PeptidesPrediction_CommonRelapseDiagnosisOnly.txt
        

import sys

#Common diagnosis relapse

Mutationi=open('./HomeQlogin/ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_weak_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_ALFE_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_strong_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_ALFE_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()
        

#Relapse Only
Mutationi=open('./HomeQlogin/ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_weak_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_ALFE_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_strong_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_ALFE_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()

#Diagnosis Only
Mutationi=open('./HomeQlogin/ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_weak_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_ALFE_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_strong_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_ALFE_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()



sh Launcher_snp2epi_Paper_Annotation.sh BESU1_BESU2_Somatic.Anno.tsv BESU1_BESU2_MHCI_NoMismatchHla HLA-A02:01,HLA-A68:01,HLA-B18:01,HLA-B56:01,HLA-C07:01,HLA-C01:02 BESU_Expression_Kallisto_NeoEpitopes.txt

       #Summary Binders Diagnosis Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A02:01\tHLA-A68:01\tHLA-B18:01\tHLA-B56:01\tHLA-C07:01\tHLA-C01:02" )
        more +2 ./RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_weak_binders_DiagOnly
#         /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/weak_binders_DiagOnly ./RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls ./RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_weak_binders_annotated_DiagOnly
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A02:01\tHLA-A68:01\tHLA-B18:01\tHLA-B56:01\tHLA-C07:01\tHLA-C01:02" )
        more +2 ./RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls |   awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_strong_binders_DiagOnly
 #        /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/strong_binders_DiagOnly ./RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls /home/fsantaniello/strong_binders_annotated_DiagOnly

        #cat /home/fsantaniello/strong_binders_annotated_DiagOnly /home/fsantaniello/weak_binders_annotated_DiagOnly > ./Summary_Binders_RESULTS_BESU1_BESU2_MHCI_NoMismatchHla__PeptidesPrediction_DiagnosisOnly.txt



        #Summary Binders Relapse Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A02:01\tHLA-A68:01\tHLA-B18:01\tHLA-B56:01\tHLA-C07:01\tHLA-C01:02" )
        more +2 ./RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_weak_binders_RelOnly
        # /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/weak_binders_RelOnly ./RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/weak_binders_annotated_RelOnly

        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A02:01\tHLA-A68:01\tHLA-B18:01\tHLA-B56:01\tHLA-C07:01\tHLA-C01:02" )
        more +2 ./RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_strong_binders_RelOnly
         #/home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/strong_binders_RelOnly ./RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/strong_binders_annotated_RelOnly
        
        #cat /home/fsantaniello/strong_binders_annotated_RelOnly /home/fsantaniello/weak_binders_annotated_RelOnly > ./Summary_Binders_RESULTS_BESU1_BESU2_MHCI_NoMismatchHla__PeptidesPrediction_RelapseOnly.txt
        
        
        #Summary Binders Common Diagnosis Relapse
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A02:01\tHLA-A68:01\tHLA-B18:01\tHLA-B56:01\tHLA-C07:01\tHLA-C01:02" )
        more +2 ./RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_weak_binders_RelDiag
  #       /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/weak_binders_RelDiag ./RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/weak_binders_annotated_RelDiag
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A02:01\tHLA-A68:01\tHLA-B18:01\tHLA-B56:01\tHLA-C07:01\tHLA-C01:02" )
        more +2 ./RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' -  | (echo -e $strong_binders; cat -) > ./RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_strong_binders_RelDiag
   #      /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/strong_binders_RelDiag ./RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/strong_binders_annotated_RelDiag
       
        #cat /home/fsantaniello/strong_binders_annotated_RelDiag /home/fsantaniello/weak_binders_annotated_RelDiag > ./Summary_Binders_RESULTS_BESU1_BESU2_MHCI_NoMismatchHla__PeptidesPrediction_CommonRelapseDiagnosisOnly.txt
        

import sys

#Common diagnosis relapse

Mutationi=open('./HomeQlogin/BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_weak_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_BESU_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_strong_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_BESU_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()
        

#Relapse Only
Mutationi=open('./HomeQlogin/BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_weak_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_BESU_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_strong_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_BESU_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()

#Diagnosis Only
Mutationi=open('./HomeQlogin/BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_weak_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_BESU_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_strong_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_BESU_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()


sh Launcher_snp2epi_Paper_Annotation.sh CALU1_CALU2_Somatic.Anno.tsv CALU1_CALU2_MHCI_NoMismatchHla HLA-A02:01,HLA-B40:02,HLA-C02:02 CALU_Expression_Kallisto_NeoEpitopes.txt


        #Summary Binders Diagnosis Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A02:01\tHLA-B40:02\tHLA-C02:02" )
        more +2 ./RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0; if ($11>50&&$11<=500) print $3,$2,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_weak_binders_DiagOnly
#         /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/weak_binders_DiagOnly ./RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls ./RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_weak_binders_annotated_DiagOnly
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A02:01\tHLA-B40:02\tHLA-C02:02" )
        more +2 ./RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls |   awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0; if ($11>0&&$11<=50) print $3,$2,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_strong_binders_DiagOnly
 #        /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/strong_binders_DiagOnly ./RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls /home/fsantaniello/strong_binders_annotated_DiagOnly

        #cat /home/fsantaniello/strong_binders_annotated_DiagOnly /home/fsantaniello/weak_binders_annotated_DiagOnly > ./Summary_Binders_RESULTS_CALU1_CALU2_MHCI_NoMismatchHla__PeptidesPrediction_DiagnosisOnly.txt



        #Summary Binders Relapse Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A02:01\tHLA-B40:02\tHLA-C02:02" )
        more +2 ./RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0; if ($11>50&&$11<=500) print $3,$2,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_weak_binders_RelOnly
        # /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/weak_binders_RelOnly ./RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/weak_binders_annotated_RelOnly

        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A02:01\tHLA-B40:02\tHLA-C02:02" )
        more +2 ./RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0; if ($11>0&&$11<=50) print $3,$2,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_strong_binders_RelOnly
         #/home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/strong_binders_RelOnly ./RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/strong_binders_annotated_RelOnly
        
        #cat /home/fsantaniello/strong_binders_annotated_RelOnly /home/fsantaniello/weak_binders_annotated_RelOnly > ./Summary_Binders_RESULTS_CALU1_CALU2_MHCI_NoMismatchHla__PeptidesPrediction_RelapseOnly.txt
        
        
        #Summary Binders Common Diagnosis Relapse
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A02:01\tHLA-B40:02\tHLA-C02:02" )
        more +2 ./RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0; if ($11>50&&$11<=500) print $3,$2,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_weak_binders_RelDiag
  #       /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/weak_binders_RelDiag ./RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/weak_binders_annotated_RelDiag
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A02:01\tHLA-B40:02\tHLA-C02:02" )
        more +2 ./RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0; if ($11>0&&$11<=50) print $3,$2,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_strong_binders_RelDiag
   #      /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/strong_binders_RelDiag ./RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/strong_binders_annotated_RelDiag
       
        #cat /home/fsantaniello/strong_binders_annotated_RelDiag /home/fsantaniello/weak_binders_annotated_RelDiag > ./Summary_Binders_RESULTS_CALU1_CALU2_MHCI_NoMismatchHla__PeptidesPrediction_CommonRelapseDiagnosisOnly.txt



import sys

#Common diagnosis relapse

Mutationi=open('./HomeQlogin/CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_weak_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_CALU_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_strong_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_CALU_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()
        

#Relapse Only
Mutationi=open('./HomeQlogin/CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_weak_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_CALU_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_strong_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_CALU_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()

#Diagnosis Only
Mutationi=open('./HomeQlogin/CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_weak_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_CALU_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_strong_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_CALU_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()


sh Launcher_snp2epi_Paper_Annotation_DEIV.sh DEIV1_DEIV2_DEIV2b_Somatic.Anno.tsv  DEIV2_DEIV2b_MHCI_NoMismatchHla HLA-A11:01,HLA-A24:02,HLA-B14:02,HLA-B35:03,HLA-C08:02,HLA-C12:03 DEIV_Expression_Kallisto_NeoEpitopes.txt 


       #Summary Binders Diagnosis Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A02:01\tHLA-A68:01\tHLA-B18:01\tHLA-B56:01\tHLA-C07:01\tHLA-C01:02" )
        more +2 ./RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_weak_binders_DiagOnly
#         /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/weak_binders_DiagOnly ./RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls ./RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_weak_binders_annotated_DiagOnly
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A02:01\tHLA-A68:01\tHLA-B18:01\tHLA-B56:01\tHLA-C07:01\tHLA-C01:02" )
        more +2 ./RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls |   awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_strong_binders_DiagOnly
 #        /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/strong_binders_DiagOnly ./RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls /home/fsantaniello/strong_binders_annotated_DiagOnly

        #cat /home/fsantaniello/strong_binders_annotated_DiagOnly /home/fsantaniello/weak_binders_annotated_DiagOnly > ./Summary_Binders_RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla__PeptidesPrediction_DiagnosisOnly.txt



        #Summary Binders Relapse Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A02:01\tHLA-A68:01\tHLA-B18:01\tHLA-B56:01\tHLA-C07:01\tHLA-C01:02" )
        more +2 ./RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_weak_binders_RelOnly
        # /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/weak_binders_RelOnly ./RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/weak_binders_annotated_RelOnly

        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A02:01\tHLA-A68:01\tHLA-B18:01\tHLA-B56:01\tHLA-C07:01\tHLA-C01:02" )
        more +2 ./RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_strong_binders_RelOnly
         #/home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/strong_binders_RelOnly ./RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/strong_binders_annotated_RelOnly
        
        #cat /home/fsantaniello/strong_binders_annotated_RelOnly /home/fsantaniello/weak_binders_annotated_RelOnly > ./Summary_Binders_RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla__PeptidesPrediction_RelapseOnly.txt
        
        
        #Summary Binders Common Diagnosis Relapse
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A02:01\tHLA-A68:01\tHLA-B18:01\tHLA-B56:01\tHLA-C07:01\tHLA-C01:02" )
        more +2 ./RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_weak_binders_RelDiag
  #       /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/weak_binders_RelDiag ./RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/weak_binders_annotated_RelDiag
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A02:01\tHLA-A68:01\tHLA-B18:01\tHLA-B56:01\tHLA-C07:01\tHLA-C01:02" )
        more +2 ./RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' -  | (echo -e $strong_binders; cat -) > ./RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_strong_binders_RelDiag
   #      /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/strong_binders_RelDiag ./RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/strong_binders_annotated_RelDiag
       
        #cat /home/fsantaniello/strong_binders_annotated_RelDiag /home/fsantaniello/weak_binders_annotated_RelDiag > ./Summary_Binders_RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla__PeptidesPrediction_CommonRelapseDiagnosisOnly.txt
        

import sys

#Common diagnosis relapse

Mutationi=open('./HomeQlogin/DEIV2_DEIV2b_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_weak_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DEIV_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/DEIV2_DEIV2b_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_strong_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DEIV_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()
        

#Relapse Only
Mutationi=open('./HomeQlogin/DEIV2_DEIV2b_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_weak_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DEIV_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/DEIV2_DEIV2b_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_strong_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DEIV_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()

#Diagnosis Only
Mutationi=open('./HomeQlogin/DEIV2_DEIV2b_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_weak_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DEIV_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/DEIV2_DEIV2b_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_strong_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DEIV2_DEIV2b_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DEIV_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()  

sh Launcher_snp2epi_Paper_Annotation.sh DEST1_DEST2_Somatic.Anno.tsv DEST1_DEST2_MHCI_NoMismatchHla  HLA-A03:01,HLA-B07:02,HLA-C07:02 DEST_Expression_Kallisto_NeoEpitopes.txt


        #Summary Binders Diagnosis Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A03:01\tHLA-B7:02\tHLA-C07:02" )
        more +2 ./RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0; if ($11>50&&$11<=500) print $3,$2,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_weak_binders_DiagOnly
#         /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/weak_binders_DiagOnly ./RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls ./RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_weak_binders_annotated_DiagOnly
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A03:01\tHLA-B7:02\tHLA-C07:02" )
        more +2 ./RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls |   awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0; if ($11>0&&$11<=50) print $3,$2,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_strong_binders_DiagOnly
 #        /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/strong_binders_DiagOnly ./RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls /home/fsantaniello/strong_binders_annotated_DiagOnly

        #cat /home/fsantaniello/strong_binders_annotated_DiagOnly /home/fsantaniello/weak_binders_annotated_DiagOnly > ./Summary_Binders_RESULTS_DEST1_DEST2_MHCI_NoMismatchHla__PeptidesPrediction_DiagnosisOnly.txt



        #Summary Binders Relapse Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A03:01\tHLA-B7:02\tHLA-C07:02" )
        more +2 ./RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0; if ($11>50&&$11<=500) print $3,$2,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_weak_binders_RelOnly
        # /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/weak_binders_RelOnly ./RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/weak_binders_annotated_RelOnly

        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A03:01\tHLA-B7:02\tHLA-C07:02" )
        more +2 ./RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0; if ($11>0&&$11<=50) print $3,$2,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_strong_binders_RelOnly
         #/home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/strong_binders_RelOnly ./RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/strong_binders_annotated_RelOnly
        
        #cat /home/fsantaniello/strong_binders_annotated_RelOnly /home/fsantaniello/weak_binders_annotated_RelOnly > ./Summary_Binders_RESULTS_DEST1_DEST2_MHCI_NoMismatchHla__PeptidesPrediction_RelapseOnly.txt
        
        
        #Summary Binders Common Diagnosis Relapse
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A03:01\tHLA-B7:02\tHLA-C07:02" )
        more +2 ./RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0; if ($11>50&&$11<=500) print $3,$2,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_weak_binders_RelDiag
  #       /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/weak_binders_RelDiag ./RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/weak_binders_annotated_RelDiag
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A03:01\tHLA-B7:02\tHLA-C07:02" )
        more +2 ./RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0; if ($11>0&&$11<=50) print $3,$2,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_strong_binders_RelDiag
   #      /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/strong_binders_RelDiag ./RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/strong_binders_annotated_RelDiag
       
        #cat /home/fsantaniello/strong_binders_annotated_RelDiag /home/fsantaniello/weak_binders_annotated_RelDiag > ./Summary_Binders_RESULTS_DEST1_DEST2_MHCI_NoMismatchHla__PeptidesPrediction_CommonRelapseDiagnosisOnly.txt
        

import sys

#Common diagnosis relapse

Mutationi=open('./HomeQlogin/DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_weak_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DEST_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_strong_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DEST_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()
        

#Relapse Only
Mutationi=open('./HomeQlogin/DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_weak_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DEST_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_strong_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DEST_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()

#Diagnosis Only
Mutationi=open('./HomeQlogin/DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_weak_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DEST_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_strong_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DEST_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()

   



sh Launcher_snp2epi_Paper_Annotation.sh DR1_Somatic_Rerun.tsv DR11_DR12_MHCI_NoMismatchHla HLA-A01:01,HLA-A03:01,HLA-B35:01,HLA-B15:01,HLA-C04:01,HLA-C03:04 DR1_Expression_Kallisto_NeoEpitopes.txt

       #Summary Binders Diagnosis Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A01:01\tHLA-A03:01\tHLA-B35:01\tHLA-B15:01\tHLA-C04:01\tHLA-C03:04" )
        more +2 ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_weak_binders_DiagOnly
#         /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/weak_binders_DiagOnly ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_weak_binders_annotated_DiagOnly
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A01:01\tHLA-A03:01\tHLA-B35:01\tHLA-B15:01\tHLA-C04:01\tHLA-C03:04" )
        more +2 ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls |   awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' -  | (echo -e $strong_binders; cat -) > ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_strong_binders_DiagOnly
 #        /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/strong_binders_DiagOnly ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls /home/fsantaniello/strong_binders_annotated_DiagOnly

        #cat /home/fsantaniello/strong_binders_annotated_DiagOnly /home/fsantaniello/weak_binders_annotated_DiagOnly > ./Summary_Binders_RESULTS_DR11_DR12_MHCI_NoMismatchHla__PeptidesPrediction_DiagnosisOnly.txt



        #Summary Binders Relapse Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A01:01\tHLA-A03:01\tHLA-B35:01\tHLA-B15:01\tHLA-C04:01\tHLA-C03:04" )
        more +2 ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_weak_binders_RelOnly
        # /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/weak_binders_RelOnly ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/weak_binders_annotated_RelOnly

        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A01:01\tHLA-A03:01\tHLA-B35:01\tHLA-B15:01\tHLA-C04:01\tHLA-C03:04" )
        more +2 ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' -  | (echo -e $strong_binders; cat -) > ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_strong_binders_RelOnly
         #/home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/strong_binders_RelOnly ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/strong_binders_annotated_RelOnly
        
        #cat /home/fsantaniello/strong_binders_annotated_RelOnly /home/fsantaniello/weak_binders_annotated_RelOnly > ./Summary_Binders_RESULTS_DR11_DR12_MHCI_NoMismatchHla__PeptidesPrediction_RelapseOnly.txt
        
        
        #Summary Binders Common Diagnosis Relapse
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A01:01\tHLA-A03:01\tHLA-B35:01\tHLA-B15:01\tHLA-C04:01\tHLA-C03:04" )
        more +2 ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_weak_binders_RelDiag
  #       /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/weak_binders_RelDiag ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/weak_binders_annotated_RelDiag
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A01:01\tHLA-A03:01\tHLA-B35:01\tHLA-B15:01\tHLA-C04:01\tHLA-C03:04" )
        more +2 ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' -  | (echo -e $strong_binders; cat -) > ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_strong_binders_RelDiag
   #      /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/strong_binders_RelDiag ./RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/strong_binders_annotated_RelDiag
       
        #cat /home/fsantaniello/strong_binders_annotated_RelDiag /home/fsantaniello/weak_binders_annotated_RelDiag > ./Summary_Binders_RESULTS_DR11_DR12_MHCI_NoMismatchHla__PeptidesPrediction_CommonRelapseDiagnosisOnly.txt


import sys

#Common diagnosis relapse

Mutationi=open('./HomeQlogin/DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR11_DR12_MHCI_NoMismatchHla_weak_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DR1_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR11_DR12_MHCI_NoMismatchHla_strong_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DR1_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()
        

#Relapse Only
Mutationi=open('./HomeQlogin/DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR11_DR12_MHCI_NoMismatchHla_weak_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DR1_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR11_DR12_MHCI_NoMismatchHla_strong_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DR1_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()

#Diagnosis Only
Mutationi=open('./HomeQlogin/DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR11_DR12_MHCI_NoMismatchHla_weak_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DR1_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR11_DR12_MHCI_NoMismatchHla_strong_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DR1_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()

 

sh Launcher_snp2epi_Paper_Annotation.sh DR41_DR42_Somatic.Anno.tsv DR41_DR42_MHCI_NoMismatchHla HLA-A29:02,HLA-A31:01,HLA-B07:02,HLA-C07:02,HLA-C03:04 DR4_Expression_Kallisto_NeoEpitopes.txt

       #Summary Binders Diagnosis Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A29:02\tHLA-A31:01\tHLA-B07:02\tHLA-C07:02\tHLA-C03:04" )
        more +2 ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_weak_binders_DiagOnly
#         /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/weak_binders_DiagOnly ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_weak_binders_annotated_DiagOnly
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A29:02\tHLA-A31:01\tHLA-B07:02\tHLA-C07:02\tHLA-C03:04" )
        more +2 ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls |   awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_strong_binders_DiagOnly
 #        /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/strong_binders_DiagOnly ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls /home/fsantaniello/strong_binders_annotated_DiagOnly

        #cat /home/fsantaniello/strong_binders_annotated_DiagOnly /home/fsantaniello/weak_binders_annotated_DiagOnly > ./Summary_Binders_RESULTS_DR41_DR42_MHCI_NoMismatchHla__PeptidesPrediction_DiagnosisOnly.txt



        #Summary Binders Relapse Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A29:02\tHLA-A31:01\tHLA-B07:02\tHLA-C07:02\tHLA-C03:04" )
        more +2 ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1}' -| (echo -e $weak_binders; cat -) > ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_weak_binders_RelOnly
        # /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/weak_binders_RelOnly ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/weak_binders_annotated_RelOnly

        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A29:02\tHLA-A31:01\tHLA-B07:02\tHLA-C07:02\tHLA-C03:04" )
        more +2 ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_strong_binders_RelOnly
         #/home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/strong_binders_RelOnly ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/strong_binders_annotated_RelOnly
        
        #cat /home/fsantaniello/strong_binders_annotated_RelOnly /home/fsantaniello/weak_binders_annotated_RelOnly > ./Summary_Binders_RESULTS_DR41_DR42_MHCI_NoMismatchHla__PeptidesPrediction_RelapseOnly.txt
        
        
        #Summary Binders Common Diagnosis Relapse
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A29:02\tHLA-A31:01\tHLA-B07:02\tHLA-C07:02\tHLA-C03:04" )
        more +2 ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_weak_binders_RelDiag
  #       /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/weak_binders_RelDiag ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/weak_binders_annotated_RelDiag
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A29:02\tHLA-A31:01\tHLA-B07:02\tHLA-C07:02\tHLA-C03:04" )
        more +2 ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_strong_binders_RelDiag
   #      /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/strong_binders_RelDiag ./RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/strong_binders_annotated_RelDiag
       
        #cat /home/fsantaniello/strong_binders_annotated_RelDiag /home/fsantaniello/weak_binders_annotated_RelDiag > ./Summary_Binders_RESULTS_DR41_DR42_MHCI_NoMismatchHla__PeptidesPrediction_CommonRelapseDiagnosisOnly.txt
        




import sys

#Common diagnosis relapse

Mutationi=open('./HomeQlogin/DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR41_DR42_MHCI_NoMismatchHla_weak_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DR4_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR41_DR42_MHCI_NoMismatchHla_strong_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DR4_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))

out.close()
        

#Relapse Only
Mutationi=open('./HomeQlogin/DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR41_DR42_MHCI_NoMismatchHla_weak_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DR4_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR41_DR42_MHCI_NoMismatchHla_strong_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DR4_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))

out.close()

#Diagnosis Only
Mutationi=open('./HomeQlogin/DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR41_DR42_MHCI_NoMismatchHla_weak_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DR4_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR41_DR42_MHCI_NoMismatchHla_strong_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DR4_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))

out.close()

 


sh Launcher_snp2epi_Paper_Annotation.sh DR51_DR52_Somatic.Anno.tsv DR51_DR52_MHCI_NoMismatchHla HLA-A01:01,HLA-A02:01,HLA-B44:03,HLA-B57:01,HLA-C16:01,HLA-C06:02 DR5_Expression_Kallisto_NeoEpitopes.txt

       #Summary Binders Diagnosis Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A01:01\tHLA-A02:01\tHLA-B44:03\tHLA-B57:01\tHLA-C16:01\tHLA-C06:02" )
        more +2 ./RESULTS_DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_DR51_DR52_MHCI_NoMismatchHla_weak_binders_DiagOnly
#         /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/weak_binders_DiagOnly ./RESULTS_DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls ./RESULTS_DR51_DR52_MHCI_NoMismatchHla_weak_binders_annotated_DiagOnly
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A01:01\tHLA-A02:01\tHLA-B44:03\tHLA-B57:01\tHLA-C16:01\tHLA-C06:02" )
        more +2 ./RESULTS_DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls |   awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' -  | (echo -e $strong_binders; cat -) > ./RESULTS_DR51_DR52_MHCI_NoMismatchHla_strong_binders_DiagOnly
 #        /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/strong_binders_DiagOnly ./RESULTS_DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls /home/fsantaniello/strong_binders_annotated_DiagOnly

        #cat /home/fsantaniello/strong_binders_annotated_DiagOnly /home/fsantaniello/weak_binders_annotated_DiagOnly > ./Summary_Binders_RESULTS_DR51_DR52_MHCI_NoMismatchHla__PeptidesPrediction_DiagnosisOnly.txt



        #Summary Binders Relapse Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A01:01\tHLA-A02:01\tHLA-B44:03\tHLA-B57:01\tHLA-C16:01\tHLA-C06:02" )
        more +2 ./RESULTS_DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_DR51_DR52_MHCI_NoMismatchHla_weak_binders_RelOnly
        # /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/weak_binders_RelOnly ./RESULTS_DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/weak_binders_annotated_RelOnly

        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A01:01\tHLA-A02:01\tHLA-B44:03\tHLA-B57:01\tHLA-C16:01\tHLA-C06:02" )
        more +2 ./RESULTS_DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' -  | (echo -e $strong_binders; cat -) > ./RESULTS_DR51_DR52_MHCI_NoMismatchHla_strong_binders_RelOnly
         #/home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/strong_binders_RelOnly ./RESULTS_DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/strong_binders_annotated_RelOnly
        
        #cat /home/fsantaniello/strong_binders_annotated_RelOnly /home/fsantaniello/weak_binders_annotated_RelOnly > ./Summary_Binders_RESULTS_DR51_DR52_MHCI_NoMismatchHla__PeptidesPrediction_RelapseOnly.txt
        
        
        #Summary Binders Common Diagnosis Relapse
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A01:01\tHLA-A02:01\tHLA-B44:03\tHLA-B57:01\tHLA-C16:01\tHLA-C06:02" )
        more +2 ./RESULTS_DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_DR51_DR52_MHCI_NoMismatchHla_weak_binders_RelDiag
  #       /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/weak_binders_RelDiag ./RESULTS_DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/weak_binders_annotated_RelDiag
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A01:01\tHLA-A02:01\tHLA-B44:03\tHLA-B57:01\tHLA-C16:01\tHLA-C06:02" )
        more +2 ./RESULTS_DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' -  | (echo -e $strong_binders; cat -) > ./RESULTS_DR51_DR52_MHCI_NoMismatchHla_strong_binders_RelDiag
   #      /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/strong_binders_RelDiag ./RESULTS_DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/strong_binders_annotated_RelDiag
       
        #cat /home/fsantaniello/strong_binders_annotated_RelDiag /home/fsantaniello/weak_binders_annotated_RelDiag > ./Summary_Binders_RESULTS_DR51_DR52_MHCI_NoMismatchHla__PeptidesPrediction_CommonRelapseDiagnosisOnly.txt
        

import sys

#Common diagnosis relapse

Mutationi=open('./HomeQlogin/DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR51_DR52_MHCI_NoMismatchHla_weak_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DR5_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR51_DR52_MHCI_NoMismatchHla_strong_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DR5_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()
        

#Relapse Only
Mutationi=open('./HomeQlogin/DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR51_DR52_MHCI_NoMismatchHla_weak_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DR5_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR51_DR52_MHCI_NoMismatchHla_strong_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DR5_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()

#Diagnosis Only
Mutationi=open('./HomeQlogin/DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR51_DR52_MHCI_NoMismatchHla_weak_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DR5_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR51_DR52_MHCI_NoMismatchHla_strong_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DR5_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()

 


sh Launcher_snp2epi_Paper_Annotation.sh FOCA1_FOCA2_Somatic.Anno.tsv FOCA1_FOCA2_MHCI_NoMismatchHla HLA-A02:01,HLA-A01:01,HLA-B51:01,HLA-B08:01,HLA-C01:02,HLA-C07:01 FOCA_Expression_Kallisto_NeoEpitopes.txt

       #Summary Binders Diagnosis Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A02:01\tHLA-A01:01\tHLA-B51:01\tHLA-B08:01\tHLA-C01:02\tHLA-C07:01" )
        more +2 ./RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_weak_binders_DiagOnly
#         /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/weak_binders_DiagOnly ./RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls ./RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_weak_binders_annotated_DiagOnly
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A02:01\tHLA-A01:01\tHLA-B51:01\tHLA-B08:01\tHLA-C01:02\tHLA-C07:01" )
        more +2 ./RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls |   awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' -  | (echo -e $strong_binders; cat -) > ./RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_strong_binders_DiagOnly
 #        /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/strong_binders_DiagOnly ./RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls /home/fsantaniello/strong_binders_annotated_DiagOnly

        #cat /home/fsantaniello/strong_binders_annotated_DiagOnly /home/fsantaniello/weak_binders_annotated_DiagOnly > ./Summary_Binders_RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla__PeptidesPrediction_DiagnosisOnly.txt



        #Summary Binders Relapse Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A02:01\tHLA-A01:01\tHLA-B51:01\tHLA-B08:01\tHLA-C01:02\tHLA-C07:01" )
        more +2 ./RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_weak_binders_RelOnly
        # /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/weak_binders_RelOnly ./RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/weak_binders_annotated_RelOnly

        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A02:01\tHLA-A01:01\tHLA-B51:01\tHLA-B08:01\tHLA-C01:02\tHLA-C07:01" )
        more +2 ./RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' -  | (echo -e $strong_binders; cat -) > ./RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_strong_binders_RelOnly
         #/home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/strong_binders_RelOnly ./RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/strong_binders_annotated_RelOnly
        
        #cat /home/fsantaniello/strong_binders_annotated_RelOnly /home/fsantaniello/weak_binders_annotated_RelOnly > ./Summary_Binders_RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla__PeptidesPrediction_RelapseOnly.txt
        
        
        #Summary Binders Common Diagnosis Relapse
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A02:01\tHLA-A01:01\tHLA-B51:01\tHLA-B08:01\tHLA-C01:02\tHLA-C07:01" )
        more +2 ./RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_weak_binders_RelDiag
  #       /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/weak_binders_RelDiag ./RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/weak_binders_annotated_RelDiag
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A02:01\tHLA-A01:01\tHLA-B51:01\tHLA-B08:01\tHLA-C01:02\tHLA-C07:01" )
        more +2 ./RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_strong_binders_RelDiag
   #      /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/strong_binders_RelDiag ./RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/strong_binders_annotated_RelDiag
       
        #cat /home/fsantaniello/strong_binders_annotated_RelDiag /home/fsantaniello/weak_binders_annotated_RelDiag > ./Summary_Binders_RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla__PeptidesPrediction_CommonRelapseDiagnosisOnly.txt
        

import sys

#Common diagnosis relapse

Mutationi=open('./HomeQlogin/FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_weak_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_FOCA_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_strong_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_FOCA_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()
        

#Relapse Only
Mutationi=open('./HomeQlogin/FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_weak_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_FOCA_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_strong_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_FOCA_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()

#Diagnosis Only
Mutationi=open('./HomeQlogin/FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_weak_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_FOCA_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_strong_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_FOCA_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()

 


sh Launcher_snp2epi_Paper_Annotation.sh GAGRA1_GAGRA2_Somatic.Anno.tsv GAGRA1_GAGRA2_MHCI_NoMismatchHla  HLA-A24:02,HLA-A26:01,HLA-B18:01,HLA-B57:01,HLA-C12:03,HLA-C06:02 GAGRA_Expression_Kallisto_NeoEpitopes.txt


       #Summary Binders Diagnosis Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A24:02\tHLA-A26:01\tHLA-B18:01\tHLA-B57:01\tHLA-C12:03\tHLA-C06:02" )
        more +2 ./RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_weak_binders_DiagOnly
#         /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/weak_binders_DiagOnly ./RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls ./RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_weak_binders_annotated_DiagOnly
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A24:02\tHLA-A26:01\tHLA-B18:01\tHLA-B57:01\tHLA-C12:03\tHLA-C06:02" )
        more +2 ./RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls |   awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_strong_binders_DiagOnly
 #        /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/strong_binders_DiagOnly ./RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls /home/fsantaniello/strong_binders_annotated_DiagOnly

        #cat /home/fsantaniello/strong_binders_annotated_DiagOnly /home/fsantaniello/weak_binders_annotated_DiagOnly > ./Summary_Binders_RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla__PeptidesPrediction_DiagnosisOnly.txt



        #Summary Binders Relapse Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A24:02\tHLA-A26:01\tHLA-B18:01\tHLA-B57:01\tHLA-C12:03\tHLA-C06:02" )
        more +2 ./RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_weak_binders_RelOnly
        # /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/weak_binders_RelOnly ./RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/weak_binders_annotated_RelOnly

        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A24:02\tHLA-A26:01\tHLA-B18:01\tHLA-B57:01\tHLA-C12:03\tHLA-C06:02" )
        more +2 ./RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' -  | (echo -e $strong_binders; cat -) > ./RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_strong_binders_RelOnly
         #/home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/strong_binders_RelOnly ./RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/strong_binders_annotated_RelOnly
        
        #cat /home/fsantaniello/strong_binders_annotated_RelOnly /home/fsantaniello/weak_binders_annotated_RelOnly > ./Summary_Binders_RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla__PeptidesPrediction_RelapseOnly.txt
        
        
        #Summary Binders Common Diagnosis Relapse
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A24:02\tHLA-A26:01\tHLA-B18:01\tHLA-B57:01\tHLA-C12:03\tHLA-C06:02" )
        more +2 ./RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_weak_binders_RelDiag
  #       /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/weak_binders_RelDiag ./RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/weak_binders_annotated_RelDiag
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A24:02\tHLA-A26:01\tHLA-B18:01\tHLA-B57:01\tHLA-C12:03\tHLA-C06:02" )
        more +2 ./RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_strong_binders_RelDiag
   #      /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/strong_binders_RelDiag ./RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/strong_binders_annotated_RelDiag
       
        #cat /home/fsantaniello/strong_binders_annotated_RelDiag /home/fsantaniello/weak_binders_annotated_RelDiag > ./Summary_Binders_RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla__PeptidesPrediction_CommonRelapseDiagnosisOnly.txt


import sys

#Common diagnosis relapse

Mutationi=open('./HomeQlogin/GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_weak_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_GAGRA_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_strong_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_GAGRA_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()
        

#Relapse Only
Mutationi=open('./HomeQlogin/GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_weak_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_GAGRA_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_strong_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_GAGRA_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()

#Diagnosis Only
Mutationi=open('./HomeQlogin/GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_weak_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_GAGRA_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_strong_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_GAGRA_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()

 
      

sh Launcher_snp2epi_Paper_Annotation.sh LUAN1_LUAN2_Somatic.Anno.tsv LUAN1_LUAN2_MHCI_NoMismatchHla HLA-A02:01,HLA-A32:01,HLA-B18:01,HLA-B52:01,HLA-C07:01 LUAN_Expression_Kallisto_NeoEpitopes.txt

       #Summary Binders Diagnosis Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A02:01\tHLA-A32:01\tHLA-B18:01\tHLA-B52:01\tHLA-C07:01" )
        more +2 ./RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_weak_binders_DiagOnly
#         /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/weak_binders_DiagOnly ./RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls ./RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_weak_binders_annotated_DiagOnly
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A02:01\tHLA-A32:01\tHLA-B18:01\tHLA-B52:01\tHLA-C07:01" )
        more +2 ./RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls |   awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_strong_binders_DiagOnly
 #        /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/strong_binders_DiagOnly ./RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls /home/fsantaniello/strong_binders_annotated_DiagOnly

        #cat /home/fsantaniello/strong_binders_annotated_DiagOnly /home/fsantaniello/weak_binders_annotated_DiagOnly > ./Summary_Binders_RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla__PeptidesPrediction_DiagnosisOnly.txt



        #Summary Binders Relapse Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A02:01\tHLA-A32:01\tHLA-B18:01\tHLA-B52:01\tHLA-C07:01" )
        more +2 ./RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1}' -| (echo -e $weak_binders; cat -) > ./RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_weak_binders_RelOnly
        # /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/weak_binders_RelOnly ./RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/weak_binders_annotated_RelOnly

        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A02:01\tHLA-A32:01\tHLA-B18:01\tHLA-B52:01\tHLA-C07:01" )
        more +2 ./RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_strong_binders_RelOnly
         #/home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/strong_binders_RelOnly ./RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/strong_binders_annotated_RelOnly
        
        #cat /home/fsantaniello/strong_binders_annotated_RelOnly /home/fsantaniello/weak_binders_annotated_RelOnly > ./Summary_Binders_RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla__PeptidesPrediction_RelapseOnly.txt
        
        
        #Summary Binders Common Diagnosis Relapse
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A02:01\tHLA-A32:01\tHLA-B18:01\tHLA-B52:01\tHLA-C07:01" )
        more +2 ./RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_weak_binders_RelDiag
  #       /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/weak_binders_RelDiag ./RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/weak_binders_annotated_RelDiag
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A02:01\tHLA-A32:01\tHLA-B18:01\tHLA-B52:01\tHLA-C07:01" )
        more +2 ./RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_strong_binders_RelDiag
   #      /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/strong_binders_RelDiag ./RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/strong_binders_annotated_RelDiag
       
        #cat /home/fsantaniello/strong_binders_annotated_RelDiag /home/fsantaniello/weak_binders_annotated_RelDiag > ./Summary_Binders_RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla__PeptidesPrediction_CommonRelapseDiagnosisOnly.txt
        


import sys

#Common diagnosis relapse

Mutationi=open('./HomeQlogin/LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_weak_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_LUAN_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_strong_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_LUAN_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))

out.close()
        

#Relapse Only
Mutationi=open('./HomeQlogin/LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_weak_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_LUAN_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_strong_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_LUAN_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))

out.close()

#Diagnosis Only
Mutationi=open('./HomeQlogin/LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_weak_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_LUAN_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_strong_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_LUAN_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))

out.close()

 


sh Launcher_snp2epi_Paper_Annotation.sh MABI1_MABI2_Somatic.Anno.tsv MABI1_MABI2_MHCI_NoMismatchHla  HLA-A23:01,HLA-B7:02,HLA-C12:03 MABI_Expression_Kallisto_NeoEpitopes.txt


        #Summary Binders Diagnosis Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A23:01\tHLA-B7:02\tHLA-C12:03" )
        more +2 ./RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0; if ($11>50&&$11<=500) print $3,$2,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_weak_binders_DiagOnly
#         /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/weak_binders_DiagOnly ./RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls ./RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_weak_binders_annotated_DiagOnly
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A23:01\tHLA-B7:02\tHLA-C12:03" )
        more +2 ./RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls |   awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0; if ($11>0&&$11<=50) print $3,$2,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_strong_binders_DiagOnly
 #        /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/strong_binders_DiagOnly ./RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls /home/fsantaniello/strong_binders_annotated_DiagOnly

        #cat /home/fsantaniello/strong_binders_annotated_DiagOnly /home/fsantaniello/weak_binders_annotated_DiagOnly > ./Summary_Binders_RESULTS_MABI1_MABI2_MHCI_NoMismatchHla__PeptidesPrediction_DiagnosisOnly.txt



        #Summary Binders Relapse Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A23:01\tHLA-B7:02\tHLA-C12:03" )
        more +2 ./RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0; if ($11>50&&$11<=500) print $3,$2,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_weak_binders_RelOnly
        # /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/weak_binders_RelOnly ./RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/weak_binders_annotated_RelOnly

        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A23:01\tHLA-B7:02\tHLA-C12:03" )
        more +2 ./RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0; if ($11>0&&$11<=50) print $3,$2,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_strong_binders_RelOnly
         #/home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/strong_binders_RelOnly ./RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/strong_binders_annotated_RelOnly
        
        #cat /home/fsantaniello/strong_binders_annotated_RelOnly /home/fsantaniello/weak_binders_annotated_RelOnly > ./Summary_Binders_RESULTS_MABI1_MABI2_MHCI_NoMismatchHla__PeptidesPrediction_RelapseOnly.txt
        
        
        #Summary Binders Common Diagnosis Relapse
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A23:01\tHLA-B7:02\tHLA-C12:03" )
        more +2 ./RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0; if ($11>50&&$11<=500) print $3,$2,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_weak_binders_RelDiag
  #       /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/weak_binders_RelDiag ./RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/weak_binders_annotated_RelDiag
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A23:01\tHLA-B7:02\tHLA-C12:03" )
        more +2 ./RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0; if ($11>0&&$11<=50) print $3,$2,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_strong_binders_RelDiag
   #      /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/strong_binders_RelDiag ./RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/strong_binders_annotated_RelDiag
       
        #cat /home/fsantaniello/strong_binders_annotated_RelDiag /home/fsantaniello/weak_binders_annotated_RelDiag > ./Summary_Binders_RESULTS_MABI1_MABI2_MHCI_NoMismatchHla__PeptidesPrediction_CommonRelapseDiagnosisOnly.txt
      
      
import sys

#Common diagnosis relapse

Mutationi=open('./HomeQlogin/MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_weak_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_MABI_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_strong_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_MABI_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()
        

#Relapse Only
Mutationi=open('./HomeQlogin/MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_weak_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_MABI_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_strong_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_MABI_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()

#Diagnosis Only
Mutationi=open('./HomeQlogin/MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_weak_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_MABI_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_strong_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_MABI_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()

     
   

sh Launcher_snp2epi_Paper_Annotation.sh MOGE1_MOGE2_Somatic.Anno.tsv MOGE1_MOGE2_MHCI_NoMismatchHla HLA-A24:02,HLA-B13:02,HLA-B44:02,HLA-C05:01,HLA-C06:02 MOGE_Expression_Kallisto_NeoEpitopes.txt

       #Summary Binders Diagnosis Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A24:02\tHLA-B13:02\tHLA-B44:02\tHLA-C05:01\tHLA-C06:02" )
        more +2 ./RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_weak_binders_DiagOnly
#         /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/weak_binders_DiagOnly ./RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls ./RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_weak_binders_annotated_DiagOnly
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A24:02\tHLA-B13:02\tHLA-B44:02\tHLA-C05:01\tHLA-C06:02" )
        more +2 ./RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls |   awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_strong_binders_DiagOnly
 #        /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/strong_binders_DiagOnly ./RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls /home/fsantaniello/strong_binders_annotated_DiagOnly

        #cat /home/fsantaniello/strong_binders_annotated_DiagOnly /home/fsantaniello/weak_binders_annotated_DiagOnly > ./Summary_Binders_RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla__PeptidesPrediction_DiagnosisOnly.txt



        #Summary Binders Relapse Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A24:02\tHLA-B13:02\tHLA-B44:02\tHLA-C05:01\tHLA-C06:02" )
        more +2 ./RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1}' -| (echo -e $weak_binders; cat -) > ./RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_weak_binders_RelOnly
        # /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/weak_binders_RelOnly ./RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/weak_binders_annotated_RelOnly

        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A24:02\tHLA-B13:02\tHLA-B44:02\tHLA-C05:01\tHLA-C06:02" )
        more +2 ./RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_strong_binders_RelOnly
         #/home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/strong_binders_RelOnly ./RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/strong_binders_annotated_RelOnly
        
        #cat /home/fsantaniello/strong_binders_annotated_RelOnly /home/fsantaniello/weak_binders_annotated_RelOnly > ./Summary_Binders_RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla__PeptidesPrediction_RelapseOnly.txt
        
        
        #Summary Binders Common Diagnosis Relapse
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A24:02\tHLA-B13:02\tHLA-B44:02\tHLA-C05:01\tHLA-C06:02" )
        more +2 ./RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_weak_binders_RelDiag
  #       /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/weak_binders_RelDiag ./RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/weak_binders_annotated_RelDiag
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A24:02\tHLA-B13:02\tHLA-B44:02\tHLA-C05:01\tHLA-C06:02" )
        more +2 ./RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_strong_binders_RelDiag
   #      /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/strong_binders_RelDiag ./RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/strong_binders_annotated_RelDiag
       
        #cat /home/fsantaniello/strong_binders_annotated_RelDiag /home/fsantaniello/weak_binders_annotated_RelDiag > ./Summary_Binders_RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla__PeptidesPrediction_CommonRelapseDiagnosisOnly.txt
        


import sys

#Common diagnosis relapse

Mutationi=open('./HomeQlogin/MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_weak_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_MOGE_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_strong_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_MOGE_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))

out.close()
        

#Relapse Only
Mutationi=open('./HomeQlogin/MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_weak_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_MOGE_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_strong_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_MOGE_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))

out.close()

#Diagnosis Only
Mutationi=open('./HomeQlogin/MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_weak_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_MOGE_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_strong_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_MOGE_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))

out.close()

 


sh Launcher_snp2epi_Paper_Annotation.sh PIAG1_PIAG2_Somatic.Anno.tsv PIAG1_PIAG2_MHCI_NoMismatchHla HLA-A03:02,HLA-A24:02,HLA-B44:03,HLA-B35:02,HLA-C04:01 PIAG_Expression_Kallisto_NeoEpitopes.txt


       #Summary Binders Diagnosis Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A03:02\tHLA-A24:02\tHLA-B44:03\tHLA-B35:02\tHLA-C04:01" )
        more +2 ./RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_weak_binders_DiagOnly
#         /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/weak_binders_DiagOnly ./RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls ./RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_weak_binders_annotated_DiagOnly
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A03:02\tHLA-A24:02\tHLA-B44:03\tHLA-B35:02\tHLA-C04:01" )
        more +2 ./RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls |   awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_strong_binders_DiagOnly
 #        /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/strong_binders_DiagOnly ./RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls /home/fsantaniello/strong_binders_annotated_DiagOnly

        #cat /home/fsantaniello/strong_binders_annotated_DiagOnly /home/fsantaniello/weak_binders_annotated_DiagOnly > ./Summary_Binders_RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla__PeptidesPrediction_DiagnosisOnly.txt



        #Summary Binders Relapse Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A03:02\tHLA-A24:02\tHLA-B44:03\tHLA-B35:02\tHLA-C04:01" )
        more +2 ./RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1}' -| (echo -e $weak_binders; cat -) > ./RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_weak_binders_RelOnly
        # /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/weak_binders_RelOnly ./RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/weak_binders_annotated_RelOnly

        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A03:02\tHLA-A24:02\tHLA-B44:03\tHLA-B35:02\tHLA-C04:01" )
        more +2 ./RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_strong_binders_RelOnly
         #/home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/strong_binders_RelOnly ./RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/strong_binders_annotated_RelOnly
        
        #cat /home/fsantaniello/strong_binders_annotated_RelOnly /home/fsantaniello/weak_binders_annotated_RelOnly > ./Summary_Binders_RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla__PeptidesPrediction_RelapseOnly.txt
        
        
        #Summary Binders Common Diagnosis Relapse
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A03:02\tHLA-A24:02\tHLA-B44:03\tHLA-B35:02\tHLA-C04:01" )
        more +2 ./RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_weak_binders_RelDiag
  #       /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/weak_binders_RelDiag ./RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/weak_binders_annotated_RelDiag
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A03:02\tHLA-A24:02\tHLA-B44:03\tHLA-B35:02\tHLA-C04:01" )
        more +2 ./RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_strong_binders_RelDiag
   #      /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/strong_binders_RelDiag ./RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/strong_binders_annotated_RelDiag
       
        #cat /home/fsantaniello/strong_binders_annotated_RelDiag /home/fsantaniello/weak_binders_annotated_RelDiag > ./Summary_Binders_RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla__PeptidesPrediction_CommonRelapseDiagnosisOnly.txt
        


import sys

#Common diagnosis relapse

Mutationi=open('./HomeQlogin/PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_weak_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_PIAG_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_strong_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_PIAG_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))

out.close()
        

#Relapse Only
Mutationi=open('./HomeQlogin/PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_weak_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_PIAG_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_strong_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_PIAG_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))

out.close()

#Diagnosis Only
Mutationi=open('./HomeQlogin/PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_weak_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_PIAG_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_strong_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_PIAG_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))

out.close()

 

sh Launcher_snp2epi_Paper_Annotation.sh PRELU1_PRELU2_Somatic.Anno.tsv PRELU1_PRELU2_MHCI_NoMismatchHla  HLA-A01:01,HLA-B44:02,HLA-C07:04 PRELU_Expression_Kallisto_NeoEpitopes.txt


        #Summary Binders Diagnosis Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A24:02\tHLA-B35:01\tHLA-C04:01" )
        more +2 ./RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0; if ($11>50&&$11<=500) print $3,$2,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_weak_binders_DiagOnly
#         /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/weak_binders_DiagOnly ./RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls ./RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_weak_binders_annotated_DiagOnly
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A24:02\tHLA-B35:01\tHLA-C04:01")
        more +2 ./RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls |   awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0; if ($11>0&&$11<=50) print $3,$2,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_strong_binders_DiagOnly
 #        /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_DiagnosisOnly.txt /home/fsantaniello/strong_binders_DiagOnly ./RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls /home/fsantaniello/strong_binders_annotated_DiagOnly

        #cat /home/fsantaniello/strong_binders_annotated_DiagOnly /home/fsantaniello/weak_binders_annotated_DiagOnly > ./Summary_Binders_RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla__PeptidesPrediction_DiagnosisOnly.txt



        #Summary Binders Relapse Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A24:02\tHLA-B35:01\tHLA-C04:01")
        more +2 ./RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0; if ($11>50&&$11<=500) print $3,$2,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_weak_binders_RelOnly
        # /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/weak_binders_RelOnly ./RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/weak_binders_annotated_RelOnly

        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A24:02\tHLA-B35:01\tHLA-C04:01" )
        more +2 ./RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0; if ($11>0&&$11<=50) print $3,$2,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_strong_binders_RelOnly
         #/home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_RelapseOnly.txt /home/fsantaniello/strong_binders_RelOnly ./RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls /home/fsantaniello/strong_binders_annotated_RelOnly
        
        #cat /home/fsantaniello/strong_binders_annotated_RelOnly /home/fsantaniello/weak_binders_annotated_RelOnly > ./Summary_Binders_RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla__PeptidesPrediction_RelapseOnly.txt
        
        
        #Summary Binders Common Diagnosis Relapse
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A24:02\tHLA-B35:01\tHLA-C04:01")
        more +2 ./RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0; if ($11>50&&$11<=500) print $3,$2,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_weak_binders_RelDiag
  #       /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/weak_binders_RelDiag ./RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/weak_binders_annotated_RelDiag
        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A24:02\tHLA-B35:01\tHLA-C04:01")
        more +2 ./RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0; if ($11>0&&$11<=50) print $3,$2,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_strong_binders_RelDiag
   #      /home/fsantaniello/Summary_Epitopi.py input_dir/name_input_ForEpitopes_CommonDiagnosisRelapse.txt /home/fsantaniello/strong_binders_RelDiag ./RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls /home/fsantaniello/strong_binders_annotated_RelDiag
       
        #cat /home/fsantaniello/strong_binders_annotated_RelDiag /home/fsantaniello/weak_binders_annotated_RelDiag > ./Summary_Binders_RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla__PeptidesPrediction_CommonRelapseDiagnosisOnly.txt
       
      import sys

#Common diagnosis relapse

Mutationi=open('./HomeQlogin/PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_weak_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_PRELU_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_strong_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_PRELU_RelDiag','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()
        

#Relapse Only
Mutationi=open('./HomeQlogin/PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_weak_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_PRELU_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_strong_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_PRELU_RelOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()

#Diagnosis Only
Mutationi=open('./HomeQlogin/PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_weak_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_PRELU_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_strong_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_PRELU_DiagOnly','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip():
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()


#***********************************************************************************************************************************************
#***********************************************************************************************************************************************
#***************************************************************Antigeni Minori classe I********************************************************
#***********************************************************************************************************************************************
#***********************************************************************************************************************************************

#ALFE

sh Launcher_snp2epi_Paper_Annotation_MinorAntigens_ClassI.sh /lustre2/scratch/fsantaniello/VagoRelapse/ALFE_Germline_vs_Donor.tsv.gz MinorAntigens_ALFE1_ALFE2_MHCI_NoMismatchHla HLA-A24:02,HLA-B35:01,HLA-C04:01 ALFE_Expression_Kallisto_NeoEpitopes.txt


        #Summary Binders Diagnosis Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A24:02\tHLA-B35:01\tHLA-C04:01" )
        more +2 ./RESULTS_MinorAntigens_ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls |  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0; if ($11>50&&$11<=500) print $3,$2,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_MinorAntigens_ALFE1_ALFE2_MHCI_NoMismatchHla_weak_binders_MinorAntigens

        strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A24:02\tHLA-B35:01\tHLA-C04:01")
        more +2 ./RESULTS_MinorAntigens_ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls |   awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0; if ($11>0&&$11<=50) print $3,$2,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_MinorAntigens_ALFE1_ALFE2_MHCI_NoMismatchHla_strong_binders_MinorAntigens



import sys

#Common diagnosis relapse

Mutationi=open('./HomeQlogin/MinorAntigens_ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_ALFE1_ALFE2_MHCI_NoMismatchHla_weak_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_ALFE_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))



out.close()

Mutationi=open('./HomeQlogin/MinorAntigens_ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_ALFE1_ALFE2_MHCI_NoMismatchHla_strong_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_ALFE1_ALFE2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_ALFE_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))



out.close()


#ALFE
 cat StrongBinders_Annotati_ALFE_MinorAntigens WeakBinders_Annotati_ALFE_MinorAntigens | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ | grep -v synonymous_variant >  Summary_ALFE_Annotato_MHCI_NoHeaders_NoMismatches_MinorAntigens.txt


ALFE1 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE1_Allfiles_sorted.bam
ALFE2 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE2_Allfiles_sorted.bam

awk '{OFS="\t"; print $8":"$9"-"$9}'   /home/fsantaniello/Summary_ALFE_Annotato_MHCI_NoHeaders_NoMismatches_MinorAntigens.txt | sort | uniq  > Summary_ALFE_Annotato_MHCI_NoHeaders_NoMismatches_MinorAntigens_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE1_Allfiles_sorted.bam >> ALFE_NeoAGs_I_MHCI_Binders_Regions_Diagnosis_NoMismatches_MinorAntigens.txt; done < Summary_ALFE_Annotato_MHCI_NoHeaders_NoMismatches_MinorAntigens_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_ALFE2_Allfiles_sorted.bam >> ALFE_NeoAGs_I_MHCI_Binders_Regions_Relapse_NoMismatches_MinorAntigens.txt; done < Summary_ALFE_Annotato_MHCI_NoHeaders_NoMismatches_MinorAntigens_Regions.txt


library(stringr)
alfe_NeoAG_mhci_RelDiag=read.table('/home/fsantaniello/Summary_ALFE_Annotato_MHCI_NoHeaders_NoMismatches_MinorAntigens.txt')
colnames(alfe_NeoAG_mhci_RelDiag)=c("Gene_ID","Peptide","HLA-A24:02","HLA-B35:01","HLA-C04:01","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A24:02","IC50_HLA-B35:01","IC50_HLA-C04:01")

strong=read.table('StrongBinders_Annotati_ALFE_MinorAntigens',skip=2,head=F)
colnames(strong)=c("Gene_ID","Peptide","HLA-A24:02","HLA-B35:01","HLA-C04:01","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A24:02","IC50_HLA-B35:01","IC50_HLA-C04:01")
strong=strong[grep('^WT_',strong$Gene_ID,invert=T),]
strong=strong[!duplicated(strong),]
strong=merge(strong,ALFE_RNK,by="genes",all.x=T)
relapse=read.table('ALFE_NeoAGs_I_MHCI_Binders_Regions_Relapse_NoMismatches_MinorAntigens.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('ALFE_NeoAGs_I_MHCI_Binders_Regions_Diagnosis_NoMismatches_MinorAntigens.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
strong=merge(strong,relapse,by=c("chr","pos"),all.x=T)
strong=merge(strong,diagnosis,by=c("chr","pos"),all.x=T)
strong$Rel_Mut=str_count(str_to_upper(strong$rel_bases), as.character(strong$mut_nt))
strong$Rel_Wt=str_count(str_to_upper(strong$rel_bases), as.character(strong$wt_nt))
strong$Dx_Mut=str_count(str_to_upper(strong$dx_bases), as.character(strong$mut_nt))
strong$Dx_Wt=str_count(str_to_upper(strong$dx_bases), as.character(strong$wt_nt))
strong=strong[strong$Mutation!='synonymous_variant',]
#112		


#Antigeni Minori strong Espressi
strong_filtered=strong[complete.cases(strong[ , c(17,18,27,29)]),]
strong_filtered=strong_filtered[strong_filtered$Rel_Mut>=3 | strong_filtered$Dx_Mut>=3, ]
#25


#Numero antigeni Minori Allele Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Dx_Mut>=3,])
#22
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3,])
#22



#Numero antigeni Minori Allele Dx_Mut espresso e Rel_Mut no
nrow(strong_filtered[strong_filtered$Dx_Mut>=3 & strong_filtered$Rel_Mut<3,])
#3
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut<3,])
#3
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut>=3,])
#19

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut > strong_filtered$Dx_Mut,])
#7
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut < strong_filtered$Dx_Mut,])
#17

#Numero Antigeni Minori Trascritto Relapse > Trascritto Diagnosi
nrow(strong_filtered[strong_filtered$ALFE1_kallisto > strong_filtered$ALFE2_kallisto,])
#16
#Numero Antigeni Minori Trascritto Relapse < Trascritto Diagnosi
nrow(strong_filtered[strong_filtered$ALFE1_kallisto < strong_filtered$ALFE2_kallisto,])
#9


weak=read.table('WeakBinders_Annotati_ALFE_MinorAntigens',skip=2,head=F)
colnames(weak)=c("Gene_ID","Peptide","HLA-A24:02","HLA-B35:01","HLA-C04:01","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A24:02","IC50_HLA-B35:01","IC50_HLA-C04:01")
weak=weak[grep('^WT_',weak$Gene_ID,invert=T),]
weak=weak[!duplicated(weak),]
weak=merge(weak,ALFE_RNK,by="genes",all.x=T)
relapse=read.table('ALFE_NeoAGs_I_MHCI_Binders_Regions_Relapse_NoMismatches_MinorAntigens.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('ALFE_NeoAGs_I_MHCI_Binders_Regions_Diagnosis_NoMismatches_MinorAntigens.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
weak=merge(weak,relapse,by=c("chr","pos"),all.x=T)
weak=merge(weak,diagnosis,by=c("chr","pos"),all.x=T)
weak$Rel_Mut=str_count(str_to_upper(weak$rel_bases), as.character(weak$mut_nt))
weak$Rel_Wt=str_count(str_to_upper(weak$rel_bases), as.character(weak$wt_nt))
weak$Dx_Mut=str_count(str_to_upper(weak$dx_bases), as.character(weak$mut_nt))
weak$Dx_Wt=str_count(str_to_upper(weak$dx_bases), as.character(weak$wt_nt))
weak=weak[weak$Mutation!='synonymous_variant',]
#305


#Antigeni Minori Espressi
weak_filtered=weak[complete.cases(weak[ , c(17,18,27,29)]),]
weak_filtered=weak_filtered[weak_filtered$Rel_Mut>=3 | weak_filtered$Dx_Mut>=3, ]
#85
#Numero antigeni Minori Allele Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Dx_Mut>=3,])
#77
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3,])
#68

#Numero antigeni Minori Allele Dx_Mut espresso e Rel_Mut no
nrow(weak_filtered[weak_filtered$Dx_Mut>=3 & weak_filtered$Rel_Mut<3,])
#17
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut<3,])
#8
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut>=3,])
#60



#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut > weak_filtered$Dx_Mut,])
#17
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut < weak_filtered$Dx_Mut,])
#66

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(weak_filtered[weak_filtered$ALFE1_kallisto > weak_filtered$ALFE2_kallisto,])
#65
#Numero Antigeni Minori Trascritto Relapse < Trascritto Diagnosi
nrow(weak_filtered[weak_filtered$ALFE1_kallisto < weak_filtered$ALFE2_kallisto,])
#20


Alfe_Binders_NoExprFilt=rbind(strong,weak)
Alfe_Binders_NoExprFilt_2=Alfe_Binders_NoExprFilt[,c(1,2,3,5,9:18,27:30)]
sel <- grepl("IC50_HLA",names(Alfe_Binders_NoExprFilt_2))
Alfe_Binders_NoExprFilt_2[sel] <- lapply(Alfe_Binders_NoExprFilt_2[sel], function(x) replace(x,x >500 , 0) )
Alfe_Binders_Ic50_NoExprFilt=data.table::as.data.table(Alfe_Binders_NoExprFilt_2)
Alfe_Binders_Ic50_NoExprFilt=data.table::melt.data.table(data=Alfe_Binders_Ic50_NoExprFilt,id=c(1:9,13:18))
Alfe_Binders_Ic50_NoExprFilt=Alfe_Binders_Ic50_NoExprFilt[Alfe_Binders_Ic50_NoExprFilt$value>0,]
#Alfe_Binders_Ic50_NoExprFilt=unique(Alfe_Binders_Ic50_NoExprFilt)
Alfe_Binders_Ic50_NoExprFilt_ord=Alfe_Binders_Ic50_NoExprFilt[order(Alfe_Binders_Ic50_NoExprFilt$value),]
Alfe_Binders_Ic50_NoExprFilt_ord=unique(Alfe_Binders_Ic50_NoExprFilt_ord)
Alfe_Binders_Ic50_NoExprFilt_Dx=Alfe_Binders_Ic50_NoExprFilt_ord[Alfe_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3,]
Alfe_Binders_Ic50_NoExprFilt_Rel=Alfe_Binders_Ic50_NoExprFilt_ord[Alfe_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3,]


#BESU

sh Launcher_snp2epi_Paper_Annotation_MinorAntigens_ClassI.sh  /lustre2/scratch/fsantaniello/VagoRelapse/BESU_Germline_vs_Donor.tsv.gz MinorAntigens_BESU1_BESU2_MHCI_NoMismatchHla HLA-A02:01,HLA-A68:01,HLA-B18:01,HLA-B56:01,HLA-C07:01,HLA-C01:02 BESU_Expression_Kallisto_NeoEpitopes.txt
	   #Summary Binders Diagnosis Only
		weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A02:01\tHLA-A68:01\tHLA-B18:01\tHLA-B56:01\tHLA-C07:01\tHLA-C01:02" )
		more +2 ./RESULTS_MinorAntigens_BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_MinorAntigens_BESU1_BESU2_MHCI_NoMismatchHla_weak_binders_MinorAntigens

		strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A02:01\tHLA-A68:01\tHLA-B18:01\tHLA-B56:01\tHLA-C07:01\tHLA-C01:02" )
		more +2 ./RESULTS_MinorAntigens_BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls |   awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_MinorAntigens_BESU1_BESU2_MHCI_NoMismatchHla_strong_binders_MinorAntigens





import sys

#MinorAntigens
Mutationi=open('./HomeQlogin/MinorAntigens_BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_BESU1_BESU2_MHCI_NoMismatchHla_weak_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_BESU_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))



out.close()


Mutationi=open('./Home2/MinorAntigens_BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./Home2/RESULTS_MinorAntigens_BESU1_BESU2_MHCI_NoMismatchHla_strong_binders_MinorAntigens').readlines()
Excel=open('./Home2/RESULTS_MinorAntigens_BESU1_BESU2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./Home2/StrongBinders_Annotati_BESU_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


	for i in binders:
		for j in Excel[2:]:
			if str(i.split('\t')[1])==str(j.split('\t')[1]):
				out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))
	
	out.close()
			
cat StrongBinders_Annotati_BESU_MinorAntigens WeakBinders_Annotati_BESU_MinorAntigens | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ | grep -v synonymous_variant >  Summary_BESU_Annotato_MinorAntigens_NoHeaders_NoMismatches.txt

BESU1 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/BESU1_Allfiles_sorted.bam
BESU2 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/BESU2_Allfiles_sorted.bam


awk '{OFS="\t"; print $11":"$12"-"$12}'   Summary_BESU_Annotato_MinorAntigens_NoHeaders_NoMismatches.txt | sort | uniq  > Summary_BESU_Annotato_MinorAntigens_NoHeaders_NoMismatches_Regions.txt
while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/BESU1_Allfiles_sorted.bam >> BESU_NeoAGs_I_MinorAntigens_NoMismatches_Binders_Regions_Diagnosis.txt; done < Summary_BESU_Annotato_MinorAntigens_NoHeaders_NoMismatches_Regions.txt
while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/BESU2_Allfiles_sorted.bam >> BESU_NeoAGs_I_MinorAntigens_NoMismatches_Binders_Regions_Relapse.txt; done < Summary_BESU_Annotato_MinorAntigens_NoHeaders_NoMismatches_Regions.txt
 


HLA-A02:01\tHLA-A68:01\tHLA-B18:01\tHLA-B56:01\tHLA-C07:01\tHLA-C01:02
besu_NeoAG_mhci_MinorAntigens=read.table('Summary_BESU_Annotato_MinorAntigens_NoHeaders_NoMismatches.txt')
colnames(besu_NeoAG_mhci_MinorAntigens)=c("Gene_ID","Peptide","HLA-A02:01","HLA-A68:01","HLA-B18:01","HLA-B56:01","HLA-C07:01","HLA-C01:02","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A02:01","IC50_HLA-A68:01","IC50_HLA-B18:01","IC50_HLA-B56:01","IC50_HLA-C07:01","IC50_HLA-C01:02")

strong=read.table('StrongBinders_Annotati_BESU_MinorAntigens',skip=2,head=F)
colnames(strong)=c("Gene_ID","Peptide","HLA-A02:01","HLA-A68:01","HLA-B18:01","HLA-B56:01","HLA-C07:01","HLA-C01:02","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A02:01","IC50_HLA-A68:01","IC50_HLA-B18:01","IC50_HLA-B56:01","IC50_HLA-C07:01","IC50_HLA-C01:02")
strong=strong[!duplicated(strong),]
strong=strong[grep('^WT_',strong$Gene_ID,invert=T),]
strong=merge(strong,BESU_RNK,by="genes",all.x=T)
relapse=read.table('BESU_NeoAGs_I_MinorAntigens_NoMismatches_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('BESU_NeoAGs_I_MinorAntigens_NoMismatches_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
strong=merge(strong,relapse,by=c("chr","pos"),all.x=T)
strong=merge(strong,diagnosis,by=c("chr","pos"),all.x=T)
strong$Rel_Mut=str_count(str_to_upper(strong$rel_bases), as.character(strong$mut_nt))
strong$Rel_Wt=str_count(str_to_upper(strong$rel_bases), as.character(strong$wt_nt))
strong$Dx_Mut=str_count(str_to_upper(strong$dx_bases), as.character(strong$mut_nt))
strong$Dx_Wt=str_count(str_to_upper(strong$dx_bases), as.character(strong$wt_nt))
strong=strong[strong$Mutation!='synonymous_variant',]
#511


#Antigeni Minori strong Espressi
strong_filtered=strong[complete.cases(strong[ , c(23,24,33,35)]),]
strong_filtered=strong_filtered[strong_filtered$Rel_Mut>=3 | strong_filtered$Dx_Mut>=3, ]
#140

#sel <- grepl("IC50_HLA",names(strong_filtered))
#strong_filtered[sel] <- lapply(strong_filtered[sel], function(x) replace(x,x >50 , 0) )
#strong_filtered$IC50_Binder=rowSums(strong_filtered[sel])
#besu_strong_filtered=strong_filtered[order(strong_filtered$IC50_Binder),]
#


#Numero antigeni Minori Allele Dx_Mut espresso e Rel_Mut no 
nrow(strong_filtered[strong_filtered$Dx_Mut>=3 & strong_filtered$Rel_Mut<3,])
#7
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut<3,])
#7
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut>=3,])
#126



#Numero antigeni Minori Allele Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Dx_Mut>=3,])
#133
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3,])
#133

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut > strong_filtered$Dx_Mut,])
#89
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut < strong_filtered$Dx_Mut,])
#45

#Numero Antigeni Minori Trascritto Relapse > Trascritto Diagnosi
nrow(strong_filtered[strong_filtered$BESU1_kallisto > strong_filtered$BESU2_kallisto,])
#41
#Numero Antigeni Minori Trascritto Relapse < Trascritto Diagnosi
nrow(strong_filtered[strong_filtered$BESU1_kallisto < strong_filtered$BESU2_kallisto,])
#99


weak=read.table('WeakBinders_Annotati_BESU_MinorAntigens',skip=2,head=F)
colnames(weak)=c("Gene_ID","Peptide","HLA-A02:01","HLA-A68:01","HLA-B18:01","HLA-B56:01","HLA-C07:01","HLA-C01:02","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A02:01","IC50_HLA-A68:01","IC50_HLA-B18:01","IC50_HLA-B56:01","IC50_HLA-C07:01","IC50_HLA-C01:02")
weak=weak[grep('^WT_',weak$Gene_ID,invert=T),]
weak=weak[!duplicated(weak),]
weak=merge(weak,BESU_RNK,by="genes",all.x=T)
relapse=read.table('BESU_NeoAGs_I_MinorAntigens_NoMismatches_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('BESU_NeoAGs_I_MinorAntigens_NoMismatches_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
weak=merge(weak,relapse,by=c("chr","pos"),all.x=T)
weak=merge(weak,diagnosis,by=c("chr","pos"),all.x=T)
weak$Rel_Mut=str_count(str_to_upper(weak$rel_bases), as.character(weak$mut_nt))
weak$Rel_Wt=str_count(str_to_upper(weak$rel_bases), as.character(weak$wt_nt))
weak$Dx_Mut=str_count(str_to_upper(weak$dx_bases), as.character(weak$mut_nt))
weak$Dx_Wt=str_count(str_to_upper(weak$dx_bases), as.character(weak$wt_nt))
weak=weak[weak$Mutation!='synonymous_variant',]
#1416


#Antigeni Minori Espressi
weak_filtered=weak[complete.cases(weak[ , c(23,24,33,35)]),]
weak_filtered=weak_filtered[weak_filtered$Rel_Mut>=3 | weak_filtered$Dx_Mut>=3, ]
#476

#Numero antigeni Minori Allele Dx_Mut espresso e Rel_Mut no 
nrow(weak_filtered[weak_filtered$Dx_Mut>=3 & weak_filtered$Rel_Mut<3,])
#35
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut<3,])
#43
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut>=3,])
#398


#Numero antigeni Minori Allele Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Dx_Mut>=3,])
#433
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3,])
#441

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut > weak_filtered$Dx_Mut,])
#276
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut < weak_filtered$Dx_Mut,])
#178

#Numero Antigeni Minori Trascritto Relapse > Trascritto Diagnosi
nrow(weak_filtered[weak_filtered$BESU1_kallisto > weak_filtered$BESU2_kallisto,])
#152
#Numero Antigeni Minori Trascritto Relapse < Trascritto Diagnosi
nrow(weak_filtered[weak_filtered$BESU1_kallisto < weak_filtered$BESU2_kallisto,])
#324

Besu_Binders_NoExprFilt=rbind(strong,weak)
Besu_Binders_NoExprFilt_2=Besu_Binders_NoExprFilt[,c(1,2,3,5,12,13,14,17:24,33:36)]
sel <- grepl("IC50_HLA",names(Besu_Binders_NoExprFilt_2))
Besu_Binders_NoExprFilt_2[sel] <- lapply(Besu_Binders_NoExprFilt_2[sel], function(x) replace(x,x >500 , 0) )
Besu_Binders_Ic50_NoExprFilt=data.table::as.data.table(Besu_Binders_NoExprFilt_2)
Besu_Binders_Ic50_NoExprFilt=data.table::melt.data.table(data=Besu_Binders_Ic50_NoExprFilt,id=c(1:7,14:19))
Besu_Binders_Ic50_NoExprFilt=Besu_Binders_Ic50_NoExprFilt[Besu_Binders_Ic50_NoExprFilt$value>0,]
Besu_Binders_Ic50_NoExprFilt_ord=Besu_Binders_Ic50_NoExprFilt[order(Besu_Binders_Ic50_NoExprFilt$value),]
Besu_Binders_Ic50_NoExprFilt_ord=unique(Besu_Binders_Ic50_NoExprFilt_ord)
Besu_Binders_Ic50_NoExprFilt_Dx=Besu_Binders_Ic50_NoExprFilt_ord[Besu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3,]
Besu_Binders_Ic50_NoExprFilt_Rel=Besu_Binders_Ic50_NoExprFilt_ord[Besu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3,]



#CALU
sh Launcher_snp2epi_Paper_Annotation_MinorAntigens_ClassI.sh  /lustre2/scratch/fsantaniello/VagoRelapse/CALU_Germline_vs_Donor.tsv.gz MinorAntigens_CALU1_CALU2_MHCI_NoMismatchHla HLA-A02:01,HLA-B40:02,HLA-C02:02 CALU_Expression_Kallisto_NeoEpitopes.txt

		#Summary Binders Diagnosis Only
		weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A02:01\tHLA-B40:02\tHLA-C02:02" )
		more +2 ./RESULTS_MinorAntigens_CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0; if ($11>50&&$11<=500) print $3,$2,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_MinorAntigens_CALU1_CALU2_MHCI_NoMismatchHla_weak_binders_MinorAntigens

		strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A02:01\tHLA-B40:02\tHLA-C02:02" )
		more +2 ./RESULTS_MinorAntigens_CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls |   awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0; if ($11>0&&$11<=50) print $3,$2,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_MinorAntigens_CALU1_CALU2_MHCI_NoMismatchHla_strong_binders_MinorAntigens




Mutationi=open('./HomeQlogin/MinorAntigens_CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_CALU1_CALU2_MHCI_NoMismatchHla_weak_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_CALU_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()



Mutationi=open('./HomeQlogin/MinorAntigens_CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_CALU1_CALU2_MHCI_NoMismatchHla_strong_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_CALU1_CALU2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_CALU_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()



#RelDiag
cat StrongBinders_Annotati_CALU_MinorAntigens WeakBinders_Annotati_CALU_MinorAntigens | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ | grep -v synonymous_variant >  Summary_CALU_Annotato_MinorAntigens_NoHeaders_NoMismatches.txt

#CALU
CALU1 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/CALU1_Allfiles_sorted.bam
CALU2 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/CALU2_Allfiles_sorted.bam



awk '{OFS="\t"; print $8":"$9"-"$9}'   Summary_CALU_Annotato_MinorAntigens_NoHeaders_NoMismatches.txt | sort | uniq  > Summary_CALU_Annotato_MinorAntigens_NoHeaders_NoMismatches_Regions.txt
while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/CALU1_Allfiles_sorted.bam >> CALU_NeoAGs_I_MinorAntigens_Binders_Regions_Diagnosis_NoMismatches.txt; done < Summary_CALU_Annotato_MinorAntigens_NoHeaders_NoMismatches_Regions.txt
while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/CALU2_Allfiles_sorted.bam >> CALU_NeoAGs_I_MinorAntigens_Binders_Regions_Relapse_NoMismatches.txt; done < Summary_CALU_Annotato_MinorAntigens_NoHeaders_NoMismatches_Regions.txt
 

calu_NeoAG_mhci_MinorAntigens=read.table('Summary_CALU_Annotato_MinorAntigens_NoHeaders_NoMismatches.txt')
colnames(calu_NeoAG_mhci_MinorAntigens)=c("Gene_ID","Peptide","HLA-A02:01","HLA-B40:02","HLA-C02:02","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A02:01","IC50_HLA-B40:02","IC50_HLA-C02:02")


strong=read.table('StrongBinders_Annotati_CALU_MinorAntigens',skip=2,head=F)
colnames(strong)=c("Gene_ID","Peptide","HLA-A02:01","HLA-B40:02","HLA-C02:02","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A02:01","IC50_HLA-B40:02","IC50_HLA-C02:02")
strong=strong[grep('^WT_',strong$Gene_ID,invert=T),]
strong=strong[!duplicated(strong),]
strong=merge(strong,CALU_RNK,by="genes",all.x=T)
relapse=read.table('CALU_NeoAGs_I_MinorAntigens_Binders_Regions_Relapse_NoMismatches.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('CALU_NeoAGs_I_MinorAntigens_Binders_Regions_Diagnosis_NoMismatches.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
strong=merge(strong,relapse,by=c("chr","pos"),all.x=T)
strong=merge(strong,diagnosis,by=c("chr","pos"),all.x=T)
strong$Rel_Mut=str_count(str_to_upper(strong$rel_bases), as.character(strong$mut_nt))
strong$Rel_Wt=str_count(str_to_upper(strong$rel_bases), as.character(strong$wt_nt))
strong$Dx_Mut=str_count(str_to_upper(strong$dx_bases), as.character(strong$mut_nt))
strong$Dx_Wt=str_count(str_to_upper(strong$dx_bases), as.character(strong$wt_nt))
strong=strong[strong$Mutation!='synonymous_variant',]
#Strong Binders Totali
#130


#Antigeni Minori strong Espressi
strong_filtered=strong[complete.cases(strong[ , c(17,18,27,29)]),]
strong_filtered=strong_filtered[strong_filtered$Rel_Mut>=3 | strong_filtered$Dx_Mut>=3, ]
nrow(strong_filtered)
#Strong binders espressi
#52
#sel <- grepl("IC50_HLA",names(strong_filtered))
#strong_filtered[sel] <- lapply(strong_filtered[sel], function(x) replace(x,x >50 , 0) )
#strong_filtered$IC50_Binder=rowSums(strong_filtered[sel])
#calu_strong_filtered=strong_filtered[order(strong_filtered$IC50_Binder),]
#


#Numero antigeni Minori Allele Dx_Mut espresso e Rel_Mut no 
nrow(strong_filtered[strong_filtered$Dx_Mut>=3 & strong_filtered$Rel_Mut<3,])
#2
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut<3,])
#2
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut>=3,])
#48




#Numero antigeni Minori Allele Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Dx_Mut>=3,])
#50
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3,])
#50

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut > strong_filtered$Dx_Mut,])
#22
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut < strong_filtered$Dx_Mut,])
#27

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(strong_filtered[strong_filtered$CALU1_kallisto > strong_filtered$CALU2_kallisto,])
#35
#Numero Antigeni Minori Trascritto Diagnosi < Trascritto Relapse
nrow(strong_filtered[strong_filtered$CALU1_kallisto < strong_filtered$CALU2_kallisto,])
#17


weak=read.table('WeakBinders_Annotati_CALU_MinorAntigens',skip=2,head=F)
colnames(weak)=c("Gene_ID","Peptide","HLA-A02:01","HLA-B40:02","HLA-C02:02","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A02:01","IC50_HLA-B40:02","IC50_HLA-C02:02")
weak=weak[grep('^WT_',weak$Gene_ID,invert=T),]
weak=weak[!duplicated(weak),]
weak=merge(weak,CALU_RNK,by="genes",all.x=T)
relapse=read.table('CALU_NeoAGs_I_MinorAntigens_Binders_Regions_Relapse_NoMismatches.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('CALU_NeoAGs_I_MinorAntigens_Binders_Regions_Diagnosis_NoMismatches.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
weak=merge(weak,relapse,by=c("chr","pos"),all.x=T)
weak=merge(weak,diagnosis,by=c("chr","pos"),all.x=T)
weak$Rel_Mut=str_count(str_to_upper(weak$rel_bases), as.character(weak$mut_nt))
weak$Rel_Wt=str_count(str_to_upper(weak$rel_bases), as.character(weak$wt_nt))
weak$Dx_Mut=str_count(str_to_upper(weak$dx_bases), as.character(weak$mut_nt))
weak$Dx_Wt=str_count(str_to_upper(weak$dx_bases), as.character(weak$wt_nt))
weak=weak[weak$Mutation!='synonymous_variant',]
#Weak Binders Totali
nrow(weak)
#[1] 365


#Antigeni Minori weak Espressi
weak_filtered=weak[complete.cases(weak[ , c(17,18,27,29)]),]
weak_filtered=weak_filtered[weak_filtered$Rel_Mut>=3 | weak_filtered$Dx_Mut>=3, ]
nrow(weak_filtered)
#Weak binders espressi
#142

#Numero antigeni Minori Allele Dx_Mut espresso e Rel_Mut no 
nrow(weak_filtered[weak_filtered$Dx_Mut>=3 & weak_filtered$Rel_Mut<3,])
#8
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut<3,])
#15
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut>=3,])
#119

#Numero antigeni Minori Allele Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Dx_Mut>=3,])
#127
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3,])
#134

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut > weak_filtered$Dx_Mut,])
#88
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut < weak_filtered$Dx_Mut,])
#51

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(weak_filtered[weak_filtered$CALU1_kallisto > weak_filtered$CALU2_kallisto,])
#65
#Numero Antigeni Minori Trascritto Diagnosi < Trascritto Relapse
nrow(weak_filtered[weak_filtered$CALU1_kallisto < weak_filtered$CALU2_kallisto,])
#77


Calu_Binders_NoExprFilt=rbind(strong,weak)

Calu_Binders_NoExprFilt_2=Calu_Binders_NoExprFilt[,c(1,2,3,5,9:11,14:18,27:30)]
sel <- grepl("IC50_HLA",names(Calu_Binders_NoExprFilt_2))
Calu_Binders_NoExprFilt_2[sel] <- lapply(Calu_Binders_NoExprFilt_2[sel], function(x) replace(x,x >500 , 0) )
Calu_Binders_Ic50_NoExprFilt=data.table::as.data.table(Calu_Binders_NoExprFilt_2)
Calu_Binders_Ic50_NoExprFilt=data.table::melt.data.table(data=Calu_Binders_Ic50_NoExprFilt,id=c(1:7,11:16))
Calu_Binders_Ic50_NoExprFilt=Calu_Binders_Ic50_NoExprFilt[Calu_Binders_Ic50_NoExprFilt$value>0,]
Calu_Binders_Ic50_NoExprFilt_ord=Calu_Binders_Ic50_NoExprFilt[order(Calu_Binders_Ic50_NoExprFilt$value),]
Calu_Binders_Ic50_NoExprFilt_ord=unique(Calu_Binders_Ic50_NoExprFilt_ord)
Calu_Binders_Ic50_NoExprFilt_ord_homo=Calu_Binders_Ic50_NoExprFilt_ord[Calu_Binders_Ic50_NoExprFilt_ord$variable=='IC50_HLA-A02:01',]
Calu_Binders_Ic50_NoExprFilt_ord_homo=rbind(Calu_Binders_Ic50_NoExprFilt_ord_homo,Calu_Binders_Ic50_NoExprFilt_ord)
Calu_Binders_Ic50_NoExprFilt_ord_homo=Calu_Binders_Ic50_NoExprFilt_ord_homo[order(Calu_Binders_Ic50_NoExprFilt_ord_homo$value),]
Calu_Binders_Ic50_NoExprFilt_Dx=Calu_Binders_Ic50_NoExprFilt_ord_homo[Calu_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3,]
Calu_Binders_Ic50_NoExprFilt_Rel=Calu_Binders_Ic50_NoExprFilt_ord_homo[Calu_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3,]


#DEIV
sh Launcher_snp2epi_Paper_Annotation_MinorAntigens_ClassI.sh  /lustre2/scratch/fsantaniello/VagoRelapse/DEIV_Germline_vs_Donor.tsv.gz MinorAntigens_DEIV1_DEIV2_MHCI_NoMismatchHla  HLA-A11:01,HLA-A24:02,HLA-B14:02,HLA-B35:03,HLA-C08:02,HLA-C12:03 DEIV_Expression_Kallisto_NeoEpitopes.txt


	   #Summary Binders Diagnosis Only
		weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A11:01\tHLA-A24:02\tHLA-B14:02\tHLA-B35:03\tHLA-C08:02\tHLA-C12:03" )
		more +2 ./RESULTS_MinorAntigens_DEIV1_DEIV2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_MinorAntigens_DEIV1_DEIV2_MHCI_NoMismatchHla_weak_binders_MinorAntigens

		strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A11:01\tHLA-A24:02\tHLA-B14:02\tHLA-B35:03\tHLA-C08:02\tHLA-C12:03" )
		more +2 ./RESULTS_MinorAntigens_DEIV1_DEIV2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls |   awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' -  | (echo -e $strong_binders; cat -) > ./RESULTS_MinorAntigens_DEIV1_DEIV2_MHCI_NoMismatchHla_strong_binders_MinorAntigens


Mutationi=open('./HomeQlogin/MinorAntigens_DEIV1_DEIV2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_DEIV1_DEIV2_MHCI_NoMismatchHla_weak_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_DEIV1_DEIV2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DEIV_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()



Mutationi=open('./HomeQlogin/MinorAntigens_DEIV1_DEIV2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_DEIV1_DEIV2_MHCI_NoMismatchHla_strong_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_DEIV1_DEIV2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DEIV_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()) :
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()
			

cat StrongBinders_Annotati_DEIV_MinorAntigens WeakBinders_Annotati_DEIV_MinorAntigens | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ | grep -v synonymous_variant >  Summary_DEIV_Annotato_MinorAntigens_NoHeaders_NoMismatches.txt



awk '{OFS="\t"; print $11":"$12"-"$12}'  Summary_DEIV_Annotato_MinorAntigens_NoHeaders_NoMismatches.txt | sort | uniq  > Summary_DEIV_Annotato_MinorAntigens_NoHeaders_NoMismatches_Regions.txt
 
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2_Allfiles_sorted.bam >> DEIV_NeoAGs_I_MinorAntigens_Binders_Regions_Diagnosis.txt; done < Summary_DEIV_Annotato_MinorAntigens_NoHeaders_NoMismatches_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEIV2b_Allfiles_sorted.bam >> DEIV_NeoAGs_I_MinorAntigens_Binders_Regions_Relapse.txt; done < Summary_DEIV_Annotato_MinorAntigens_NoHeaders_NoMismatches_Regions.txt
 

cat StrongBinders_Annotati_DEIV_MinorAntigens WeakBinders_Annotati_DEIV_MinorAntigens | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ | grep -v synonymous_variant >  Summary_DEIV_Annotato_MinorAntigens_NoHeaders_NoMismatches.txt



strong=read.table('StrongBinders_Annotati_DEIV_MinorAntigens',skip=2,head=F)
colnames(strong)=c("Gene_ID","Peptide","HLA-A11:01","HLA-A24:02","HLA-B14:02","HLA-B35:03","HLA-C08:02","HLA-C12:03","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A11:01","IC50_HLA-A24:02","IC50_HLA-B14:02","IC50_HLA-B35:03","IC50_HLA-C08:02","IC50_HLA-C12:03")
strong=strong[grep('^WT_',strong$Gene_ID,invert=T),]
strong=strong[!duplicated(strong),]
strong=merge(strong,DEIV_RNK,by="genes",all.x=T)
relapse=read.table('DEIV_NeoAGs_I_MinorAntigens_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DEIV_NeoAGs_I_MinorAntigens_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
strong=merge(strong,relapse,by=c("chr","pos"),all.x=T)
strong=merge(strong,diagnosis,by=c("chr","pos"),all.x=T)
strong$Rel_Mut=str_count(str_to_upper(strong$rel_bases), as.character(strong$mut_nt))
strong$Rel_Wt=str_count(str_to_upper(strong$rel_bases), as.character(strong$wt_nt))
strong$Dx_Mut=str_count(str_to_upper(strong$dx_bases), as.character(strong$mut_nt))
strong$Dx_Wt=str_count(str_to_upper(strong$dx_bases), as.character(strong$wt_nt))
strong=strong[strong$Mutation!='synonymous_variant',]

nrow(strong)
#Antigeni minori Strong totali
#335

#Antigeni Minori Espressi
strong_filtered=strong[complete.cases(strong[ , c(23,24,33,35)]),]
strong_filtered=strong_filtered[strong_filtered$Rel_Mut>=3 | strong_filtered$Dx_Mut>=3, ]
#98

#sel <- grepl("IC50_HLA",names(strong_filtered))
#strong_filtered[sel] <- lapply(strong_filtered[sel], function(x) replace(x,x >50 , 0) )
#strong_filtered$IC50_Binder=rowSums(strong_filtered[sel])
#deiv_strong_filtered=strong_filtered[order(strong_filtered$IC50_Binder),]


#Numero antigeni Minori Allele Dx_Mut espresso e Rel_Mut no 
nrow(strong_filtered[strong_filtered$Dx_Mut>=3 & strong_filtered$Rel_Mut<3,])
#8
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut<3,])
#13
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut>=3,])
#77


#Numero antigeni Minori Allele Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Dx_Mut>=3,])
#85
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3,])
#90

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut > strong_filtered$Dx_Mut,])
#65
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut < strong_filtered$Dx_Mut,])
#28

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(strong_filtered[strong_filtered$DEIV2_kallisto > strong_filtered$DEIV2b_kallisto,])
#23
#Numero Antigeni Minori Trascritto Diagnosi < Trascritto Relapse
nrow(strong_filtered[strong_filtered$DEIV2_kallisto < strong_filtered$DEIV2b_kallisto,])
#75




weak=read.table('WeakBinders_Annotati_DEIV_MinorAntigens',skip=2,head=F)
colnames(weak)=c("Gene_ID","Peptide","HLA-A11:01","HLA-A24:02","HLA-B14:02","HLA-B35:03","HLA-C08:02","HLA-C12:03","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A11:01","IC50_HLA-A24:02","IC50_HLA-B14:02","IC50_HLA-B35:03","IC50_HLA-C08:02","IC50_HLA-C12:03")
weak=weak[grep('^WT_',weak$Gene_ID,invert=T),]
weak=weak[!duplicated(weak),]
weak=merge(weak,DEIV_RNK,by="genes",all.x=T)
relapse=read.table('DEIV_NeoAGs_I_MinorAntigens_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DEIV_NeoAGs_I_MinorAntigens_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
weak=merge(weak,relapse,by=c("chr","pos"),all.x=T)
weak=merge(weak,diagnosis,by=c("chr","pos"),all.x=T)
weak$Rel_Mut=str_count(str_to_upper(weak$rel_bases), as.character(weak$mut_nt))
weak$Rel_Wt=str_count(str_to_upper(weak$rel_bases), as.character(weak$wt_nt))
weak$Dx_Mut=str_count(str_to_upper(weak$dx_bases), as.character(weak$mut_nt))
weak$Dx_Wt=str_count(str_to_upper(weak$dx_bases), as.character(weak$wt_nt))
weak=weak[weak$Mutation!='synonymous_variant',]
nrow(weak)
#Antigeni Minori weak totali
#1500

weak_filtered=weak[complete.cases(weak[ , c(23,24,33,35)]),]
weak_filtered=weak_filtered[weak_filtered$Rel_Mut>=3 | weak_filtered$Dx_Mut>=3, ]
#Antigeni minori weak espressi
#445


#Numero antigeni Minori Allele Dx_Mut espresso e Rel_Mut no 
nrow(weak_filtered[weak_filtered$Dx_Mut>=3 & weak_filtered$Rel_Mut<3,])
#35
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut<3,])
#54
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut>=3,])
#356


#Numero antigeni Minori Allele Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Dx_Mut>=3,])
#391
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3,])
#410

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut > weak_filtered$Dx_Mut,])
#289
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut < weak_filtered$Dx_Mut,])
#137

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(weak_filtered[weak_filtered$DEIV2_kallisto > weak_filtered$DEIV2b_kallisto,])
#112
#Numero Antigeni Minori Trascritto Diagnosi < Trascritto Relapse
nrow(weak_filtered[weak_filtered$DEIV2_kallisto < weak_filtered$DEIV2b_kallisto,])
#333


Deiv_Binders_NoExprFilt=rbind(strong,weak)

Deiv_Binders_NoExprFilt_2=Deiv_Binders_NoExprFilt[,c(1,2,3,5,12:14,17:24,33:36)]
sel <- grepl("IC50_HLA",names(Deiv_Binders_NoExprFilt_2))
Deiv_Binders_NoExprFilt_2[sel] <- lapply(Deiv_Binders_NoExprFilt_2[sel], function(x) replace(x,x >500 , 0) )
Deiv_Binders_Ic50_NoExprFilt=data.table::as.data.table(Deiv_Binders_NoExprFilt_2)
Deiv_Binders_Ic50_NoExprFilt=data.table::melt.data.table(data=Deiv_Binders_Ic50_NoExprFilt,id=c(1:7,14:19))
Deiv_Binders_Ic50_NoExprFilt=Deiv_Binders_Ic50_NoExprFilt[Deiv_Binders_Ic50_NoExprFilt$value>0,]
Deiv_Binders_Ic50_NoExprFilt_ord=Deiv_Binders_Ic50_NoExprFilt[order(Deiv_Binders_Ic50_NoExprFilt$value),]
Deiv_Binders_Ic50_NoExprFilt_ord=unique(Deiv_Binders_Ic50_NoExprFilt_ord)
Deiv_Binders_Ic50_NoExprFilt_Dx=Deiv_Binders_Ic50_NoExprFilt_ord[Deiv_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3,]
Deiv_Binders_Ic50_NoExprFilt_Rel=Deiv_Binders_Ic50_NoExprFilt_ord[Deiv_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3,]



#DEST 
sh Launcher_snp2epi_Paper_Annotation_MinorAntigens_ClassI.sh  /lustre2/scratch/fsantaniello/VagoRelapse/DEST_Germline_vs_Donor.tsv.gz MinorAntigens_DEST1_DEST2_MHCI_NoMismatchHla  HLA-A03:01,HLA-B7:02,HLA-C07:02 DEST_Expression_Kallisto_NeoEpitopes.txt
		

		weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A03:01\tHLA-B7:02\tHLA-C07:02" )
		more +2 ./RESULTS_MinorAntigens_DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0; if ($11>50&&$11<=500) print $3,$2,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_MinorAntigens_DEST1_DEST2_MHCI_NoMismatchHla_weak_binders_MinorAntigens

		strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A03:01\tHLA-B7:02\tHLA-C07:02" )
		more +2 ./RESULTS_MinorAntigens_DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0; if ($11>0&&$11<=50) print $3,$2,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_MinorAntigens_DEST1_DEST2_MHCI_NoMismatchHla_strong_binders_MinorAntigens

import sys

Mutationi=open('./HomeQlogin/MinorAntigens_DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_DEST1_DEST2_MHCI_NoMismatchHla_weak_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DEST_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/MinorAntigens_DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_DEST1_DEST2_MHCI_NoMismatchHla_strong_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_DEST1_DEST2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DEST_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()
		

cat StrongBinders_Annotati_DEST_MinorAntigens WeakBinders_Annotati_DEST_MinorAntigens | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ | grep -v synonymous_variant >  Summary_DEST_Annotato_MinorAntigens_NoHeaders_NoMismatches.txt



awk '{OFS="\t"; print $8":"$9"-"$9}'   Summary_DEST_Annotato_MinorAntigens_NoHeaders_NoMismatches.txt | sort | uniq  > Summary_DEST_Annotato_NeoEpitopi_MinorAntigens_MHCI_Regions.txt

DEST1 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST1_Allfiles_sorted.bam
DEST2 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST2_Allfiles_sorted.bam


 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST1_Allfiles_sorted.bam >> DEST_NeoAGs_I_MinorAntigens_Binders_Regions_Diagnosis.txt; done < Summary_DEST_Annotato_NeoEpitopi_MinorAntigens_MHCI_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DEST2_Allfiles_sorted.bam >> DEST_NeoAGs_I_MinorAntigens_Binders_Regions_Relapse.txt; done < Summary_DEST_Annotato_NeoEpitopi_MinorAntigens_MHCI_Regions.txt
 


dest_NeoAG_mhci_MinorAntigens=read.table('Summary_DEST_Annotato_MinorAntigens_NoHeaders_NoMismatches.txt')
colnames(dest_NeoAG_mhci_MinorAntigens)=c("Gene_ID","Peptide","HLA-A03:01","HLA-B07:02","HLA-C07:02","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A03:01","IC50_HLA-B07:02","IC50_HLA-C07:02")

strong=read.table('StrongBinders_Annotati_DEST_MinorAntigens',skip=2,head=F)
colnames(strong)=c("Gene_ID","Peptide","HLA-A03:01","HLA-B07:02","HLA-C07:02","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A03:01","IC50_HLA-B07:02","IC50_HLA-C07:02")
strong=strong[grep('^WT_',strong$Gene_ID,invert=T),]
strong=strong[!duplicated(strong),]
strong=merge(strong,DEST_RNK,by="genes",all.x=T)
relapse=read.table('DEST_NeoAGs_I_MinorAntigens_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DEST_NeoAGs_I_MinorAntigens_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
strong=merge(strong,relapse,by=c("chr","pos"),all.x=T)
strong=merge(strong,diagnosis,by=c("chr","pos"),all.x=T)
strong$Rel_Mut=str_count(str_to_upper(strong$rel_bases), as.character(strong$mut_nt))
strong$Rel_Wt=str_count(str_to_upper(strong$rel_bases), as.character(strong$wt_nt))
strong$Dx_Mut=str_count(str_to_upper(strong$dx_bases), as.character(strong$mut_nt))
strong$Dx_Wt=str_count(str_to_upper(strong$dx_bases), as.character(strong$wt_nt))
strong=strong[strong$Mutation!='synonymous_variant',]
#Strong Binders Totali
#121


#Antigeni Minori strong Espressi
strong_filtered=strong[complete.cases(strong[ , c(17,18,27,29)]),]
strong_filtered=strong_filtered[strong_filtered$Rel_Mut>=3 | strong_filtered$Dx_Mut>=3, ]
nrow(strong_filtered)
#Strong binders espressi
#35

#sel <- grepl("IC50_HLA",names(strong_filtered))
#strong_filtered[sel] <- lapply(strong_filtered[sel], function(x) replace(x,x >50 , 0) )
#strong_filtered$IC50_Binder=rowSums(strong_filtered[sel])
#dest_strong_filtered=strong_filtered[order(strong_filtered$IC50_Binder),]
#


#Numero antigeni Minori Allele Dx_Mut espresso e Rel_Mut no 
nrow(strong_filtered[strong_filtered$Dx_Mut>=3 & strong_filtered$Rel_Mut<3,])
#1
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut<3,])
#4
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut>=3,])
#30



#Numero antigeni Minori Allele Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Dx_Mut>=3,])
#31
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3,])
#34

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut > strong_filtered$Dx_Mut,])
#23
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut < strong_filtered$Dx_Mut,])
#11

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(strong_filtered[strong_filtered$DEST1_kallisto > strong_filtered$DEST2_kallisto,])
#14
#Numero Antigeni Minori Trascritto Diagnosi < Trascritto Relapse
nrow(strong_filtered[strong_filtered$DEST1_kallisto < strong_filtered$DEST2_kallisto,])
#21


weak=read.table('WeakBinders_Annotati_DEST_MinorAntigens',skip=2,head=F)
colnames(weak)=c("Gene_ID","Peptide","HLA-A03:01","HLA-B07:02","HLA-C07:02","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A03:01","IC50_HLA-B07:02","IC50_HLA-C07:02")
weak=weak[grep('^WT_',weak$Gene_ID,invert=T),]
weak=weak[!duplicated(weak),]
weak=merge(weak,DEST_RNK,by="genes",all.x=T)
relapse=read.table('DEST_NeoAGs_I_MinorAntigens_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DEST_NeoAGs_I_MinorAntigens_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
weak=merge(weak,relapse,by=c("chr","pos"),all.x=T)
weak=merge(weak,diagnosis,by=c("chr","pos"),all.x=T)
weak$Rel_Mut=str_count(str_to_upper(weak$rel_bases), as.character(weak$mut_nt))
weak$Rel_Wt=str_count(str_to_upper(weak$rel_bases), as.character(weak$wt_nt))
weak$Dx_Mut=str_count(str_to_upper(weak$dx_bases), as.character(weak$mut_nt))
weak$Dx_Wt=str_count(str_to_upper(weak$dx_bases), as.character(weak$wt_nt))
weak=weak[weak$Mutation!='synonymous_variant',]
nrow(weak)
#Weak Binders Totali
#447



#Antigeni Minori weak Espressi
weak_filtered=weak[complete.cases(weak[ , c(17,18,27,29)]),]
weak_filtered=weak_filtered[weak_filtered$Rel_Mut>=3 | weak_filtered$Dx_Mut>=3, ]
nrow(weak_filtered)
#Weak binders espressi
#145

#Numero antigeni Minori Allele Dx_Mut espresso e Rel_Mut no 
nrow(weak_filtered[weak_filtered$Dx_Mut>=3 & weak_filtered$Rel_Mut<3,])
#10
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut<3,])
#27
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut>=3,])
#108



#Numero antigeni Minori Allele Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Dx_Mut>=3,])
#118
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3,])
#135

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut > weak_filtered$Dx_Mut,])
#89
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut < weak_filtered$Dx_Mut,])
#49

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(weak_filtered[weak_filtered$DEST1_kallisto > weak_filtered$DEST2_kallisto,])
#83
#Numero Antigeni Minori Trascritto Diagnosi < Trascritto Relapse
nrow(weak_filtered[weak_filtered$DEST1_kallisto < weak_filtered$DEST2_kallisto,])
#62

Dest_Binders_NoExprFilt=rbind(strong,weak)

Dest_Binders_NoExprFilt_2=Dest_Binders_NoExprFilt[,c(1,2,3,5,9:11,14:18,27:30)]
sel <- grepl("IC50_HLA",names(Dest_Binders_NoExprFilt_2))
Dest_Binders_NoExprFilt_2[sel] <- lapply(Dest_Binders_NoExprFilt_2[sel], function(x) replace(x,x >500 , 0) )
Dest_Binders_Ic50_NoExprFilt=data.table::as.data.table(Dest_Binders_NoExprFilt_2)
Dest_Binders_Ic50_NoExprFilt=data.table::melt.data.table(data=Dest_Binders_Ic50_NoExprFilt,id=c(1:7,11:16))
Dest_Binders_Ic50_NoExprFilt=Dest_Binders_Ic50_NoExprFilt[Dest_Binders_Ic50_NoExprFilt$value>0,]
Dest_Binders_Ic50_NoExprFilt_ord=Dest_Binders_Ic50_NoExprFilt[order(Dest_Binders_Ic50_NoExprFilt$value),]
Dest_Binders_Ic50_NoExprFilt_ord=unique(Dest_Binders_Ic50_NoExprFilt_ord)
Dest_Binders_Ic50_NoExprFilt_Dx=Dest_Binders_Ic50_NoExprFilt_ord[Dest_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3,]
Dest_Binders_Ic50_NoExprFilt_Rel=Dest_Binders_Ic50_NoExprFilt[Dest_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3,]



#DR1
sh Launcher_snp2epi_Paper_Annotation_MinorAntigens_ClassI.sh  /lustre2/scratch/fsantaniello/VagoRelapse/DR1_Germline_vs_Donor.tsv.gz MinorAntigens_DR11_DR12_MHCI_NoMismatchHla HLA-A01:01,HLA-A03:01,HLA-B35:01,HLA-B15:01,HLA-C04:01,HLA-C03:04 DR1_Expression_Kallisto_NeoEpitopes.txt

	   #Summary Binders Diagnosis Only
		weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A01:01\tHLA-A03:01\tHLA-B35:01\tHLA-B15:01\tHLA-C04:01\tHLA-C03:04" )
		more +2 ./RESULTS_MinorAntigens_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_MinorAntigens_DR11_DR12_MHCI_NoMismatchHla_weak_binders_MinorAntigens

		strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A01:01\tHLA-A03:01\tHLA-B35:01\tHLA-B15:01\tHLA-C04:01\tHLA-C03:04" )
		more +2 ./RESULTS_MinorAntigens_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls |   awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' -  | (echo -e $strong_binders; cat -) > ./RESULTS_MinorAntigens_DR11_DR12_MHCI_NoMismatchHla_strong_binders_MinorAntigens


Mutationi=open('./HomeQlogin/MinorAntigens_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_DR11_DR12_MHCI_NoMismatchHla_weak_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DR1_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()



Mutationi=open('./HomeQlogin/MinorAntigens_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_DR11_DR12_MHCI_NoMismatchHla_strong_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_DR11_DR12_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DR1_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()


cat StrongBinders_Annotati_DR1_MinorAntigens WeakBinders_Annotati_DR1_MinorAntigens | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ | grep -v synonymous_variant >  Summary_DR1_Annotato_MinorAntigens_NoHeaders_NoMismatches.txt

awk '{OFS="\t"; print $11":"$12"-"$12}'   Summary_DR1_Annotato_MinorAntigens_NoHeaders_NoMismatches.txt | sort | uniq  > Summary_DR1_Annotato_NeoEpitopi_MinorAntigens_MHCI_Regions.txt


DR1_1 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DR1_1_Allfiles_sorted.bam
DR1_2 =  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DR1_2_Allfiles_sorted.bam


 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DR1_1_Allfiles_sorted.bam >> DR1_MHCI_MinorAntigens_Binders_Regions_Diagnosis.txt; done < Summary_DR1_Annotato_NeoEpitopi_MinorAntigens_MHCI_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_DR1_2_Allfiles_sorted.bam >> DR1_MHCI_MinorAntigens_Binders_Regions_Relapse.txt; done < Summary_DR1_Annotato_NeoEpitopi_MinorAntigens_MHCI_Regions.txt
 

HLA-A01:01\tHLA-A03:01\tHLA-B35:01\tHLA-B15:01\tHLA-C04:01\tHLA-C03:04

dr1_NeoAG_mhci_MinorAntigens=read.table('Summary_DR1_Annotato_MinorAntigens_NoHeaders_NoMismatches.txt')
colnames(dr1_NeoAG_mhci_MinorAntigens)=c("Gene_ID","Peptide","HLA-A01:01","HLA-A03:01","HLA-B35:01","HLA-B15:01","HLA-C04:01","HLA-C03:04","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A01:01","IC50_HLA-A03:01","IC50_HLA-B35:01","IC50_HLA-B15:01","IC50_HLA-C04:01","IC50_HLA-C03:04")

strong=read.table('StrongBinders_Annotati_DR1_MinorAntigens',skip=2,head=F)
colnames(strong)=c("Gene_ID","Peptide","HLA-A01:01","HLA-A03:01","HLA-B35:01","HLA-B15:01","HLA-C04:01","HLA-C03:04","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A01:01","IC50_HLA-A03:01","IC50_HLA-B35:01","IC50_HLA-B15:01","IC50_HLA-C04:01","IC50_HLA-C03:04")
strong=strong[grep('^WT_',strong$Gene_ID,invert=T),]
strong=strong[!duplicated(strong),]
strong=merge(strong,DR1_RNK,by="genes",all.x=T)
relapse=read.table('DR1_MHCI_MinorAntigens_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DR1_MHCI_MinorAntigens_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
strong=merge(strong,relapse,by=c("chr","pos"),all.x=T)
strong=merge(strong,diagnosis,by=c("chr","pos"),all.x=T)
strong$Rel_Mut=str_count(str_to_upper(strong$rel_bases), as.character(strong$mut_nt))
strong$Rel_Wt=str_count(str_to_upper(strong$rel_bases), as.character(strong$wt_nt))
strong$Dx_Mut=str_count(str_to_upper(strong$dx_bases), as.character(strong$mut_nt))
strong$Dx_Wt=str_count(str_to_upper(strong$dx_bases), as.character(strong$wt_nt))
strong=strong[strong$Mutation!='synonymous_variant',]

nrow(strong)
#Antigeni minori Strong totali
#437

#Antigeni Minori Espressi
strong_filtered=strong[complete.cases(strong[ , c(23,24,33,35)]),]
strong_filtered=strong_filtered[strong_filtered$Rel_Mut>=3 | strong_filtered$Dx_Mut>=3, ]
#122

#sel <- grepl("IC50_HLA",names(strong_filtered))
#strong_filtered[sel] <- lapply(strong_filtered[sel], function(x) replace(x,x >50 , 0) )
#strong_filtered$IC50_Binder=rowSums(strong_filtered[sel])
#dr1_strong_filtered=strong_filtered[order(strong_filtered$IC50_Binder),]
#


#Numero antigeni Minori Allele Dx_Mut espresso e Rel_Mut no 
nrow(strong_filtered[strong_filtered$Dx_Mut>=3 & strong_filtered$Rel_Mut<3,])
#13
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut<3,])
#13
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut>=3,])
#96


#Numero antigeni Minori Allele Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Dx_Mut>=3,])
#109
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3,])
#109

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut > strong_filtered$Dx_Mut,])
#54
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut < strong_filtered$Dx_Mut,])
#60

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(strong_filtered[strong_filtered$DR1_1_kallisto > strong_filtered$DR1_2_kallisto,])
#79
#Numero Antigeni Minori Trascritto Diagnosi < Trascritto Relapse
nrow(strong_filtered[strong_filtered$DR1_1_kallisto < strong_filtered$DR1_2_kallisto,])
#43


weak=read.table('WeakBinders_Annotati_DR1_MinorAntigens',skip=2,head=F)
colnames(weak)=c("Gene_ID","Peptide","HLA-A01:01","HLA-A03:01","HLA-B35:01","HLA-B15:01","HLA-C04:01","HLA-C03:04","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A01:01","IC50_HLA-A03:01","IC50_HLA-B35:01","IC50_HLA-B15:01","IC50_HLA-C04:01","IC50_HLA-C03:04")
weak=weak[grep('^WT_',weak$Gene_ID,invert=T),]
weak=weak[!duplicated(weak),]
weak=merge(weak,DR1_RNK,by="genes",all.x=T)
relapse=read.table('DR1_MHCI_MinorAntigens_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DR1_MHCI_MinorAntigens_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
weak=merge(weak,relapse,by=c("chr","pos"),all.x=T)
weak=merge(weak,diagnosis,by=c("chr","pos"),all.x=T)
weak$Rel_Mut=str_count(str_to_upper(weak$rel_bases), as.character(weak$mut_nt))
weak$Rel_Wt=str_count(str_to_upper(weak$rel_bases), as.character(weak$wt_nt))
weak$Dx_Mut=str_count(str_to_upper(weak$dx_bases), as.character(weak$mut_nt))
weak$Dx_Wt=str_count(str_to_upper(weak$dx_bases), as.character(weak$wt_nt))
weak=weak[weak$Mutation!='synonymous_variant',]
nrow(weak)
#Antigeni minori Weak totali
#1572



#Antigeni Minori Espressi
weak_filtered=weak[complete.cases(weak[ , c(23,24,33,35)]),]
weak_filtered=weak_filtered[weak_filtered$Rel_Mut>=3 | weak_filtered$Dx_Mut>=3, ]
nrow(weak_filtered)
#453




#Numero antigeni Minori Allele Dx_Mut espresso e Rel_Mut no 
nrow(weak_filtered[weak_filtered$Dx_Mut>=3 & weak_filtered$Rel_Mut<3,])
#53
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut<3,])
#63
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut>=3,])
#336




#Numero antigeni Minori Allele Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Dx_Mut>=3,])
#389
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3,])
#400

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut > weak_filtered$Dx_Mut,])
#205
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut < weak_filtered$Dx_Mut,])
#218

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(weak_filtered[weak_filtered$DR1_1_kallisto > weak_filtered$DR1_2_kallisto,])
#287
#Numero Antigeni Minori Trascritto Diagnosi < Trascritto Relapse
nrow(weak_filtered[weak_filtered$DR1_1_kallisto < weak_filtered$DR1_2_kallisto,])
#166

Dr1_Binders_NoExprFilt=rbind(strong,weak)

Dr1_Binders_NoExprFilt_2=Dr1_Binders_NoExprFilt[,c(1,2,3,5,12:14,17:24,33:36)]
sel <- grepl("IC50_HLA",names(Dr1_Binders_NoExprFilt_2))
Dr1_Binders_NoExprFilt_2[sel] <- lapply(Dr1_Binders_NoExprFilt_2[sel], function(x) replace(x,x >500 , 0) )
Dr1_Binders_Ic50_NoExprFilt=data.table::as.data.table(Dr1_Binders_NoExprFilt_2)
Dr1_Binders_Ic50_NoExprFilt=data.table::melt.data.table(data=Dr1_Binders_Ic50_NoExprFilt,id=c(1:7,14:19))
#Dr1_Binders_Ic50_NoExprFilt=unique(Dr1_Binders_Ic50_NoExprFilt)
Dr1_Binders_Ic50_NoExprFilt=Dr1_Binders_Ic50_NoExprFilt[Dr1_Binders_Ic50_NoExprFilt$value>0,]
Dr1_Binders_Ic50_NoExprFilt_ord=Dr1_Binders_Ic50_NoExprFilt[order(Dr1_Binders_Ic50_NoExprFilt$value),]
Dr1_Binders_Ic50_NoExprFilt_ord=unique(Dr1_Binders_Ic50_NoExprFilt_ord)
Dr1_Binders_Ic50_NoExprFilt_Dx=Dr1_Binders_Ic50_NoExprFilt_ord[Dr1_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3,]
Dr1_Binders_Ic50_NoExprFilt_Rel=Dr1_Binders_Ic50_NoExprFilt_ord[Dr1_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3,]


#DR4
sh Launcher_snp2epi_Paper_Annotation_MinorAntigens_ClassI.sh  /lustre2/scratch/fsantaniello/VagoRelapse/DR4_Germline_vs_Donor.tsv.gz MinorAntigens_DR41_DR42_MHCI_NoMismatchHla HLA-A29:02,HLA-A31:01,HLA-B07:02,HLA-C07:02,HLA-C03:04 DR4_Expression_Kallisto_NeoEpitopes.txt

	   #Summary Binders Diagnosis Only
		weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A29:02\tHLA-A31:01\tHLA-B07:02\tHLA-C07:02\tHLA-C03:04" )
		more +2 ./RESULTS_MinorAntigens_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_MinorAntigens_DR41_DR42_MHCI_NoMismatchHla_weak_binders_MinorAntigens

		strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A29:02\tHLA-A31:01\tHLA-B07:02\tHLA-C07:02\tHLA-C03:04" )
		more +2 ./RESULTS_MinorAntigens_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls |   awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_MinorAntigens_DR41_DR42_MHCI_NoMismatchHla_strong_binders_MinorAntigens



import sys

Mutationi=open('./HomeQlogin/MinorAntigens_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_DR41_DR42_MHCI_NoMismatchHla_weak_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DR4_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/MinorAntigens_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_DR41_DR42_MHCI_NoMismatchHla_strong_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_DR41_DR42_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DR4_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))

out.close()
		


#DR4 RelDiag
cat StrongBinders_Annotati_DR4_MinorAntigens WeakBinders_Annotati_DR4_MinorAntigens | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ | grep -v synonymous_variant >  Summary_DR4_Annotato_MinorAntigens_NoHeaders_NoMismatches.txt


DR41 =  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR41_Allfiles_sorted.bam
DR42 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR42_Allfiles_sorted.bam

 
awk '{OFS="\t"; print $10":"$11"-"$11}'   Summary_DR4_Annotato_MinorAntigens_NoHeaders_NoMismatches.txt | sort | uniq  > Summary_DR4_Annotato_NeoEpitopi_MinorAntigens_MHCI_Regions.txt

while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR41_Allfiles_sorted.bam >> DR4_MHC_I_MinorAntigens_Binders_Regions_Diagnosis.txt; done < Summary_DR4_Annotato_NeoEpitopi_MinorAntigens_MHCI_Regions.txt
while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR42_Allfiles_sorted.bam >> DR4_MHC_I_MinorAntigens_Binders_Regions_Relapse.txt; done < Summary_DR4_Annotato_NeoEpitopi_MinorAntigens_MHCI_Regions.txt


dr4_NeoAG_mhci_MinorAntigens=read.table('Summary_DR4_Annotato_MinorAntigens_NoHeaders_NoMismatches.txt')
colnames(dr4_NeoAG_mhci_MinorAntigens)=c("Gene_ID","Peptide","HLA-A29:02","HLA-A31:01","HLA-B07:02","HLA-C07:02","HLA-C03:04","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A29:02","IC50_HLA-A31:01","IC50_HLA-B07:02","IC50_HLA-C07:02","IC50_HLA-C03:04")

strong=read.table('StrongBinders_Annotati_DR4_MinorAntigens',skip=2,head=F)
colnames(strong)=c("Gene_ID","Peptide","HLA-A29:02","HLA-A31:01","HLA-B07:02","HLA-C07:02","HLA-C03:04","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A29:02","IC50_HLA-A31:01","IC50_HLA-B07:02","IC50_HLA-C07:02","IC50_HLA-C03:04")
strong=strong[grep('^WT_',strong$Gene_ID,invert=T),]
strong=strong[!duplicated(strong),]
strong=merge(strong,DR4_RNK,by="genes",all.x=T)
relapse=read.table('DR4_MHC_I_MinorAntigens_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DR4_MHC_I_MinorAntigens_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
strong=merge(strong,relapse,by=c("chr","pos"),all.x=T)
strong=merge(strong,diagnosis,by=c("chr","pos"),all.x=T)
strong$Rel_Mut=str_count(str_to_upper(strong$rel_bases), as.character(strong$mut_nt))
strong$Rel_Wt=str_count(str_to_upper(strong$rel_bases), as.character(strong$wt_nt))
strong$Dx_Mut=str_count(str_to_upper(strong$dx_bases), as.character(strong$mut_nt))
strong$Dx_Wt=str_count(str_to_upper(strong$dx_bases), as.character(strong$wt_nt))
strong=strong[strong$Mutation!='synonymous_variant',]

nrow(strong)
#Antigeni minori Strong totali
#280

#Antigeni Minori Espressi
strong_filtered=strong[complete.cases(strong[ , c(21,22,31,33)]),]
strong_filtered=strong_filtered[strong_filtered$Rel_Mut>=3 | strong_filtered$Dx_Mut>=3, ]
nrow(strong_filtered)
#73

#sel <- grepl("IC50_HLA",names(strong_filtered))
#strong_filtered[sel] <- lapply(strong_filtered[sel], function(x) replace(x,x >50 , 0) )
#strong_filtered$IC50_Binder=rowSums(strong_filtered[sel])
#dr4_strong_filtered=strong_filtered[order(strong_filtered$IC50_Binder),]


#Numero antigeni Minori Allele Dx_Mut espresso e Rel_Mut no 
nrow(strong_filtered[strong_filtered$Dx_Mut>=3 & strong_filtered$Rel_Mut<3,])
#12
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut<3,])
#9
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut>=3,])
#52

#Numero antigeni Minori Allele Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Dx_Mut>=3,])
#64
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3,])
#61

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut > strong_filtered$Dx_Mut,])
#29
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut < strong_filtered$Dx_Mut,])
#43

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(strong_filtered[strong_filtered$DR41_kallisto > strong_filtered$DR42_kallisto,])
#46
#Numero Antigeni Minori Trascritto Diagnosi < Trascritto Relapse
nrow(strong_filtered[strong_filtered$DR41_kallisto < strong_filtered$DR42_kallisto,])
#27


weak=read.table('WeakBinders_Annotati_DR4_MinorAntigens',skip=2,head=F)
colnames(weak)=c("Gene_ID","Peptide","HLA-A29:02","HLA-A31:01","HLA-B07:02","HLA-C07:02","HLA-C03:04","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A29:02","IC50_HLA-A31:01","IC50_HLA-B07:02","IC50_HLA-C07:02","IC50_HLA-C03:04")
weak=weak[grep('^WT_',weak$Gene_ID,invert=T),]
weak=weak[!duplicated(weak),]
weak=merge(weak,DR4_RNK,by="genes",all.x=T)
relapse=read.table('DR4_MHC_I_MinorAntigens_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DR4_MHC_I_MinorAntigens_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
weak=merge(weak,relapse,by=c("chr","pos"),all.x=T)
weak=merge(weak,diagnosis,by=c("chr","pos"),all.x=T)
weak$Rel_Mut=str_count(str_to_upper(weak$rel_bases), as.character(weak$mut_nt))
weak$Rel_Wt=str_count(str_to_upper(weak$rel_bases), as.character(weak$wt_nt))
weak$Dx_Mut=str_count(str_to_upper(weak$dx_bases), as.character(weak$mut_nt))
weak$Dx_Wt=str_count(str_to_upper(weak$dx_bases), as.character(weak$wt_nt))
weak=weak[weak$Mutation!='synonymous_variant',]
nrow(weak)
#Antigeni minori Weak totali
#865

#Antigeni Minori Espressi
weak_filtered=weak[complete.cases(weak[ , c(21,22,31,33)]),]
weak_filtered=weak_filtered[weak_filtered$Rel_Mut>=3 | weak_filtered$Dx_Mut>=3, ]
nrow(weak_filtered)
#245

#Numero antigeni Minori Allele Dx_Mut espresso e Rel_Mut no 
nrow(weak_filtered[weak_filtered$Dx_Mut>=3 & weak_filtered$Rel_Mut<3,])
#40
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut<3,])
#28
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut>=3,])
#177



#Numero antigeni Minori Allele Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Dx_Mut>=3,])
#217
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3,])
#205

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut > weak_filtered$Dx_Mut,])
#110
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut < weak_filtered$Dx_Mut,])
#126

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(weak_filtered[weak_filtered$DR41_kallisto > weak_filtered$DR42_kallisto,])
#157
#Numero Antigeni Minori Trascritto Diagnosi < Trascritto Relapse
nrow(weak_filtered[weak_filtered$DR41_kallisto < weak_filtered$DR42_kallisto,])
#88


Dr4_Binders_NoExprFilt=rbind(strong,weak)

Dr4_Binders_NoExprFilt_2=Dr4_Binders_NoExprFilt[,c(1,2,3,5,11:13,16:22,31:34)]
sel <- grepl("IC50_HLA",names(Dr4_Binders_NoExprFilt_2))
Dr4_Binders_NoExprFilt_2[sel] <- lapply(Dr4_Binders_NoExprFilt_2[sel], function(x) replace(x,x >500 , 0) )
Dr4_Binders_Ic50_NoExprFilt=data.table::as.data.table(Dr4_Binders_NoExprFilt_2)
Dr4_Binders_Ic50_NoExprFilt=data.table::melt.data.table(data=Dr4_Binders_Ic50_NoExprFilt,id=c(1:7,13:18))
Dr4_Binders_Ic50_NoExprFilt=Dr4_Binders_Ic50_NoExprFilt[Dr4_Binders_Ic50_NoExprFilt$value>0,]
Dr4_Binders_Ic50_NoExprFilt_ord=Dr4_Binders_Ic50_NoExprFilt[order(Dr4_Binders_Ic50_NoExprFilt$value),]
Dr4_Binders_Ic50_NoExprFilt_ord=unique(Dr4_Binders_Ic50_NoExprFilt_ord)
Dr4_Binders_Ic50_NoExprFilt_Dx=Dr4_Binders_Ic50_NoExprFilt_ord[Dr4_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3,]
Dr4_Binders_Ic50_NoExprFilt_Rel=Dr4_Binders_Ic50_NoExprFilt_ord[Dr4_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3,]



#DR5
sh Launcher_snp2epi_Paper_Annotation_MinorAntigens_ClassI.sh  /lustre2/scratch/fsantaniello/VagoRelapse/DR5_Germline_vs_Donor.tsv.gz MinorAntigens_DR51_DR52_MHCI_NoMismatchHla HLA-A01:01,HLA-A02:01,HLA-B44:03,HLA-B57:01,HLA-C16:01,HLA-C06:02 DR5_Expression_Kallisto_NeoEpitopes.txt



		weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A01:01\tHLA-A02:01\tHLA-B44:03\tHLA-B57:01\tHLA-C16:01\tHLA-C06:02" )
		more +2 ./RESULTS_MinorAntigens_DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_MinorAntigens_DR51_DR52_MHCI_NoMismatchHla_weak_binders_MinorAntigens

		strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A01:01\tHLA-A02:01\tHLA-B44:03\tHLA-B57:01\tHLA-C16:01\tHLA-C06:02" )
		more +2 ./RESULTS_MinorAntigens_DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' -  | (echo -e $strong_binders; cat -) > ./RESULTS_MinorAntigens_DR51_DR52_MHCI_NoMismatchHla_strong_binders_MinorAntigens
		

import sys

Mutationi=open('./HomeQlogin/MinorAntigens_DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_DR51_DR52_MHCI_NoMismatchHla_weak_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DR5_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/MinorAntigens_DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_DR51_DR52_MHCI_NoMismatchHla_strong_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_DR51_DR52_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DR5_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()

cat StrongBinders_Annotati_DR5_MinorAntigens WeakBinders_Annotati_DR5_MinorAntigens | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ | grep -v synonymous_variant >  Summary_DR5_Annotato_MinorAntigens_MCHI_NoHeaders_NoMismatches.txt



DR51 =  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR51_Allfiles_sorted.bam
DR52 =  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR52_Allfiles_sorted.bam
 
awk '{OFS="\t"; print $11":"$12"-"$12}'   Summary_DR5_Annotato_MinorAntigens_MCHI_NoHeaders_NoMismatches.txt | sort | uniq  > Summary_DR5_Annotato_MinorAntigens_MCHI_NoHeaders_NoMismatches_Regions.txt

while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR51_Allfiles_sorted.bam >> DR5_MHCI_MinorAntigens_Binders_Regions_Diagnosis.txt; done < Summary_DR5_Annotato_MinorAntigens_MCHI_NoHeaders_NoMismatches_Regions.txt
while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/DR52_Allfiles_sorted.bam >> DR5_MHCI_MinorAntigens_Binders_Regions_Relapse.txt; done < Summary_DR5_Annotato_MinorAntigens_MCHI_NoHeaders_NoMismatches_Regions.txt
 

HLA-A01:01\tHLA-A02:01\tHLA-B44:03\tHLA-B57:01\tHLA-C16:01\tHLA-C06:02

dr5_NeoAG_mhci_MinorAntigens=read.table('Summary_DR5_Annotato_MinorAntigens_MCHI_NoHeaders_NoMismatches.txt')
colnames(dr5_NeoAG_mhci_MinorAntigens)=c("Gene_ID","Peptide","HLA-A01:01","HLA-A02:01","HLA-B44:03","HLA-B57:01","HLA-C16:01","HLA-C06:02","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A01:01","IC50_HLA-A02:01","IC50_HLA-B44:03","IC50_HLA-B57:01","IC50_HLA-C16:01","IC50_HLA-C06:02")


strong=read.table('StrongBinders_Annotati_DR5_MinorAntigens',skip=2,head=F)
colnames(strong)=c("Gene_ID","Peptide","HLA-A01:01","HLA-A02:01","HLA-B44:03","HLA-B57:01","HLA-C16:01","HLA-C06:02","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A01:01","IC50_HLA-A02:01","IC50_HLA-B44:03","IC50_HLA-B57:01","IC50_HLA-C16:01","IC50_HLA-C06:02")
strong=strong[grep('^WT_',strong$Gene_ID,invert=T),]
strong=strong[!duplicated(strong),]
strong=merge(strong,DR5_RNK,by="genes",all.x=T)
relapse=read.table('DR5_MHCI_MinorAntigens_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DR1_MHCI_MinorAntigens_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
strong=merge(strong,relapse,by=c("chr","pos"),all.x=T)
strong=merge(strong,diagnosis,by=c("chr","pos"),all.x=T)
strong$Rel_Mut=str_count(str_to_upper(strong$rel_bases), as.character(strong$mut_nt))
strong$Rel_Wt=str_count(str_to_upper(strong$rel_bases), as.character(strong$wt_nt))
strong$Dx_Mut=str_count(str_to_upper(strong$dx_bases), as.character(strong$mut_nt))
strong$Dx_Wt=str_count(str_to_upper(strong$dx_bases), as.character(strong$wt_nt))
strong=strong[strong$Mutation!='synonymous_variant',]

nrow(strong)
#Antigeni minori Strong totali
#185

#Antigeni Minori Espressi
strong_filtered=strong[complete.cases(strong[ , c(23,24,33,35)]),]
strong_filtered=strong_filtered[strong_filtered$Rel_Mut>=3 | strong_filtered$Dx_Mut>=3, ]
#5


#sel <- grepl("IC50_HLA",names(strong_filtered))
#strong_filtered[sel] <- lapply(strong_filtered[sel], function(x) replace(x,x >50 , 0) )
#strong_filtered$IC50_Binder=rowSums(strong_filtered[sel])
#dr5_strong_filtered=strong_filtered[order(strong_filtered$IC50_Binder),]
#


#Numero antigeni Minori Allele Dx_Mut espresso e Rel_Mut no 
nrow(strong_filtered[strong_filtered$Dx_Mut>=3 & strong_filtered$Rel_Mut<3,])
#1
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut<3,])
#0
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut>=3,])
#4



#Numero antigeni Minori Allele Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Dx_Mut>=3,])
#5
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3,])
#4

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut > strong_filtered$Dx_Mut,])
#2
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut < strong_filtered$Dx_Mut,])
#4

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(strong_filtered[strong_filtered$DR51_kallisto > strong_filtered$DR52_kallisto,])
#3
#Numero Antigeni Minori Trascritto Diagnosi < Trascritto Relapse
nrow(strong_filtered[strong_filtered$DR51_kallisto < strong_filtered$DR52_kallisto,])
#2



weak=read.table('WeakBinders_Annotati_DR5_MinorAntigens',skip=2,head=F)
colnames(weak)=c("Gene_ID","Peptide","HLA-A01:01","HLA-A02:01","HLA-B44:03","HLA-B57:01","HLA-C16:01","HLA-C06:02","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A01:01","IC50_HLA-A02:01","IC50_HLA-B44:03","IC50_HLA-B57:01","IC50_HLA-C16:01","IC50_HLA-C06:02")
weak=weak[grep('^WT_',weak$Gene_ID,invert=T),]
weak=weak[!duplicated(weak),]
weak=merge(weak,DR5_RNK,by="genes",all.x=T)
relapse=read.table('DR5_MHCI_MinorAntigens_Binders_Regions_Relapse.txt')
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('DR1_MHCI_MinorAntigens_Binders_Regions_Diagnosis.txt')
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
weak=merge(weak,relapse,by=c("chr","pos"),all.x=T)
weak=merge(weak,diagnosis,by=c("chr","pos"),all.x=T)
weak$Rel_Mut=str_count(str_to_upper(weak$rel_bases), as.character(weak$mut_nt))
weak$Rel_Wt=str_count(str_to_upper(weak$rel_bases), as.character(weak$wt_nt))
weak$Dx_Mut=str_count(str_to_upper(weak$dx_bases), as.character(weak$mut_nt))
weak$Dx_Wt=str_count(str_to_upper(weak$dx_bases), as.character(weak$wt_nt))
weak=weak[weak$Mutation!='synonymous_variant',]
nrow(weak)
#Antigeni minori Weak totali
#628

#Antigeni Minori Espressi
weak_filtered=weak[complete.cases(weak[ , c(23,24,33,35)]),]
weak_filtered=weak_filtered[weak_filtered$Rel_Mut>=3 | weak_filtered$Dx_Mut>=3, ]
#14


#Numero antigeni Minori Allele Dx_Mut espresso e Rel_Mut no 
nrow(weak_filtered[weak_filtered$Dx_Mut>=3 & weak_filtered$Rel_Mut<3,])
#1
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut<3,])
#0
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut>=3,])
#13




#Numero antigeni Minori Allele Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Dx_Mut>=3,])
#14
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3,])
#13

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut > weak_filtered$Dx_Mut,])
#7
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut < weak_filtered$Dx_Mut,])
#6

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(weak_filtered[weak_filtered$DR51_kallisto > weak_filtered$DR52_kallisto,])
#9
#Numero Antigeni Minori Trascritto Diagnosi < Trascritto Relapse
nrow(weak_filtered[weak_filtered$DR51_kallisto < weak_filtered$DR52_kallisto,])
#5


Dr5_Binders_NoExprFilt=rbind(strong,weak)

Dr5_Binders_NoExprFilt_2=Dr5_Binders_NoExprFilt[,c(1,2,3,5,12:14,17:24,33:36)]
sel <- grepl("IC50_HLA",names(Dr5_Binders_NoExprFilt_2))
Dr5_Binders_NoExprFilt_2[sel] <- lapply(Dr5_Binders_NoExprFilt_2[sel], function(x) replace(x,x >500 , 0) )
Dr5_Binders_Ic50_NoExprFilt=data.table::as.data.table(Dr5_Binders_NoExprFilt_2)
Dr5_Binders_Ic50_NoExprFilt=data.table::melt.data.table(data=Dr5_Binders_Ic50_NoExprFilt,id=c(1:7,14:19))
Dr5_Binders_Ic50_NoExprFilt=Dr5_Binders_Ic50_NoExprFilt[Dr5_Binders_Ic50_NoExprFilt$value>0,]
Dr5_Binders_Ic50_NoExprFilt_ord=Dr5_Binders_Ic50_NoExprFilt[order(Dr5_Binders_Ic50_NoExprFilt$value),]
Dr5_Binders_Ic50_NoExprFilt_ord=unique(Dr5_Binders_Ic50_NoExprFilt_ord)
Dr5_Binders_Ic50_NoExprFilt_Dx=Dr5_Binders_Ic50_NoExprFilt_ord[Dr5_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3,]
Dr5_Binders_Ic50_NoExprFilt_Rel=Dr5_Binders_Ic50_NoExprFilt_ord[Dr5_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3,]



#FOCA
sh Launcher_snp2epi_Paper_Annotation_MinorAntigens_ClassI.sh /lustre2/scratch/fsantaniello/VagoRelapse/FOCA_Germline_vs_Donor.tsv.gz MinorAntigens_FOCA1_FOCA2_MHCI_NoMismatchHla HLA-A02:01,HLA-A01:01,HLA-B51:01,HLA-B08:01,HLA-C01:02,HLA-C07:01 FOCA_Expression_Kallisto_NeoEpitopes.txt


		#Summary Binders Common Diagnosis Relapse
		weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A02:01\tHLA-A01:01\tHLA-B51:01\tHLA-B08:01\tHLA-C01:02\tHLA-C07:01" )
		more +2 ./RESULTS_MinorAntigens_FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_MinorAntigens_FOCA1_FOCA2_MHCI_NoMismatchHla_weak_binders_MinorAntigens

		strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A02:01\tHLA-A01:01\tHLA-B51:01\tHLA-B08:01\tHLA-C01:02\tHLA-C07:01" )
		more +2 ./RESULTS_MinorAntigens_FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_MinorAntigens_FOCA1_FOCA2_MHCI_NoMismatchHla_strong_binders_MinorAntigens

		

import sys

#Common diagnosis relapse

Mutationi=open('./HomeQlogin/MinorAntigens_FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_FOCA1_FOCA2_MHCI_NoMismatchHla_weak_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_FOCA_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/MinorAntigens_FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_FOCA1_FOCA2_MHCI_NoMismatchHla_strong_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_FOCA1_FOCA2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_FOCA_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()
		
cat StrongBinders_Annotati_FOCA_MinorAntigens WeakBinders_Annotati_FOCA_MinorAntigens | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ | grep -v synonymous_variant >  Summary_FOCA_Annotato_MinorAntigens_NoHeaders_NoMismatches.txt



FOCA1 =  /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/FOCA1_Allfiles_sorted.bam
FOCA2 =  /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/FOCA2_Allfiles_sorted.bam

awk '{OFS="\t"; print $11":"$12"-"$12}' Summary_FOCA_Annotato_MinorAntigens_NoHeaders_NoMismatches.txt | sort | uniq  > Summary_FOCA_Annotato_MinorAntigens_NoHeaders_NoMismatches_Regions.txt
while read region;  do samtools mpileup -r $region    /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/FOCA1_Allfiles_sorted.bam >> FOCA_MHC_I_MinorAntigens_Binders_Regions_Diagnosis.txt; done < Summary_FOCA_Annotato_MinorAntigens_NoHeaders_NoMismatches_Regions.txt
while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/FOCA2_Allfiles_sorted.bam >> FOCA_MHC_I_MinorAntigens_Binders_Regions_Relapse.txt; done < Summary_FOCA_Annotato_MinorAntigens_NoHeaders_NoMismatches_Regions.txt

HLA-A02:01\tHLA-A01:01\tHLA-B51:01\tHLA-B08:01\tHLA-C01:02\tHLA-C07:01

foca_NeoAG_mhci_MinorAntigens=read.table('Summary_FOCA_Annotato_MinorAntigens_NoHeaders_NoMismatches.txt')
colnames(foca_NeoAG_mhci_MinorAntigens)=c("Gene_ID","Peptide","HLA-A02:01","HLA-A01:01","HLA-B51:01","HLA-B09:01","HLA-C01:02","HLA-C07:01","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A02:01","IC50_HLA-A01:01","IC50_HLA-B51:01","IC50_HLA-B09:01","IC50_HLA-C01:02","IC50_HLA-C07:01")

strong=read.table('StrongBinders_Annotati_FOCA_MinorAntigens',skip=2,head=F)
colnames(strong)=c("Gene_ID","Peptide","HLA-A02:01","HLA-A01:01","HLA-B51:01","HLA-B09:01","HLA-C01:02","HLA-C07:01","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A02:01","IC50_HLA-A01:01","IC50_HLA-B51:01","IC50_HLA-B09:01","IC50_HLA-C01:02","IC50_HLA-C07:01")
strong=strong[grep('^WT_',strong$Gene_ID,invert=T),]
strong=strong[!duplicated(strong),]
strong=merge(strong,FOCA_RNK,by="genes",all.x=T)
relapse=read.table('FOCA_MHC_I_MinorAntigens_Binders_Regions_Relapse.txt',fill=T)
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('FOCA_MHC_I_MinorAntigens_Binders_Regions_Diagnosis.txt',fill=T)
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
strong=merge(strong,relapse,by=c("chr","pos"),all.x=T)
strong=merge(strong,diagnosis,by=c("chr","pos"),all.x=T)
strong$Rel_Mut=str_count(str_to_upper(strong$rel_bases), as.character(strong$mut_nt))
strong$Rel_Wt=str_count(str_to_upper(strong$rel_bases), as.character(strong$wt_nt))
strong$Dx_Mut=str_count(str_to_upper(strong$dx_bases), as.character(strong$mut_nt))
strong$Dx_Wt=str_count(str_to_upper(strong$dx_bases), as.character(strong$wt_nt))
strong=strong[strong$Mutation!='synonymous_variant',]
nrow(strong)
#Antigeni minori Strong totali
#286

#Antigeni Minori Espressi
strong_filtered=strong[complete.cases(strong[ , c(23,24,33,35)]),]
strong_filtered=strong_filtered[strong_filtered$Rel_Mut>=3 | strong_filtered$Dx_Mut>=3, ]
#97


#sel <- grepl("IC50_HLA",names(strong_filtered))
#strong_filtered[sel] <- lapply(strong_filtered[sel], function(x) replace(x,x >50 , 0) )
#strong_filtered$IC50_Binder=rowSums(strong_filtered[sel])
#foca_strong_filtered=strong_filtered[order(strong_filtered$IC50_Binder),]
#


#Numero antigeni Minori Allele Dx_Mut espresso e Rel_Mut no 
nrow(strong_filtered[strong_filtered$Dx_Mut>=3 & strong_filtered$Rel_Mut<3,])
#10
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut<3,])
#6
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut>=3,])
#81



#Numero antigeni Minori Allele Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Dx_Mut>=3,])
#91
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3,])
#87

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut > strong_filtered$Dx_Mut,])
#36
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut < strong_filtered$Dx_Mut,])
#5

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(strong_filtered[strong_filtered$FOCA1_kallisto > strong_filtered$FOCA2_kallisto,])
#53
#Numero Antigeni Minori Trascritto Diagnosi < Trascritto Relapse
nrow(strong_filtered[strong_filtered$FOCA1_kallisto < strong_filtered$FOCA2_kallisto,])
#44



weak=read.table('WeakBinders_Annotati_FOCA_MinorAntigens',skip=2,head=F)
colnames(weak)=c("Gene_ID","Peptide","HLA-A02:01","HLA-A01:01","HLA-B51:01","HLA-B09:01","HLA-C01:02","HLA-C07:01","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A02:01","IC50_HLA-A01:01","IC50_HLA-B51:01","IC50_HLA-B09:01","IC50_HLA-C01:02","IC50_HLA-C07:01")
weak=weak[grep('^WT_',weak$Gene_ID,invert=T),]
weak=weak[!duplicated(weak),]
weak=merge(weak,FOCA_RNK,by="genes",all.x=T)
relapse=read.table('FOCA_MHC_I_MinorAntigens_Binders_Regions_Relapse.txt',fill=T)
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('FOCA_MHC_I_MinorAntigens_Binders_Regions_Diagnosis.txt',fill=T)
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
weak=merge(weak,relapse,by=c("chr","pos"),all.x=T)
weak=merge(weak,diagnosis,by=c("chr","pos"),all.x=T)
weak$Rel_Mut=str_count(str_to_upper(weak$rel_bases), as.character(weak$mut_nt))
weak$Rel_Wt=str_count(str_to_upper(weak$rel_bases), as.character(weak$wt_nt))
weak$Dx_Mut=str_count(str_to_upper(weak$dx_bases), as.character(weak$mut_nt))
weak$Dx_Wt=str_count(str_to_upper(weak$dx_bases), as.character(weak$wt_nt))
weak=weak[weak$Mutation!='synonymous_variant',]

nrow(weak)
#Antigeni minori Weak totali
#845

#Antigeni Minori Espressi
weak_filtered=weak[complete.cases(weak[ , c(23,24,33,35)]),]
weak_filtered=weak_filtered[weak_filtered$Rel_Mut>=3 | weak_filtered$Dx_Mut>=3, ]
#281

#Numero antigeni Minori Allele Dx_Mut espresso e Rel_Mut no 
nrow(weak_filtered[weak_filtered$Dx_Mut>=3 & weak_filtered$Rel_Mut<3,])
#14
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut<3,])
#24
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut>=3,])
#243



#Numero antigeni Minori Allele Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Dx_Mut>=3,])
#257
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3,])
#267

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut > weak_filtered$Dx_Mut,])
#139
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut < weak_filtered$Dx_Mut,])
#131

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(weak_filtered[weak_filtered$FOCA1_kallisto > weak_filtered$FOCA2_kallisto,])
#154
#Numero Antigeni Minori Trascritto Diagnosi < Trascritto Relapse
nrow(weak_filtered[weak_filtered$FOCA1_kallisto < weak_filtered$FOCA2_kallisto,])
#127


Foca_Binders_NoExprFilt=rbind(strong,weak)

Foca_Binders_NoExprFilt_2=Foca_Binders_NoExprFilt[,c(1,2,3,5,12:14,17:24,33:36)]
sel <- grepl("IC50_HLA",names(Foca_Binders_NoExprFilt_2))
Foca_Binders_NoExprFilt_2[sel] <- lapply(Foca_Binders_NoExprFilt_2[sel], function(x) replace(x,x >500 , 0) )
Foca_Binders_Ic50_NoExprFilt=data.table::as.data.table(Foca_Binders_NoExprFilt_2)
Foca_Binders_Ic50_NoExprFilt=data.table::melt.data.table(data=Foca_Binders_Ic50_NoExprFilt,id=c(1:7,14:19))
Foca_Binders_Ic50_NoExprFilt=Foca_Binders_Ic50_NoExprFilt[Foca_Binders_Ic50_NoExprFilt$value>0,]
Foca_Binders_Ic50_NoExprFilt_ord=Foca_Binders_Ic50_NoExprFilt[order(Foca_Binders_Ic50_NoExprFilt$value),]
Foca_Binders_Ic50_NoExprFilt_ord=unique(Foca_Binders_Ic50_NoExprFilt_ord)
Foca_Binders_Ic50_NoExprFilt_Dx=Foca_Binders_Ic50_NoExprFilt_ord[Foca_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3,]
Foca_Binders_Ic50_NoExprFilt_Rel=Foca_Binders_Ic50_NoExprFilt_ord[Foca_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3,]



#GAGRA
sh Launcher_snp2epi_Paper_Annotation_MinorAntigens_ClassI.sh /lustre2/scratch/fsantaniello/VagoRelapse/GAGRA_Germline_vs_Donor.tsv.gz MinorAntigens_GAGRA1_GAGRA2_MHCI_NoMismatchHla  HLA-A24:02,HLA-A26:01,HLA-B18:01,HLA-B57:01,HLA-C12:03,HLA-C06:02 GAGRA_Expression_Kallisto_NeoEpitopes.txt

	   
		weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A24:02\tHLA-A26:01\tHLA-B18:01\tHLA-B57:01\tHLA-C12:03\tHLA-C06:02" )
		more +2 ./RESULTS_MinorAntigens_GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1,0; if ($20>50&&$20<=500) print $3,$2,0,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_MinorAntigens_GAGRA1_GAGRA2_MHCI_NoMismatchHla_weak_binders_MinorAntigens
  
		strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A24:02\tHLA-A26:01\tHLA-B18:01\tHLA-B57:01\tHLA-C12:03\tHLA-C06:02" )
		more +2 ./RESULTS_MinorAntigens_GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1,0; if ($20>0&&$20<=50) print $3,$2,0,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_MinorAntigens_GAGRA1_GAGRA2_MHCI_NoMismatchHla_strong_binders_MinorAntigens
   

import sys

Mutationi=open('./HomeQlogin/MinorAntigens_GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_GAGRA1_GAGRA2_MHCI_NoMismatchHla_weak_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_GAGRA_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/MinorAntigens_GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_GAGRA1_GAGRA2_MHCI_NoMismatchHla_strong_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_GAGRA1_GAGRA2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_GAGRA_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\t'+str(j.split('\t')[19])+'\n'))

out.close()
		
cat StrongBinders_Annotati_GAGRA_MinorAntigens WeakBinders_Annotati_GAGRA_MinorAntigens | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ | grep -v synonymous_variant >  Summary_GAGRA_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches.txt


GAGRA1 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/GAGRA1_Allfiles_sorted.bam
GAGRA2 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/GAGRA2_Allfiles_sorted.bam

awk '{OFS="\t"; print $11":"$12"-"$12}' Summary_GAGRA_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches.txt | sort | uniq  > Summary_GAGRA_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches_Regions.txt
while read region;  do samtools mpileup -r $region    /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/GAGRA1_Allfiles_sorted.bam >> GAGRA_MHC_I_MinorAntigens_Binders_Regions_Diagnosis.txt; done < Summary_GAGRA_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches_Regions.txt
while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/GAGRA2_Allfiles_sorted.bam >> GAGRA_MHC_I_MinorAntigens_Binders_Regions_Relapse.txt; done < Summary_GAGRA_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches_Regions.txt



gagra_NeoAG_mhci_MinorAntigens=read.table('Summary_GAGRA_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches.txt')
colnames(gagra_NeoAG_mhci_MinorAntigens)=c("Gene_ID","Peptide","HLA-A24:02","HLA-A26:01","HLA-B18:01","HLA-B57:01","HLA-C12:03","HLA-C06:02","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A24:02","IC50_HLA-A26:01","IC50_HLA-B18:01","IC50_HLA-B57:01","IC50_HLA-C12:03","IC50_HLA-C06:02")

strong=read.table('StrongBinders_Annotati_GAGRA_MinorAntigens',skip=2,head=F)
colnames(strong)=c("Gene_ID","Peptide","HLA-A24:02","HLA-A26:01","HLA-B18:01","HLA-B57:01","HLA-C12:03","HLA-C06:02","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A24:02","IC50_HLA-A26:01","IC50_HLA-B18:01","IC50_HLA-B57:01","IC50_HLA-C12:03","IC50_HLA-C06:02")
strong=strong[grep('^WT_',strong$Gene_ID,invert=T),]
strong=strong[!duplicated(strong),]
strong=merge(strong,GAGRA_RNK,by="genes",all.x=T)
relapse=read.table('GAGRA_MHC_I_MinorAntigens_Binders_Regions_Relapse.txt',fill=T)
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('GAGRA_MHC_I_MinorAntigens_Binders_Regions_Diagnosis.txt',fill=T)
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
strong=merge(strong,relapse,by=c("chr","pos"),all.x=T)
strong=merge(strong,diagnosis,by=c("chr","pos"),all.x=T)
strong$Rel_Mut=str_count(str_to_upper(strong$rel_bases), as.character(strong$mut_nt))
strong$Rel_Wt=str_count(str_to_upper(strong$rel_bases), as.character(strong$wt_nt))
strong$Dx_Mut=str_count(str_to_upper(strong$dx_bases), as.character(strong$mut_nt))
strong$Dx_Wt=str_count(str_to_upper(strong$dx_bases), as.character(strong$wt_nt))
strong=strong[strong$Mutation!='synonymous_variant',]

nrow(strong)
#Antigeni minori Strong totali
#211

#Antigeni Minori Espressi
strong_filtered=strong[complete.cases(strong[ , c(23,24,33,35)]),]
strong_filtered=strong_filtered[strong_filtered$Rel_Mut>=3 | strong_filtered$Dx_Mut>=3, ]
#75

#sel <- grepl("IC50_HLA",names(strong_filtered))
#strong_filtered[sel] <- lapply(strong_filtered[sel], function(x) replace(x,x >50 , 0) )
#strong_filtered$IC50_Binder=rowSums(strong_filtered[sel])
#gagra_strong_filtered=strong_filtered[order(strong_filtered$IC50_Binder),]


#Numero antigeni Minori Allele Dx_Mut espresso e Rel_Mut no 
nrow(strong_filtered[strong_filtered$Dx_Mut>=3 & strong_filtered$Rel_Mut<3,])
#11
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut<3,])
#2
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut>=3,])
#62


#Numero antigeni Minori Allele Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Dx_Mut>=3,])
#73
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3,])
#64

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut > strong_filtered$Dx_Mut,])
#27
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut < strong_filtered$Dx_Mut,])
#43

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(strong_filtered[strong_filtered$GAGRA1_kallisto > strong_filtered$GAGRA2_kallisto,])
#35
#Numero Antigeni Minori Trascritto Diagnosi < Trascritto Relapse
nrow(strong_filtered[strong_filtered$GAGRA1_kallisto < strong_filtered$GAGRA2_kallisto,])
#40


weak=read.table('WeakBinders_Annotati_GAGRA_MinorAntigens',skip=2,head=F)
colnames(weak)=c("Gene_ID","Peptide","HLA-A24:02","HLA-A26:01","HLA-B18:01","HLA-B57:01","HLA-C12:03","HLA-C06:02","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A24:02","IC50_HLA-A26:01","IC50_HLA-B18:01","IC50_HLA-B57:01","IC50_HLA-C12:03","IC50_HLA-C06:02")
weak=weak[grep('^WT_',weak$Gene_ID,invert=T),]
weak=weak[!duplicated(weak),]
weak=merge(weak,GAGRA_RNK,by="genes",all.x=T)
relapse=read.table('GAGRA_MHC_I_MinorAntigens_Binders_Regions_Relapse.txt',fill=T)
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('GAGRA_MHC_I_MinorAntigens_Binders_Regions_Diagnosis.txt',fill=T)
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
weak=merge(weak,relapse,by=c("chr","pos"),all.x=T)
weak=merge(weak,diagnosis,by=c("chr","pos"),all.x=T)
weak$Rel_Mut=str_count(str_to_upper(weak$rel_bases), as.character(weak$mut_nt))
weak$Rel_Wt=str_count(str_to_upper(weak$rel_bases), as.character(weak$wt_nt))
weak$Dx_Mut=str_count(str_to_upper(weak$dx_bases), as.character(weak$mut_nt))
weak$Dx_Wt=str_count(str_to_upper(weak$dx_bases), as.character(weak$wt_nt))
weak=weak[weak$Mutation!='synonymous_variant',]

nrow(weak)
#Antigeni minori Weak totali
#1080

#Antigeni Minori Espressi
weak_filtered=weak[complete.cases(weak[ , c(23,24,33,35)]),]
weak_filtered=weak_filtered[weak_filtered$Rel_Mut>=3 | weak_filtered$Dx_Mut>=3, ]
#420


#Numero antigeni Minori Allele Dx_Mut espresso e Rel_Mut no 
nrow(weak_filtered[weak_filtered$Dx_Mut>=3 & weak_filtered$Rel_Mut<3,])
#41
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut<3,])
#24
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut>=3,])
#355


#Numero antigeni Minori Allele Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Dx_Mut>=3,])
#396
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3,])
#379

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut > weak_filtered$Dx_Mut,])
#157
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut < weak_filtered$Dx_Mut,])
#246

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(weak_filtered[weak_filtered$GAGRA1_kallisto > weak_filtered$GAGRA2_kallisto,])
#204
#Numero Antigeni Minori Trascritto Diagnosi < Trascritto Relapse
nrow(weak_filtered[weak_filtered$GAGRA1_kallisto < weak_filtered$GAGRA2_kallisto,])
#216



Gagra_Binders_NoExprFilt=rbind(strong,weak)

Gagra_Binders_NoExprFilt_2=Gagra_Binders_NoExprFilt[,c(1,2,3,5,12:14,17:24,33:36)]
sel <- grepl("IC50_HLA",names(Gagra_Binders_NoExprFilt_2))
Gagra_Binders_NoExprFilt_2[sel] <- lapply(Gagra_Binders_NoExprFilt_2[sel], function(x) replace(x,x >500 , 0) )
Gagra_Binders_Ic50_NoExprFilt=data.table::as.data.table(Gagra_Binders_NoExprFilt_2)
Gagra_Binders_Ic50_NoExprFilt=data.table::melt.data.table(data=Gagra_Binders_Ic50_NoExprFilt,id=c(1:7,14:19))
Gagra_Binders_Ic50_NoExprFilt=Gagra_Binders_Ic50_NoExprFilt[Gagra_Binders_Ic50_NoExprFilt$value>0,]
Gagra_Binders_Ic50_NoExprFilt_ord=Gagra_Binders_Ic50_NoExprFilt[order(Gagra_Binders_Ic50_NoExprFilt$value),]
Gagra_Binders_Ic50_NoExprFilt_ord=unique(Gagra_Binders_Ic50_NoExprFilt_ord)
Gagra_Binders_Ic50_NoExprFilt_Dx=Gagra_Binders_Ic50_NoExprFilt_ord[Gagra_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3,]
Gagra_Binders_Ic50_NoExprFilt_Rel=Gagra_Binders_Ic50_NoExprFilt_ord[Gagra_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3,]



#LUAN
sh Launcher_snp2epi_Paper_Annotation_MinorAntigens_ClassI.sh /lustre2/scratch/fsantaniello/VagoRelapse/LUAN_Germline_vs_Donor.tsv.gz MinorAntigens_LUAN1_LUAN2_MHCI_NoMismatchHla HLA-A02:01,HLA-A32:01,HLA-B18:01,HLA-B52:01,HLA-C07:01 LUAN_Expression_Kallisto_NeoEpitopes.txt

		weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A02:01\tHLA-A32:01\tHLA-B18:01\tHLA-B52:01\tHLA-C07:01" )
		more +2 ./RESULTS_MinorAntigens_LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_MinorAntigens_LUAN1_LUAN2_MHCI_NoMismatchHla_weak_binders_MinorAntigens
  
		strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A02:01\tHLA-A32:01\tHLA-B18:01\tHLA-B52:01\tHLA-C07:01" )
		more +2 ./RESULTS_MinorAntigens_LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_MinorAntigens_LUAN1_LUAN2_MHCI_NoMismatchHla_strong_binders_MinorAntigens
		


import sys

Mutationi=open('./HomeQlogin/MinorAntigens_LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_LUAN1_LUAN2_MHCI_NoMismatchHla_weak_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_LUAN_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/MinorAntigens_LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_LUAN1_LUAN2_MHCI_NoMismatchHla_strong_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_LUAN1_LUAN2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_LUAN_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))

out.close()
		

cat StrongBinders_Annotati_LUAN_MinorAntigens WeakBinders_Annotati_LUAN_MinorAntigens | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ | grep -v synonymous_variant >  Summary_LUAN_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches.txt

LUAN1 = /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN1_Allfiles_sorted.bam
LUAN2 = /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN2_Allfiles_sorted.bam

awk '{OFS="\t"; print $10":"$11"-"$11}' Summary_LUAN_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches.txt | sort | uniq  >  Summary_LUAN_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches_Regions.txt
while read region;  do samtools mpileup -r $region    /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN1_Allfiles_sorted.bam >> LUAN_MHC_I_MinorAntigens_Binders_Regions_Diagnosis.txt; done < Summary_LUAN_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches_Regions.txt
while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/170706_SN859_0455_AHMT32BCXY/Project_Vago_161_Relapsing_Leukemia/Sample_LUAN2_Allfiles_sorted.bam >> LUAN_MHC_I_MinorAntigens_Binders_Regions_Relapse.txt; done < Summary_LUAN_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches_Regions.txt


luan_NeoAG_mhci_MinorAntigens=read.table('Summary_LUAN_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches.txt')
colnames(luan_NeoAG_mhci_MinorAntigens)=c("Gene_ID","Peptide","HLA-A02:01","HLA-A32:01","HLA-B18:01","HLA-B52:01","HLA-C07:01","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A02:01","IC50_HLA-A32:01","IC50_HLA-B18:01","IC50_HLA-B52:01","IC50_HLA-C07:01")

strong=read.table('StrongBinders_Annotati_LUAN_MinorAntigens',skip=2,head=F)
colnames(strong)=c("Gene_ID","Peptide","HLA-A02:01","HLA-A32:01","HLA-B18:01","HLA-B52:01","HLA-C07:01","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A02:01","IC50_HLA-A32:01","IC50_HLA-B18:01","IC50_HLA-B52:01","IC50_HLA-C07:01")
strong=strong[grep('^WT_',strong$Gene_ID,invert=T),]
strong=strong[!duplicated(strong),]
strong=merge(strong,LUAN_RNK,by="genes",all.x=T)
relapse=read.table('LUAN_MHC_I_MinorAntigens_Binders_Regions_Relapse.txt',fill=T)
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('LUAN_MHC_I_MinorAntigens_Binders_Regions_Diagnosis.txt',fill=T)
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
strong=merge(strong,relapse,by=c("chr","pos"),all.x=T)
strong=merge(strong,diagnosis,by=c("chr","pos"),all.x=T)
strong$Rel_Mut=str_count(str_to_upper(strong$rel_bases), as.character(strong$mut_nt))
strong$Rel_Wt=str_count(str_to_upper(strong$rel_bases), as.character(strong$wt_nt))
strong$Dx_Mut=str_count(str_to_upper(strong$dx_bases), as.character(strong$mut_nt))
strong$Dx_Wt=str_count(str_to_upper(strong$dx_bases), as.character(strong$wt_nt))
strong=strong[strong$Mutation!='synonymous_variant',]
nrow(strong)
#Antigeni minori Strong totali
#174

#Antigeni Minori Espressi
strong_filtered=strong[complete.cases(strong[ , c(21,22,31,33)]),]
strong_filtered=strong_filtered[strong_filtered$Rel_Mut>=3 | strong_filtered$Dx_Mut>=3, ]
#53

#sel <- grepl("IC50_HLA",names(strong_filtered))
#strong_filtered[sel] <- lapply(strong_filtered[sel], function(x) replace(x,x >50 , 0) )
#strong_filtered$IC50_Binder=rowSums(strong_filtered[sel])
#luan_strong_filtered=strong_filtered[order(strong_filtered$IC50_Binder),]
#


#Numero antigeni Minori Allele Dx_Mut espresso e Rel_Mut no 
nrow(strong_filtered[strong_filtered$Dx_Mut>=3 & strong_filtered$Rel_Mut<3,])
#5
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut<3,])
#8
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut>=3,])
#40



#Numero antigeni Minori Allele Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Dx_Mut>=3,])
#45
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3,])
#48

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut > strong_filtered$Dx_Mut,])
#27
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut < strong_filtered$Dx_Mut,])
#26

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(strong_filtered[strong_filtered$LUAN1_kallisto > strong_filtered$LUAN2_kallisto,])
#31
#Numero Antigeni Minori Trascritto Diagnosi < Trascritto Relapse
nrow(strong_filtered[strong_filtered$LUAN1_kallisto < strong_filtered$LUAN2_kallisto,])
#22


weak=read.table('WeakBinders_Annotati_LUAN_MinorAntigens',skip=2,head=F)
colnames(weak)=c("Gene_ID","Peptide","HLA-A02:01","HLA-A32:01","HLA-B18:01","HLA-B52:01","HLA-C07:01","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A02:01","IC50_HLA-A32:01","IC50_HLA-B18:01","IC50_HLA-B52:01","IC50_HLA-C07:01")
weak=weak[grep('^WT_',weak$Gene_ID,invert=T),]
weak=weak[!duplicated(weak),]
weak=merge(weak,LUAN_RNK,by="genes",all.x=T)
relapse=read.table('LUAN_MHC_I_MinorAntigens_Binders_Regions_Relapse.txt',fill=T)
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('LUAN_MHC_I_MinorAntigens_Binders_Regions_Diagnosis.txt',fill=T)
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
weak=merge(weak,relapse,by=c("chr","pos"),all.x=T)
weak=merge(weak,diagnosis,by=c("chr","pos"),all.x=T)
weak$Rel_Mut=str_count(str_to_upper(weak$rel_bases), as.character(weak$mut_nt))
weak$Rel_Wt=str_count(str_to_upper(weak$rel_bases), as.character(weak$wt_nt))
weak$Dx_Mut=str_count(str_to_upper(weak$dx_bases), as.character(weak$mut_nt))
weak$Dx_Wt=str_count(str_to_upper(weak$dx_bases), as.character(weak$wt_nt))
weak=weak[weak$Mutation!='synonymous_variant',]

nrow(weak)
#Antigeni minori Weak totali
#490

#Antigeni Minori Espressi
weak_filtered=weak[complete.cases(weak[ , c(21,22,31,33)]),]
weak_filtered=weak_filtered[weak_filtered$Rel_Mut>=3 | weak_filtered$Dx_Mut>=3, ]
#139

#Numero antigeni Minori Allele Dx_Mut espresso e Rel_Mut no 
nrow(weak_filtered[weak_filtered$Dx_Mut>=3 & weak_filtered$Rel_Mut<3,])
#9
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut<3,])
#19
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut>=3,])
#111




#Numero antigeni Minori Allele Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Dx_Mut>=3,])
#120
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3,])
#130

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut > weak_filtered$Dx_Mut,])
#62
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut < weak_filtered$Dx_Mut,])
#69

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(weak_filtered[weak_filtered$LUAN1_kallisto > weak_filtered$LUAN2_kallisto,])
#95
#Numero Antigeni Minori Trascritto Diagnosi < Trascritto Relapse
nrow(weak_filtered[weak_filtered$LUAN1_kallisto < weak_filtered$LUAN2_kallisto,])
#44


Luan_Binders_NoExprFilt=rbind(strong,weak)

Luan_Binders_NoExprFilt_2=Luan_Binders_NoExprFilt[,c(1,2,3,5,11:13,16:22,31:34)]
sel <- grepl("IC50_HLA",names(Luan_Binders_NoExprFilt_2))
Luan_Binders_NoExprFilt_2[sel] <- lapply(Luan_Binders_NoExprFilt_2[sel], function(x) replace(x,x >500 , 0) )
Luan_Binders_Ic50_NoExprFilt=data.table::as.data.table(Luan_Binders_NoExprFilt_2)
Luan_Binders_Ic50_NoExprFilt=data.table::melt.data.table(data=Luan_Binders_Ic50_NoExprFilt,id=c(1:7,13:18))
Luan_Binders_Ic50_NoExprFilt=Luan_Binders_Ic50_NoExprFilt[Luan_Binders_Ic50_NoExprFilt$value>0,]
Luan_Binders_Ic50_NoExprFilt_ord=Luan_Binders_Ic50_NoExprFilt[order(Luan_Binders_Ic50_NoExprFilt$value),]
Luan_Binders_Ic50_NoExprFilt_ord=unique(Luan_Binders_Ic50_NoExprFilt_ord)
Luan_Binders_Ic50_NoExprFilt_ord_homo=Luan_Binders_Ic50_NoExprFilt_ord[Luan_Binders_Ic50_NoExprFilt_ord$variable=='IC50_HLA-C07:01',]
Luan_Binders_Ic50_NoExprFilt_ord_homo=rbind(Luan_Binders_Ic50_NoExprFilt_ord_homo,Luan_Binders_Ic50_NoExprFilt_ord)
Luan_Binders_Ic50_NoExprFilt_ord_homo=Luan_Binders_Ic50_NoExprFilt_ord_homo[order(Luan_Binders_Ic50_NoExprFilt_ord_homo$value),]
Luan_Binders_Ic50_NoExprFilt_Dx=Luan_Binders_Ic50_NoExprFilt_ord_homo[Luan_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3,]
Luan_Binders_Ic50_NoExprFilt_Rel=Luan_Binders_Ic50_NoExprFilt_ord_homo[Luan_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3,]


#MABI
sh Launcher_snp2epi_Paper_Annotation_MinorAntigens_ClassI.sh /lustre2/scratch/fsantaniello/VagoRelapse/MABI_Germline_vs_Donor.tsv.gz MinorAntigens_MABI1_MABI2_MHCI_NoMismatchHla  HLA-A23:01,HLA-B7:02,HLA-C12:03 MABI_Expression_Kallisto_NeoEpitopes.txt

		weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A23:01\tHLA-B7:02\tHLA-C12:03" )
		more +2 ./RESULTS_MinorAntigens_MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0; if ($11>50&&$11<=500) print $3,$2,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_MinorAntigens_MABI1_MABI2_MHCI_NoMismatchHla_weak_binders_MinorAntigens

		strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A23:01\tHLA-B7:02\tHLA-C12:03" )
		more +2 ./RESULTS_MinorAntigens_MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0; if ($11>0&&$11<=50) print $3,$2,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_MinorAntigens_MABI1_MABI2_MHCI_NoMismatchHla_strong_binders_MinorAntigens

	  
import sys

Mutationi=open('./HomeQlogin/MinorAntigens_MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_MABI1_MABI2_MHCI_NoMismatchHla_weak_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_MABI_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/MinorAntigens_MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_MABI1_MABI2_MHCI_NoMismatchHla_strong_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_MABI1_MABI2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_MABI_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()


cat StrongBinders_Annotati_MABI_MinorAntigens WeakBinders_Annotati_MABI_MinorAntigens | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ | grep -v synonymous_variant >  Summary_MABI_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches.txt


MABI1 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/MABI1_Allfiles_sorted.bam
MABI2 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/MABI2_Allfiles_sorted.bam

awk '{OFS="\t"; print $8":"$9"-"$9}' Summary_MABI_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches.txt |  sort | uniq  > Summary_MABI_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches_Regions.txt
 while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/MABI1_Allfiles_sorted.bam >> MABI_MHC_I_MinorAntigens_Binders_Regions_Diagnosis.txt; done < Summary_MABI_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches_Regions.txt
 while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/MABI2_Allfiles_sorted.bam >> MABI_MHC_I_MinorAntigens_Binders_Regions_Relapse.txt; done < Summary_MABI_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches_Regions.txt


mabi_NeoAG_mhci_MinorAntigens=read.table('Summary_MABI_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches.txt')
colnames(mabi_NeoAG_mhci_MinorAntigens)=c("Gene_ID","Peptide","HLA-A23:01","HLA-B7:02","HLA-C12:03","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A23:01","IC50_HLA-B7:02","IC50_HLA-C12:03")


strong=read.table('StrongBinders_Annotati_MABI_MinorAntigens',skip=2,head=F)
colnames(strong)=c("Gene_ID","Peptide","HLA-A23:01","HLA-B7:02","HLA-C12:03","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A23:01","IC50_HLA-B7:02","IC50_HLA-C12:03")
strong=strong[grep('^WT_',strong$Gene_ID,invert=T),]
strong=strong[!duplicated(strong),]
strong=merge(strong,MABI_RNK,by="genes",all.x=T)
relapse=read.table('MABI_MHC_I_MinorAntigens_Binders_Regions_Relapse.txt',fill=T)
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('MABI_MHC_I_MinorAntigens_Binders_Regions_Diagnosis.txt',fill=T)
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
strong=merge(strong,relapse,by=c("chr","pos"),all.x=T)
strong=merge(strong,diagnosis,by=c("chr","pos"),all.x=T)
strong$Rel_Mut=str_count(str_to_upper(strong$rel_bases), as.character(strong$mut_nt))
strong$Rel_Wt=str_count(str_to_upper(strong$rel_bases), as.character(strong$wt_nt))
strong$Dx_Mut=str_count(str_to_upper(strong$dx_bases), as.character(strong$mut_nt))
strong$Dx_Wt=str_count(str_to_upper(strong$dx_bases), as.character(strong$wt_nt))
strong=strong[strong$Mutation!='synonymous_variant',]
nrow(strong)
#Antigeni minori Strong totali
#170

#Antigeni Minori Espressi
strong_filtered=strong[complete.cases(strong[ , c(17,18,27,29)]),]
strong_filtered=strong_filtered[strong_filtered$Rel_Mut>=3 | strong_filtered$Dx_Mut>=3, ]
#51

#sel <- grepl("IC50_HLA",names(strong_filtered))
#strong_filtered[sel] <- lapply(strong_filtered[sel], function(x) replace(x,x >50 , 0) )
#strong_filtered$IC50_Binder=rowSums(strong_filtered[sel])
#mabi_strong_filtered=strong_filtered[order(strong_filtered$IC50_Binder),]




#Numero antigeni Minori Allele Dx_Mut espresso e Rel_Mut no 
nrow(strong_filtered[strong_filtered$Dx_Mut>=3 & strong_filtered$Rel_Mut<3,])
#9
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut<3,])
#5
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut>=3,])
#37



#Numero antigeni Minori Allele Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Dx_Mut>=3,])
#46
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3,])
#42

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut > strong_filtered$Dx_Mut,])
#27
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut < strong_filtered$Dx_Mut,])
#24

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(strong_filtered[strong_filtered$MABI1_kallisto > strong_filtered$MABI2_kallisto,])
#27
#Numero Antigeni Minori Trascritto Diagnosi < Trascritto Relapse
nrow(strong_filtered[strong_filtered$MABI1_kallisto < strong_filtered$MABI2_kallisto,])
#24


weak=read.table('WeakBinders_Annotati_MABI_MinorAntigens',skip=2,head=F)
colnames(weak)=c("Gene_ID","Peptide","HLA-A23:01","HLA-B7:02","HLA-C12:03","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A23:01","IC50_HLA-B7:02","IC50_HLA-C12:03")
weak=weak[grep('^WT_',weak$Gene_ID,invert=T),]
weak=weak[!duplicated(weak),]
weak=merge(weak,MABI_RNK,by="genes",all.x=T)
relapse=read.table('MABI_MHC_I_MinorAntigens_Binders_Regions_Relapse.txt',fill=T)
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('MABI_MHC_I_MinorAntigens_Binders_Regions_Diagnosis.txt',fill=T)
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
weak=merge(weak,relapse,by=c("chr","pos"),all.x=T)
weak=merge(weak,diagnosis,by=c("chr","pos"),all.x=T)
weak$Rel_Mut=str_count(str_to_upper(weak$rel_bases), as.character(weak$mut_nt))
weak$Rel_Wt=str_count(str_to_upper(weak$rel_bases), as.character(weak$wt_nt))
weak$Dx_Mut=str_count(str_to_upper(weak$dx_bases), as.character(weak$mut_nt))
weak$Dx_Wt=str_count(str_to_upper(weak$dx_bases), as.character(weak$wt_nt))
weak=weak[weak$Mutation!='synonymous_variant',]
nrow(weak)
#Antigeni minori Weak totali
#597

#Antigeni Minori Espressi
weak_filtered=weak[complete.cases(weak[ , c(17,18,27,29)]),]
weak_filtered=weak_filtered[weak_filtered$Rel_Mut>=3 | weak_filtered$Dx_Mut>=3, ]
#201

#Numero antigeni Minori Allele Dx_Mut espresso e Rel_Mut no 
nrow(weak_filtered[weak_filtered$Dx_Mut>=3 & weak_filtered$Rel_Mut<3,])
#17
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut<3,])
#21
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut>=3,])
#163



#Numero antigeni Minori Allele Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Dx_Mut>=3,])
#180
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3,])
#180

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut > weak_filtered$Dx_Mut,])
#95
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut < weak_filtered$Dx_Mut,])
#94

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(weak_filtered[weak_filtered$MABI1_kallisto > weak_filtered$MABI2_kallisto,])
#114
#Numero Antigeni Minori Trascritto Diagnosi < Trascritto Relapse
nrow(weak_filtered[weak_filtered$MABI1_kallisto < weak_filtered$MABI2_kallisto,])
#87

Mabi_Binders_NoExprFilt=rbind(strong,weak)

Mabi_Binders_NoExprFilt_2=Mabi_Binders_NoExprFilt[,c(1,2,3,5,9:11,14:18,27:30)]
sel <- grepl("IC50_HLA",names(Mabi_Binders_NoExprFilt_2))
Mabi_Binders_NoExprFilt_2[sel] <- lapply(Mabi_Binders_NoExprFilt_2[sel], function(x) replace(x,x >500 , 0) )
Mabi_Binders_Ic50_NoExprFilt=data.table::as.data.table(Mabi_Binders_NoExprFilt_2)
Mabi_Binders_Ic50_NoExprFilt=data.table::melt.data.table(data=Mabi_Binders_Ic50_NoExprFilt,id=c(1:7,11:16))
Mabi_Binders_Ic50_NoExprFilt=Mabi_Binders_Ic50_NoExprFilt[Mabi_Binders_Ic50_NoExprFilt$value>0,]
Mabi_Binders_Ic50_NoExprFilt_ord=Mabi_Binders_Ic50_NoExprFilt[order(Mabi_Binders_Ic50_NoExprFilt$value),]
Mabi_Binders_Ic50_NoExprFilt_ord=unique(Mabi_Binders_Ic50_NoExprFilt_ord)
Mabi_Binders_Ic50_NoExprFilt_Dx=Mabi_Binders_Ic50_NoExprFilt_ord[Mabi_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3,]
Mabi_Binders_Ic50_NoExprFilt_Rel=Mabi_Binders_Ic50_NoExprFilt_ord[Mabi_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3,]


#MOGE
sh Launcher_snp2epi_Paper_Annotation_MinorAntigens_ClassI.sh /lustre2/scratch/fsantaniello/VagoRelapse/MOGE_Germline_vs_Donor.tsv.gz MinorAntigens_MOGE1_MOGE2_MHCI_NoMismatchHla HLA-A24:02,HLA-B13:02,HLA-B44:02,HLA-C05:01,HLA-C06:02 MOGE_Expression_Kallisto_NeoEpitopes.txt


		weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A24:02\tHLA-B13:02\tHLA-B44:02\tHLA-C05:01\tHLA-C06:02" )
		more +2 ./RESULTS_MinorAntigens_MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_MinorAntigens_MOGE1_MOGE2_MHCI_NoMismatchHla_weak_binders_MinorAntigens

		strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A24:02\tHLA-B13:02\tHLA-B44:02\tHLA-C05:01\tHLA-C06:02" )
		more +2 ./RESULTS_MinorAntigens_MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_MinorAntigens_MOGE1_MOGE2_MHCI_NoMismatchHla_strong_binders_MinorAntigens

		
import sys

Mutationi=open('./HomeQlogin/MinorAntigens_MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_MOGE1_MOGE2_MHCI_NoMismatchHla_weak_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_MOGE_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))

out.close()


Mutationi=open('./HomeQlogin/	').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_MOGE1_MOGE2_MHCI_NoMismatchHla_strong_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_MOGE1_MOGE2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_MOGE_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))

out.close()



cat StrongBinders_Annotati_MOGE_MinorAntigens WeakBinders_Annotati_MOGE_MinorAntigens | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ | grep -v synonymous_variant >  Summary_MOGE_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches.txt


MOGE1 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE1_Allfiles_sorted.bam
MOGE2 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE2_Allfiles_sorted.bam

awk '{OFS="\t"; print $10":"$11"-"$11}'  Summary_MOGE_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches.txt  |  sort | uniq  > Summary_MOGE_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches_Regions.txt
while read region;  do samtools mpileup -r $region   /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE1_Allfiles_sorted.bam >> MOGE_MHC_I_MinorAntigens_Binders_Regions_Diagnosis.txt; done < Summary_MOGE_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches_Regions.txt
while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/MOGE2_Allfiles_sorted.bam >> MOGE_MHC_I_MinorAntigens_Binders_Regions_Relapse.txt; done < Summary_MOGE_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches_Regions.txt

moge_NeoAG_mhci_MinorAntigens=read.table('Summary_MOGE_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches.txt')
colnames(moge_NeoAG_mhci_MinorAntigens)=c("Gene_ID","Peptide","HLA-A24:02","HLA-B13:02","HLA-B44:02","HLA-C05:01","HLA-C06:02","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A24:02","IC50_HLA-B13:02","IC50_HLA-B44:02","IC50_HLA-C05:01","IC50_HLA-C06:02")


strong=read.table('StrongBinders_Annotati_MOGE_MinorAntigens',skip=2,head=F)
colnames(strong)=c("Gene_ID","Peptide","HLA-A24:02","HLA-B13:02","HLA-B44:02","HLA-C05:01","HLA-C06:02","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A24:02","IC50_HLA-B13:02","IC50_HLA-B44:02","IC50_HLA-C05:01","IC50_HLA-C06:02")
strong=strong[grep('^WT_',strong$Gene_ID,invert=T),]
strong=strong[!duplicated(strong),]
strong=merge(strong,MOGE_RNK,by="genes",all.x=T)
relapse=read.table('MOGE_MHC_I_MinorAntigens_Binders_Regions_Relapse.txt',fill=T)
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('MOGE_MHC_I_MinorAntigens_Binders_Regions_Diagnosis.txt',fill=T)
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
strong=merge(strong,relapse,by=c("chr","pos"),all.x=T)
strong=merge(strong,diagnosis,by=c("chr","pos"),all.x=T)
strong$Rel_Mut=str_count(str_to_upper(strong$rel_bases), as.character(strong$mut_nt))
strong$Rel_Wt=str_count(str_to_upper(strong$rel_bases), as.character(strong$wt_nt))
strong$Dx_Mut=str_count(str_to_upper(strong$dx_bases), as.character(strong$mut_nt))
strong$Dx_Wt=str_count(str_to_upper(strong$dx_bases), as.character(strong$wt_nt))
strong=strong[strong$Mutation!='synonymous_variant',]
nrow(strong)
#Antigeni minori Strong totali
#86

#Antigeni Minori Espressi
strong_filtered=strong[complete.cases(strong[ , c(21,22,31,33)]),]
strong_filtered=strong_filtered[strong_filtered$Rel_Mut>=3 | strong_filtered$Dx_Mut>=3, ]
#28

#sel <- grepl("IC50_HLA",names(strong_filtered))
#strong_filtered[sel] <- lapply(strong_filtered[sel], function(x) replace(x,x >50 , 0) )
#strong_filtered$IC50_Binder=rowSums(strong_filtered[sel])
#moge_strong_filtered=strong_filtered[order(strong_filtered$IC50_Binder),]

#Numero antigeni Minori Allele Dx_Mut espresso e Rel_Mut no 
nrow(strong_filtered[strong_filtered$Dx_Mut>=3 & strong_filtered$Rel_Mut<3,])
#4
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut<3,])
#2
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut>=3,])
#22


#sel <- grepl("IC50_HLA",names(strong_filtered))
#strong_filtered[sel] <- lapply(strong_filtered[sel], function(x) replace(x,x >50 , 0) )
#strong_filtered$IC50_Binder=rowSums(strong_filtered[sel])
#moge_strong_filtered=strong_filtered[order(strong_filtered$IC50_Binder),]




#Numero antigeni Minori Allele Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Dx_Mut>=3,])
#26
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3,])
#24

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut > strong_filtered$Dx_Mut,])
#15
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut < strong_filtered$Dx_Mut,])
#12

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(strong_filtered[strong_filtered$MOGE1_kallisto > strong_filtered$MOGE2_kallisto,])
#15
#Numero Antigeni Minori Trascritto Diagnosi < Trascritto Relapse
nrow(strong_filtered[strong_filtered$MOGE1_kallisto < strong_filtered$MOGE2_kallisto,])
#13


weak=read.table('WeakBinders_Annotati_MOGE_MinorAntigens',skip=2,head=F)
colnames(weak)=c("Gene_ID","Peptide","HLA-A24:02","HLA-B13:02","HLA-B44:02","HLA-C05:01","HLA-C06:02","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A24:02","IC50_HLA-B13:02","IC50_HLA-B44:02","IC50_HLA-C05:01","IC50_HLA-C06:02")
weak=weak[grep('^WT_',weak$Gene_ID,invert=T),]
weak=weak[!duplicated(weak),]
weak=merge(weak,MOGE_RNK,by="genes",all.x=T)
relapse=read.table('MOGE_MHC_I_MinorAntigens_Binders_Regions_Relapse.txt',fill=T)
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('MOGE_MHC_I_MinorAntigens_Binders_Regions_Diagnosis.txt',fill=T)
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
weak=merge(weak,relapse,by=c("chr","pos"),all.x=T)
weak=merge(weak,diagnosis,by=c("chr","pos"),all.x=T)
weak$Rel_Mut=str_count(str_to_upper(weak$rel_bases), as.character(weak$mut_nt))
weak$Rel_Wt=str_count(str_to_upper(weak$rel_bases), as.character(weak$wt_nt))
weak$Dx_Mut=str_count(str_to_upper(weak$dx_bases), as.character(weak$mut_nt))
weak$Dx_Wt=str_count(str_to_upper(weak$dx_bases), as.character(weak$wt_nt))
weak=weak[weak$Mutation!='synonymous_variant',]

nrow(weak)
#Antigeni minori Weak totali
#526

#Antigeni Minori Espressi
weak_filtered=weak[complete.cases(weak[ , c(21,22,31,33)]),]
weak_filtered=weak_filtered[weak_filtered$Rel_Mut>=3 | weak_filtered$Dx_Mut>=3, ]
#168

#Numero antigeni Minori Allele Dx_Mut espresso e Rel_Mut no 
nrow(weak_filtered[weak_filtered$Dx_Mut>=3 & weak_filtered$Rel_Mut<3,])
#16
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut<3,])
#14
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut>=3,])
#138

#Numero antigeni Minori Allele Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Dx_Mut>=3,])
#154
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3,])
#152

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut > weak_filtered$Dx_Mut,])
#89
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut < weak_filtered$Dx_Mut,])
#71

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(weak_filtered[weak_filtered$MOGE1_kallisto > weak_filtered$MOGE2_kallisto,])
#89
#Numero Antigeni Minori Trascritto Diagnosi < Trascritto Relapse
nrow(weak_filtered[weak_filtered$MOGE1_kallisto < weak_filtered$MOGE2_kallisto,])
#79

Moge_Binders_NoExprFilt=rbind(strong,weak)

Moge_Binders_NoExprFilt_2=Moge_Binders_NoExprFilt[,c(1,2,3,5,11:13,16:22,31:34)]
sel <- grepl("IC50_HLA",names(Moge_Binders_NoExprFilt_2))
Moge_Binders_NoExprFilt_2[sel] <- lapply(Moge_Binders_NoExprFilt_2[sel], function(x) replace(x,x >500 , 0) )
Moge_Binders_Ic50_NoExprFilt=data.table::as.data.table(Moge_Binders_NoExprFilt_2)
Moge_Binders_Ic50_NoExprFilt=data.table::melt.data.table(data=Moge_Binders_Ic50_NoExprFilt,id=c(1:7,13:18))
Moge_Binders_Ic50_NoExprFilt=Moge_Binders_Ic50_NoExprFilt[Moge_Binders_Ic50_NoExprFilt$value>0,]
Moge_Binders_Ic50_NoExprFilt_ord=Moge_Binders_Ic50_NoExprFilt[order(Moge_Binders_Ic50_NoExprFilt$value),]
Moge_Binders_Ic50_NoExprFilt_ord=unique(Moge_Binders_Ic50_NoExprFilt_ord)
Moge_Binders_Ic50_NoExprFilt_ord_homo=Moge_Binders_Ic50_NoExprFilt_ord[Moge_Binders_Ic50_NoExprFilt_ord$variable=='IC50_HLA-A24:02',]
Moge_Binders_Ic50_NoExprFilt_ord_homo=rbind(Moge_Binders_Ic50_NoExprFilt_ord_homo,Moge_Binders_Ic50_NoExprFilt_ord)
Moge_Binders_Ic50_NoExprFilt_ord_homo=Moge_Binders_Ic50_NoExprFilt_ord_homo[order(Moge_Binders_Ic50_NoExprFilt_ord_homo$value),]
Moge_Binders_Ic50_NoExprFilt_Dx=Moge_Binders_Ic50_NoExprFilt_ord_homo[Moge_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3,]
Moge_Binders_Ic50_NoExprFilt_Rel=Moge_Binders_Ic50_NoExprFilt_ord_homo[Moge_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3,]



#PIAG
sh Launcher_snp2epi_Paper_Annotation_MinorAntigens_ClassI.sh /lustre2/scratch/fsantaniello/VagoRelapse/PIAG_Germline_vs_Donor.tsv.gz MinorAntigens_PIAG1_PIAG2_MHCI_NoMismatchHla HLA-A03:02,HLA-A24:02,HLA-B44:03,HLA-B35:02,HLA-C04:01 PIAG_Expression_Kallisto_NeoEpitopes.txt

	
		weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A03:02\tHLA-A24:02\tHLA-B44:03\tHLA-B35:02\tHLA-C04:01" )
		more +2 ./RESULTS_MinorAntigens_PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0,0,0; if ($11>50&&$11<=500) print $3,$2,0,0,1,0,0; if ($14>50&&$14<=500) print $3,$2,0,0,0,1,0; if ($17>50&&$17<=500) print $3,$2,0,0,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_MinorAntigens_PIAG1_PIAG2_MHCI_NoMismatchHla_weak_binders_MinorAntigens

		strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A03:02\tHLA-A24:02\tHLA-B44:03\tHLA-B35:02\tHLA-C04:01" )
		more +2 ./RESULTS_MinorAntigens_PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0,0,0; if ($11>0&&$11<=50) print $3,$2,0,0,1,0,0; if ($14>0&&$14<=50) print $3,$2,0,0,0,1,0; if ($17>0&&$17<=50) print $3,$2,0,0,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_MinorAntigens_PIAG1_PIAG2_MHCI_NoMismatchHla_strong_binders_MinorAntigens

import sys

Mutationi=open('./HomeQlogin/MinorAntigens_PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_PIAG1_PIAG2_MHCI_NoMismatchHla_weak_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_PIAG_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))


out.close()


Mutationi=open('./HomeQlogin/MinorAntigens_PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_PIAG1_PIAG2_MHCI_NoMismatchHla_strong_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_PIAG1_PIAG2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_PIAG_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\t'+str(j.split('\t')[13])+'\t'+str(j.split('\t')[16])+'\n'))

out.close()
		


cat StrongBinders_Annotati_PIAG_MinorAntigens WeakBinders_Annotati_PIAG_MinorAntigens | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ | grep -v synonymous_variant  > Summary_PIAG_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches.txt

PIAG1 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/PIAG1_Allfiles_sorted.bam
PIAG2 = /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/PIAG2_Allfiles_sorted.bam


awk '{OFS="\t"; print $10":"$11"-"$11}'  Summary_PIAG_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches.txt  |  sort | uniq  > Summary_PIAG_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches_Regions.txt
while read region;  do samtools mpileup -r $region   /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/PIAG1_Allfiles_sorted.bam >> PIAG_MHC_I_MinorAntigens_Binders_Regions_Diagnosis.txt; done < Summary_PIAG_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches_Regions.txt
while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/170531_SN859_0446_BHHM2GBCXY/Project_Vago_161_Relapsing_Leukemia/PIAG2_Allfiles_sorted.bam >> PIAG_MHC_I_MinorAntigens_Binders_Regions_Relapse.txt;     done < Summary_PIAG_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches_Regions.txt

piag_NeoAG_mhci_MinorAntigens=read.table('Summary_PIAG_Annotato_MinorAntigens_MHCI_NoHeaders_NoMismatches.txt')
colnames(piag_NeoAG_mhci_MinorAntigens)=c("Gene_ID","Peptide","HLA-A03:02","HLA-A24:02","HLA-B44:03","HLA-B35:02","HLA-C04:01","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A03:02","IC50_HLA-A24:02","IC50_HLA-B44:03","IC50_HLA-B35:02","IC50_HLA-C04:01")


strong=read.table('StrongBinders_Annotati_PIAG_MinorAntigens',skip=2,head=F)
colnames(strong)=c("Gene_ID","Peptide","HLA-A03:02","HLA-A24:02","HLA-B44:03","HLA-B35:02","HLA-C04:01","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A03:02","IC50_HLA-A24:02","IC50_HLA-B44:03","IC50_HLA-B35:02","IC50_HLA-C04:01")
strong=strong[grep('^WT_',strong$Gene_ID,invert=T),]
strong=strong[!duplicated(strong),]
strong=merge(strong,PIAG_RNK,by="genes",all.x=T)
relapse=read.table('PIAG_MHC_I_MinorAntigens_Binders_Regions_Relapse.txt',fill=T)
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('PIAG_MHC_I_MinorAntigens_Binders_Regions_Diagnosis.txt',fill=T)
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
strong=merge(strong,relapse,by=c("chr","pos"),all.x=T)
strong=merge(strong,diagnosis,by=c("chr","pos"),all.x=T)
strong$Rel_Mut=str_count(str_to_upper(strong$rel_bases), as.character(strong$mut_nt))
strong$Rel_Wt=str_count(str_to_upper(strong$rel_bases), as.character(strong$wt_nt))
strong$Dx_Mut=str_count(str_to_upper(strong$dx_bases), as.character(strong$mut_nt))
strong$Dx_Wt=str_count(str_to_upper(strong$dx_bases), as.character(strong$wt_nt))
strong=strong[strong$Mutation!='synonymous_variant',]
nrow(strong)
#Antigeni minori Strong totali
#72

#Antigeni Minori Espressi
strong_filtered=strong[complete.cases(strong[ , c(21,22,31,33)]),]
strong_filtered=strong_filtered[strong_filtered$Rel_Mut>=3 | strong_filtered$Dx_Mut>=3, ]
#19

#sel <- grepl("IC50_HLA",names(strong_filtered))
#strong_filtered[sel] <- lapply(strong_filtered[sel], function(x) replace(x,x >50 , 0) )
#strong_filtered$IC50_Binder=rowSums(strong_filtered[sel])
#piag_strong_filtered=strong_filtered[order(strong_filtered$IC50_Binder),]
#


nrow(strong_filtered[strong_filtered$Dx_Mut>=3 & strong_filtered$Rel_Mut<3,])
#0
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut<3,])
#1
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut>=3,])
#18


#Numero antigeni Minori Allele Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Dx_Mut>=3,])
#18
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3,])
#19

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut > strong_filtered$Dx_Mut,])
#10
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut < strong_filtered$Dx_Mut,])
#9

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(strong_filtered[strong_filtered$PIAG1_kallisto > strong_filtered$PIAG2_kallisto,])
#9
#Numero Antigeni Minori Trascritto Diagnosi < Trascritto Relapse
nrow(strong_filtered[strong_filtered$PIAG1_kallisto < strong_filtered$PIAG2_kallisto,])
#10



weak=read.table('WeakBinders_Annotati_PIAG_MinorAntigens',skip=2,head=F)
colnames(weak)=c("Gene_ID","Peptide","HLA-A03:02","HLA-A24:02","HLA-B44:03","HLA-B35:02","HLA-C04:01","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A03:02","IC50_HLA-A24:02","IC50_HLA-B44:03","IC50_HLA-B35:02","IC50_HLA-C04:01")
weak=weak[grep('^WT_',weak$Gene_ID,invert=T),]
weak=weak[!duplicated(weak),]
weak=merge(weak,PIAG_RNK,by="genes",all.x=T)
relapse=read.table('PIAG_MHC_I_MinorAntigens_Binders_Regions_Relapse.txt',fill=T)
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('PIAG_MHC_I_MinorAntigens_Binders_Regions_Diagnosis.txt',fill=T)
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
weak=merge(weak,relapse,by=c("chr","pos"),all.x=T)
weak=merge(weak,diagnosis,by=c("chr","pos"),all.x=T)
weak$Rel_Mut=str_count(str_to_upper(weak$rel_bases), as.character(weak$mut_nt))
weak$Rel_Wt=str_count(str_to_upper(weak$rel_bases), as.character(weak$wt_nt))
weak$Dx_Mut=str_count(str_to_upper(weak$dx_bases), as.character(weak$mut_nt))
weak$Dx_Wt=str_count(str_to_upper(weak$dx_bases), as.character(weak$wt_nt))
weak=weak[weak$Mutation!='synonymous_variant',]
nrow(weak)
#Antigeni minori Weak totali
#502

#Antigeni Minori Espressi
weak_filtered=weak[complete.cases(weak[ , c(21,22,31,33)]),]
weak_filtered=weak_filtered[weak_filtered$Rel_Mut>=3 | weak_filtered$Dx_Mut>=3, ]
#157

nrow(weak_filtered[weak_filtered$Dx_Mut>=3 & weak_filtered$Rel_Mut<3,])
#20
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut<3,])
#13
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut>=3,])
#124




#Numero antigeni Minori Allele Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Dx_Mut>=3,])
#144
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3,])
#137

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut > weak_filtered$Dx_Mut,])
#51
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut < weak_filtered$Dx_Mut,])
#103

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(weak_filtered[weak_filtered$PIAG1_kallisto > weak_filtered$PIAG2_kallisto,])
#93
#Numero Antigeni Minori Trascritto Diagnosi < Trascritto Relapse
nrow(weak_filtered[weak_filtered$PIAG1_kallisto < weak_filtered$PIAG2_kallisto,])
#64

Piag_Binders_NoExprFilt=rbind(strong,weak)

Piag_Binders_NoExprFilt_2=Piag_Binders_NoExprFilt[,c(1,2,3,5,11:13,16:22,31:34)]
sel <- grepl("IC50_HLA",names(Piag_Binders_NoExprFilt_2))
Piag_Binders_NoExprFilt_2[sel] <- lapply(Piag_Binders_NoExprFilt_2[sel], function(x) replace(x,x >500 , 0) )
Piag_Binders_Ic50_NoExprFilt=data.table::as.data.table(Piag_Binders_NoExprFilt_2)
Piag_Binders_Ic50_NoExprFilt=data.table::melt.data.table(data=Piag_Binders_Ic50_NoExprFilt,id=c(1:7,13:18))
Piag_Binders_Ic50_NoExprFilt=Piag_Binders_Ic50_NoExprFilt[Piag_Binders_Ic50_NoExprFilt$value>0,]
Piag_Binders_Ic50_NoExprFilt_ord=Piag_Binders_Ic50_NoExprFilt[order(Piag_Binders_Ic50_NoExprFilt$value),]
Piag_Binders_Ic50_NoExprFilt_ord=unique(Piag_Binders_Ic50_NoExprFilt_ord)
Piag_Binders_Ic50_NoExprFilt_ord_homo=Piag_Binders_Ic50_NoExprFilt_ord[Piag_Binders_Ic50_NoExprFilt_ord$variable=='IC50_HLA-C04:01',]
Piag_Binders_Ic50_NoExprFilt_ord_homo=rbind(Piag_Binders_Ic50_NoExprFilt_ord_homo,Piag_Binders_Ic50_NoExprFilt_ord)
Piag_Binders_Ic50_NoExprFilt_ord_homo=Piag_Binders_Ic50_NoExprFilt_ord_homo[order(Piag_Binders_Ic50_NoExprFilt_ord_homo$value),]
Piag_Binders_Ic50_NoExprFilt_Dx=Piag_Binders_Ic50_NoExprFilt_ord_homo[Piag_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3,]
Piag_Binders_Ic50_NoExprFilt_Rel=Piag_Binders_Ic50_NoExprFilt_ord_homo[Piag_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3,]



#PRELU
sh Launcher_snp2epi_Paper_Annotation_MinorAntigens_ClassI.sh /lustre2/scratch/fsantaniello/VagoRelapse/PRELU_Germline_vs_Donor.tsv.gz MinorAntigens_PRELU1_PRELU2_MHCI_NoMismatchHla  HLA-A01:01,HLA-B44:02,HLA-C07:04 PRELU_Expression_Kallisto_NeoEpitopes.txt

		weak_binders=$(echo "Weak binders\nGene\tPeptide\tHLA-A24:02\tHLA-B35:01\tHLA-C04:01")
		more +2 ./RESULTS_MinorAntigens_PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls|  awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0,0;  if ($8>50&&$8<=500) print $3,$2,0,1,0; if ($11>50&&$11<=500) print $3,$2,0,0,1}' - | (echo -e $weak_binders; cat -) > ./RESULTS_MinorAntigens_PRELU1_PRELU2_MHCI_NoMismatchHla_weak_binders_MinorAntigens

		strong_binders=$(echo "Strong binders\nGene\tPeptide\tHLA-A24:02\tHLA-B35:01\tHLA-C04:01")
		more +2 ./RESULTS_MinorAntigens_PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls |  awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0,0;  if ($8>0&&$8<=50) print $3,$2,0,1,0; if ($11>0&&$11<=50) print $3,$2,0,0,1}' - | (echo -e $strong_binders; cat -) > ./RESULTS_MinorAntigens_PRELU1_PRELU2_MHCI_NoMismatchHla_strong_binders_MinorAntigens

	   
	  import sys

Mutationi=open('./HomeQlogin/MinorAntigens_PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_PRELU1_PRELU2_MHCI_NoMismatchHla_weak_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_PRELU_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/MinorAntigens_PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_MinorAntigens_PRELU1_PRELU2_MHCI_NoMismatchHla_strong_binders_MinorAntigens').readlines()
Excel=open('./HomeQlogin/RESULTS_MinorAntigens_PRELU1_PRELU2_MHCI_NoMismatchHla_PeptidesPrediction_MinorAntigens.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_PRELU_MinorAntigens','w')
header=Summary[0:2]

out.write(''.join(header).replace(',HLA','\tHLA').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('HLA','IC50_HLA').rstrip()+'\n')


binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\t'+str(j.split('\t')[10])+'\n'))

out.close()
		


cat StrongBinders_Annotati_PRELU_MinorAntigens WeakBinders_Annotati_PRELU_MinorAntigens | sed '/Strong/,+1 d' - | sed '/Weak/,+1 d' - | grep -v ^WT_ | grep -v synonymous_variant >  Summary_PRELU_Annotato_MinorAntigens_NoHeaders_NoMismatches.txt

PRELU1 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/PRELU1_Allfiles_sorted.bam
PRELU2 = /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/PRELU2_Allfiles_sorted.bam



awk '{OFS="\t"; print $8":"$9"-"$9}'   Summary_PRELU_Annotato_MinorAntigens_NoHeaders_NoMismatches.txt | sort | uniq  > Summary_PRELU_Annotato_MinorAntigens_NoHeaders_NoMismatches_Regions.txt
while read region;  do samtools mpileup -r $region  /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/PRELU1_Allfiles_sorted.bam >> PRELU_MHC_I_MinorAntigens_Binders_Regions_Diagnosis_NoMismatches.txt; done < Summary_PRELU_Annotato_MinorAntigens_NoHeaders_NoMismatches_Regions.txt
while read region;  do samtools mpileup -r $region /lustre1/workspace/Ciceri/161_Leukemia/150727_SN859_0224_AHKYMJADXX/Project_Vago_161_Leukemia/PRELU2_Allfiles_sorted.bam  >> PRELU_MHC_I_MinorAntigens_Binders_Regions_Relapse_NoMismatches.txt; done < Summary_PRELU_Annotato_MinorAntigens_NoHeaders_NoMismatches_Regions.txt
 


prelu_NeoAG_mhci_MinorAntigens=read.table('Summary_PRELU_Annotato_MinorAntigens_NoHeaders_NoMismatches.txt')
colnames(prelu_NeoAG_mhci_MinorAntigens)=c("Gene_ID","Peptide","HLA-A01:01","HLA-B44:02","HLA-C07:04","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A01:01","IC50_HLA-B44:02","IC50_HLA-C07:04")

strong=read.table('StrongBinders_Annotati_PRELU_MinorAntigens',skip=2,head=F)
colnames(strong)=c("Gene_ID","Peptide","HLA-A01:01","HLA-B44:02","HLA-C07:04","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A01:01","IC50_HLA-B44:02","IC50_HLA-C07:04")
strong=strong[grep('^WT_',strong$Gene_ID,invert=T),]
strong=strong[!duplicated(strong),]
strong=merge(strong,PRELU_RNK,by="genes",all.x=T)
relapse=read.table('PRELU_MHC_I_MinorAntigens_Binders_Regions_Relapse_NoMismatches.txt',fill=T)
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('PRELU_MHC_I_MinorAntigens_Binders_Regions_Diagnosis_NoMismatches.txt',fill=T)
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
strong=merge(strong,relapse,by=c("chr","pos"),all.x=T)
strong=merge(strong,diagnosis,by=c("chr","pos"),all.x=T)
strong$Rel_Mut=str_count(str_to_upper(strong$rel_bases), as.character(strong$mut_nt))
strong$Rel_Wt=str_count(str_to_upper(strong$rel_bases), as.character(strong$wt_nt))
strong$Dx_Mut=str_count(str_to_upper(strong$dx_bases), as.character(strong$mut_nt))
strong$Dx_Wt=str_count(str_to_upper(strong$dx_bases), as.character(strong$wt_nt))
strong=strong[strong$Mutation!='synonymous_variant',]
nrow(strong)
#Antigeni minori Strong totali
#19

#Antigeni Minori Espressi
strong_filtered=strong[complete.cases(strong[ , c(17,18,27,29)]),]
strong_filtered=strong_filtered[strong_filtered$Rel_Mut>=3 | strong_filtered$Dx_Mut>=3, ]
#7

#sel <- grepl("IC50_HLA",names(strong_filtered))
#strong_filtered[sel] <- lapply(strong_filtered[sel], function(x) replace(x,x >50 , 0) )
#strong_filtered$IC50_Binder=rowSums(strong_filtered[sel])
#prelu_strong_filtered=strong_filtered[order(strong_filtered$IC50_Binder),]



nrow(strong_filtered[strong_filtered$Dx_Mut>=3 & strong_filtered$Rel_Mut<3,])
#1
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut<3,])
#2
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3 & strong_filtered$Dx_Mut>=3,])
#4




#Numero antigeni Minori Allele Dx_Mut espresso
nrow(strong_filtered[strong_filtered$Dx_Mut>=3,])
#5
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(strong_filtered[strong_filtered$Rel_Mut>=3,])
#6

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut > strong_filtered$Dx_Mut,])
#5
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(strong_filtered[strong_filtered$Rel_Mut < strong_filtered$Dx_Mut,])
#2

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(strong_filtered[strong_filtered$PRELU1_kallisto > strong_filtered$PRELU2_kallisto,])
#2
#Numero Antigeni Minori Trascritto Diagnosi < Trascritto Relapse
nrow(strong_filtered[strong_filtered$PRELU1_kallisto < strong_filtered$PRELU2_kallisto,])
#5


weak=read.table('WeakBinders_Annotati_PRELU_MinorAntigens',skip=2,head=F)
colnames(weak)=c("Gene_ID","Peptide","HLA-A01:01","HLA-B44:02","HLA-C07:04","Mutation","genes","chr","pos" ,"wt_nt","mut_nt","Expr_Relapse","Expr_Diag","IC50_HLA-A01:01","IC50_HLA-B44:02","IC50_HLA-C07:04")
weak=weak[grep('^WT_',weak$Gene_ID,invert=T),]
weak=weak[!duplicated(weak),]
weak=merge(weak,PRELU_RNK,by="genes",all.x=T)
relapse=read.table('PRELU_MHC_I_MinorAntigens_Binders_Regions_Relapse_NoMismatches.txt',fill=T)
colnames(relapse)=c("chr","pos","N","rel_cov","rel_bases","rel_qual")
diagnosis=read.table('PRELU_MHC_I_MinorAntigens_Binders_Regions_Diagnosis_NoMismatches.txt',fill=T)
colnames(diagnosis)=c("chr","pos","N","dx_cov","dx_bases","dx_qual")
weak=merge(weak,relapse,by=c("chr","pos"),all.x=T)
weak=merge(weak,diagnosis,by=c("chr","pos"),all.x=T)
weak$Rel_Mut=str_count(str_to_upper(weak$rel_bases), as.character(weak$mut_nt))
weak$Rel_Wt=str_count(str_to_upper(weak$rel_bases), as.character(weak$wt_nt))
weak$Dx_Mut=str_count(str_to_upper(weak$dx_bases), as.character(weak$mut_nt))
weak$Dx_Wt=str_count(str_to_upper(weak$dx_bases), as.character(weak$wt_nt))
weak=weak[weak$Mutation!='synonymous_variant',]
nrow(weak)
#Antigeni minori Weak totali
#140

#Antigeni Minori Espressi
weak_filtered=weak[complete.cases(weak[ , c(17,18,27,29)]),]
weak_filtered=weak_filtered[weak_filtered$Rel_Mut>=3 | weak_filtered$Dx_Mut>=3, ]
#48


nrow(weak_filtered[weak_filtered$Dx_Mut>=3 & weak_filtered$Rel_Mut<3,])
#2
#Numero antigeni Minori Allele Rel_Mut espresso e Dx_Mut no 
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut<3,])
#7
#Numero antigeni Minori Allele Entrambi Rel_Mut espresso e Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3 & weak_filtered$Dx_Mut>=3,])
#39


#Numero antigeni Minori Allele Dx_Mut espresso
nrow(weak_filtered[weak_filtered$Dx_Mut>=3,])
#41
#Numero antigeni Minori Allele Rel_Mut espresso
nrow(weak_filtered[weak_filtered$Rel_Mut>=3,])
#46

#Numero Antigeni minori Allele Relapse > Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut > weak_filtered$Dx_Mut,])
#31
#Numero Antigeni minori Allele Relapse < Allele Diagnosi
nrow(weak_filtered[weak_filtered$Rel_Mut < weak_filtered$Dx_Mut,])
#15

#Numero Antigeni Minori Trascritto Diagnosi > Trascritto Relapse
nrow(weak_filtered[weak_filtered$PRELU1_kallisto > weak_filtered$PRELU2_kallisto,])
#23
#Numero Antigeni Minori Trascritto Diagnosi < Trascritto Relapse
nrow(weak_filtered[weak_filtered$PRELU1_kallisto < weak_filtered$PRELU2_kallisto,])
#25


Prelu_Binders_NoExprFilt=rbind(strong,weak)

Prelu_Binders_NoExprFilt_2=Prelu_Binders_NoExprFilt[,c(1,2,3,5,9:11,14:18,27:30)]
sel <- grepl("IC50_HLA",names(Prelu_Binders_NoExprFilt_2))
Prelu_Binders_NoExprFilt_2[sel] <- lapply(Prelu_Binders_NoExprFilt_2[sel], function(x) replace(x,x >500 , 0) )
Prelu_Binders_Ic50_NoExprFilt=data.table::as.data.table(Prelu_Binders_NoExprFilt_2)
Prelu_Binders_Ic50_NoExprFilt=data.table::melt.data.table(data=Prelu_Binders_Ic50_NoExprFilt,id=c(1:7,11:16))
Prelu_Binders_Ic50_NoExprFilt=Prelu_Binders_Ic50_NoExprFilt[Prelu_Binders_Ic50_NoExprFilt$value>0,]
Prelu_Binders_Ic50_NoExprFilt_ord=Prelu_Binders_Ic50_NoExprFilt[order(Prelu_Binders_Ic50_NoExprFilt$value),]
Prelu_Binders_Ic50_NoExprFilt_ord=unique(Prelu_Binders_Ic50_NoExprFilt_ord)
Prelu_Binders_Ic50_NoExprFilt_Dx=Prelu_Binders_Ic50_NoExprFilt_ord[Prelu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3,]
Prelu_Binders_Ic50_NoExprFilt_Rel=Prelu_Binders_Ic50_NoExprFilt_ord[Prelu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3,]




alfe_auc=auc(1:length(Alfe_Binders_Ic50_NoExprFilt_ord$value), y=Alfe_Binders_Ic50_NoExprFilt_ord$value)
besu_auc=auc(1:length(Besu_Binders_Ic50_NoExprFilt_ord$value), y=Besu_Binders_Ic50_NoExprFilt_ord$value)
calu_auc=auc(1:length(Calu_Binders_Ic50_NoExprFilt_ord_homo$value), y=Calu_Binders_Ic50_NoExprFilt_ord_homo$value)
deiv_auc=auc(1:length(Deiv_Binders_Ic50_NoExprFilt_ord$value), y=Deiv_Binders_Ic50_NoExprFilt_ord$value)
dest_auc=auc(1:length(Dest_Binders_Ic50_NoExprFilt_ord$value), y=Dest_Binders_Ic50_NoExprFilt_ord$value)
dr1_auc=auc(1:length(Dr1_Binders_Ic50_NoExprFilt_ord$value), y=Dr1_Binders_Ic50_NoExprFilt_ord$value)
dr4_auc=auc(1:length(Dr4_Binders_Ic50_NoExprFilt_ord$value), y=Dr4_Binders_Ic50_NoExprFilt_ord$value)
dr5_auc=auc(1:length(Dr5_Binders_Ic50_NoExprFilt_ord$value), y=Dr5_Binders_Ic50_NoExprFilt_ord$value)
foca_auc=auc(1:length(Foca_Binders_Ic50_NoExprFilt_ord$value), y=Foca_Binders_Ic50_NoExprFilt_ord$value)
gagra_auc=auc(1:length(Gagra_Binders_Ic50_NoExprFilt_ord$value), y=Gagra_Binders_Ic50_NoExprFilt_ord$value)
luan_auc=auc(1:length(Luan_Binders_Ic50_NoExprFilt_ord_homo$value), y=Luan_Binders_Ic50_NoExprFilt_ord_homo$value)
mabi_auc=auc(1:length(Mabi_Binders_Ic50_NoExprFilt_ord$value), y=Mabi_Binders_Ic50_NoExprFilt_ord$value)
moge_auc=auc(1:length(Moge_Binders_Ic50_NoExprFilt_ord_homo$value), y=Moge_Binders_Ic50_NoExprFilt_ord_homo$value)
piag_auc=auc(1:length(Piag_Binders_Ic50_NoExprFilt_ord_homo$value), y=Piag_Binders_Ic50_NoExprFilt_ord_homo$value)
prelu_auc=auc(1:length(Prelu_Binders_Ic50_NoExprFilt_ord$value), y=Prelu_Binders_Ic50_NoExprFilt_ord$value)


all_aucs=rbind(alfe_auc,besu_auc,calu_auc,deiv_auc,dest_auc,dr1_auc,dr4_auc,dr5_auc,foca_auc,gagra_auc,luan_auc,mabi_auc,moge_auc,piag_auc,prelu_auc)
all_aucs=as.data.frame(all_aucs)
all_aucs$genes=rownames(all_aucs)
all_aucs$RelDon=c(1,0,1,0,1,0,1,1,0,0,1,1,0,0,1)
#all_aucs$RelDon=c(0,1,0,1,0,1,0,0,1,1,0,0,1,1,0)
all_aucs=all_aucs[order(all_aucs$V1),]
colnames(all_aucs)=c('aucs','patients','RelDon')
all_aucs$Donors=ifelse(all_aucs$RelDon==1,'Related','Unrelated')

write.table(all_aucs,'../../Documents/Progetto_LuciaGabri/MhAGs_AUC_Pazienti_NoExpressionFilter.txt',sep='\t',col.names=T,row.names=F,quote=F)


#Per il summary
#Alfe
#N° MhAGs Strong
nrow(Alfe_Binders_Ic50_NoExprFilt_ord[Alfe_Binders_Ic50_NoExprFilt_ord$value<=50,])
#N° MhAGs Weak
nrow(Alfe_Binders_Ic50_NoExprFilt_ord[Alfe_Binders_Ic50_NoExprFilt_ord$value>50,])
#N° Expr. MhAGs  Strong
nrow(Alfe_Binders_Ic50_NoExprFilt_ord[(Alfe_Binders_Ic50_NoExprFilt_ord$ALFE2_kallisto >= 0.58 | Alfe_Binders_Ic50_NoExprFilt_ord$ALFE1_kallisto >=0.58) & Alfe_Binders_Ic50_NoExprFilt_ord$value<=50 & (Alfe_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Alfe_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#N° Expr. MhAGs  Weak
nrow(Alfe_Binders_Ic50_NoExprFilt_ord[(Alfe_Binders_Ic50_NoExprFilt_ord$ALFE2_kallisto >= 0.58 | Alfe_Binders_Ic50_NoExprFilt_ord$ALFE1_kallisto >=0.58) & Alfe_Binders_Ic50_NoExprFilt_ord$value>50 & (Alfe_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Alfe_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele ONLY
nrow(Alfe_Binders_Ic50_NoExprFilt_ord[(Alfe_Binders_Ic50_NoExprFilt_ord$ALFE2_kallisto >= 0.58 | Alfe_Binders_Ic50_NoExprFilt_ord$ALFE1_kallisto >=0.58) & Alfe_Binders_Ic50_NoExprFilt_ord$value<=50  & Alfe_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3 & Alfe_Binders_Ic50_NoExprFilt_ord$Rel_Mut<3, ])
#STRONG  Expressed REL Allele Only
nrow(Alfe_Binders_Ic50_NoExprFilt_ord[(Alfe_Binders_Ic50_NoExprFilt_ord$ALFE2_kallisto >= 0.58 | Alfe_Binders_Ic50_NoExprFilt_ord$ALFE1_kallisto >=0.58) & Alfe_Binders_Ic50_NoExprFilt_ord$value<=50  & Alfe_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 & Alfe_Binders_Ic50_NoExprFilt_ord$Dx_Mut<3 , ])
#STRONG  Expressed Common DIAG REL
nrow(Alfe_Binders_Ic50_NoExprFilt_ord[(Alfe_Binders_Ic50_NoExprFilt_ord$ALFE2_kallisto >= 0.58 | Alfe_Binders_Ic50_NoExprFilt_ord$ALFE1_kallisto >=0.58) & Alfe_Binders_Ic50_NoExprFilt_ord$value<=50 & (Alfe_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 & Alfe_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele Diagnosis
nrow(Alfe_Binders_Ic50_NoExprFilt_ord[(Alfe_Binders_Ic50_NoExprFilt_ord$ALFE2_kallisto >= 0.58 | Alfe_Binders_Ic50_NoExprFilt_ord$ALFE1_kallisto >=0.58) & Alfe_Binders_Ic50_NoExprFilt_ord$value<=50  & Alfe_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3 , ])
#STRONG  Expressed REL Allele Relapse
nrow(Alfe_Binders_Ic50_NoExprFilt_ord[(Alfe_Binders_Ic50_NoExprFilt_ord$ALFE2_kallisto >= 0.58 | Alfe_Binders_Ic50_NoExprFilt_ord$ALFE1_kallisto >=0.58) & Alfe_Binders_Ic50_NoExprFilt_ord$value<=50  & Alfe_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3, ])
#WEAK Expressed DIAG Allele ONLY
nrow(Alfe_Binders_Ic50_NoExprFilt_ord[(Alfe_Binders_Ic50_NoExprFilt_ord$ALFE2_kallisto >= 0.58 | Alfe_Binders_Ic50_NoExprFilt_ord$ALFE1_kallisto >=0.58) & Alfe_Binders_Ic50_NoExprFilt_ord$value>50  & Alfe_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3 & Alfe_Binders_Ic50_NoExprFilt_ord$Rel_Mut < 3, ])
#WEAK  Expressed REL Allele Only
nrow(Alfe_Binders_Ic50_NoExprFilt_ord[(Alfe_Binders_Ic50_NoExprFilt_ord$ALFE2_kallisto >= 0.58 | Alfe_Binders_Ic50_NoExprFilt_ord$ALFE1_kallisto >=0.58) & Alfe_Binders_Ic50_NoExprFilt_ord$value>50  & Alfe_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3 & Alfe_Binders_Ic50_NoExprFilt_ord$Dx_Mut < 3 , ])
#WEAK  Expressed Common DIAG REL
nrow(Alfe_Binders_Ic50_NoExprFilt_ord[(Alfe_Binders_Ic50_NoExprFilt_ord$ALFE2_kallisto >= 0.58 | Alfe_Binders_Ic50_NoExprFilt_ord$ALFE1_kallisto >=0.58) & Alfe_Binders_Ic50_NoExprFilt_ord$value>50 & (Alfe_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3 & Alfe_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3), ])
#WEAK Expressed DIAG Allele Diagnosis
nrow(Alfe_Binders_Ic50_NoExprFilt_ord[(Alfe_Binders_Ic50_NoExprFilt_ord$ALFE2_kallisto >= 0.58 | Alfe_Binders_Ic50_NoExprFilt_ord$ALFE1_kallisto >=0.58) & Alfe_Binders_Ic50_NoExprFilt_ord$value>50  & Alfe_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3 , ])
#WEAK  Expressed REL Allele Relapse
nrow(Alfe_Binders_Ic50_NoExprFilt_ord[(Alfe_Binders_Ic50_NoExprFilt_ord$ALFE2_kallisto >= 0.58 | Alfe_Binders_Ic50_NoExprFilt_ord$ALFE1_kallisto >=0.58) & Alfe_Binders_Ic50_NoExprFilt_ord$value>50  & Alfe_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3, ])
#STRONG Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Alfe_Binders_Ic50_NoExprFilt_ord[(Alfe_Binders_Ic50_NoExprFilt_ord$ALFE2_kallisto >= 0.58 | Alfe_Binders_Ic50_NoExprFilt_ord$ALFE1_kallisto >=0.58) & Alfe_Binders_Ic50_NoExprFilt_ord$value<=50  & (Alfe_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Alfe_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Alfe_Binders_Ic50_NoExprFilt_ord$Rel_Mut > Alfe_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
nrow(Alfe_Binders_Ic50_NoExprFilt_ord[(Alfe_Binders_Ic50_NoExprFilt_ord$ALFE2_kallisto >= 0.58 | Alfe_Binders_Ic50_NoExprFilt_ord$ALFE1_kallisto >=0.58) & Alfe_Binders_Ic50_NoExprFilt_ord$value<=50  & (Alfe_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Alfe_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Alfe_Binders_Ic50_NoExprFilt_ord$Rel_Mut < Alfe_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
#WEAK Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Alfe_Binders_Ic50_NoExprFilt_ord[(Alfe_Binders_Ic50_NoExprFilt_ord$ALFE2_kallisto >= 0.58 | Alfe_Binders_Ic50_NoExprFilt_ord$ALFE1_kallisto >=0.58) & Alfe_Binders_Ic50_NoExprFilt_ord$value>50  & (Alfe_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Alfe_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Alfe_Binders_Ic50_NoExprFilt_ord$Rel_Mut > Alfe_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
nrow(Alfe_Binders_Ic50_NoExprFilt_ord[(Alfe_Binders_Ic50_NoExprFilt_ord$ALFE2_kallisto >= 0.58 | Alfe_Binders_Ic50_NoExprFilt_ord$ALFE1_kallisto >=0.58) & Alfe_Binders_Ic50_NoExprFilt_ord$value>50  & (Alfe_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Alfe_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Alfe_Binders_Ic50_NoExprFilt_ord$Rel_Mut < Alfe_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
#STRONG   Trascript  Expr_Relapse VS Expr_Diagnosis	
nrow(Alfe_Binders_Ic50_NoExprFilt_ord[(Alfe_Binders_Ic50_NoExprFilt_ord$ALFE2_kallisto >= 0.58 | Alfe_Binders_Ic50_NoExprFilt_ord$ALFE1_kallisto >=0.58) & Alfe_Binders_Ic50_NoExprFilt_ord$value<=50  & (Alfe_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Alfe_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Alfe_Binders_Ic50_NoExprFilt_ord$ALFE2_kallisto > Alfe_Binders_Ic50_NoExprFilt_ord$ALFE1_kallisto , ])
nrow(Alfe_Binders_Ic50_NoExprFilt_ord[(Alfe_Binders_Ic50_NoExprFilt_ord$ALFE2_kallisto >= 0.58 | Alfe_Binders_Ic50_NoExprFilt_ord$ALFE1_kallisto >=0.58) & Alfe_Binders_Ic50_NoExprFilt_ord$value<=50  & (Alfe_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Alfe_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Alfe_Binders_Ic50_NoExprFilt_ord$ALFE2_kallisto < Alfe_Binders_Ic50_NoExprFilt_ord$ALFE1_kallisto , ])
#WEAK. Transcript  Expr_Relapse VS Expr_Diagnosis	
nrow(Alfe_Binders_Ic50_NoExprFilt_ord[(Alfe_Binders_Ic50_NoExprFilt_ord$ALFE2_kallisto >= 0.58 | Alfe_Binders_Ic50_NoExprFilt_ord$ALFE1_kallisto >=0.58) & Alfe_Binders_Ic50_NoExprFilt_ord$value>50  & (Alfe_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Alfe_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Alfe_Binders_Ic50_NoExprFilt_ord$ALFE2_kallisto > Alfe_Binders_Ic50_NoExprFilt_ord$ALFE1_kallisto , ])
nrow(Alfe_Binders_Ic50_NoExprFilt_ord[(Alfe_Binders_Ic50_NoExprFilt_ord$ALFE2_kallisto >= 0.58 | Alfe_Binders_Ic50_NoExprFilt_ord$ALFE1_kallisto >=0.58) & Alfe_Binders_Ic50_NoExprFilt_ord$value>50  & (Alfe_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Alfe_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Alfe_Binders_Ic50_NoExprFilt_ord$ALFE2_kallisto < Alfe_Binders_Ic50_NoExprFilt_ord$ALFE1_kallisto , ])


#BESU
#N° MhAGs Strong
nrow(Besu_Binders_Ic50_NoExprFilt_ord[Besu_Binders_Ic50_NoExprFilt_ord$value<=50,])
#N° MhAGs Weak
nrow(Besu_Binders_Ic50_NoExprFilt_ord[Besu_Binders_Ic50_NoExprFilt_ord$value>50,])
#N° Expr. MhAGs  Strong
nrow(Besu_Binders_Ic50_NoExprFilt_ord[(Besu_Binders_Ic50_NoExprFilt_ord$BESU2_kallisto >= 0.58 | Besu_Binders_Ic50_NoExprFilt_ord$BESU1_kallisto >= 0.58) & Besu_Binders_Ic50_NoExprFilt_ord$value<=50 & (Besu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Besu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#N° Expr. MhAGs  Weak
nrow(Besu_Binders_Ic50_NoExprFilt_ord[(Besu_Binders_Ic50_NoExprFilt_ord$BESU2_kallisto >= 0.58 | Besu_Binders_Ic50_NoExprFilt_ord$BESU1_kallisto >= 0.58) & Besu_Binders_Ic50_NoExprFilt_ord$value>50 & (Besu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Besu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele ONLY
nrow(Besu_Binders_Ic50_NoExprFilt_ord[(Besu_Binders_Ic50_NoExprFilt_ord$BESU2_kallisto >= 0.58 | Besu_Binders_Ic50_NoExprFilt_ord$BESU1_kallisto >= 0.58) & Besu_Binders_Ic50_NoExprFilt_ord$value<=50  & Besu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3 & Besu_Binders_Ic50_NoExprFilt_ord$Rel_Mut<3, ])
#STRONG  Expressed REL Allele Only
nrow(Besu_Binders_Ic50_NoExprFilt_ord[(Besu_Binders_Ic50_NoExprFilt_ord$BESU2_kallisto >= 0.58 | Besu_Binders_Ic50_NoExprFilt_ord$BESU1_kallisto >= 0.58) & Besu_Binders_Ic50_NoExprFilt_ord$value<=50  & Besu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 & Besu_Binders_Ic50_NoExprFilt_ord$Dx_Mut<3 , ])
#STRONG  Expressed Common DIAG REL
nrow(Besu_Binders_Ic50_NoExprFilt_ord[(Besu_Binders_Ic50_NoExprFilt_ord$BESU2_kallisto >= 0.58 | Besu_Binders_Ic50_NoExprFilt_ord$BESU1_kallisto >= 0.58) & Besu_Binders_Ic50_NoExprFilt_ord$value<=50 & (Besu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 & Besu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele Diagnosis
nrow(Besu_Binders_Ic50_NoExprFilt_ord[(Besu_Binders_Ic50_NoExprFilt_ord$BESU2_kallisto >= 0.58 | Besu_Binders_Ic50_NoExprFilt_ord$BESU1_kallisto >= 0.58) & Besu_Binders_Ic50_NoExprFilt_ord$value<=50  & Besu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3 , ])
#STRONG  Expressed REL Allele Relapse
nrow(Besu_Binders_Ic50_NoExprFilt_ord[(Besu_Binders_Ic50_NoExprFilt_ord$BESU2_kallisto >= 0.58 | Besu_Binders_Ic50_NoExprFilt_ord$BESU1_kallisto >= 0.58) & Besu_Binders_Ic50_NoExprFilt_ord$value<=50  & Besu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3, ])
#WEAK Expressed DIAG Allele ONLY
nrow(Besu_Binders_Ic50_NoExprFilt_ord[(Besu_Binders_Ic50_NoExprFilt_ord$BESU2_kallisto >= 0.58 | Besu_Binders_Ic50_NoExprFilt_ord$BESU1_kallisto >= 0.58) & Besu_Binders_Ic50_NoExprFilt_ord$value>50  & Besu_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3 & Besu_Binders_Ic50_NoExprFilt_ord$Rel_Mut < 3, ])
#WEAK  Expressed REL Allele Only
nrow(Besu_Binders_Ic50_NoExprFilt_ord[(Besu_Binders_Ic50_NoExprFilt_ord$BESU2_kallisto >= 0.58 | Besu_Binders_Ic50_NoExprFilt_ord$BESU1_kallisto >= 0.58) & Besu_Binders_Ic50_NoExprFilt_ord$value>50  & Besu_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3 & Besu_Binders_Ic50_NoExprFilt_ord$Dx_Mut < 3 , ])
#WEAK  Expressed Common DIAG REL
nrow(Besu_Binders_Ic50_NoExprFilt_ord[(Besu_Binders_Ic50_NoExprFilt_ord$BESU2_kallisto >= 0.58 | Besu_Binders_Ic50_NoExprFilt_ord$BESU1_kallisto >= 0.58) & Besu_Binders_Ic50_NoExprFilt_ord$value>50 & (Besu_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3 & Besu_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3), ])
#WEAK Expressed DIAG Allele Diagnosis
nrow(Besu_Binders_Ic50_NoExprFilt_ord[(Besu_Binders_Ic50_NoExprFilt_ord$BESU2_kallisto >= 0.58 | Besu_Binders_Ic50_NoExprFilt_ord$BESU1_kallisto >= 0.58) & Besu_Binders_Ic50_NoExprFilt_ord$value>50  & Besu_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3 , ])
#WEAK  Expressed REL Allele Relapse
nrow(Besu_Binders_Ic50_NoExprFilt_ord[(Besu_Binders_Ic50_NoExprFilt_ord$BESU2_kallisto >= 0.58 | Besu_Binders_Ic50_NoExprFilt_ord$BESU1_kallisto >= 0.58) & Besu_Binders_Ic50_NoExprFilt_ord$value>50  & Besu_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3, ])
#STRONG Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Besu_Binders_Ic50_NoExprFilt_ord[(Besu_Binders_Ic50_NoExprFilt_ord$BESU2_kallisto >= 0.58 | Besu_Binders_Ic50_NoExprFilt_ord$BESU1_kallisto >= 0.58) & Besu_Binders_Ic50_NoExprFilt_ord$value<=50  & (Besu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Besu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Besu_Binders_Ic50_NoExprFilt_ord$Rel_Mut > Besu_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
nrow(Besu_Binders_Ic50_NoExprFilt_ord[(Besu_Binders_Ic50_NoExprFilt_ord$BESU2_kallisto >= 0.58 | Besu_Binders_Ic50_NoExprFilt_ord$BESU1_kallisto >= 0.58) & Besu_Binders_Ic50_NoExprFilt_ord$value<=50  & (Besu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Besu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Besu_Binders_Ic50_NoExprFilt_ord$Rel_Mut < Besu_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
#WEAK Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Besu_Binders_Ic50_NoExprFilt_ord[(Besu_Binders_Ic50_NoExprFilt_ord$BESU2_kallisto >= 0.58 | Besu_Binders_Ic50_NoExprFilt_ord$BESU1_kallisto >= 0.58) & Besu_Binders_Ic50_NoExprFilt_ord$value>50  & (Besu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Besu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Besu_Binders_Ic50_NoExprFilt_ord$Rel_Mut > Besu_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
nrow(Besu_Binders_Ic50_NoExprFilt_ord[(Besu_Binders_Ic50_NoExprFilt_ord$BESU2_kallisto >= 0.58 | Besu_Binders_Ic50_NoExprFilt_ord$BESU1_kallisto >= 0.58) & Besu_Binders_Ic50_NoExprFilt_ord$value>50  & (Besu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Besu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Besu_Binders_Ic50_NoExprFilt_ord$Rel_Mut < Besu_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
#STRONG   Trascript  Expr_Relapse VS Expr_Diagnosis	
nrow(Besu_Binders_Ic50_NoExprFilt_ord[(Besu_Binders_Ic50_NoExprFilt_ord$BESU2_kallisto >= 0.58 | Besu_Binders_Ic50_NoExprFilt_ord$BESU1_kallisto >= 0.58) & Besu_Binders_Ic50_NoExprFilt_ord$value<=50  & (Besu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Besu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Besu_Binders_Ic50_NoExprFilt_ord$BESU2_kallisto > Besu_Binders_Ic50_NoExprFilt_ord$BESU1_kallisto , ])
nrow(Besu_Binders_Ic50_NoExprFilt_ord[(Besu_Binders_Ic50_NoExprFilt_ord$BESU2_kallisto >= 0.58 | Besu_Binders_Ic50_NoExprFilt_ord$BESU1_kallisto >= 0.58) & Besu_Binders_Ic50_NoExprFilt_ord$value<=50  & (Besu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Besu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Besu_Binders_Ic50_NoExprFilt_ord$BESU2_kallisto < Besu_Binders_Ic50_NoExprFilt_ord$BESU1_kallisto , ])
#WEAK. Transcript  Expr_Relapse VS Expr_Diagnosis	
nrow(Besu_Binders_Ic50_NoExprFilt_ord[(Besu_Binders_Ic50_NoExprFilt_ord$BESU2_kallisto >= 0.58 | Besu_Binders_Ic50_NoExprFilt_ord$BESU1_kallisto >= 0.58) & Besu_Binders_Ic50_NoExprFilt_ord$value>50  & (Besu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Besu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Besu_Binders_Ic50_NoExprFilt_ord$BESU2_kallisto > Besu_Binders_Ic50_NoExprFilt_ord$BESU1_kallisto , ])
nrow(Besu_Binders_Ic50_NoExprFilt_ord[(Besu_Binders_Ic50_NoExprFilt_ord$BESU2_kallisto >= 0.58 | Besu_Binders_Ic50_NoExprFilt_ord$BESU1_kallisto >= 0.58) & Besu_Binders_Ic50_NoExprFilt_ord$value>50  & (Besu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Besu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Besu_Binders_Ic50_NoExprFilt_ord$BESU2_kallisto < Besu_Binders_Ic50_NoExprFilt_ord$BESU1_kallisto , ])


#CALU
#N° MhAGs Strong
nrow(Calu_Binders_Ic50_NoExprFilt_ord_homo[Calu_Binders_Ic50_NoExprFilt_ord_homo$value<=50,])
#N° MhAGs Weak
nrow(Calu_Binders_Ic50_NoExprFilt_ord_homo[Calu_Binders_Ic50_NoExprFilt_ord_homo$value>50,])
#N° Expr. MhAGs  Strong
nrow(Calu_Binders_Ic50_NoExprFilt_ord_homo[ (Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU2_kallisto>=0.58 | Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU1_kallisto>=0.58) & Calu_Binders_Ic50_NoExprFilt_ord_homo$value<=50 & (Calu_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Calu_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3), ])
#N° Expr. MhAGs  Weak
nrow(Calu_Binders_Ic50_NoExprFilt_ord_homo[ (Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU2_kallisto>=0.58 | Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU1_kallisto>=0.58) & Calu_Binders_Ic50_NoExprFilt_ord_homo$value>50 & (Calu_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Calu_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele ONLY
nrow(Calu_Binders_Ic50_NoExprFilt_ord_homo[ (Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU2_kallisto>=0.58 | Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU1_kallisto>=0.58) & Calu_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & Calu_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3 & Calu_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut<3, ])
#STRONG  Expressed REL Allele Only
nrow(Calu_Binders_Ic50_NoExprFilt_ord_homo[ (Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU2_kallisto>=0.58 | Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU1_kallisto>=0.58) & Calu_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & Calu_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 & Calu_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut<3 , ])
#STRONG  Expressed Common DIAG REL
nrow(Calu_Binders_Ic50_NoExprFilt_ord_homo[ (Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU2_kallisto>=0.58 | Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU1_kallisto>=0.58) & Calu_Binders_Ic50_NoExprFilt_ord_homo$value<=50 & (Calu_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 & Calu_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele Diagnosis
nrow(Calu_Binders_Ic50_NoExprFilt_ord_homo[ (Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU2_kallisto>=0.58 | Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU1_kallisto>=0.58) & Calu_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & Calu_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3 , ])
#STRONG  Expressed REL Allele Relapse
nrow(Calu_Binders_Ic50_NoExprFilt_ord_homo[ (Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU2_kallisto>=0.58 | Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU1_kallisto>=0.58) & Calu_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & Calu_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3, ])
#WEAK Expressed DIAG Allele ONLY
nrow(Calu_Binders_Ic50_NoExprFilt_ord_homo[ (Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU2_kallisto>=0.58 | Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU1_kallisto>=0.58) & Calu_Binders_Ic50_NoExprFilt_ord_homo$value>50  & Calu_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut >= 3 & Calu_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut < 3, ])
#WEAK  Expressed REL Allele Only
nrow(Calu_Binders_Ic50_NoExprFilt_ord_homo[ (Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU2_kallisto>=0.58 | Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU1_kallisto>=0.58) & Calu_Binders_Ic50_NoExprFilt_ord_homo$value>50  & Calu_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut >= 3 & Calu_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut < 3 , ])
#WEAK  Expressed Common DIAG REL
nrow(Calu_Binders_Ic50_NoExprFilt_ord_homo[ (Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU2_kallisto>=0.58 | Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU1_kallisto>=0.58) & Calu_Binders_Ic50_NoExprFilt_ord_homo$value>50 & (Calu_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut >= 3 & Calu_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut >= 3), ])
#WEAK Expressed DIAG Allele Diagnosis
nrow(Calu_Binders_Ic50_NoExprFilt_ord_homo[ (Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU2_kallisto>=0.58 | Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU1_kallisto>=0.58) & Calu_Binders_Ic50_NoExprFilt_ord_homo$value>50  & Calu_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut >= 3 , ])
#WEAK  Expressed REL Allele Relapse
nrow(Calu_Binders_Ic50_NoExprFilt_ord_homo[ (Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU2_kallisto>=0.58 | Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU1_kallisto>=0.58) & Calu_Binders_Ic50_NoExprFilt_ord_homo$value>50  & Calu_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut >= 3, ])
#STRONG Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Calu_Binders_Ic50_NoExprFilt_ord_homo[ (Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU2_kallisto>=0.58 | Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU1_kallisto>=0.58) & Calu_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & (Calu_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Calu_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3) & Calu_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut > Calu_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut , ])
nrow(Calu_Binders_Ic50_NoExprFilt_ord_homo[ (Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU2_kallisto>=0.58 | Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU1_kallisto>=0.58) & Calu_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & (Calu_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Calu_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3) & Calu_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut < Calu_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut , ])
#WEAK Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Calu_Binders_Ic50_NoExprFilt_ord_homo[ (Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU2_kallisto>=0.58 | Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU1_kallisto>=0.58) & Calu_Binders_Ic50_NoExprFilt_ord_homo$value>50  & (Calu_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Calu_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3) & Calu_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut > Calu_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut , ])
nrow(Calu_Binders_Ic50_NoExprFilt_ord_homo[ (Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU2_kallisto>=0.58 | Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU1_kallisto>=0.58) & Calu_Binders_Ic50_NoExprFilt_ord_homo$value>50  & (Calu_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Calu_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3) & Calu_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut < Calu_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut , ])
#STRONG   Trascript  Expr_Relapse VS Expr_Diagnosis	
nrow(Calu_Binders_Ic50_NoExprFilt_ord_homo[ (Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU2_kallisto>=0.58 | Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU1_kallisto>=0.58) & Calu_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & (Calu_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Calu_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3)  & Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU2_kallisto > Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU1_kallisto , ])
nrow(Calu_Binders_Ic50_NoExprFilt_ord_homo[ (Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU2_kallisto>=0.58 | Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU1_kallisto>=0.58) & Calu_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & (Calu_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Calu_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3)  & Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU2_kallisto < Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU1_kallisto , ])
#WEAK. Transcript  Expr_Relapse VS Expr_Diagnosis	
nrow(Calu_Binders_Ic50_NoExprFilt_ord_homo[ (Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU2_kallisto>=0.58 | Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU1_kallisto>=0.58) & Calu_Binders_Ic50_NoExprFilt_ord_homo$value>50  & (Calu_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Calu_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3)  & Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU2_kallisto > Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU1_kallisto , ])
nrow(Calu_Binders_Ic50_NoExprFilt_ord_homo[ (Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU2_kallisto>=0.58 | Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU1_kallisto>=0.58) & Calu_Binders_Ic50_NoExprFilt_ord_homo$value>50  & (Calu_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Calu_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3)  & Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU2_kallisto < Calu_Binders_Ic50_NoExprFilt_ord_homo$CALU1_kallisto , ])



#DEIV
#N° MhAGs Strong
nrow(Deiv_Binders_Ic50_NoExprFilt_ord[Deiv_Binders_Ic50_NoExprFilt_ord$value<=50,])
#N° MhAGs Weak
nrow(Deiv_Binders_Ic50_NoExprFilt_ord[Deiv_Binders_Ic50_NoExprFilt_ord$value>50,])
#N° Expr. MhAGs  Strong
nrow(Deiv_Binders_Ic50_NoExprFilt_ord[(Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2b_kallisto>=0.58| Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2_kallisto>=0.58) & Deiv_Binders_Ic50_NoExprFilt_ord$value<=50 & (Deiv_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Deiv_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#N° Expr. MhAGs  Weak
nrow(Deiv_Binders_Ic50_NoExprFilt_ord[(Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2b_kallisto>=0.58| Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2_kallisto>=0.58) & Deiv_Binders_Ic50_NoExprFilt_ord$value>50 & (Deiv_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Deiv_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele ONLY
nrow(Deiv_Binders_Ic50_NoExprFilt_ord[(Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2b_kallisto>=0.58| Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2_kallisto>=0.58) & Deiv_Binders_Ic50_NoExprFilt_ord$value<=50  & Deiv_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3 & Deiv_Binders_Ic50_NoExprFilt_ord$Rel_Mut<3, ])
#STRONG  Expressed REL Allele Only
nrow(Deiv_Binders_Ic50_NoExprFilt_ord[(Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2b_kallisto>=0.58| Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2_kallisto>=0.58) & Deiv_Binders_Ic50_NoExprFilt_ord$value<=50  & Deiv_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 & Deiv_Binders_Ic50_NoExprFilt_ord$Dx_Mut<3 , ])
#STRONG  Expressed Common DIAG REL
nrow(Deiv_Binders_Ic50_NoExprFilt_ord[(Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2b_kallisto>=0.58| Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2_kallisto>=0.58) & Deiv_Binders_Ic50_NoExprFilt_ord$value<=50 & (Deiv_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 & Deiv_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele Diagnosis
nrow(Deiv_Binders_Ic50_NoExprFilt_ord[(Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2b_kallisto>=0.58| Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2_kallisto>=0.58) & Deiv_Binders_Ic50_NoExprFilt_ord$value<=50  & Deiv_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3 , ])
#STRONG  Expressed REL Allele Relapse
nrow(Deiv_Binders_Ic50_NoExprFilt_ord[(Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2b_kallisto>=0.58| Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2_kallisto>=0.58) & Deiv_Binders_Ic50_NoExprFilt_ord$value<=50  & Deiv_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3, ])
#WEAK Expressed DIAG Allele ONLY
nrow(Deiv_Binders_Ic50_NoExprFilt_ord[(Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2b_kallisto>=0.58| Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2_kallisto>=0.58) & Deiv_Binders_Ic50_NoExprFilt_ord$value>50  & Deiv_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3 & Deiv_Binders_Ic50_NoExprFilt_ord$Rel_Mut < 3, ])
#WEAK  Expressed REL Allele Only
nrow(Deiv_Binders_Ic50_NoExprFilt_ord[(Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2b_kallisto>=0.58| Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2_kallisto>=0.58) & Deiv_Binders_Ic50_NoExprFilt_ord$value>50  & Deiv_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3 & Deiv_Binders_Ic50_NoExprFilt_ord$Dx_Mut < 3 , ])
#WEAK  Expressed Common DIAG REL
nrow(Deiv_Binders_Ic50_NoExprFilt_ord[(Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2b_kallisto>=0.58| Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2_kallisto>=0.58) & Deiv_Binders_Ic50_NoExprFilt_ord$value>50 & (Deiv_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3 & Deiv_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3), ])
#WEAK Expressed DIAG Allele Diagnosis
nrow(Deiv_Binders_Ic50_NoExprFilt_ord[(Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2b_kallisto>=0.58| Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2_kallisto>=0.58) & Deiv_Binders_Ic50_NoExprFilt_ord$value>50  & Deiv_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3 , ])
#WEAK  Expressed REL Allele Relapse
nrow(Deiv_Binders_Ic50_NoExprFilt_ord[(Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2b_kallisto>=0.58| Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2_kallisto>=0.58) & Deiv_Binders_Ic50_NoExprFilt_ord$value>50  & Deiv_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3, ])
#STRONG Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Deiv_Binders_Ic50_NoExprFilt_ord[(Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2b_kallisto>=0.58| Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2_kallisto>=0.58) & Deiv_Binders_Ic50_NoExprFilt_ord$value<=50  & (Deiv_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Deiv_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Deiv_Binders_Ic50_NoExprFilt_ord$Rel_Mut > Deiv_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
nrow(Deiv_Binders_Ic50_NoExprFilt_ord[(Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2b_kallisto>=0.58| Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2_kallisto>=0.58) & Deiv_Binders_Ic50_NoExprFilt_ord$value<=50  & (Deiv_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Deiv_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Deiv_Binders_Ic50_NoExprFilt_ord$Rel_Mut < Deiv_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
#WEAK Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Deiv_Binders_Ic50_NoExprFilt_ord[(Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2b_kallisto>=0.58| Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2_kallisto>=0.58) & Deiv_Binders_Ic50_NoExprFilt_ord$value>50  & (Deiv_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Deiv_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Deiv_Binders_Ic50_NoExprFilt_ord$Rel_Mut > Deiv_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
nrow(Deiv_Binders_Ic50_NoExprFilt_ord[(Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2b_kallisto>=0.58| Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2_kallisto>=0.58) & Deiv_Binders_Ic50_NoExprFilt_ord$value>50  & (Deiv_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Deiv_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Deiv_Binders_Ic50_NoExprFilt_ord$Rel_Mut < Deiv_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
#STRONG   Trascript  Expr_Relapse VS Expr_Diagnosis	
nrow(Deiv_Binders_Ic50_NoExprFilt_ord[(Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2b_kallisto>=0.58| Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2_kallisto>=0.58) & Deiv_Binders_Ic50_NoExprFilt_ord$value<=50  & (Deiv_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Deiv_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2b_kallisto  > Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2_kallisto , ])
nrow(Deiv_Binders_Ic50_NoExprFilt_ord[(Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2b_kallisto>=0.58| Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2_kallisto>=0.58) & Deiv_Binders_Ic50_NoExprFilt_ord$value<=50  & (Deiv_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Deiv_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2b_kallisto < Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2_kallisto , ])
#WEAK. Transcript  Expr_Relapse VS Expr_Diagnosis	
nrow(Deiv_Binders_Ic50_NoExprFilt_ord[(Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2b_kallisto>=0.58| Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2_kallisto>=0.58) & Deiv_Binders_Ic50_NoExprFilt_ord$value>50  & (Deiv_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Deiv_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2b_kallisto > Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2_kallisto , ])
nrow(Deiv_Binders_Ic50_NoExprFilt_ord[(Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2b_kallisto>=0.58| Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2_kallisto>=0.58) & Deiv_Binders_Ic50_NoExprFilt_ord$value>50  & (Deiv_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Deiv_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2b_kallisto < Deiv_Binders_Ic50_NoExprFilt_ord$DEIV2_kallisto , ])





#DEST
#N° MhAGs Strong
nrow(Dest_Binders_Ic50_NoExprFilt_ord[Dest_Binders_Ic50_NoExprFilt_ord$value<=50,])
#N° MhAGs Weak
nrow(Dest_Binders_Ic50_NoExprFilt_ord[Dest_Binders_Ic50_NoExprFilt_ord$value>50,])
#N° Expr. MhAGs  Strong
nrow(Dest_Binders_Ic50_NoExprFilt_ord[(Dest_Binders_Ic50_NoExprFilt_ord$DEST1_kallisto>=0.58  | Dest_Binders_Ic50_NoExprFilt_ord$DEST2_kallisto>=0.58) & Dest_Binders_Ic50_NoExprFilt_ord$value<=50 & (Dest_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dest_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#N° Expr. MhAGs  Weak
nrow(Dest_Binders_Ic50_NoExprFilt_ord[(Dest_Binders_Ic50_NoExprFilt_ord$DEST1_kallisto>=0.58  | Dest_Binders_Ic50_NoExprFilt_ord$DEST2_kallisto>=0.58) & Dest_Binders_Ic50_NoExprFilt_ord$value>50 & (Dest_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dest_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele ONLY
nrow(Dest_Binders_Ic50_NoExprFilt_ord[(Dest_Binders_Ic50_NoExprFilt_ord$DEST1_kallisto>=0.58  | Dest_Binders_Ic50_NoExprFilt_ord$DEST2_kallisto>=0.58) & Dest_Binders_Ic50_NoExprFilt_ord$value<=50  & Dest_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3 & Dest_Binders_Ic50_NoExprFilt_ord$Rel_Mut<3, ])
#STRONG  Expressed REL Allele Only
nrow(Dest_Binders_Ic50_NoExprFilt_ord[(Dest_Binders_Ic50_NoExprFilt_ord$DEST1_kallisto>=0.58  | Dest_Binders_Ic50_NoExprFilt_ord$DEST2_kallisto>=0.58) & Dest_Binders_Ic50_NoExprFilt_ord$value<=50  & Dest_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 & Dest_Binders_Ic50_NoExprFilt_ord$Dx_Mut<3 , ])
#STRONG  Expressed Common DIAG REL
nrow(Dest_Binders_Ic50_NoExprFilt_ord[(Dest_Binders_Ic50_NoExprFilt_ord$DEST1_kallisto>=0.58  | Dest_Binders_Ic50_NoExprFilt_ord$DEST2_kallisto>=0.58) & Dest_Binders_Ic50_NoExprFilt_ord$value<=50 & (Dest_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 & Dest_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele Diagnosis
nrow(Dest_Binders_Ic50_NoExprFilt_ord[(Dest_Binders_Ic50_NoExprFilt_ord$DEST1_kallisto>=0.58  | Dest_Binders_Ic50_NoExprFilt_ord$DEST2_kallisto>=0.58) & Dest_Binders_Ic50_NoExprFilt_ord$value<=50  & Dest_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3 , ])
#STRONG  Expressed REL Allele Relapse
nrow(Dest_Binders_Ic50_NoExprFilt_ord[(Dest_Binders_Ic50_NoExprFilt_ord$DEST1_kallisto>=0.58  | Dest_Binders_Ic50_NoExprFilt_ord$DEST2_kallisto>=0.58) & Dest_Binders_Ic50_NoExprFilt_ord$value<=50  & Dest_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3, ])
#WEAK Expressed DIAG Allele ONLY
nrow(Dest_Binders_Ic50_NoExprFilt_ord[(Dest_Binders_Ic50_NoExprFilt_ord$DEST1_kallisto>=0.58  | Dest_Binders_Ic50_NoExprFilt_ord$DEST2_kallisto>=0.58) & Dest_Binders_Ic50_NoExprFilt_ord$value>50  & Dest_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3 & Dest_Binders_Ic50_NoExprFilt_ord$Rel_Mut < 3, ])
#WEAK  Expressed REL Allele Only
nrow(Dest_Binders_Ic50_NoExprFilt_ord[(Dest_Binders_Ic50_NoExprFilt_ord$DEST1_kallisto>=0.58  | Dest_Binders_Ic50_NoExprFilt_ord$DEST2_kallisto>=0.58) & Dest_Binders_Ic50_NoExprFilt_ord$value>50  & Dest_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3 & Dest_Binders_Ic50_NoExprFilt_ord$Dx_Mut < 3 , ])
#WEAK  Expressed Common DIAG REL
nrow(Dest_Binders_Ic50_NoExprFilt_ord[(Dest_Binders_Ic50_NoExprFilt_ord$DEST1_kallisto>=0.58  | Dest_Binders_Ic50_NoExprFilt_ord$DEST2_kallisto>=0.58) & Dest_Binders_Ic50_NoExprFilt_ord$value>50 & (Dest_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3 & Dest_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3), ])
#WEAK Expressed DIAG Allele Diagnosis
nrow(Dest_Binders_Ic50_NoExprFilt_ord[(Dest_Binders_Ic50_NoExprFilt_ord$DEST1_kallisto>=0.58  | Dest_Binders_Ic50_NoExprFilt_ord$DEST2_kallisto>=0.58) & Dest_Binders_Ic50_NoExprFilt_ord$value>50  & Dest_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3 , ])
#WEAK  Expressed REL Allele Relapse
nrow(Dest_Binders_Ic50_NoExprFilt_ord[(Dest_Binders_Ic50_NoExprFilt_ord$DEST1_kallisto>=0.58  | Dest_Binders_Ic50_NoExprFilt_ord$DEST2_kallisto>=0.58) & Dest_Binders_Ic50_NoExprFilt_ord$value>50  & Dest_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3, ])
#STRONG Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Dest_Binders_Ic50_NoExprFilt_ord[(Dest_Binders_Ic50_NoExprFilt_ord$DEST1_kallisto>=0.58  | Dest_Binders_Ic50_NoExprFilt_ord$DEST2_kallisto>=0.58) & Dest_Binders_Ic50_NoExprFilt_ord$value<=50  & (Dest_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dest_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Dest_Binders_Ic50_NoExprFilt_ord$Rel_Mut > Dest_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
nrow(Dest_Binders_Ic50_NoExprFilt_ord[(Dest_Binders_Ic50_NoExprFilt_ord$DEST1_kallisto>=0.58  | Dest_Binders_Ic50_NoExprFilt_ord$DEST2_kallisto>=0.58) & Dest_Binders_Ic50_NoExprFilt_ord$value<=50  & (Dest_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dest_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Dest_Binders_Ic50_NoExprFilt_ord$Rel_Mut < Dest_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
#WEAK Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Dest_Binders_Ic50_NoExprFilt_ord[(Dest_Binders_Ic50_NoExprFilt_ord$DEST1_kallisto>=0.58  | Dest_Binders_Ic50_NoExprFilt_ord$DEST2_kallisto>=0.58) & Dest_Binders_Ic50_NoExprFilt_ord$value>50  & (Dest_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dest_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Dest_Binders_Ic50_NoExprFilt_ord$Rel_Mut > Dest_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
nrow(Dest_Binders_Ic50_NoExprFilt_ord[(Dest_Binders_Ic50_NoExprFilt_ord$DEST1_kallisto>=0.58  | Dest_Binders_Ic50_NoExprFilt_ord$DEST2_kallisto>=0.58) & Dest_Binders_Ic50_NoExprFilt_ord$value>50  & (Dest_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dest_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Dest_Binders_Ic50_NoExprFilt_ord$Rel_Mut < Dest_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
#STRONG   Trascript  Expr_Relapse VS Expr_Diagnosis	
nrow(Dest_Binders_Ic50_NoExprFilt_ord[(Dest_Binders_Ic50_NoExprFilt_ord$DEST1_kallisto>=0.58  | Dest_Binders_Ic50_NoExprFilt_ord$DEST2_kallisto>=0.58) & Dest_Binders_Ic50_NoExprFilt_ord$value<=50  & (Dest_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dest_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Dest_Binders_Ic50_NoExprFilt_ord$DEST2_kallisto  > Dest_Binders_Ic50_NoExprFilt_ord$DEST1_kallisto , ])
nrow(Dest_Binders_Ic50_NoExprFilt_ord[(Dest_Binders_Ic50_NoExprFilt_ord$DEST1_kallisto>=0.58  | Dest_Binders_Ic50_NoExprFilt_ord$DEST2_kallisto>=0.58) & Dest_Binders_Ic50_NoExprFilt_ord$value<=50  & (Dest_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dest_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Dest_Binders_Ic50_NoExprFilt_ord$DEST2_kallisto < Dest_Binders_Ic50_NoExprFilt_ord$DEST1_kallisto , ])
#WEAK. Transcript  Expr_Relapse VS Expr_Diagnosis	
nrow(Dest_Binders_Ic50_NoExprFilt_ord[(Dest_Binders_Ic50_NoExprFilt_ord$DEST1_kallisto>=0.58  | Dest_Binders_Ic50_NoExprFilt_ord$DEST2_kallisto>=0.58) & Dest_Binders_Ic50_NoExprFilt_ord$value>50  & (Dest_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dest_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Dest_Binders_Ic50_NoExprFilt_ord$DEST2_kallisto > Dest_Binders_Ic50_NoExprFilt_ord$DEST1_kallisto , ])
nrow(Dest_Binders_Ic50_NoExprFilt_ord[(Dest_Binders_Ic50_NoExprFilt_ord$DEST1_kallisto>=0.58  | Dest_Binders_Ic50_NoExprFilt_ord$DEST2_kallisto>=0.58) & Dest_Binders_Ic50_NoExprFilt_ord$value>50  & (Dest_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dest_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Dest_Binders_Ic50_NoExprFilt_ord$DEST2_kallisto < Dest_Binders_Ic50_NoExprFilt_ord$DEST1_kallisto , ])




#DR1
#N° MhAGs Strong
nrow(Dr1_Binders_Ic50_NoExprFilt_ord[Dr1_Binders_Ic50_NoExprFilt_ord$value<=50,])
#N° MhAGs Weak
nrow(Dr1_Binders_Ic50_NoExprFilt_ord[Dr1_Binders_Ic50_NoExprFilt_ord$value>50,])
#N° Expr. MhAGs  Strong
nrow(Dr1_Binders_Ic50_NoExprFilt_ord[(Dr1_Binders_Ic50_NoExprFilt_ord$DR1_1_kallisto>=0.58 | Dr1_Binders_Ic50_NoExprFilt_ord$DR1_2_kallisto>=0.58) & Dr1_Binders_Ic50_NoExprFilt_ord$value<=50 & (Dr1_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr1_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#N° Expr. MhAGs  Weak
nrow(Dr1_Binders_Ic50_NoExprFilt_ord[(Dr1_Binders_Ic50_NoExprFilt_ord$DR1_1_kallisto>=0.58 | Dr1_Binders_Ic50_NoExprFilt_ord$DR1_2_kallisto>=0.58) & Dr1_Binders_Ic50_NoExprFilt_ord$value>50 & (Dr1_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr1_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele ONLY
nrow(Dr1_Binders_Ic50_NoExprFilt_ord[(Dr1_Binders_Ic50_NoExprFilt_ord$DR1_1_kallisto>=0.58 | Dr1_Binders_Ic50_NoExprFilt_ord$DR1_2_kallisto>=0.58) & Dr1_Binders_Ic50_NoExprFilt_ord$value<=50  & Dr1_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3 & Dr1_Binders_Ic50_NoExprFilt_ord$Rel_Mut<3, ])
#STRONG  Expressed REL Allele Only
nrow(Dr1_Binders_Ic50_NoExprFilt_ord[(Dr1_Binders_Ic50_NoExprFilt_ord$DR1_1_kallisto>=0.58 | Dr1_Binders_Ic50_NoExprFilt_ord$DR1_2_kallisto>=0.58) & Dr1_Binders_Ic50_NoExprFilt_ord$value<=50  & Dr1_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 & Dr1_Binders_Ic50_NoExprFilt_ord$Dx_Mut<3 , ])
#STRONG  Expressed Common DIAG REL
nrow(Dr1_Binders_Ic50_NoExprFilt_ord[(Dr1_Binders_Ic50_NoExprFilt_ord$DR1_1_kallisto>=0.58 | Dr1_Binders_Ic50_NoExprFilt_ord$DR1_2_kallisto>=0.58) & Dr1_Binders_Ic50_NoExprFilt_ord$value<=50 & (Dr1_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 & Dr1_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele Diagnosis
nrow(Dr1_Binders_Ic50_NoExprFilt_ord[(Dr1_Binders_Ic50_NoExprFilt_ord$DR1_1_kallisto>=0.58 | Dr1_Binders_Ic50_NoExprFilt_ord$DR1_2_kallisto>=0.58) & Dr1_Binders_Ic50_NoExprFilt_ord$value<=50  & Dr1_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3 , ])
#STRONG  Expressed REL Allele Relapse
nrow(Dr1_Binders_Ic50_NoExprFilt_ord[(Dr1_Binders_Ic50_NoExprFilt_ord$DR1_1_kallisto>=0.58 | Dr1_Binders_Ic50_NoExprFilt_ord$DR1_2_kallisto>=0.58) & Dr1_Binders_Ic50_NoExprFilt_ord$value<=50  & Dr1_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3, ])
#WEAK Expressed DIAG Allele ONLY
nrow(Dr1_Binders_Ic50_NoExprFilt_ord[(Dr1_Binders_Ic50_NoExprFilt_ord$DR1_1_kallisto>=0.58 | Dr1_Binders_Ic50_NoExprFilt_ord$DR1_2_kallisto>=0.58) & Dr1_Binders_Ic50_NoExprFilt_ord$value>50  & Dr1_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3 & Dr1_Binders_Ic50_NoExprFilt_ord$Rel_Mut < 3, ])
#WEAK  Expressed REL Allele Only
nrow(Dr1_Binders_Ic50_NoExprFilt_ord[(Dr1_Binders_Ic50_NoExprFilt_ord$DR1_1_kallisto>=0.58 | Dr1_Binders_Ic50_NoExprFilt_ord$DR1_2_kallisto>=0.58) & Dr1_Binders_Ic50_NoExprFilt_ord$value>50  & Dr1_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3 & Dr1_Binders_Ic50_NoExprFilt_ord$Dx_Mut < 3 , ])
#WEAK  Expressed Common DIAG REL
nrow(Dr1_Binders_Ic50_NoExprFilt_ord[(Dr1_Binders_Ic50_NoExprFilt_ord$DR1_1_kallisto>=0.58 | Dr1_Binders_Ic50_NoExprFilt_ord$DR1_2_kallisto>=0.58) & Dr1_Binders_Ic50_NoExprFilt_ord$value>50 & (Dr1_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3 & Dr1_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3), ])
#WEAK Expressed DIAG Allele Diagnosis
nrow(Dr1_Binders_Ic50_NoExprFilt_ord[(Dr1_Binders_Ic50_NoExprFilt_ord$DR1_1_kallisto>=0.58 | Dr1_Binders_Ic50_NoExprFilt_ord$DR1_2_kallisto>=0.58) & Dr1_Binders_Ic50_NoExprFilt_ord$value>50  & Dr1_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3 , ])
#WEAK  Expressed REL Allele Relapse
nrow(Dr1_Binders_Ic50_NoExprFilt_ord[(Dr1_Binders_Ic50_NoExprFilt_ord$DR1_1_kallisto>=0.58 | Dr1_Binders_Ic50_NoExprFilt_ord$DR1_2_kallisto>=0.58) & Dr1_Binders_Ic50_NoExprFilt_ord$value>50  & Dr1_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3, ])
#STRONG Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Dr1_Binders_Ic50_NoExprFilt_ord[(Dr1_Binders_Ic50_NoExprFilt_ord$DR1_1_kallisto>=0.58 | Dr1_Binders_Ic50_NoExprFilt_ord$DR1_2_kallisto>=0.58) & Dr1_Binders_Ic50_NoExprFilt_ord$value<=50  & (Dr1_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr1_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Dr1_Binders_Ic50_NoExprFilt_ord$Rel_Mut > Dr1_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
nrow(Dr1_Binders_Ic50_NoExprFilt_ord[(Dr1_Binders_Ic50_NoExprFilt_ord$DR1_1_kallisto>=0.58 | Dr1_Binders_Ic50_NoExprFilt_ord$DR1_2_kallisto>=0.58) & Dr1_Binders_Ic50_NoExprFilt_ord$value<=50  & (Dr1_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr1_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Dr1_Binders_Ic50_NoExprFilt_ord$Rel_Mut < Dr1_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
#WEAK Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Dr1_Binders_Ic50_NoExprFilt_ord[(Dr1_Binders_Ic50_NoExprFilt_ord$DR1_1_kallisto>=0.58 | Dr1_Binders_Ic50_NoExprFilt_ord$DR1_2_kallisto>=0.58) & Dr1_Binders_Ic50_NoExprFilt_ord$value>50  & (Dr1_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr1_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Dr1_Binders_Ic50_NoExprFilt_ord$Rel_Mut > Dr1_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
nrow(Dr1_Binders_Ic50_NoExprFilt_ord[(Dr1_Binders_Ic50_NoExprFilt_ord$DR1_1_kallisto>=0.58 | Dr1_Binders_Ic50_NoExprFilt_ord$DR1_2_kallisto>=0.58) & Dr1_Binders_Ic50_NoExprFilt_ord$value>50  & (Dr1_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr1_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Dr1_Binders_Ic50_NoExprFilt_ord$Rel_Mut < Dr1_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
#STRONG   Trascript  Expr_Relapse VS Expr_Diagnosis	
nrow(Dr1_Binders_Ic50_NoExprFilt_ord[(Dr1_Binders_Ic50_NoExprFilt_ord$DR1_1_kallisto>=0.58 | Dr1_Binders_Ic50_NoExprFilt_ord$DR1_2_kallisto>=0.58) & Dr1_Binders_Ic50_NoExprFilt_ord$value<=50  & (Dr1_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr1_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Dr1_Binders_Ic50_NoExprFilt_ord$DR1_2_kallisto  > Dr1_Binders_Ic50_NoExprFilt_ord$DR1_1_kallisto , ])
nrow(Dr1_Binders_Ic50_NoExprFilt_ord[(Dr1_Binders_Ic50_NoExprFilt_ord$DR1_1_kallisto>=0.58 | Dr1_Binders_Ic50_NoExprFilt_ord$DR1_2_kallisto>=0.58) & Dr1_Binders_Ic50_NoExprFilt_ord$value<=50  & (Dr1_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr1_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Dr1_Binders_Ic50_NoExprFilt_ord$DR1_2_kallisto < Dr1_Binders_Ic50_NoExprFilt_ord$DR1_1_kallisto , ])
#WEAK. Transcript  Expr_Relapse VS Expr_Diagnosis	
nrow(Dr1_Binders_Ic50_NoExprFilt_ord[(Dr1_Binders_Ic50_NoExprFilt_ord$DR1_1_kallisto>=0.58 | Dr1_Binders_Ic50_NoExprFilt_ord$DR1_2_kallisto>=0.58) & Dr1_Binders_Ic50_NoExprFilt_ord$value>50  & (Dr1_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr1_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Dr1_Binders_Ic50_NoExprFilt_ord$DR1_2_kallisto > Dr1_Binders_Ic50_NoExprFilt_ord$DR1_1_kallisto , ])
nrow(Dr1_Binders_Ic50_NoExprFilt_ord[(Dr1_Binders_Ic50_NoExprFilt_ord$DR1_1_kallisto>=0.58 | Dr1_Binders_Ic50_NoExprFilt_ord$DR1_2_kallisto>=0.58) & Dr1_Binders_Ic50_NoExprFilt_ord$value>50  & (Dr1_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr1_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Dr1_Binders_Ic50_NoExprFilt_ord$DR1_2_kallisto < Dr1_Binders_Ic50_NoExprFilt_ord$DR1_1_kallisto , ])




#DR4
#N° MhAGs Strong
nrow(Dr4_Binders_Ic50_NoExprFilt_ord[Dr4_Binders_Ic50_NoExprFilt_ord$value<=50,])
#N° MhAGs Weak
nrow(Dr4_Binders_Ic50_NoExprFilt_ord[Dr4_Binders_Ic50_NoExprFilt_ord$value>50,])
#N° Expr. MhAGs  Strong
nrow(Dr4_Binders_Ic50_NoExprFilt_ord[(Dr4_Binders_Ic50_NoExprFilt_ord$DR42_kallisto>=0.58 | Dr4_Binders_Ic50_NoExprFilt_ord$DR41_kallisto>=0.58) & Dr4_Binders_Ic50_NoExprFilt_ord$value<=50 & (Dr4_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr4_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#N° Expr. MhAGs  Weak
nrow(Dr4_Binders_Ic50_NoExprFilt_ord[(Dr4_Binders_Ic50_NoExprFilt_ord$DR42_kallisto>=0.58 | Dr4_Binders_Ic50_NoExprFilt_ord$DR41_kallisto>=0.58) & Dr4_Binders_Ic50_NoExprFilt_ord$value>50 & (Dr4_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr4_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele ONLY
nrow(Dr4_Binders_Ic50_NoExprFilt_ord[(Dr4_Binders_Ic50_NoExprFilt_ord$DR42_kallisto>=0.58 | Dr4_Binders_Ic50_NoExprFilt_ord$DR41_kallisto>=0.58) & Dr4_Binders_Ic50_NoExprFilt_ord$value<=50  & Dr4_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3 & Dr4_Binders_Ic50_NoExprFilt_ord$Rel_Mut<3, ])
#STRONG  Expressed REL Allele Only
nrow(Dr4_Binders_Ic50_NoExprFilt_ord[(Dr4_Binders_Ic50_NoExprFilt_ord$DR42_kallisto>=0.58 | Dr4_Binders_Ic50_NoExprFilt_ord$DR41_kallisto>=0.58) & Dr4_Binders_Ic50_NoExprFilt_ord$value<=50  & Dr4_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 & Dr4_Binders_Ic50_NoExprFilt_ord$Dx_Mut<3 , ])
#STRONG  Expressed Common DIAG REL
nrow(Dr4_Binders_Ic50_NoExprFilt_ord[(Dr4_Binders_Ic50_NoExprFilt_ord$DR42_kallisto>=0.58 | Dr4_Binders_Ic50_NoExprFilt_ord$DR41_kallisto>=0.58) & Dr4_Binders_Ic50_NoExprFilt_ord$value<=50 & (Dr4_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 & Dr4_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele Diagnosis
nrow(Dr4_Binders_Ic50_NoExprFilt_ord[(Dr4_Binders_Ic50_NoExprFilt_ord$DR42_kallisto>=0.58 | Dr4_Binders_Ic50_NoExprFilt_ord$DR41_kallisto>=0.58) & Dr4_Binders_Ic50_NoExprFilt_ord$value<=50  & Dr4_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3 , ])
#STRONG  Expressed REL Allele Relapse
nrow(Dr4_Binders_Ic50_NoExprFilt_ord[(Dr4_Binders_Ic50_NoExprFilt_ord$DR42_kallisto>=0.58 | Dr4_Binders_Ic50_NoExprFilt_ord$DR41_kallisto>=0.58) & Dr4_Binders_Ic50_NoExprFilt_ord$value<=50  & Dr4_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3, ])
#WEAK Expressed DIAG Allele ONLY
nrow(Dr4_Binders_Ic50_NoExprFilt_ord[(Dr4_Binders_Ic50_NoExprFilt_ord$DR42_kallisto>=0.58 | Dr4_Binders_Ic50_NoExprFilt_ord$DR41_kallisto>=0.58) & Dr4_Binders_Ic50_NoExprFilt_ord$value>50  & Dr4_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3 & Dr4_Binders_Ic50_NoExprFilt_ord$Rel_Mut < 3, ])
#WEAK  Expressed REL Allele Only
nrow(Dr4_Binders_Ic50_NoExprFilt_ord[(Dr4_Binders_Ic50_NoExprFilt_ord$DR42_kallisto>=0.58 | Dr4_Binders_Ic50_NoExprFilt_ord$DR41_kallisto>=0.58) & Dr4_Binders_Ic50_NoExprFilt_ord$value>50  & Dr4_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3 & Dr4_Binders_Ic50_NoExprFilt_ord$Dx_Mut < 3 , ])
#WEAK  Expressed Common DIAG REL
nrow(Dr4_Binders_Ic50_NoExprFilt_ord[(Dr4_Binders_Ic50_NoExprFilt_ord$DR42_kallisto>=0.58 | Dr4_Binders_Ic50_NoExprFilt_ord$DR41_kallisto>=0.58) & Dr4_Binders_Ic50_NoExprFilt_ord$value>50 & (Dr4_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3 & Dr4_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3), ])
#WEAK Expressed DIAG Allele Diagnosis
nrow(Dr4_Binders_Ic50_NoExprFilt_ord[(Dr4_Binders_Ic50_NoExprFilt_ord$DR42_kallisto>=0.58 | Dr4_Binders_Ic50_NoExprFilt_ord$DR41_kallisto>=0.58) & Dr4_Binders_Ic50_NoExprFilt_ord$value>50  & Dr4_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3 , ])
#WEAK  Expressed REL Allele Relapse
nrow(Dr4_Binders_Ic50_NoExprFilt_ord[(Dr4_Binders_Ic50_NoExprFilt_ord$DR42_kallisto>=0.58 | Dr4_Binders_Ic50_NoExprFilt_ord$DR41_kallisto>=0.58) & Dr4_Binders_Ic50_NoExprFilt_ord$value>50  & Dr4_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3, ])
#STRONG Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Dr4_Binders_Ic50_NoExprFilt_ord[(Dr4_Binders_Ic50_NoExprFilt_ord$DR42_kallisto>=0.58 | Dr4_Binders_Ic50_NoExprFilt_ord$DR41_kallisto>=0.58) & Dr4_Binders_Ic50_NoExprFilt_ord$value<=50  & (Dr4_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr4_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Dr4_Binders_Ic50_NoExprFilt_ord$Rel_Mut > Dr4_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
nrow(Dr4_Binders_Ic50_NoExprFilt_ord[(Dr4_Binders_Ic50_NoExprFilt_ord$DR42_kallisto>=0.58 | Dr4_Binders_Ic50_NoExprFilt_ord$DR41_kallisto>=0.58) & Dr4_Binders_Ic50_NoExprFilt_ord$value<=50  & (Dr4_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr4_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Dr4_Binders_Ic50_NoExprFilt_ord$Rel_Mut < Dr4_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
#WEAK Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Dr4_Binders_Ic50_NoExprFilt_ord[(Dr4_Binders_Ic50_NoExprFilt_ord$DR42_kallisto>=0.58 | Dr4_Binders_Ic50_NoExprFilt_ord$DR41_kallisto>=0.58) & Dr4_Binders_Ic50_NoExprFilt_ord$value>50  & (Dr4_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr4_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Dr4_Binders_Ic50_NoExprFilt_ord$Rel_Mut > Dr4_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
nrow(Dr4_Binders_Ic50_NoExprFilt_ord[(Dr4_Binders_Ic50_NoExprFilt_ord$DR42_kallisto>=0.58 | Dr4_Binders_Ic50_NoExprFilt_ord$DR41_kallisto>=0.58) & Dr4_Binders_Ic50_NoExprFilt_ord$value>50  & (Dr4_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr4_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Dr4_Binders_Ic50_NoExprFilt_ord$Rel_Mut < Dr4_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
#STRONG   Trascript  Expr_Relapse VS Expr_Diagnosis	
nrow(Dr4_Binders_Ic50_NoExprFilt_ord[(Dr4_Binders_Ic50_NoExprFilt_ord$DR42_kallisto>=0.58 | Dr4_Binders_Ic50_NoExprFilt_ord$DR41_kallisto>=0.58) & Dr4_Binders_Ic50_NoExprFilt_ord$value<=50  & (Dr4_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr4_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Dr4_Binders_Ic50_NoExprFilt_ord$DR42_kallisto  > Dr4_Binders_Ic50_NoExprFilt_ord$DR41_kallisto , ])
nrow(Dr4_Binders_Ic50_NoExprFilt_ord[(Dr4_Binders_Ic50_NoExprFilt_ord$DR42_kallisto>=0.58 | Dr4_Binders_Ic50_NoExprFilt_ord$DR41_kallisto>=0.58) & Dr4_Binders_Ic50_NoExprFilt_ord$value<=50  & (Dr4_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr4_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Dr4_Binders_Ic50_NoExprFilt_ord$DR42_kallisto < Dr4_Binders_Ic50_NoExprFilt_ord$DR41_kallisto , ])
#WEAK. Transcript  Expr_Relapse VS Expr_Diagnosis	
nrow(Dr4_Binders_Ic50_NoExprFilt_ord[(Dr4_Binders_Ic50_NoExprFilt_ord$DR42_kallisto>=0.58 | Dr4_Binders_Ic50_NoExprFilt_ord$DR41_kallisto>=0.58) & Dr4_Binders_Ic50_NoExprFilt_ord$value>50  & (Dr4_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr4_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Dr4_Binders_Ic50_NoExprFilt_ord$DR42_kallisto > Dr4_Binders_Ic50_NoExprFilt_ord$DR41_kallisto , ])
nrow(Dr4_Binders_Ic50_NoExprFilt_ord[(Dr4_Binders_Ic50_NoExprFilt_ord$DR42_kallisto>=0.58 | Dr4_Binders_Ic50_NoExprFilt_ord$DR41_kallisto>=0.58) & Dr4_Binders_Ic50_NoExprFilt_ord$value>50  & (Dr4_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr4_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Dr4_Binders_Ic50_NoExprFilt_ord$DR42_kallisto < Dr4_Binders_Ic50_NoExprFilt_ord$DR41_kallisto , ])




#DR5
#N° MhAGs Strong
nrow(Dr5_Binders_Ic50_NoExprFilt_ord[Dr5_Binders_Ic50_NoExprFilt_ord$value<=50,])
#N° MhAGs Weak
nrow(Dr5_Binders_Ic50_NoExprFilt_ord[Dr5_Binders_Ic50_NoExprFilt_ord$value>50,])
#N° Expr. MhAGs  Strong
nrow(Dr5_Binders_Ic50_NoExprFilt_ord[(Dr5_Binders_Ic50_NoExprFilt_ord$DR52_kallisto>=0.58 | Dr5_Binders_Ic50_NoExprFilt_ord$DR51_kallisto>=0.58 ) & Dr5_Binders_Ic50_NoExprFilt_ord$value<=50 & (Dr5_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr5_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#N° Expr. MhAGs  Weak
nrow(Dr5_Binders_Ic50_NoExprFilt_ord[(Dr5_Binders_Ic50_NoExprFilt_ord$DR52_kallisto>=0.58 | Dr5_Binders_Ic50_NoExprFilt_ord$DR51_kallisto>=0.58 ) & Dr5_Binders_Ic50_NoExprFilt_ord$value>50 & (Dr5_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr5_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele ONLY
nrow(Dr5_Binders_Ic50_NoExprFilt_ord[(Dr5_Binders_Ic50_NoExprFilt_ord$DR52_kallisto>=0.58 | Dr5_Binders_Ic50_NoExprFilt_ord$DR51_kallisto>=0.58 ) & Dr5_Binders_Ic50_NoExprFilt_ord$value<=50  & Dr5_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3 & Dr5_Binders_Ic50_NoExprFilt_ord$Rel_Mut<3, ])
#STRONG  Expressed REL Allele Only
nrow(Dr5_Binders_Ic50_NoExprFilt_ord[(Dr5_Binders_Ic50_NoExprFilt_ord$DR52_kallisto>=0.58 | Dr5_Binders_Ic50_NoExprFilt_ord$DR51_kallisto>=0.58 ) & Dr5_Binders_Ic50_NoExprFilt_ord$value<=50  & Dr5_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 & Dr5_Binders_Ic50_NoExprFilt_ord$Dx_Mut<3 , ])
#STRONG  Expressed Common DIAG REL
nrow(Dr5_Binders_Ic50_NoExprFilt_ord[(Dr5_Binders_Ic50_NoExprFilt_ord$DR52_kallisto>=0.58 | Dr5_Binders_Ic50_NoExprFilt_ord$DR51_kallisto>=0.58 ) & Dr5_Binders_Ic50_NoExprFilt_ord$value<=50 & (Dr5_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 & Dr5_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele Diagnosis
nrow(Dr5_Binders_Ic50_NoExprFilt_ord[(Dr5_Binders_Ic50_NoExprFilt_ord$DR52_kallisto>=0.58 | Dr5_Binders_Ic50_NoExprFilt_ord$DR51_kallisto>=0.58 ) & Dr5_Binders_Ic50_NoExprFilt_ord$value<=50  & Dr5_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3 , ])
#STRONG  Expressed REL Allele Relapse
nrow(Dr5_Binders_Ic50_NoExprFilt_ord[(Dr5_Binders_Ic50_NoExprFilt_ord$DR52_kallisto>=0.58 | Dr5_Binders_Ic50_NoExprFilt_ord$DR51_kallisto>=0.58 ) & Dr5_Binders_Ic50_NoExprFilt_ord$value<=50  & Dr5_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3, ])
#WEAK Expressed DIAG Allele ONLY
nrow(Dr5_Binders_Ic50_NoExprFilt_ord[(Dr5_Binders_Ic50_NoExprFilt_ord$DR52_kallisto>=0.58 | Dr5_Binders_Ic50_NoExprFilt_ord$DR51_kallisto>=0.58 ) & Dr5_Binders_Ic50_NoExprFilt_ord$value>50  & Dr5_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3 & Dr5_Binders_Ic50_NoExprFilt_ord$Rel_Mut < 3, ])
#WEAK  Expressed REL Allele Only
nrow(Dr5_Binders_Ic50_NoExprFilt_ord[(Dr5_Binders_Ic50_NoExprFilt_ord$DR52_kallisto>=0.58 | Dr5_Binders_Ic50_NoExprFilt_ord$DR51_kallisto>=0.58 ) & Dr5_Binders_Ic50_NoExprFilt_ord$value>50  & Dr5_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3 & Dr5_Binders_Ic50_NoExprFilt_ord$Dx_Mut < 3 , ])
#WEAK  Expressed Common DIAG REL
nrow(Dr5_Binders_Ic50_NoExprFilt_ord[(Dr5_Binders_Ic50_NoExprFilt_ord$DR52_kallisto>=0.58 | Dr5_Binders_Ic50_NoExprFilt_ord$DR51_kallisto>=0.58 ) & Dr5_Binders_Ic50_NoExprFilt_ord$value>50 & (Dr5_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3 & Dr5_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3), ])
#WEAK Expressed DIAG Allele Diagnosis
nrow(Dr5_Binders_Ic50_NoExprFilt_ord[(Dr5_Binders_Ic50_NoExprFilt_ord$DR52_kallisto>=0.58 | Dr5_Binders_Ic50_NoExprFilt_ord$DR51_kallisto>=0.58 ) & Dr5_Binders_Ic50_NoExprFilt_ord$value>50  & Dr5_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3 , ])
#WEAK  Expressed REL Allele Relapse
nrow(Dr5_Binders_Ic50_NoExprFilt_ord[(Dr5_Binders_Ic50_NoExprFilt_ord$DR52_kallisto>=0.58 | Dr5_Binders_Ic50_NoExprFilt_ord$DR51_kallisto>=0.58 ) & Dr5_Binders_Ic50_NoExprFilt_ord$value>50  & Dr5_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3, ])
#STRONG Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Dr5_Binders_Ic50_NoExprFilt_ord[(Dr5_Binders_Ic50_NoExprFilt_ord$DR52_kallisto>=0.58 | Dr5_Binders_Ic50_NoExprFilt_ord$DR51_kallisto>=0.58 ) & Dr5_Binders_Ic50_NoExprFilt_ord$value<=50  & (Dr5_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr5_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Dr5_Binders_Ic50_NoExprFilt_ord$Rel_Mut > Dr5_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
nrow(Dr5_Binders_Ic50_NoExprFilt_ord[(Dr5_Binders_Ic50_NoExprFilt_ord$DR52_kallisto>=0.58 | Dr5_Binders_Ic50_NoExprFilt_ord$DR51_kallisto>=0.58 ) & Dr5_Binders_Ic50_NoExprFilt_ord$value<=50  & (Dr5_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr5_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Dr5_Binders_Ic50_NoExprFilt_ord$Rel_Mut < Dr5_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
#WEAK Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Dr5_Binders_Ic50_NoExprFilt_ord[(Dr5_Binders_Ic50_NoExprFilt_ord$DR52_kallisto>=0.58 | Dr5_Binders_Ic50_NoExprFilt_ord$DR51_kallisto>=0.58 ) & Dr5_Binders_Ic50_NoExprFilt_ord$value>50  & (Dr5_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr5_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Dr5_Binders_Ic50_NoExprFilt_ord$Rel_Mut > Dr5_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
nrow(Dr5_Binders_Ic50_NoExprFilt_ord[(Dr5_Binders_Ic50_NoExprFilt_ord$DR52_kallisto>=0.58 | Dr5_Binders_Ic50_NoExprFilt_ord$DR51_kallisto>=0.58 ) & Dr5_Binders_Ic50_NoExprFilt_ord$value>50  & (Dr5_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr5_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Dr5_Binders_Ic50_NoExprFilt_ord$Rel_Mut < Dr5_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
#STRONG   Trascript  Expr_Relapse VS Expr_Diagnosis	
nrow(Dr5_Binders_Ic50_NoExprFilt_ord[(Dr5_Binders_Ic50_NoExprFilt_ord$DR52_kallisto>=0.58 | Dr5_Binders_Ic50_NoExprFilt_ord$DR51_kallisto>=0.58 ) & Dr5_Binders_Ic50_NoExprFilt_ord$value<=50  & (Dr5_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr5_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Dr5_Binders_Ic50_NoExprFilt_ord$DR52_kallisto  > Dr5_Binders_Ic50_NoExprFilt_ord$DR51_kallisto , ])
nrow(Dr5_Binders_Ic50_NoExprFilt_ord[(Dr5_Binders_Ic50_NoExprFilt_ord$DR52_kallisto>=0.58 | Dr5_Binders_Ic50_NoExprFilt_ord$DR51_kallisto>=0.58 ) & Dr5_Binders_Ic50_NoExprFilt_ord$value<=50  & (Dr5_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr5_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Dr5_Binders_Ic50_NoExprFilt_ord$DR52_kallisto < Dr5_Binders_Ic50_NoExprFilt_ord$DR51_kallisto , ])
#WEAK. Transcript  Expr_Relapse VS Expr_Diagnosis	
nrow(Dr5_Binders_Ic50_NoExprFilt_ord[(Dr5_Binders_Ic50_NoExprFilt_ord$DR52_kallisto>=0.58 | Dr5_Binders_Ic50_NoExprFilt_ord$DR51_kallisto>=0.58 ) & Dr5_Binders_Ic50_NoExprFilt_ord$value>50  & (Dr5_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr5_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Dr5_Binders_Ic50_NoExprFilt_ord$DR52_kallisto > Dr5_Binders_Ic50_NoExprFilt_ord$DR51_kallisto , ])
nrow(Dr5_Binders_Ic50_NoExprFilt_ord[(Dr5_Binders_Ic50_NoExprFilt_ord$DR52_kallisto>=0.58 | Dr5_Binders_Ic50_NoExprFilt_ord$DR51_kallisto>=0.58 ) & Dr5_Binders_Ic50_NoExprFilt_ord$value>50  & (Dr5_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Dr5_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Dr5_Binders_Ic50_NoExprFilt_ord$DR52_kallisto < Dr5_Binders_Ic50_NoExprFilt_ord$DR51_kallisto , ])




#FOCA
#N° MhAGs Strong
nrow(Foca_Binders_Ic50_NoExprFilt_ord[Foca_Binders_Ic50_NoExprFilt_ord$value<=50,])
#N° MhAGs Weak
nrow(Foca_Binders_Ic50_NoExprFilt_ord[Foca_Binders_Ic50_NoExprFilt_ord$value>50,])
#N° Expr. MhAGs  Strong
nrow(Foca_Binders_Ic50_NoExprFilt_ord[(Foca_Binders_Ic50_NoExprFilt_ord$FOCA2_kallisto>=0.58 | Foca_Binders_Ic50_NoExprFilt_ord$FOCA1_kallisto>=0.58 ) & Foca_Binders_Ic50_NoExprFilt_ord$value<=50 & (Foca_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Foca_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#N° Expr. MhAGs  Weak
nrow(Foca_Binders_Ic50_NoExprFilt_ord[(Foca_Binders_Ic50_NoExprFilt_ord$FOCA2_kallisto>=0.58 | Foca_Binders_Ic50_NoExprFilt_ord$FOCA1_kallisto>=0.58 ) & Foca_Binders_Ic50_NoExprFilt_ord$value>50 & (Foca_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Foca_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele ONLY
nrow(Foca_Binders_Ic50_NoExprFilt_ord[(Foca_Binders_Ic50_NoExprFilt_ord$FOCA2_kallisto>=0.58 | Foca_Binders_Ic50_NoExprFilt_ord$FOCA1_kallisto>=0.58 ) & Foca_Binders_Ic50_NoExprFilt_ord$value<=50  & Foca_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3 & Foca_Binders_Ic50_NoExprFilt_ord$Rel_Mut<3, ])
#STRONG  Expressed REL Allele Only
nrow(Foca_Binders_Ic50_NoExprFilt_ord[(Foca_Binders_Ic50_NoExprFilt_ord$FOCA2_kallisto>=0.58 | Foca_Binders_Ic50_NoExprFilt_ord$FOCA1_kallisto>=0.58 ) & Foca_Binders_Ic50_NoExprFilt_ord$value<=50  & Foca_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 & Foca_Binders_Ic50_NoExprFilt_ord$Dx_Mut<3 , ])
#STRONG  Expressed Common DIAG REL
nrow(Foca_Binders_Ic50_NoExprFilt_ord[(Foca_Binders_Ic50_NoExprFilt_ord$FOCA2_kallisto>=0.58 | Foca_Binders_Ic50_NoExprFilt_ord$FOCA1_kallisto>=0.58 ) & Foca_Binders_Ic50_NoExprFilt_ord$value<=50 & (Foca_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 & Foca_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele Diagnosis
nrow(Foca_Binders_Ic50_NoExprFilt_ord[(Foca_Binders_Ic50_NoExprFilt_ord$FOCA2_kallisto>=0.58 | Foca_Binders_Ic50_NoExprFilt_ord$FOCA1_kallisto>=0.58 ) & Foca_Binders_Ic50_NoExprFilt_ord$value<=50  & Foca_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3 , ])
#STRONG  Expressed REL Allele Relapse
nrow(Foca_Binders_Ic50_NoExprFilt_ord[(Foca_Binders_Ic50_NoExprFilt_ord$FOCA2_kallisto>=0.58 | Foca_Binders_Ic50_NoExprFilt_ord$FOCA1_kallisto>=0.58 ) & Foca_Binders_Ic50_NoExprFilt_ord$value<=50  & Foca_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3, ])
#WEAK Expressed DIAG Allele ONLY
nrow(Foca_Binders_Ic50_NoExprFilt_ord[(Foca_Binders_Ic50_NoExprFilt_ord$FOCA2_kallisto>=0.58 | Foca_Binders_Ic50_NoExprFilt_ord$FOCA1_kallisto>=0.58 ) & Foca_Binders_Ic50_NoExprFilt_ord$value>50  & Foca_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3 & Foca_Binders_Ic50_NoExprFilt_ord$Rel_Mut < 3, ])
#WEAK  Expressed REL Allele Only
nrow(Foca_Binders_Ic50_NoExprFilt_ord[(Foca_Binders_Ic50_NoExprFilt_ord$FOCA2_kallisto>=0.58 | Foca_Binders_Ic50_NoExprFilt_ord$FOCA1_kallisto>=0.58 ) & Foca_Binders_Ic50_NoExprFilt_ord$value>50  & Foca_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3 & Foca_Binders_Ic50_NoExprFilt_ord$Dx_Mut < 3 , ])
#WEAK  Expressed Common DIAG REL
nrow(Foca_Binders_Ic50_NoExprFilt_ord[(Foca_Binders_Ic50_NoExprFilt_ord$FOCA2_kallisto>=0.58 | Foca_Binders_Ic50_NoExprFilt_ord$FOCA1_kallisto>=0.58 ) & Foca_Binders_Ic50_NoExprFilt_ord$value>50 & (Foca_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3 & Foca_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3), ])
#WEAK Expressed DIAG Allele Diagnosis
nrow(Foca_Binders_Ic50_NoExprFilt_ord[(Foca_Binders_Ic50_NoExprFilt_ord$FOCA2_kallisto>=0.58 | Foca_Binders_Ic50_NoExprFilt_ord$FOCA1_kallisto>=0.58 ) & Foca_Binders_Ic50_NoExprFilt_ord$value>50  & Foca_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3 , ])
#WEAK  Expressed REL Allele Relapse
nrow(Foca_Binders_Ic50_NoExprFilt_ord[(Foca_Binders_Ic50_NoExprFilt_ord$FOCA2_kallisto>=0.58 | Foca_Binders_Ic50_NoExprFilt_ord$FOCA1_kallisto>=0.58 ) & Foca_Binders_Ic50_NoExprFilt_ord$value>50  & Foca_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3, ])
#STRONG Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Foca_Binders_Ic50_NoExprFilt_ord[(Foca_Binders_Ic50_NoExprFilt_ord$FOCA2_kallisto>=0.58 | Foca_Binders_Ic50_NoExprFilt_ord$FOCA1_kallisto>=0.58 ) & Foca_Binders_Ic50_NoExprFilt_ord$value<=50  & (Foca_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Foca_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Foca_Binders_Ic50_NoExprFilt_ord$Rel_Mut > Foca_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
nrow(Foca_Binders_Ic50_NoExprFilt_ord[(Foca_Binders_Ic50_NoExprFilt_ord$FOCA2_kallisto>=0.58 | Foca_Binders_Ic50_NoExprFilt_ord$FOCA1_kallisto>=0.58 ) & Foca_Binders_Ic50_NoExprFilt_ord$value<=50  & (Foca_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Foca_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Foca_Binders_Ic50_NoExprFilt_ord$Rel_Mut < Foca_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
#WEAK Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Foca_Binders_Ic50_NoExprFilt_ord[(Foca_Binders_Ic50_NoExprFilt_ord$FOCA2_kallisto>=0.58 | Foca_Binders_Ic50_NoExprFilt_ord$FOCA1_kallisto>=0.58 ) & Foca_Binders_Ic50_NoExprFilt_ord$value>50  & (Foca_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Foca_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Foca_Binders_Ic50_NoExprFilt_ord$Rel_Mut > Foca_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
nrow(Foca_Binders_Ic50_NoExprFilt_ord[(Foca_Binders_Ic50_NoExprFilt_ord$FOCA2_kallisto>=0.58 | Foca_Binders_Ic50_NoExprFilt_ord$FOCA1_kallisto>=0.58 ) & Foca_Binders_Ic50_NoExprFilt_ord$value>50  & (Foca_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Foca_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Foca_Binders_Ic50_NoExprFilt_ord$Rel_Mut < Foca_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
#STRONG   Trascript  Expr_Relapse VS Expr_Diagnosis	
nrow(Foca_Binders_Ic50_NoExprFilt_ord[(Foca_Binders_Ic50_NoExprFilt_ord$FOCA2_kallisto>=0.58 | Foca_Binders_Ic50_NoExprFilt_ord$FOCA1_kallisto>=0.58 ) & Foca_Binders_Ic50_NoExprFilt_ord$value<=50  & (Foca_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Foca_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Foca_Binders_Ic50_NoExprFilt_ord$FOCA2_kallisto  > Foca_Binders_Ic50_NoExprFilt_ord$FOCA1_kallisto , ])
nrow(Foca_Binders_Ic50_NoExprFilt_ord[(Foca_Binders_Ic50_NoExprFilt_ord$FOCA2_kallisto>=0.58 | Foca_Binders_Ic50_NoExprFilt_ord$FOCA1_kallisto>=0.58 ) & Foca_Binders_Ic50_NoExprFilt_ord$value<=50  & (Foca_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Foca_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Foca_Binders_Ic50_NoExprFilt_ord$FOCA2_kallisto < Foca_Binders_Ic50_NoExprFilt_ord$FOCA1_kallisto , ])
#WEAK. Transcript  Expr_Relapse VS Expr_Diagnosis	
nrow(Foca_Binders_Ic50_NoExprFilt_ord[(Foca_Binders_Ic50_NoExprFilt_ord$FOCA2_kallisto>=0.58 | Foca_Binders_Ic50_NoExprFilt_ord$FOCA1_kallisto>=0.58 ) & Foca_Binders_Ic50_NoExprFilt_ord$value>50  & (Foca_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Foca_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Foca_Binders_Ic50_NoExprFilt_ord$FOCA2_kallisto > Foca_Binders_Ic50_NoExprFilt_ord$FOCA1_kallisto , ])
nrow(Foca_Binders_Ic50_NoExprFilt_ord[(Foca_Binders_Ic50_NoExprFilt_ord$FOCA2_kallisto>=0.58 | Foca_Binders_Ic50_NoExprFilt_ord$FOCA1_kallisto>=0.58 ) & Foca_Binders_Ic50_NoExprFilt_ord$value>50  & (Foca_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Foca_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Foca_Binders_Ic50_NoExprFilt_ord$FOCA2_kallisto < Foca_Binders_Ic50_NoExprFilt_ord$FOCA1_kallisto , ])





#GAGRA
#N° MhAGs Strong
nrow(Gagra_Binders_Ic50_NoExprFilt_ord[Gagra_Binders_Ic50_NoExprFilt_ord$value<=50,])
#N° MhAGs Weak
nrow(Gagra_Binders_Ic50_NoExprFilt_ord[Gagra_Binders_Ic50_NoExprFilt_ord$value>50,])
#N° Expr. MhAGs  Strong
nrow(Gagra_Binders_Ic50_NoExprFilt_ord[(Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA2_kallisto >=0.58 | Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA1_kallisto> 0.58) & Gagra_Binders_Ic50_NoExprFilt_ord$value<=50 & (Gagra_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Gagra_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#N° Expr. MhAGs  Weak
nrow(Gagra_Binders_Ic50_NoExprFilt_ord[(Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA2_kallisto >=0.58 | Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA1_kallisto> 0.58) & Gagra_Binders_Ic50_NoExprFilt_ord$value>50 & (Gagra_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Gagra_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele ONLY
nrow(Gagra_Binders_Ic50_NoExprFilt_ord[(Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA2_kallisto >=0.58 | Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA1_kallisto> 0.58) & Gagra_Binders_Ic50_NoExprFilt_ord$value<=50  & Gagra_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3 & Gagra_Binders_Ic50_NoExprFilt_ord$Rel_Mut<3, ])
#STRONG  Expressed REL Allele Only
nrow(Gagra_Binders_Ic50_NoExprFilt_ord[(Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA2_kallisto >=0.58 | Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA1_kallisto> 0.58) & Gagra_Binders_Ic50_NoExprFilt_ord$value<=50  & Gagra_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 & Gagra_Binders_Ic50_NoExprFilt_ord$Dx_Mut<3 , ])
#STRONG  Expressed Common DIAG REL
nrow(Gagra_Binders_Ic50_NoExprFilt_ord[(Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA2_kallisto >=0.58 | Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA1_kallisto> 0.58) & Gagra_Binders_Ic50_NoExprFilt_ord$value<=50 & (Gagra_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 & Gagra_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele Diagnosis
nrow(Gagra_Binders_Ic50_NoExprFilt_ord[(Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA2_kallisto >=0.58 | Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA1_kallisto> 0.58) & Gagra_Binders_Ic50_NoExprFilt_ord$value<=50  & Gagra_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3 , ])
#STRONG  Expressed REL Allele Relapse
nrow(Gagra_Binders_Ic50_NoExprFilt_ord[(Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA2_kallisto >=0.58 | Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA1_kallisto> 0.58) & Gagra_Binders_Ic50_NoExprFilt_ord$value<=50  & Gagra_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3, ])
#WEAK Expressed DIAG Allele ONLY
nrow(Gagra_Binders_Ic50_NoExprFilt_ord[(Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA2_kallisto >=0.58 | Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA1_kallisto> 0.58) & Gagra_Binders_Ic50_NoExprFilt_ord$value>50  & Gagra_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3 & Gagra_Binders_Ic50_NoExprFilt_ord$Rel_Mut < 3, ])
#WEAK  Expressed REL Allele Only
nrow(Gagra_Binders_Ic50_NoExprFilt_ord[(Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA2_kallisto >=0.58 | Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA1_kallisto> 0.58) & Gagra_Binders_Ic50_NoExprFilt_ord$value>50  & Gagra_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3 & Gagra_Binders_Ic50_NoExprFilt_ord$Dx_Mut < 3 , ])
#WEAK  Expressed Common DIAG REL
nrow(Gagra_Binders_Ic50_NoExprFilt_ord[(Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA2_kallisto >=0.58 | Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA1_kallisto> 0.58) & Gagra_Binders_Ic50_NoExprFilt_ord$value>50 & (Gagra_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3 & Gagra_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3), ])
#WEAK Expressed DIAG Allele Diagnosis
nrow(Gagra_Binders_Ic50_NoExprFilt_ord[(Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA2_kallisto >=0.58 | Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA1_kallisto> 0.58) & Gagra_Binders_Ic50_NoExprFilt_ord$value>50  & Gagra_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3 , ])
#WEAK  Expressed REL Allele Relapse
nrow(Gagra_Binders_Ic50_NoExprFilt_ord[(Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA2_kallisto >=0.58 | Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA1_kallisto> 0.58) & Gagra_Binders_Ic50_NoExprFilt_ord$value>50  & Gagra_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3, ])
#STRONG Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Gagra_Binders_Ic50_NoExprFilt_ord[(Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA2_kallisto >=0.58 | Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA1_kallisto> 0.58) & Gagra_Binders_Ic50_NoExprFilt_ord$value<=50  & (Gagra_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Gagra_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Gagra_Binders_Ic50_NoExprFilt_ord$Rel_Mut > Gagra_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
nrow(Gagra_Binders_Ic50_NoExprFilt_ord[(Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA2_kallisto >=0.58 | Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA1_kallisto> 0.58) & Gagra_Binders_Ic50_NoExprFilt_ord$value<=50  & (Gagra_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Gagra_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Gagra_Binders_Ic50_NoExprFilt_ord$Rel_Mut < Gagra_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
#WEAK Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Gagra_Binders_Ic50_NoExprFilt_ord[(Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA2_kallisto >=0.58 | Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA1_kallisto> 0.58) & Gagra_Binders_Ic50_NoExprFilt_ord$value>50  & (Gagra_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Gagra_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Gagra_Binders_Ic50_NoExprFilt_ord$Rel_Mut > Gagra_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
nrow(Gagra_Binders_Ic50_NoExprFilt_ord[(Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA2_kallisto >=0.58 | Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA1_kallisto> 0.58) & Gagra_Binders_Ic50_NoExprFilt_ord$value>50  & (Gagra_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Gagra_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Gagra_Binders_Ic50_NoExprFilt_ord$Rel_Mut < Gagra_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
#STRONG   Trascript  Expr_Relapse VS Expr_Diagnosis	
nrow(Gagra_Binders_Ic50_NoExprFilt_ord[(Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA2_kallisto >=0.58 | Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA1_kallisto> 0.58) & Gagra_Binders_Ic50_NoExprFilt_ord$value<=50  & (Gagra_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Gagra_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA2_kallisto  > Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA1_kallisto , ])
nrow(Gagra_Binders_Ic50_NoExprFilt_ord[(Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA2_kallisto >=0.58 | Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA1_kallisto> 0.58) & Gagra_Binders_Ic50_NoExprFilt_ord$value<=50  & (Gagra_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Gagra_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA2_kallisto < Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA1_kallisto , ])
#WEAK. Transcript  Expr_Relapse VS Expr_Diagnosis	
nrow(Gagra_Binders_Ic50_NoExprFilt_ord[(Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA2_kallisto >=0.58 | Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA1_kallisto> 0.58) & Gagra_Binders_Ic50_NoExprFilt_ord$value>50  & (Gagra_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Gagra_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA2_kallisto > Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA1_kallisto , ])
nrow(Gagra_Binders_Ic50_NoExprFilt_ord[(Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA2_kallisto >=0.58 | Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA1_kallisto> 0.58) & Gagra_Binders_Ic50_NoExprFilt_ord$value>50  & (Gagra_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Gagra_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA2_kallisto < Gagra_Binders_Ic50_NoExprFilt_ord$GAGRA1_kallisto , ])






#LUAN
#N° MhAGs Strong
nrow(Luan_Binders_Ic50_NoExprFilt_ord_homo[Luan_Binders_Ic50_NoExprFilt_ord_homo$value<=50,])
#N° MhAGs Weak
nrow(Luan_Binders_Ic50_NoExprFilt_ord_homo[Luan_Binders_Ic50_NoExprFilt_ord_homo$value>50,])
#N° Expr. MhAGs  Strong
nrow(Luan_Binders_Ic50_NoExprFilt_ord_homo[ (Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN2_kallisto>=0.58 | Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN1_kallisto >= 0.58) & Luan_Binders_Ic50_NoExprFilt_ord_homo$value<=50 & (Luan_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Luan_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3), ])
#N° Expr. MhAGs  Weak
nrow(Luan_Binders_Ic50_NoExprFilt_ord_homo[ (Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN2_kallisto>=0.58 | Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN1_kallisto >= 0.58) & Luan_Binders_Ic50_NoExprFilt_ord_homo$value>50 & (Luan_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Luan_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele ONLY
nrow(Luan_Binders_Ic50_NoExprFilt_ord_homo[ (Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN2_kallisto>=0.58 | Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN1_kallisto >= 0.58) & Luan_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & Luan_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3 & Luan_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut<3, ])
#STRONG  Expressed REL Allele Only
nrow(Luan_Binders_Ic50_NoExprFilt_ord_homo[ (Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN2_kallisto>=0.58 | Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN1_kallisto >= 0.58) & Luan_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & Luan_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 & Luan_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut<3 , ])
#STRONG  Expressed Common DIAG REL
nrow(Luan_Binders_Ic50_NoExprFilt_ord_homo[ (Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN2_kallisto>=0.58 | Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN1_kallisto >= 0.58) & Luan_Binders_Ic50_NoExprFilt_ord_homo$value<=50 & (Luan_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 & Luan_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele Diagnosis
nrow(Luan_Binders_Ic50_NoExprFilt_ord_homo[ (Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN2_kallisto>=0.58 | Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN1_kallisto >= 0.58) & Luan_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & Luan_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3 , ])
#STRONG  Expressed REL Allele Relapse
nrow(Luan_Binders_Ic50_NoExprFilt_ord_homo[ (Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN2_kallisto>=0.58 | Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN1_kallisto >= 0.58) & Luan_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & Luan_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3, ])
#WEAK Expressed DIAG Allele ONLY
nrow(Luan_Binders_Ic50_NoExprFilt_ord_homo[ (Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN2_kallisto>=0.58 | Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN1_kallisto >= 0.58) & Luan_Binders_Ic50_NoExprFilt_ord_homo$value>50  & Luan_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut >= 3 & Luan_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut < 3, ])
#WEAK  Expressed REL Allele Only
nrow(Luan_Binders_Ic50_NoExprFilt_ord_homo[ (Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN2_kallisto>=0.58 | Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN1_kallisto >= 0.58) & Luan_Binders_Ic50_NoExprFilt_ord_homo$value>50  & Luan_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut >= 3 & Luan_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut < 3 , ])
#WEAK  Expressed Common DIAG REL
nrow(Luan_Binders_Ic50_NoExprFilt_ord_homo[ (Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN2_kallisto>=0.58 | Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN1_kallisto >= 0.58) & Luan_Binders_Ic50_NoExprFilt_ord_homo$value>50 & (Luan_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut >= 3 & Luan_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut >= 3), ])
#WEAK Expressed DIAG Allele Diagnosis
nrow(Luan_Binders_Ic50_NoExprFilt_ord_homo[ (Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN2_kallisto>=0.58 | Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN1_kallisto >= 0.58) & Luan_Binders_Ic50_NoExprFilt_ord_homo$value>50  & Luan_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut >= 3 , ])
#WEAK  Expressed REL Allele Relapse
nrow(Luan_Binders_Ic50_NoExprFilt_ord_homo[ (Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN2_kallisto>=0.58 | Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN1_kallisto >= 0.58) & Luan_Binders_Ic50_NoExprFilt_ord_homo$value>50  & Luan_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut >= 3, ])
#STRONG Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Luan_Binders_Ic50_NoExprFilt_ord_homo[ (Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN2_kallisto>=0.58 | Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN1_kallisto >= 0.58) & Luan_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & (Luan_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Luan_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3) & Luan_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut > Luan_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut , ])
nrow(Luan_Binders_Ic50_NoExprFilt_ord_homo[ (Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN2_kallisto>=0.58 | Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN1_kallisto >= 0.58) & Luan_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & (Luan_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Luan_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3) & Luan_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut < Luan_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut , ])
#WEAK Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Luan_Binders_Ic50_NoExprFilt_ord_homo[ (Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN2_kallisto>=0.58 | Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN1_kallisto >= 0.58) & Luan_Binders_Ic50_NoExprFilt_ord_homo$value>50  & (Luan_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Luan_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3) & Luan_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut > Luan_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut , ])
nrow(Luan_Binders_Ic50_NoExprFilt_ord_homo[ (Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN2_kallisto>=0.58 | Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN1_kallisto >= 0.58) & Luan_Binders_Ic50_NoExprFilt_ord_homo$value>50  & (Luan_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Luan_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3) & Luan_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut < Luan_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut , ])
#STRONG   Trascript  Expr_Relapse VS Expr_Diagnosis	
nrow(Luan_Binders_Ic50_NoExprFilt_ord_homo[ (Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN2_kallisto>=0.58 | Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN1_kallisto >= 0.58) & Luan_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & (Luan_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Luan_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3)  & Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN2_kallisto  > Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN1_kallisto , ])
nrow(Luan_Binders_Ic50_NoExprFilt_ord_homo[ (Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN2_kallisto>=0.58 | Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN1_kallisto >= 0.58) & Luan_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & (Luan_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Luan_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3)  & Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN2_kallisto < Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN1_kallisto , ])
#WEAK. Transcript  Expr_Relapse VS Expr_Diagnosis	
nrow(Luan_Binders_Ic50_NoExprFilt_ord_homo[ (Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN2_kallisto>=0.58 | Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN1_kallisto >= 0.58) & Luan_Binders_Ic50_NoExprFilt_ord_homo$value>50  & (Luan_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Luan_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3)  & Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN2_kallisto > Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN1_kallisto , ])
nrow(Luan_Binders_Ic50_NoExprFilt_ord_homo[ (Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN2_kallisto>=0.58 | Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN1_kallisto >= 0.58) & Luan_Binders_Ic50_NoExprFilt_ord_homo$value>50  & (Luan_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Luan_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3)  & Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN2_kallisto < Luan_Binders_Ic50_NoExprFilt_ord_homo$LUAN1_kallisto , ])





#MABI
#N° MhAGs Strong
nrow(Mabi_Binders_Ic50_NoExprFilt_ord[Mabi_Binders_Ic50_NoExprFilt_ord$value<=50,])
#N° MhAGs Weak
nrow(Mabi_Binders_Ic50_NoExprFilt_ord[Mabi_Binders_Ic50_NoExprFilt_ord$value>50,])
#N° Expr. MhAGs  Strong
nrow(Mabi_Binders_Ic50_NoExprFilt_ord[(Mabi_Binders_Ic50_NoExprFilt_ord$MABI2_kallisto>=0.58 | Mabi_Binders_Ic50_NoExprFilt_ord$MABI1_kallisto>=0.58) & Mabi_Binders_Ic50_NoExprFilt_ord$value<=50 & (Mabi_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Mabi_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#N° Expr. MhAGs  Weak
nrow(Mabi_Binders_Ic50_NoExprFilt_ord[(Mabi_Binders_Ic50_NoExprFilt_ord$MABI2_kallisto>=0.58 | Mabi_Binders_Ic50_NoExprFilt_ord$MABI1_kallisto>=0.58) & Mabi_Binders_Ic50_NoExprFilt_ord$value>50 & (Mabi_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Mabi_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele ONLY
nrow(Mabi_Binders_Ic50_NoExprFilt_ord[(Mabi_Binders_Ic50_NoExprFilt_ord$MABI2_kallisto>=0.58 | Mabi_Binders_Ic50_NoExprFilt_ord$MABI1_kallisto>=0.58) & Mabi_Binders_Ic50_NoExprFilt_ord$value<=50  & Mabi_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3 & Mabi_Binders_Ic50_NoExprFilt_ord$Rel_Mut<3, ])
#STRONG  Expressed REL Allele Only
nrow(Mabi_Binders_Ic50_NoExprFilt_ord[(Mabi_Binders_Ic50_NoExprFilt_ord$MABI2_kallisto>=0.58 | Mabi_Binders_Ic50_NoExprFilt_ord$MABI1_kallisto>=0.58) & Mabi_Binders_Ic50_NoExprFilt_ord$value<=50  & Mabi_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 & Mabi_Binders_Ic50_NoExprFilt_ord$Dx_Mut<3 , ])
#STRONG  Expressed Common DIAG REL
nrow(Mabi_Binders_Ic50_NoExprFilt_ord[(Mabi_Binders_Ic50_NoExprFilt_ord$MABI2_kallisto>=0.58 | Mabi_Binders_Ic50_NoExprFilt_ord$MABI1_kallisto>=0.58) & Mabi_Binders_Ic50_NoExprFilt_ord$value<=50 & (Mabi_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 & Mabi_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele Diagnosis
nrow(Mabi_Binders_Ic50_NoExprFilt_ord[(Mabi_Binders_Ic50_NoExprFilt_ord$MABI2_kallisto>=0.58 | Mabi_Binders_Ic50_NoExprFilt_ord$MABI1_kallisto>=0.58) & Mabi_Binders_Ic50_NoExprFilt_ord$value<=50  & Mabi_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3 , ])
#STRONG  Expressed REL Allele Relapse
nrow(Mabi_Binders_Ic50_NoExprFilt_ord[(Mabi_Binders_Ic50_NoExprFilt_ord$MABI2_kallisto>=0.58 | Mabi_Binders_Ic50_NoExprFilt_ord$MABI1_kallisto>=0.58) & Mabi_Binders_Ic50_NoExprFilt_ord$value<=50  & Mabi_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3, ])
#WEAK Expressed DIAG Allele ONLY
nrow(Mabi_Binders_Ic50_NoExprFilt_ord[(Mabi_Binders_Ic50_NoExprFilt_ord$MABI2_kallisto>=0.58 | Mabi_Binders_Ic50_NoExprFilt_ord$MABI1_kallisto>=0.58) & Mabi_Binders_Ic50_NoExprFilt_ord$value>50  & Mabi_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3 & Mabi_Binders_Ic50_NoExprFilt_ord$Rel_Mut < 3, ])
#WEAK  Expressed REL Allele Only
nrow(Mabi_Binders_Ic50_NoExprFilt_ord[(Mabi_Binders_Ic50_NoExprFilt_ord$MABI2_kallisto>=0.58 | Mabi_Binders_Ic50_NoExprFilt_ord$MABI1_kallisto>=0.58) & Mabi_Binders_Ic50_NoExprFilt_ord$value>50  & Mabi_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3 & Mabi_Binders_Ic50_NoExprFilt_ord$Dx_Mut < 3 , ])
#WEAK  Expressed Common DIAG REL
nrow(Mabi_Binders_Ic50_NoExprFilt_ord[(Mabi_Binders_Ic50_NoExprFilt_ord$MABI2_kallisto>=0.58 | Mabi_Binders_Ic50_NoExprFilt_ord$MABI1_kallisto>=0.58) & Mabi_Binders_Ic50_NoExprFilt_ord$value>50 & (Mabi_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3 & Mabi_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3), ])
#WEAK Expressed DIAG Allele Diagnosis
nrow(Mabi_Binders_Ic50_NoExprFilt_ord[(Mabi_Binders_Ic50_NoExprFilt_ord$MABI2_kallisto>=0.58 | Mabi_Binders_Ic50_NoExprFilt_ord$MABI1_kallisto>=0.58) & Mabi_Binders_Ic50_NoExprFilt_ord$value>50  & Mabi_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3 , ])
#WEAK  Expressed REL Allele Relapse
nrow(Mabi_Binders_Ic50_NoExprFilt_ord[(Mabi_Binders_Ic50_NoExprFilt_ord$MABI2_kallisto>=0.58 | Mabi_Binders_Ic50_NoExprFilt_ord$MABI1_kallisto>=0.58) & Mabi_Binders_Ic50_NoExprFilt_ord$value>50  & Mabi_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3, ])
#STRONG Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Mabi_Binders_Ic50_NoExprFilt_ord[(Mabi_Binders_Ic50_NoExprFilt_ord$MABI2_kallisto>=0.58 | Mabi_Binders_Ic50_NoExprFilt_ord$MABI1_kallisto>=0.58) & Mabi_Binders_Ic50_NoExprFilt_ord$value<=50  & (Mabi_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Mabi_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Mabi_Binders_Ic50_NoExprFilt_ord$Rel_Mut > Mabi_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
nrow(Mabi_Binders_Ic50_NoExprFilt_ord[(Mabi_Binders_Ic50_NoExprFilt_ord$MABI2_kallisto>=0.58 | Mabi_Binders_Ic50_NoExprFilt_ord$MABI1_kallisto>=0.58) & Mabi_Binders_Ic50_NoExprFilt_ord$value<=50  & (Mabi_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Mabi_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Mabi_Binders_Ic50_NoExprFilt_ord$Rel_Mut < Mabi_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
#WEAK Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Mabi_Binders_Ic50_NoExprFilt_ord[(Mabi_Binders_Ic50_NoExprFilt_ord$MABI2_kallisto>=0.58 | Mabi_Binders_Ic50_NoExprFilt_ord$MABI1_kallisto>=0.58) & Mabi_Binders_Ic50_NoExprFilt_ord$value>50  & (Mabi_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Mabi_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Mabi_Binders_Ic50_NoExprFilt_ord$Rel_Mut > Mabi_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
nrow(Mabi_Binders_Ic50_NoExprFilt_ord[(Mabi_Binders_Ic50_NoExprFilt_ord$MABI2_kallisto>=0.58 | Mabi_Binders_Ic50_NoExprFilt_ord$MABI1_kallisto>=0.58) & Mabi_Binders_Ic50_NoExprFilt_ord$value>50  & (Mabi_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Mabi_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Mabi_Binders_Ic50_NoExprFilt_ord$Rel_Mut < Mabi_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
#STRONG   Trascript  Expr_Relapse VS Expr_Diagnosis	
nrow(Mabi_Binders_Ic50_NoExprFilt_ord[(Mabi_Binders_Ic50_NoExprFilt_ord$MABI2_kallisto>=0.58 | Mabi_Binders_Ic50_NoExprFilt_ord$MABI1_kallisto>=0.58) & Mabi_Binders_Ic50_NoExprFilt_ord$value<=50  & (Mabi_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Mabi_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Mabi_Binders_Ic50_NoExprFilt_ord$MABI2_kallisto  > Mabi_Binders_Ic50_NoExprFilt_ord$MABI1_kallisto , ])
nrow(Mabi_Binders_Ic50_NoExprFilt_ord[(Mabi_Binders_Ic50_NoExprFilt_ord$MABI2_kallisto>=0.58 | Mabi_Binders_Ic50_NoExprFilt_ord$MABI1_kallisto>=0.58) & Mabi_Binders_Ic50_NoExprFilt_ord$value<=50  & (Mabi_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Mabi_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Mabi_Binders_Ic50_NoExprFilt_ord$MABI2_kallisto < Mabi_Binders_Ic50_NoExprFilt_ord$MABI1_kallisto , ])
#WEAK. Transcript  Expr_Relapse VS Expr_Diagnosis	
nrow(Mabi_Binders_Ic50_NoExprFilt_ord[(Mabi_Binders_Ic50_NoExprFilt_ord$MABI2_kallisto>=0.58 | Mabi_Binders_Ic50_NoExprFilt_ord$MABI1_kallisto>=0.58) & Mabi_Binders_Ic50_NoExprFilt_ord$value>50  & (Mabi_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Mabi_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Mabi_Binders_Ic50_NoExprFilt_ord$MABI2_kallisto > Mabi_Binders_Ic50_NoExprFilt_ord$MABI1_kallisto , ])
nrow(Mabi_Binders_Ic50_NoExprFilt_ord[(Mabi_Binders_Ic50_NoExprFilt_ord$MABI2_kallisto>=0.58 | Mabi_Binders_Ic50_NoExprFilt_ord$MABI1_kallisto>=0.58) & Mabi_Binders_Ic50_NoExprFilt_ord$value>50  & (Mabi_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Mabi_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Mabi_Binders_Ic50_NoExprFilt_ord$MABI2_kallisto < Mabi_Binders_Ic50_NoExprFilt_ord$MABI1_kallisto , ])






#MOGE
#N° MhAGs Strong
nrow(Moge_Binders_Ic50_NoExprFilt_ord_homo[Moge_Binders_Ic50_NoExprFilt_ord_homo$value<=50,])
#N° MhAGs Weak
nrow(Moge_Binders_Ic50_NoExprFilt_ord_homo[Moge_Binders_Ic50_NoExprFilt_ord_homo$value>50,])
#N° Expr. MhAGs  Strong
nrow(Moge_Binders_Ic50_NoExprFilt_ord_homo[(Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE2_kallisto >=0.58 | Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE1_kallisto >= 0.58) &  Moge_Binders_Ic50_NoExprFilt_ord_homo$value<=50 & (Moge_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Moge_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3), ])
#N° Expr. MhAGs  Weak
nrow(Moge_Binders_Ic50_NoExprFilt_ord_homo[(Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE2_kallisto >=0.58 | Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE1_kallisto >= 0.58) &  Moge_Binders_Ic50_NoExprFilt_ord_homo$value>50 & (Moge_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Moge_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele ONLY
nrow(Moge_Binders_Ic50_NoExprFilt_ord_homo[(Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE2_kallisto >=0.58 | Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE1_kallisto >= 0.58) &  Moge_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & Moge_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3 & Moge_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut<3, ])
#STRONG  Expressed REL Allele Only
nrow(Moge_Binders_Ic50_NoExprFilt_ord_homo[(Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE2_kallisto >=0.58 | Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE1_kallisto >= 0.58) &  Moge_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & Moge_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 & Moge_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut<3 , ])
#STRONG  Expressed Common DIAG REL
nrow(Moge_Binders_Ic50_NoExprFilt_ord_homo[(Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE2_kallisto >=0.58 | Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE1_kallisto >= 0.58) &  Moge_Binders_Ic50_NoExprFilt_ord_homo$value<=50 & (Moge_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 & Moge_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele Diagnosis
nrow(Moge_Binders_Ic50_NoExprFilt_ord_homo[(Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE2_kallisto >=0.58 | Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE1_kallisto >= 0.58) &  Moge_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & Moge_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3 , ])
#STRONG  Expressed REL Allele Relapse
nrow(Moge_Binders_Ic50_NoExprFilt_ord_homo[(Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE2_kallisto >=0.58 | Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE1_kallisto >= 0.58) &  Moge_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & Moge_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3, ])
#WEAK Expressed DIAG Allele ONLY
nrow(Moge_Binders_Ic50_NoExprFilt_ord_homo[(Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE2_kallisto >=0.58 | Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE1_kallisto >= 0.58) &  Moge_Binders_Ic50_NoExprFilt_ord_homo$value>50  & Moge_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut >= 3 & Moge_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut < 3, ])
#WEAK  Expressed REL Allele Only
nrow(Moge_Binders_Ic50_NoExprFilt_ord_homo[(Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE2_kallisto >=0.58 | Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE1_kallisto >= 0.58) &  Moge_Binders_Ic50_NoExprFilt_ord_homo$value>50  & Moge_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut >= 3 & Moge_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut < 3 , ])
#WEAK  Expressed Common DIAG REL
nrow(Moge_Binders_Ic50_NoExprFilt_ord_homo[(Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE2_kallisto >=0.58 | Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE1_kallisto >= 0.58) &  Moge_Binders_Ic50_NoExprFilt_ord_homo$value>50 & (Moge_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut >= 3 & Moge_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut >= 3), ])
#WEAK Expressed DIAG Allele Diagnosis
nrow(Moge_Binders_Ic50_NoExprFilt_ord_homo[(Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE2_kallisto >=0.58 | Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE1_kallisto >= 0.58) &  Moge_Binders_Ic50_NoExprFilt_ord_homo$value>50  & Moge_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut >= 3 , ])
#WEAK  Expressed REL Allele Relapse
nrow(Moge_Binders_Ic50_NoExprFilt_ord_homo[(Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE2_kallisto >=0.58 | Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE1_kallisto >= 0.58) &  Moge_Binders_Ic50_NoExprFilt_ord_homo$value>50  & Moge_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut >= 3, ])
#STRONG Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Moge_Binders_Ic50_NoExprFilt_ord_homo[(Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE2_kallisto >=0.58 | Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE1_kallisto >= 0.58) &  Moge_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & (Moge_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Moge_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3) & Moge_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut > Moge_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut , ])
nrow(Moge_Binders_Ic50_NoExprFilt_ord_homo[(Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE2_kallisto >=0.58 | Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE1_kallisto >= 0.58) &  Moge_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & (Moge_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Moge_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3) & Moge_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut < Moge_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut , ])
#WEAK Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Moge_Binders_Ic50_NoExprFilt_ord_homo[(Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE2_kallisto >=0.58 | Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE1_kallisto >= 0.58) &  Moge_Binders_Ic50_NoExprFilt_ord_homo$value>50  & (Moge_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Moge_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3) & Moge_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut > Moge_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut , ])
nrow(Moge_Binders_Ic50_NoExprFilt_ord_homo[(Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE2_kallisto >=0.58 | Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE1_kallisto >= 0.58) &  Moge_Binders_Ic50_NoExprFilt_ord_homo$value>50  & (Moge_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Moge_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3) & Moge_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut < Moge_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut , ])
#STRONG   Trascript  Expr_Relapse VS Expr_Diagnosis	
nrow(Moge_Binders_Ic50_NoExprFilt_ord_homo[(Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE2_kallisto >=0.58 | Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE1_kallisto >= 0.58) &  Moge_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & (Moge_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Moge_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3)  & Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE2_kallisto  > Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE1_kallisto , ])
nrow(Moge_Binders_Ic50_NoExprFilt_ord_homo[(Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE2_kallisto >=0.58 | Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE1_kallisto >= 0.58) &  Moge_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & (Moge_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Moge_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3)  & Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE2_kallisto < Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE1_kallisto , ])
#WEAK. Transcript  Expr_Relapse VS Expr_Diagnosis	
nrow(Moge_Binders_Ic50_NoExprFilt_ord_homo[(Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE2_kallisto >=0.58 | Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE1_kallisto >= 0.58) &  Moge_Binders_Ic50_NoExprFilt_ord_homo$value>50  & (Moge_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Moge_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3)  & Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE2_kallisto > Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE1_kallisto , ])
nrow(Moge_Binders_Ic50_NoExprFilt_ord_homo[(Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE2_kallisto >=0.58 | Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE1_kallisto >= 0.58) &  Moge_Binders_Ic50_NoExprFilt_ord_homo$value>50  & (Moge_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Moge_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3)  & Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE2_kallisto < Moge_Binders_Ic50_NoExprFilt_ord_homo$MOGE1_kallisto , ])



#PIAG
#N° MhAGs Strong
nrow(Piag_Binders_Ic50_NoExprFilt_ord_homo[Piag_Binders_Ic50_NoExprFilt_ord_homo$value<=50,])
#N° MhAGs Weak
nrow(Piag_Binders_Ic50_NoExprFilt_ord_homo[Piag_Binders_Ic50_NoExprFilt_ord_homo$value>50,])
#N° Expr. MhAGs  Strong
nrow(Piag_Binders_Ic50_NoExprFilt_ord_homo[(Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG2_kallisto >=0.58 | Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG1_kallisto >= 0.58) &  Piag_Binders_Ic50_NoExprFilt_ord_homo$value<=50 & (Piag_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Piag_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3), ])
#N° Expr. MhAGs  Weak
nrow(Piag_Binders_Ic50_NoExprFilt_ord_homo[(Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG2_kallisto >=0.58 | Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG1_kallisto >= 0.58) &  Piag_Binders_Ic50_NoExprFilt_ord_homo$value>50 & (Piag_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Piag_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele ONLY
nrow(Piag_Binders_Ic50_NoExprFilt_ord_homo[(Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG2_kallisto >=0.58 | Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG1_kallisto >= 0.58) &  Piag_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & Piag_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3 & Piag_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut<3, ])
#STRONG  Expressed REL Allele Only
nrow(Piag_Binders_Ic50_NoExprFilt_ord_homo[(Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG2_kallisto >=0.58 | Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG1_kallisto >= 0.58) &  Piag_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & Piag_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 & Piag_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut<3 , ])
#STRONG  Expressed Common DIAG REL
nrow(Piag_Binders_Ic50_NoExprFilt_ord_homo[(Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG2_kallisto >=0.58 | Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG1_kallisto >= 0.58) &  Piag_Binders_Ic50_NoExprFilt_ord_homo$value<=50 & (Piag_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 & Piag_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele Diagnosis
nrow(Piag_Binders_Ic50_NoExprFilt_ord_homo[(Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG2_kallisto >=0.58 | Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG1_kallisto >= 0.58) &  Piag_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & Piag_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3 , ])
#STRONG  Expressed REL Allele Relapse
nrow(Piag_Binders_Ic50_NoExprFilt_ord_homo[(Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG2_kallisto >=0.58 | Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG1_kallisto >= 0.58) &  Piag_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & Piag_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3, ])
#WEAK Expressed DIAG Allele ONLY
nrow(Piag_Binders_Ic50_NoExprFilt_ord_homo[(Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG2_kallisto >=0.58 | Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG1_kallisto >= 0.58) &  Piag_Binders_Ic50_NoExprFilt_ord_homo$value>50  & Piag_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut >= 3 & Piag_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut < 3, ])
#WEAK  Expressed REL Allele Only
nrow(Piag_Binders_Ic50_NoExprFilt_ord_homo[(Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG2_kallisto >=0.58 | Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG1_kallisto >= 0.58) &  Piag_Binders_Ic50_NoExprFilt_ord_homo$value>50  & Piag_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut >= 3 & Piag_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut < 3 , ])
#WEAK  Expressed Common DIAG REL
nrow(Piag_Binders_Ic50_NoExprFilt_ord_homo[(Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG2_kallisto >=0.58 | Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG1_kallisto >= 0.58) &  Piag_Binders_Ic50_NoExprFilt_ord_homo$value>50 & (Piag_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut >= 3 & Piag_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut >= 3), ])
#WEAK Expressed DIAG Allele Diagnosis
nrow(Piag_Binders_Ic50_NoExprFilt_ord_homo[(Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG2_kallisto >=0.58 | Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG1_kallisto >= 0.58) &  Piag_Binders_Ic50_NoExprFilt_ord_homo$value>50  & Piag_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut >= 3 , ])
#WEAK  Expressed REL Allele Relapse
nrow(Piag_Binders_Ic50_NoExprFilt_ord_homo[(Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG2_kallisto >=0.58 | Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG1_kallisto >= 0.58) &  Piag_Binders_Ic50_NoExprFilt_ord_homo$value>50  & Piag_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut >= 3, ])
#STRONG Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Piag_Binders_Ic50_NoExprFilt_ord_homo[(Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG2_kallisto >=0.58 | Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG1_kallisto >= 0.58) &  Piag_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & (Piag_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Piag_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3) & Piag_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut > Piag_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut , ])
nrow(Piag_Binders_Ic50_NoExprFilt_ord_homo[(Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG2_kallisto >=0.58 | Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG1_kallisto >= 0.58) &  Piag_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & (Piag_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Piag_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3) & Piag_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut < Piag_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut , ])
#WEAK Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Piag_Binders_Ic50_NoExprFilt_ord_homo[(Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG2_kallisto >=0.58 | Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG1_kallisto >= 0.58) &  Piag_Binders_Ic50_NoExprFilt_ord_homo$value>50  & (Piag_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Piag_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3) & Piag_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut > Piag_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut , ])
nrow(Piag_Binders_Ic50_NoExprFilt_ord_homo[(Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG2_kallisto >=0.58 | Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG1_kallisto >= 0.58) &  Piag_Binders_Ic50_NoExprFilt_ord_homo$value>50  & (Piag_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Piag_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3) & Piag_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut < Piag_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut , ])
#STRONG   Trascript  Expr_Relapse VS Expr_Diagnosis	
nrow(Piag_Binders_Ic50_NoExprFilt_ord_homo[(Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG2_kallisto >=0.58 | Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG1_kallisto >= 0.58) &  Piag_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & (Piag_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Piag_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3)  & Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG2_kallisto  > Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG1_kallisto , ])
nrow(Piag_Binders_Ic50_NoExprFilt_ord_homo[(Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG2_kallisto >=0.58 | Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG1_kallisto >= 0.58) &  Piag_Binders_Ic50_NoExprFilt_ord_homo$value<=50  & (Piag_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Piag_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3)  & Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG2_kallisto < Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG1_kallisto , ])
#WEAK. Transcript  Expr_Relapse VS Expr_Diagnosis	
nrow(Piag_Binders_Ic50_NoExprFilt_ord_homo[(Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG2_kallisto >=0.58 | Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG1_kallisto >= 0.58) &  Piag_Binders_Ic50_NoExprFilt_ord_homo$value>50  & (Piag_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Piag_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3)  & Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG2_kallisto > Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG1_kallisto , ])
nrow(Piag_Binders_Ic50_NoExprFilt_ord_homo[(Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG2_kallisto >=0.58 | Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG1_kallisto >= 0.58) &  Piag_Binders_Ic50_NoExprFilt_ord_homo$value>50  & (Piag_Binders_Ic50_NoExprFilt_ord_homo$Rel_Mut>=3 | Piag_Binders_Ic50_NoExprFilt_ord_homo$Dx_Mut>=3)  & Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG2_kallisto < Piag_Binders_Ic50_NoExprFilt_ord_homo$PIAG1_kallisto , ])






#PRELU
#N° MhAGs Strong
nrow(Prelu_Binders_Ic50_NoExprFilt_ord[Prelu_Binders_Ic50_NoExprFilt_ord$value<=50,])
#N° MhAGs Weak
nrow(Prelu_Binders_Ic50_NoExprFilt_ord[Prelu_Binders_Ic50_NoExprFilt_ord$value>50,])
#N° Expr. MhAGs  Strong
nrow(Prelu_Binders_Ic50_NoExprFilt_ord[(Prelu_Binders_Ic50_NoExprFilt_ord$PRELU2_kallisto >=0.58 |  Prelu_Binders_Ic50_NoExprFilt_ord$PRELU1_kallisto>=0.58) &  Prelu_Binders_Ic50_NoExprFilt_ord$value<=50 & (Prelu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Prelu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#N° Expr. MhAGs  Weak
nrow(Prelu_Binders_Ic50_NoExprFilt_ord[(Prelu_Binders_Ic50_NoExprFilt_ord$PRELU2_kallisto >=0.58 |  Prelu_Binders_Ic50_NoExprFilt_ord$PRELU1_kallisto>=0.58) &  Prelu_Binders_Ic50_NoExprFilt_ord$value>50 & (Prelu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Prelu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele ONLY
nrow(Prelu_Binders_Ic50_NoExprFilt_ord[(Prelu_Binders_Ic50_NoExprFilt_ord$PRELU2_kallisto >=0.58 |  Prelu_Binders_Ic50_NoExprFilt_ord$PRELU1_kallisto>=0.58) &  Prelu_Binders_Ic50_NoExprFilt_ord$value<=50  & Prelu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3 & Prelu_Binders_Ic50_NoExprFilt_ord$Rel_Mut<3, ])
#STRONG  Expressed REL Allele Only
nrow(Prelu_Binders_Ic50_NoExprFilt_ord[(Prelu_Binders_Ic50_NoExprFilt_ord$PRELU2_kallisto >=0.58 |  Prelu_Binders_Ic50_NoExprFilt_ord$PRELU1_kallisto>=0.58) &  Prelu_Binders_Ic50_NoExprFilt_ord$value<=50  & Prelu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 & Prelu_Binders_Ic50_NoExprFilt_ord$Dx_Mut<3 , ])
#STRONG  Expressed Common DIAG REL
nrow(Prelu_Binders_Ic50_NoExprFilt_ord[(Prelu_Binders_Ic50_NoExprFilt_ord$PRELU2_kallisto >=0.58 |  Prelu_Binders_Ic50_NoExprFilt_ord$PRELU1_kallisto>=0.58) &  Prelu_Binders_Ic50_NoExprFilt_ord$value<=50 & (Prelu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 & Prelu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3), ])
#STRONG Expressed DIAG Allele Diagnosis
nrow(Prelu_Binders_Ic50_NoExprFilt_ord[(Prelu_Binders_Ic50_NoExprFilt_ord$PRELU2_kallisto >=0.58 |  Prelu_Binders_Ic50_NoExprFilt_ord$PRELU1_kallisto>=0.58) &  Prelu_Binders_Ic50_NoExprFilt_ord$value<=50  & Prelu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3 , ])
#STRONG  Expressed REL Allele Relapse
nrow(Prelu_Binders_Ic50_NoExprFilt_ord[(Prelu_Binders_Ic50_NoExprFilt_ord$PRELU2_kallisto >=0.58 |  Prelu_Binders_Ic50_NoExprFilt_ord$PRELU1_kallisto>=0.58) &  Prelu_Binders_Ic50_NoExprFilt_ord$value<=50  & Prelu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3, ])
#WEAK Expressed DIAG Allele ONLY
nrow(Prelu_Binders_Ic50_NoExprFilt_ord[(Prelu_Binders_Ic50_NoExprFilt_ord$PRELU2_kallisto >=0.58 |  Prelu_Binders_Ic50_NoExprFilt_ord$PRELU1_kallisto>=0.58) &  Prelu_Binders_Ic50_NoExprFilt_ord$value>50  & Prelu_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3 & Prelu_Binders_Ic50_NoExprFilt_ord$Rel_Mut < 3, ])
#WEAK  Expressed REL Allele Only
nrow(Prelu_Binders_Ic50_NoExprFilt_ord[(Prelu_Binders_Ic50_NoExprFilt_ord$PRELU2_kallisto >=0.58 |  Prelu_Binders_Ic50_NoExprFilt_ord$PRELU1_kallisto>=0.58) &  Prelu_Binders_Ic50_NoExprFilt_ord$value>50  & Prelu_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3 & Prelu_Binders_Ic50_NoExprFilt_ord$Dx_Mut < 3 , ])
#WEAK  Expressed Common DIAG REL
nrow(Prelu_Binders_Ic50_NoExprFilt_ord[(Prelu_Binders_Ic50_NoExprFilt_ord$PRELU2_kallisto >=0.58 |  Prelu_Binders_Ic50_NoExprFilt_ord$PRELU1_kallisto>=0.58) &  Prelu_Binders_Ic50_NoExprFilt_ord$value>50 & (Prelu_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3 & Prelu_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3), ])
#WEAK Expressed DIAG Allele Diagnosis
nrow(Prelu_Binders_Ic50_NoExprFilt_ord[(Prelu_Binders_Ic50_NoExprFilt_ord$PRELU2_kallisto >=0.58 |  Prelu_Binders_Ic50_NoExprFilt_ord$PRELU1_kallisto>=0.58) &  Prelu_Binders_Ic50_NoExprFilt_ord$value>50  & Prelu_Binders_Ic50_NoExprFilt_ord$Dx_Mut >= 3 , ])
#WEAK  Expressed REL Allele Relapse
nrow(Prelu_Binders_Ic50_NoExprFilt_ord[(Prelu_Binders_Ic50_NoExprFilt_ord$PRELU2_kallisto >=0.58 |  Prelu_Binders_Ic50_NoExprFilt_ord$PRELU1_kallisto>=0.58) &  Prelu_Binders_Ic50_NoExprFilt_ord$value>50  & Prelu_Binders_Ic50_NoExprFilt_ord$Rel_Mut >= 3, ])
#STRONG Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Prelu_Binders_Ic50_NoExprFilt_ord[(Prelu_Binders_Ic50_NoExprFilt_ord$PRELU2_kallisto >=0.58 |  Prelu_Binders_Ic50_NoExprFilt_ord$PRELU1_kallisto>=0.58) &  Prelu_Binders_Ic50_NoExprFilt_ord$value<=50  & (Prelu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Prelu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Prelu_Binders_Ic50_NoExprFilt_ord$Rel_Mut > Prelu_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
nrow(Prelu_Binders_Ic50_NoExprFilt_ord[(Prelu_Binders_Ic50_NoExprFilt_ord$PRELU2_kallisto >=0.58 |  Prelu_Binders_Ic50_NoExprFilt_ord$PRELU1_kallisto>=0.58) &  Prelu_Binders_Ic50_NoExprFilt_ord$value<=50  & (Prelu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Prelu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Prelu_Binders_Ic50_NoExprFilt_ord$Rel_Mut < Prelu_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
#WEAK Allele  Expr_Relapse VS Expr_Diagnosis	
nrow(Prelu_Binders_Ic50_NoExprFilt_ord[(Prelu_Binders_Ic50_NoExprFilt_ord$PRELU2_kallisto >=0.58 |  Prelu_Binders_Ic50_NoExprFilt_ord$PRELU1_kallisto>=0.58) &  Prelu_Binders_Ic50_NoExprFilt_ord$value>50  & (Prelu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Prelu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Prelu_Binders_Ic50_NoExprFilt_ord$Rel_Mut > Prelu_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
nrow(Prelu_Binders_Ic50_NoExprFilt_ord[(Prelu_Binders_Ic50_NoExprFilt_ord$PRELU2_kallisto >=0.58 |  Prelu_Binders_Ic50_NoExprFilt_ord$PRELU1_kallisto>=0.58) &  Prelu_Binders_Ic50_NoExprFilt_ord$value>50  & (Prelu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Prelu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3) & Prelu_Binders_Ic50_NoExprFilt_ord$Rel_Mut < Prelu_Binders_Ic50_NoExprFilt_ord$Dx_Mut , ])
#STRONG   Trascript  Expr_Relapse VS Expr_Diagnosis	
nrow(Prelu_Binders_Ic50_NoExprFilt_ord[(Prelu_Binders_Ic50_NoExprFilt_ord$PRELU2_kallisto >=0.58 |  Prelu_Binders_Ic50_NoExprFilt_ord$PRELU1_kallisto>=0.58) &  Prelu_Binders_Ic50_NoExprFilt_ord$value<=50  & (Prelu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Prelu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Prelu_Binders_Ic50_NoExprFilt_ord$PRELU2_kallisto  > Prelu_Binders_Ic50_NoExprFilt_ord$PRELU1_kallisto , ])
nrow(Prelu_Binders_Ic50_NoExprFilt_ord[(Prelu_Binders_Ic50_NoExprFilt_ord$PRELU2_kallisto >=0.58 |  Prelu_Binders_Ic50_NoExprFilt_ord$PRELU1_kallisto>=0.58) &  Prelu_Binders_Ic50_NoExprFilt_ord$value<=50  & (Prelu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Prelu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Prelu_Binders_Ic50_NoExprFilt_ord$PRELU2_kallisto < Prelu_Binders_Ic50_NoExprFilt_ord$PRELU1_kallisto , ])
#WEAK. Transcript  Expr_Relapse VS Expr_Diagnosis	
nrow(Prelu_Binders_Ic50_NoExprFilt_ord[(Prelu_Binders_Ic50_NoExprFilt_ord$PRELU2_kallisto >=0.58 |  Prelu_Binders_Ic50_NoExprFilt_ord$PRELU1_kallisto>=0.58) &  Prelu_Binders_Ic50_NoExprFilt_ord$value>50  & (Prelu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Prelu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Prelu_Binders_Ic50_NoExprFilt_ord$PRELU2_kallisto > Prelu_Binders_Ic50_NoExprFilt_ord$PRELU1_kallisto , ])
nrow(Prelu_Binders_Ic50_NoExprFilt_ord[(Prelu_Binders_Ic50_NoExprFilt_ord$PRELU2_kallisto >=0.58 |  Prelu_Binders_Ic50_NoExprFilt_ord$PRELU1_kallisto>=0.58) &  Prelu_Binders_Ic50_NoExprFilt_ord$value>50  & (Prelu_Binders_Ic50_NoExprFilt_ord$Rel_Mut>=3 | Prelu_Binders_Ic50_NoExprFilt_ord$Dx_Mut>=3)  & Prelu_Binders_Ic50_NoExprFilt_ord$PRELU2_kallisto < Prelu_Binders_Ic50_NoExprFilt_ord$PRELU1_kallisto , ])








pdf('../../Documents/Progetto_LuciaGabri/Plot_AUCS_Patients_NoExprFilt.pdf')

ggplot(all_aucs,aes(x=seq(1,length(all_aucs$aucs)),y=aucs,col=Donors)) +
geom_point()  +
geom_rug() +
scale_color_manual(values=c( Related= "#E69F00", Unrelated="darkorchid2")) +
geom_smooth( method = "lm",mapping = aes(fill=all_aucs$Donors)) +
scale_fill_manual(values=c(Related= "#E69F00", Unrelated="darkorchid2")) +
ggtitle('MRD vs MUD AUCs Binders\n linear regression') +
xlab('Patients')


summary(glm(factor(all_aucs$Donors) ~ all_aucs$aucs,family = binomial(link = "probit")))

```
Call:
glm(formula = factor(all_aucs$Donors) ~ all_aucs$aucs, family = binomial(link = "probit"))

Deviance Residuals: 
	Min       1Q   Median       3Q      Max  
-1.5010  -0.6655  -0.1185   0.2412   1.8148  

Coefficients:
				Estimate Std. Error z value Pr(>|z|)  
(Intercept)   -3.056e+00  1.510e+00  -2.024    0.043 *
all_aucs$aucs  1.673e-05  8.897e-06   1.881    0.060 .
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

	Null deviance: 20.728  on 14  degrees of freedom
Residual deviance: 11.140  on 13  degrees of freedom
AIC: 15.14

Number of Fisher Scoring iterations: 7
```
plot.new()
color=c25[c(1:11,13,18,22,24)]
plot(lines(loess.smooth(1:length(Alfe_Binders_Ic50_NoExprFilt_ord$value), y=Alfe_Binders_Ic50_NoExprFilt_ord$value)),xlim=c(0,3000),ylim=c(0,500),col=color[1],xlab='Number of Predicted mhAGs Binders',ylab='IC50 (nM)')
lines(loess.smooth(1:length(Alfe_Binders_Ic50_NoExprFilt_ord$value), y=Alfe_Binders_Ic50_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[1],xlab='',ylab='IC50 (nM)',lwd=3)
#polygon(1:length(Besu_Binders_ord$value), y=Besu_Binders_ord$value,col='skyblue')
par(new=TRUE)
lines(loess.smooth(1:length(Calu_Binders_Ic50_NoExprFilt_ord_homo$value), y=Calu_Binders_Ic50_NoExprFilt_ord_homo$value,degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[2],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Calu_Binders_Ic50_NoExprFilt_ord_homo$value), y=Calu_Binders_Ic50_NoExprFilt_ord_homo$value,degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[3],xlab='',ylab='',lwd=3)
par(new=TRUE)
lines(loess.smooth(1:length(Deiv_Binders_Ic50_NoExprFilt_ord$value), y=Deiv_Binders_Ic50_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[4],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Dest_Binders_Ic50_NoExprFilt_ord$value), y=Dest_Binders_Ic50_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[5],xlab='',ylab='',lwd=3)
par(new=TRUE)
lines(loess.smooth(1:length(Dr1_Binders_Ic50_NoExprFilt_ord$value), y=Dr1_Binders_Ic50_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[6],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Dr4_Binders_Ic50_NoExprFilt_ord$value), y=Dr4_Binders_Ic50_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[7],xlab='',ylab='',lwd=3,)
par(new=TRUE)
lines(loess.smooth(1:length(Dr5_Binders_Ic50_NoExprFilt_ord$value), y=Dr5_Binders_Ic50_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[8],xlab='',ylab='',lwd=3,)
par(new=TRUE)
lines(loess.smooth(1:length(Foca_Binders_Ic50_NoExprFilt_ord$value), y=Foca_Binders_Ic50_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[9],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Gagra_Binders_Ic50_NoExprFilt_ord$value), y=Gagra_Binders_Ic50_NoExprFilt_ord$value, degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[10],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Luan_Binders_Ic50_NoExprFilt_ord_homo$value), y=Luan_Binders_Ic50_NoExprFilt_ord_homo$value,  degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[11],xlab='',ylab='',lwd=3,)
par(new=TRUE)
lines(loess.smooth(1:length(Mabi_Binders_Ic50_NoExprFilt_ord$value), y=Mabi_Binders_Ic50_NoExprFilt_ord$value,  degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[12],xlab='',ylab='',lwd=3)
par(new=TRUE)
lines(loess.smooth(1:length(Moge_Binders_Ic50_NoExprFilt_ord_homo$value), y=Moge_Binders_Ic50_NoExprFilt_ord_homo$value,  degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[13],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Piag_Binders_Ic50_NoExprFilt_ord_homo$value), y=Piag_Binders_Ic50_NoExprFilt_ord_homo$value,  degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[14],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Prelu_Binders_Ic50_NoExprFilt_ord$value), y=Prelu_Binders_Ic50_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[15],xlab='',ylab='',lwd=3)
legend("topright", inset=c(-0,0), c('alfe','besu','calu','deiv','dest','dr1','dr4','dr5','foca','gagra','luan','mabi','moge','piag','prelu'),col=color[1:15], title="Patients",cex=0.7,lwd=rep(2,15),lty=c(1,2,1,2,1,2,1,1,2,2,1,1,2,2,1))
dev.off()



write.table(Alfe_Binders_Ic50_NoExprFilt_ord,'../../Documents/Progetto_LuciaGabri/Alfe_Binders_allInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)
write.table(Calu_Binders_Ic50_NoExprFilt_ord_homo,'../../Documents/Progetto_LuciaGabri/Besu_Binders_allInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)
write.table(Calu_Binders_Ic50_NoExprFilt_ord_homo,'../../Documents/Progetto_LuciaGabri/Calu_Binders_allInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)
write.table(Deiv_Binders_Ic50_NoExprFilt_ord,'../../Documents/Progetto_LuciaGabri/Deiv_Binders_allInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)
write.table(Dest_Binders_Ic50_NoExprFilt_ord,'../../Documents/Progetto_LuciaGabri/Dest_Binders_allInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)
write.table(Dr1_Binders_Ic50_NoExprFilt_ord,'../../Documents/Progetto_LuciaGabri/Dr1_Binders_allInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)
write.table(Dr4_Binders_Ic50_NoExprFilt_ord,'../../Documents/Progetto_LuciaGabri/Dr4_Binders_allInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)
write.table(Dr5_Binders_Ic50_NoExprFilt_ord,'../../Documents/Progetto_LuciaGabri/Dr5_Binders_allInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)
write.table(Foca_Binders_Ic50_NoExprFilt_ord,'../../Documents/Progetto_LuciaGabri/Foca_Binders_allInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)
write.table(Gagra_Binders_Ic50_NoExprFilt_ord,'../../Documents/Progetto_LuciaGabri/Gagra_Binders_allInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)
write.table(Luan_Binders_Ic50_NoExprFilt_ord_homo,'../../Documents/Progetto_LuciaGabri/Luan_Binders_allInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)
write.table(Mabi_Binders_Ic50_NoExprFilt_ord,'../../Documents/Progetto_LuciaGabri/Mabi_Binders_allInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)
write.table(Moge_Binders_Ic50_NoExprFilt_ord_homo,'../../Documents/Progetto_LuciaGabri/Moge_Binders_allInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)
write.table(Piag_Binders_Ic50_NoExprFilt_ord_homo,'../../Documents/Progetto_LuciaGabri/Piag_Binders_allInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)
write.table(Prelu_Binders_Ic50_NoExprFilt_ord,'../../Documents/Progetto_LuciaGabri/Prelu_Binders_allInfos.txt',sep='\t',col.names=T,row.names=F,quote=F)








alfe_auc_dx=auc(1:length(Alfe_Binders_Ic50_NoExprFilt_Dx$value), y=Alfe_Binders_Ic50_NoExprFilt_Dx$value)
besu_auc_dx=auc(1:length(Besu_Binders_Ic50_NoExprFilt_Dx$value), y=Besu_Binders_Ic50_NoExprFilt_Dx$value)
calu_auc_dx=auc(1:length(Calu_Binders_Ic50_NoExprFilt_Dx$value), y=Calu_Binders_Ic50_NoExprFilt_Dx$value)
deiv_auc_dx=auc(1:length(Deiv_Binders_Ic50_NoExprFilt_Dx$value), y=Deiv_Binders_Ic50_NoExprFilt_Dx$value)
dest_auc_dx=auc(1:length(Dest_Binders_Ic50_NoExprFilt_Dx$value), y=Dest_Binders_Ic50_NoExprFilt_Dx$value)
dr1_auc_dx=auc(1:length(Dr1_Binders_Ic50_NoExprFilt_Dx$value), y=Dr1_Binders_Ic50_NoExprFilt_Dx$value)
dr4_auc_dx=auc(1:length(Dr4_Binders_Ic50_NoExprFilt_Dx$value), y=Dr4_Binders_Ic50_NoExprFilt_Dx$value)
dr5_auc_dx=auc(1:length(Dr5_Binders_Ic50_NoExprFilt_Dx$value), y=Dr5_Binders_Ic50_NoExprFilt_Dx$value)
foca_auc_dx=auc(1:length(Foca_Binders_Ic50_NoExprFilt_Dx$value), y=Foca_Binders_Ic50_NoExprFilt_Dx$value)
gagra_auc_dx=auc(1:length(Gagra_Binders_Ic50_NoExprFilt_Dx$value), y=Gagra_Binders_Ic50_NoExprFilt_Dx$value)
luan_auc_dx=auc(1:length(Luan_Binders_Ic50_NoExprFilt_Dx$value), y=Luan_Binders_Ic50_NoExprFilt_Dx$value)
mabi_auc_dx=auc(1:length(Mabi_Binders_Ic50_NoExprFilt_Dx$value), y=Mabi_Binders_Ic50_NoExprFilt_Dx$value)
moge_auc_dx=auc(1:length(Moge_Binders_Ic50_NoExprFilt_Dx$value), y=Moge_Binders_Ic50_NoExprFilt_Dx$value)
piag_auc_dx=auc(1:length(Piag_Binders_Ic50_NoExprFilt_Dx$value), y=Piag_Binders_Ic50_NoExprFilt_Dx$value)
prelu_auc_dx=auc(1:length(Prelu_Binders_Ic50_NoExprFilt_Dx$value), y=Prelu_Binders_Ic50_NoExprFilt_Dx$value)
			

all_aucs_Dx=rbind(alfe_auc_Dx,besu_auc_Dx,calu_auc_Dx,deiv_auc_Dx,dest_auc_Dx,dr1_auc_Dx,dr4_auc_Dx,dr5_auc_Dx,foca_auc_Dx,gagra_auc_Dx,luan_auc_Dx,mabi_auc_Dx,moge_auc_Dx,piag_auc_Dx,prelu_auc)
all_aucs_Dx=as.data.frame(all_aucs_Dx)
all_aucs_Dx$genes=rownames(all_aucs_Dx)
all_aucs_Dx$RelDon=c(1,0,1,0,1,0,1,1,0,0,1,1,0,0,1)
#all_aucs_Dx$RelDon=c(0,1,0,1,0,1,0,0,1,1,0,0,1,1,0)
all_aucs_Dx=all_aucs_Dx[order(all_aucs_Dx$V1),]
colnames(all_aucs_Dx)=c('aucs_Dx','patients','RelDon')
all_aucs_Dx$Donors=ifelse(all_aucs_Dx$RelDon==1,'Related','Unrelated')

write.table(all_aucs_Dx,'../../Documents/Progetto_LuciaGabri/MhAGs_AUC_Pazienti_DiagnosisExpressionFiltered.txt',sep='\t',col.names=T,row.names=F,quote=F)

pdf('../../Documents/Progetto_LuciaGabri/Plot_AUCS_Patients_DiagnosisExprFilt.pdf')

ggplot(all_aucs_Dx,aes(x=seq(1,length(all_aucs_Dx$aucs_Dx)),y=aucs_Dx,col=Donors)) +
geom_point()  +
geom_rug() +
scale_color_manual(values=c( Related= "#E69F00", Unrelated="darkorchid2")) +
geom_smooth( method = "lm",mapping = aes(fill=all_aucs_Dx$Donors)) +
scale_fill_manual(values=c(Related= "#E69F00", Unrelated="darkorchid2")) +
ggtitle('MRD vs MUD AUCs Binders\n linear regression') +
xlab('Patients')


plot(log(all_aucs_Dx$aucs_Dx),col=as.factor(all_aucs_Dx$Donor),lwd=3)
text(log(all_aucs_Dx$aucs_Dx),labels=all_aucs_Dx$patients, cex= 0.7,adj=1,pos=3,col=color[1:15])

summary(glm(factor(all_aucs_Dx$Donors) ~ all_aucs_Dx$aucs_Dx,family = binomial(link = "probit")))
			
Call:
glm(formula = factor(all_aucs_Dx$Donors) ~ all_aucs_Dx$aucs_Dx, 
    family = binomial(link = "probit"))

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-1.3628  -0.6074  -0.1389   0.1006   1.8166  

Coefficients:
                      Estimate Std. Error z value Pr(>|z|)  
(Intercept)         -3.030e+00  1.501e+00  -2.018   0.0435 *
all_aucs_Dx$aucs_Dx  5.822e-05  3.189e-05   1.825   0.0679 .
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 20.728  on 14  degrees of freedom
Residual deviance: 10.417  on 13  degrees of freedom
AIC: 14.417

Number of Fisher Scoring iterations: 7
plot.new()
color=c25[c(1:11,13,18,22,24)]
plot(lines(loess.smooth(1:length(Alfe_Binders_Ic50_NoExprFilt_Dx$value), y=Alfe_Binders_Ic50_NoExprFilt_Dx$value)),xlim=c(0,800),ylim=c(0,500),col=color[1],xlab='Number of Predicted mhAGs Binders',ylab='IC50 (nM)')
lines(loess.smooth(1:length(Alfe_Binders_Ic50_NoExprFilt_Dx$value), y=Alfe_Binders_Ic50_NoExprFilt_Dx$value,degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[1],xlab='',ylab='IC50 (nM)',lwd=3)
#polygon(1:length(Besu_Binders_ord$value), y=Besu_Binders_ord$value,col='skyblue')
par(new=TRUE)
lines(loess.smooth(1:length(Besu_Binders_Ic50_NoExprFilt_Dx$value), y=Besu_Binders_Ic50_NoExprFilt_Dx$value,degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[2],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Calu_Binders_Ic50_NoExprFilt_Dx$value), y=Calu_Binders_Ic50_NoExprFilt_Dx$value,degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[3],xlab='',ylab='',lwd=3)
par(new=TRUE)
lines(loess.smooth(1:length(Deiv_Binders_Ic50_NoExprFilt_Dx$value), y=Deiv_Binders_Ic50_NoExprFilt_Dx$value,degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[4],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Dest_Binders_Ic50_NoExprFilt_Dx$value), y=Dest_Binders_Ic50_NoExprFilt_Dx$value,degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[5],xlab='',ylab='',lwd=3)
par(new=TRUE)
lines(loess.smooth(1:length(Dr1_Binders_Ic50_NoExprFilt_Dx$value), y=Dr1_Binders_Ic50_NoExprFilt_Dx$value,degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[6],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Dr4_Binders_Ic50_NoExprFilt_Dx$value), y=Dr4_Binders_Ic50_NoExprFilt_Dx$value,degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[7],xlab='',ylab='',lwd=3,)
par(new=TRUE)
lines(loess.smooth(1:length(Dr5_Binders_Ic50_NoExprFilt_Dx$value), y=Dr5_Binders_Ic50_NoExprFilt_Dx$value,degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[8],xlab='',ylab='',lwd=3,)
par(new=TRUE)
lines(loess.smooth(1:length(Foca_Binders_Ic50_NoExprFilt_Dx$value), y=Foca_Binders_Ic50_NoExprFilt_Dx$value,degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[9],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Gagra_Binders_Ic50_NoExprFilt_Dx$value), y=Gagra_Binders_Ic50_NoExprFilt_Dx$value, degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[10],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Luan_Binders_Ic50_NoExprFilt_Dx$value), y=Luan_Binders_Ic50_NoExprFilt_Dx$value,  degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[11],xlab='',ylab='',lwd=3,)
par(new=TRUE)
lines(loess.smooth(1:length(Mabi_Binders_Ic50_NoExprFilt_Dx$value), y=Mabi_Binders_Ic50_NoExprFilt_Dx$value,  degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[12],xlab='',ylab='',lwd=3)
par(new=TRUE)
lines(loess.smooth(1:length(Moge_Binders_Ic50_NoExprFilt_Dx$value), y=Moge_Binders_Ic50_NoExprFilt_Dx$value,  degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[13],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Piag_Binders_Ic50_NoExprFilt_Dx$value), y=Piag_Binders_Ic50_NoExprFilt_Dx$value,  degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[14],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Prelu_Binders_Ic50_NoExprFilt_Dx$value), y=Prelu_Binders_Ic50_NoExprFilt_Dx$value,degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[15],xlab='',ylab='',lwd=3)
legend("topright", inset=c(-0,0), c('alfe','besu','calu','deiv','dest','dr1','dr4','dr5','foca','gagra','luan','mabi','moge','piag','prelu'),col=color[1:15], title="Patients",cex=0.7,lwd=rep(2,15),lty=c(1,2,1,2,1,2,1,1,2,2,1,1,2,2,1))






alfe_auc_rel=auc(1:length(Alfe_Binders_Ic50_NoExprFilt_Rel$value), y=Alfe_Binders_Ic50_NoExprFilt_Rel$value)
besu_auc_rel=auc(1:length(Besu_Binders_Ic50_NoExprFilt_Rel$value), y=Besu_Binders_Ic50_NoExprFilt_Rel$value)
calu_auc_rel=auc(1:length(Calu_Binders_Ic50_NoExprFilt_Rel$value), y=Calu_Binders_Ic50_NoExprFilt_Rel$value)
deiv_auc_rel=auc(1:length(Deiv_Binders_Ic50_NoExprFilt_Rel$value), y=Deiv_Binders_Ic50_NoExprFilt_Rel$value)
dest_auc_rel=auc(1:length(Dest_Binders_Ic50_NoExprFilt_Rel$value), y=Dest_Binders_Ic50_NoExprFilt_Rel$value)
dr1_auc_rel=auc(1:length(Dr1_Binders_Ic50_NoExprFilt_Rel$value), y=Dr1_Binders_Ic50_NoExprFilt_Rel$value)
dr4_auc_rel=auc(1:length(Dr4_Binders_Ic50_NoExprFilt_Rel$value), y=Dr4_Binders_Ic50_NoExprFilt_Rel$value)
dr5_auc_rel=auc(1:length(Dr5_Binders_Ic50_NoExprFilt_Rel$value), y=Dr5_Binders_Ic50_NoExprFilt_Rel$value)
foca_auc_rel=auc(1:length(Foca_Binders_Ic50_NoExprFilt_Rel$value), y=Foca_Binders_Ic50_NoExprFilt_Rel$value)
gagra_auc_rel=auc(1:length(Gagra_Binders_Ic50_NoExprFilt_Rel$value), y=Gagra_Binders_Ic50_NoExprFilt_Rel$value)
luan_auc_rel=auc(1:length(Luan_Binders_Ic50_NoExprFilt_Rel$value), y=Luan_Binders_Ic50_NoExprFilt_Rel$value)
mabi_auc_rel=auc(1:length(Mabi_Binders_Ic50_NoExprFilt_Rel$value), y=Mabi_Binders_Ic50_NoExprFilt_Rel$value)
moge_auc_rel=auc(1:length(Moge_Binders_Ic50_NoExprFilt_Rel$value), y=Moge_Binders_Ic50_NoExprFilt_Rel$value)
piag_auc_rel=auc(1:length(Piag_Binders_Ic50_NoExprFilt_Rel$value), y=Piag_Binders_Ic50_NoExprFilt_Rel$value)
prelu_auc_rel=auc(1:length(Prelu_Binders_Ic50_NoExprFilt_Rel$value), y=Prelu_Binders_Ic50_NoExprFilt_Rel$value)
			

all_aucs_Rel=rbind(alfe_auc_Rel,besu_auc_Rel,calu_auc_Rel,deiv_auc_Rel,dest_auc_Rel,dr1_auc_Rel,dr4_auc_Rel,dr5_auc_Rel,foca_auc_Rel,gagra_auc_Rel,luan_auc_Rel,mabi_auc_Rel,moge_auc_Rel,piag_auc_Rel,prelu_auc)
all_aucs_Rel=as.data.frame(all_aucs_Rel)
all_aucs_Rel$genes=rownames(all_aucs_Rel)
all_aucs_Rel$RelDon=c(1,0,1,0,1,0,1,1,0,0,1,1,0,0,1)
#all_aucs_Rel$RelDon=c(0,1,0,1,0,1,0,0,1,1,0,0,1,1,0)
all_aucs_Rel=all_aucs_Rel[order(all_aucs_Rel$V1),]
colnames(all_aucs_Rel)=c('aucs_Rel','patients','RelDon')
all_aucs_Rel$Donors=ifelse(all_aucs_Rel$RelDon==1,'Related','Unrelated')

pdf('../../Documents/Progetto_LuciaGabri/Plot_AUCS_Patients_RelapseExprFilt.pdf')

ggplot(all_aucs_Rel,aes(x=seq(1,length(all_aucs_Rel$aucs_Rel)),y=aucs_Rel,col=Donors)) +
geom_point()  +
geom_rug() +
scale_color_manual(values=c( Related= "#E69F00", Unrelated="darkorchid2")) +
geom_smooth( method = "lm",mapping = aes(fill=all_aucs_Rel$Donors)) +
scale_fill_manual(values=c(Related= "#E69F00", Unrelated="darkorchid2")) +
ggtitle('MRD vs MUD AUCs Binders\n linear regression') +
xlab('Patients')


plot(log(all_aucs_Rel$aucs_Rel),col=as.factor(all_aucs_Rel$Donor),lwd=3)
text(log(all_aucs_Rel$aucs_Rel),labels=all_aucs_Rel$patients, cex= 0.7,adj=1,pos=3,col=color[1:15])

summary(glm(factor(all_aucs_Rel$Donors) ~ all_aucs_Rel$aucs_Rel,family = binomial(link = "probit")))

write.table(all_aucs_Rel,'../../Documents/Progetto_LuciaGabri/MhAGs_AUC_Pazienti_RelapseExpressionFiltered.txt',sep='\t',col.names=T,row.names=F,quote=F)

plot.new()
color=c25[c(1:11,13,18,22,24)]
plot(lines(loess.smooth(1:length(Alfe_Binders_Ic50_NoExprFilt_Rel$value), y=Alfe_Binders_Ic50_NoExprFilt_Rel$value)),xlim=c(0,800),ylim=c(0,500),col=color[1],xlab='Number of Predicted mhAGs Binders',ylab='IC50 (nM)')
lines(loess.smooth(1:length(Alfe_Binders_Ic50_NoExprFilt_Rel$value), y=Alfe_Binders_Ic50_NoExprFilt_Rel$value,degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[1],xlab='',ylab='IC50 (nM)',lwd=3)
#polygon(1:length(Besu_Binders_ord$value), y=Besu_Binders_ord$value,col='skyblue')
par(new=TRUE)
lines(loess.smooth(1:length(Besu_Binders_Ic50_NoExprFilt_Rel$value), y=Besu_Binders_Ic50_NoExprFilt_Rel$value,degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[2],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Calu_Binders_Ic50_NoExprFilt_Rel$value), y=Calu_Binders_Ic50_NoExprFilt_Rel$value,degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[3],xlab='',ylab='',lwd=3)
par(new=TRUE)
lines(loess.smooth(1:length(Deiv_Binders_Ic50_NoExprFilt_Rel$value), y=Deiv_Binders_Ic50_NoExprFilt_Rel$value,degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[4],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Dest_Binders_Ic50_NoExprFilt_Rel$value), y=Dest_Binders_Ic50_NoExprFilt_Rel$value,degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[5],xlab='',ylab='',lwd=3)
par(new=TRUE)
lines(loess.smooth(1:length(Dr1_Binders_Ic50_NoExprFilt_Rel$value), y=Dr1_Binders_Ic50_NoExprFilt_Rel$value,degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[6],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Dr4_Binders_Ic50_NoExprFilt_Rel$value), y=Dr4_Binders_Ic50_NoExprFilt_Rel$value,degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[7],xlab='',ylab='',lwd=3,)
par(new=TRUE)
lines(loess.smooth(1:length(Dr5_Binders_Ic50_NoExprFilt_Rel$value), y=Dr5_Binders_Ic50_NoExprFilt_Rel$value,degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[8],xlab='',ylab='',lwd=3,)
par(new=TRUE)
lines(loess.smooth(1:length(Foca_Binders_Ic50_NoExprFilt_Rel$value), y=Foca_Binders_Ic50_NoExprFilt_Rel$value,degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[9],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Gagra_Binders_Ic50_NoExprFilt_Rel$value), y=Gagra_Binders_Ic50_NoExprFilt_Rel$value, degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[10],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Luan_Binders_Ic50_NoExprFilt_Rel$value), y=Luan_Binders_Ic50_NoExprFilt_Rel$value,  degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[11],xlab='',ylab='',lwd=3,)
par(new=TRUE)
lines(loess.smooth(1:length(Mabi_Binders_Ic50_NoExprFilt_Rel$value), y=Mabi_Binders_Ic50_NoExprFilt_Rel$value,  degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[12],xlab='',ylab='',lwd=3)
par(new=TRUE)
lines(loess.smooth(1:length(Moge_Binders_Ic50_NoExprFilt_Rel$value), y=Moge_Binders_Ic50_NoExprFilt_Rel$value,  degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[13],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Piag_Binders_Ic50_NoExprFilt_Rel$value), y=Piag_Binders_Ic50_NoExprFilt_Rel$value,  degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[14],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Prelu_Binders_Ic50_NoExprFilt_Rel$value), y=Prelu_Binders_Ic50_NoExprFilt_Rel$value,degree=1,span=0.8),xlim=c(0,800),ylim=c(0,500),col=color[15],xlab='',ylab='',lwd=3)
legend("topright", inset=c(-0,0), c('alfe','besu','calu','deiv','dest','dr1','dr4','dr5','foca','gagra','luan','mabi','moge','piag','prelu'),col=color[1:15], title="Patients",cex=0.7,lwd=rep(2,15),lty=c(1,2,1,2,1,2,1,1,2,2,1,1,2,2,1))









#*************************************************************************************************************************************************
#*************************************************************************************************************************************************
#***************************************************************HLA classe II**********************************************************************
#*************************************************************************************************************************************************
#*************************************************************************************************************************************************

sh Launcher_snp2epi_Paper_Annotation_NewForMhcII.sh ALFE_Somatic_Rerun.tsv ALFE1_ALFE2_NoMismatchHla_MHCII DRB1_0801 ALFE_Expression_Kallisto_NeoEpitopes.txt 

sh Launcher_snp2epi_Paper_Annotation_NewForMhcII.sh BESU1_BESU2_Somatic.Anno.tsv BESU1_BESU2_NoMismatchHla_MHCII DRB1_0104,DRB1_1104 BESU_Expression_Kallisto_NeoEpitopes.txt

sh Launcher_snp2epi_Paper_Annotation_NewForMhcII.sh CALU1_CALU2_Somatic.Anno.tsv CALU1_CALU2_NoMismatchHla_MHCII DRB1_1101 CALU_Expression_Kallisto_NeoEpitopes.txt

sh Launcher_snp2epi_Paper_Annotation_MHCII_DEIV.sh  DEIV1_DEIV2_DEIV2b_Somatic.Anno.tsv  DEIV2_DEIV2b_NoMismatchHla_MHCII DRB1_0102,DRB1_1501 DEIV_Expression_Kallisto_NeoEpitopes.txt

sh Launcher_snp2epi_Paper_Annotation_NewForMhcII.sh DEST1_DEST2_Somatic.Anno.tsv DEST1_DEST2_NoMismatchHla_MHCII  DRB1_1101,DRB1_0301 DEST_Expression_Kallisto_NeoEpitopes.txt

sh Launcher_snp2epi_Paper_Annotation_NewForMhcII.sh DR1_Somatic_Rerun.tsv DR11_DR12_NoMismatchHla_MHCII DRB1_0401,DRB1_0103 DR1_Expression_Kallisto_NeoEpitopes.txt

sh Launcher_snp2epi_Paper_Annotation_NewForMhcII.sh DR41_DR42_Somatic.Anno.tsv DR41_DR42_NoMismatchHla_MHCII DRB1_1501,DRB1_0401 DR4_Expression_Kallisto_NeoEpitopes.txt

sh Launcher_snp2epi_Paper_Annotation_NewForMhcII.sh DR51_DR52_Somatic.Anno.tsv DR51_DR52_NoMismatchHla_MHCII DRB1_0701 DR5_Expression_Kallisto_NeoEpitopes.txt

sh Launcher_snp2epi_Paper_Annotation_NewForMhcII.sh FOCA1_FOCA2_Somatic.Anno.tsv FOCA1_FOCA2_NoMismatchHla_MHCII DRB1_0801,DRB1_0301 FOCA_Expression_Kallisto_NeoEpitopes.txt

sh Launcher_snp2epi_Paper_Annotation_NewForMhcII.sh GAGRA1_GAGRA2_Somatic.Anno.tsv GAGRA1_GAGRA2_NoMismatchHla_MHCII  DRB1_1104,DRB1_0701 GAGRA_Expression_Kallisto_NeoEpitopes.txt

sh Launcher_snp2epi_Paper_Annotation_NewForMhcII.sh LUAN1_LUAN2_Somatic.Anno.tsv LUAN1_LUAN2_NoMismatchHla_MHCII DRB1_1454,DRB1_1601 LUAN_Expression_Kallisto_NeoEpitopes.txt

sh Launcher_snp2epi_Paper_Annotation_NewForMhcII.sh MABI1_MABI2_Somatic.Anno.tsv MABI1_MABI2_NoMismatchHla_MHCII DRB1_1501 MABI_Expression_Kallisto_NeoEpitopes.txt

sh Launcher_snp2epi_Paper_Annotation_NewForMhcII.sh MOGE1_MOGE2_Somatic.Anno.tsv MOGE1_MOGE2_NoMismatchHla_MHCII DRB1_0402,DRB1_0701 MOGE_Expression_Kallisto_NeoEpitopes.txt

sh Launcher_snp2epi_Paper_Annotation_NewForMhcII.sh PIAG1_PIAG2_Somatic.Anno.tsv PIAG1_PIAG2_NoMismatchHla_MHCII DRB1_1104 PIAG_Expression_Kallisto_NeoEpitopes.txt

sh Launcher_snp2epi_Paper_Annotation_NewForMhcII.sh PRELU1_PRELU2_Somatic.Anno.tsv PRELU1_PRELU2_NoMismatchHla_MHCII  DRB1_1104 PRELU_Expression_Kallisto_NeoEpitopes.txt





#ALFE
#Common diagnosis relapse
Mutationi=open('./HomeQlogin/ALFE1_ALFE2_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_ALFE1_ALFE2_NoMismatchHla_MHCII_weak_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_ALFE1_ALFE2_NoMismatchHla_MHCII_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_ALFE_RelDiag_MHCII','w')
header=Summary[0:2]

#out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
#binders=[]
#for i in Summary[2:]:
#	for j in Mutationi:
#		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip():
#			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))
#
#
#for i in binders:
#	for j in Excel[2:]:
#		if str(i.split('\t')[1])==str(j.split('\t')[1]):
#			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\n'))
#
#
#
#out.close()

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/ALFE1_ALFE2_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_ALFE1_ALFE2_NoMismatchHla_MHCII_strong_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_ALFE1_ALFE2_NoMismatchHla_MHCII_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_ALFE_RelDiag_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\n'))



out.close()


#Relapse Only
Mutationi=open('./HomeQlogin/ALFE1_ALFE2_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_ALFE1_ALFE2_NoMismatchHla_MHCII_weak_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_ALFE1_ALFE2_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_ALFE_RelOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/ALFE1_ALFE2_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_ALFE1_ALFE2_NoMismatchHla_MHCII_strong_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_ALFE1_ALFE2_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_ALFE_RelOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\n'))



out.close()




#BESU
#Common diagnosis relapse
Mutationi=open('./HomeQlogin/BESU1_BESU2_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_BESU1_BESU2_NoMismatchHla_MHCII_weak_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_BESU1_BESU2_NoMismatchHla_MHCII_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_BESU_RelDiag_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/BESU1_BESU2_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_BESU1_BESU2_NoMismatchHla_MHCII_strong_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_BESU1_BESU2_NoMismatchHla_MHCII_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_BESU_RelDiag_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()


#Diagnosis Only
Mutationi=open('./HomeQlogin/BESU1_BESU2_NoMismatchHla_MHCII_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_BESU1_BESU2_NoMismatchHla_MHCII_weak_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_BESU1_BESU2_NoMismatchHla_MHCII_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_BESU_DiagOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/BESU1_BESU2_NoMismatchHla_MHCII_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_BESU1_BESU2_NoMismatchHla_MHCII_strong_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_BESU1_BESU2_NoMismatchHla_MHCII_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_BESU_DiagOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()





#CALU
#Common diagnosis relapse
Mutationi=open('./HomeQlogin/CALU1_CALU2_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_CALU1_CALU2_NoMismatchHla_MHCII_weak_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_CALU1_CALU2_NoMismatchHla_MHCII_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_CALU_RelDiag_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\n'))



out.close()



Mutationi=open('./HomeQlogin/CALU1_CALU2_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_CALU1_CALU2_NoMismatchHla_MHCII_strong_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_CALU1_CALU2_NoMismatchHla_MHCII_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_CALU_RelDiag_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\n'))



out.close()



#Relapse Only
Mutationi=open('./HomeQlogin/CALU1_CALU2_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_CALU1_CALU2_NoMismatchHla_MHCII_weak_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_CALU1_CALU2_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_CALU_RelOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\n'))



out.close()



Mutationi=open('./HomeQlogin/CALU1_CALU2_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_CALU1_CALU2_NoMismatchHla_MHCII_strong_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_CALU1_CALU2_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_CALU_RelOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\n'))



out.close()



#DEIV
#Relapse Only
Mutationi=open('./HomeQlogin/DEIV2_DEIV2b_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DEIV2_DEIV2b_NoMismatchHla_MHCII_weak_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DEIV2_DEIV2b_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DEIV_RelOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()


#Relapse Only
Mutationi=open('./HomeQlogin/DEIV2_DEIV2b_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DEIV2_DEIV2b_NoMismatchHla_MHCII_strong_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DEIV2_DEIV2b_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DEIV_RelOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()


#Diagnosis Only
Mutationi=open('./HomeQlogin/DEIV2_DEIV2b_NoMismatchHla_MHCII_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DEIV2_DEIV2b_NoMismatchHla_MHCII_weak_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DEIV2_DEIV2b_NoMismatchHla_MHCII_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DEIV_DiagOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+str(j.split('\t')[7])+'\n'))



out.close()



Mutationi=open('./HomeQlogin/DEIV2_DEIV2b_NoMismatchHla_MHCII_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DEIV2_DEIV2b_NoMismatchHla_MHCII_strong_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DEIV2_DEIV2b_NoMismatchHla_MHCII_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DEIV_DiagOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()



#DEST
#Common diagnosis relapse
Mutationi=open('./HomeQlogin/DEST1_DEST2_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DEST1_DEST2_NoMismatchHla_MHCII_weak_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_DEST1_DEST2_NoMismatchHla_MHCII_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DEST_RelDiag_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()

Mutationi=open('./HomeQlogin/DEST1_DEST2_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DEST1_DEST2_NoMismatchHla_MHCII_strong_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_DEST1_DEST2_NoMismatchHla_MHCII_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DEST_RelDiag_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()






#DR1
#Common diagnosis relapse
Mutationi=open('./HomeQlogin/DR11_DR12_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR11_DR12_NoMismatchHla_MHCII_weak_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_DR11_DR12_NoMismatchHla_MHCII_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DR1_RelDiag_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()

Mutationi=open('./HomeQlogin/DR11_DR12_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR11_DR12_NoMismatchHla_MHCII_strong_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_DR11_DR12_NoMismatchHla_MHCII_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DR1_RelDiag_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()




#DR4
#Common diagnosis relapse
Mutationi=open('./HomeQlogin/DR41_DR42_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR41_DR42_NoMismatchHla_MHCII_weak_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_DR41_DR42_NoMismatchHla_MHCII_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DR4_RelDiag_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/DR41_DR42_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR41_DR42_NoMismatchHla_MHCII_strong_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_DR41_DR42_NoMismatchHla_MHCII_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DR4_RelDiag_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()




#Diagnosis Only
Mutationi=open('./HomeQlogin/DR41_DR42_NoMismatchHla_MHCII_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR41_DR42_NoMismatchHla_MHCII_weak_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DR41_DR42_NoMismatchHla_MHCII_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DR4_DiagOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/DR41_DR42_NoMismatchHla_MHCII_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR41_DR42_NoMismatchHla_MHCII_strong_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DR41_DR42_NoMismatchHla_MHCII_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DR4_DiagOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()


#Relapse Only
Mutationi=open('./HomeQlogin/DR41_DR42_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR41_DR42_NoMismatchHla_MHCII_weak_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DR41_DR42_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DR4_RelOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Relapse'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/DR41_DR42_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR41_DR42_NoMismatchHla_MHCII_strong_binders_RelOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DR41_DR42_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DR4_RelOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Relapse'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()




#DR5
#Common diagnosis relapse
Mutationi=open('./HomeQlogin/DR51_DR52_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR51_DR52_NoMismatchHla_MHCII_weak_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_DR51_DR52_NoMismatchHla_MHCII_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DR5_RelDiag_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/DR51_DR52_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR51_DR52_NoMismatchHla_MHCII_strong_binders_RelDiag').readlines()
Excel=open('./HomeQlogin/RESULTS_DR51_DR52_NoMismatchHla_MHCII_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DR5_RelDiag_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()




#Diagnosis Only
Mutationi=open('./HomeQlogin/DR51_DR52_NoMismatchHla_MHCII_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR51_DR52_NoMismatchHla_MHCII_weak_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DR51_DR52_NoMismatchHla_MHCII_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/WeakBinders_Annotati_DR5_DiagOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()


Mutationi=open('./HomeQlogin/DR51_DR52_NoMismatchHla_MHCII_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./HomeQlogin/RESULTS_DR51_DR52_NoMismatchHla_MHCII_strong_binders_DiagOnly').readlines()
Excel=open('./HomeQlogin/RESULTS_DR51_DR52_NoMismatchHla_MHCII_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./HomeQlogin/StrongBinders_Annotati_DR5_DiagOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()


#FOCA
#Common diagnosis relapse
Mutationi=open('./Home2/FOCA1_FOCA2_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./Home2/RESULTS_FOCA1_FOCA2_NoMismatchHla_MHCII_weak_binders_RelDiag').readlines()
Excel=open('./Home2/RESULTS_FOCA1_FOCA2_NoMismatchHla_MHCII_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./Home2/WeakBinders_Annotati_FOCA_RelDiag_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()


Mutationi=open('./Home2/FOCA1_FOCA2_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./Home2/RESULTS_FOCA1_FOCA2_NoMismatchHla_MHCII_strong_binders_RelDiag').readlines()
Excel=open('./Home2/RESULTS_FOCA1_FOCA2_NoMismatchHla_MHCII_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./Home2/StrongBinders_Annotati_FOCA_RelDiag_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()



#GAGRA
#Common diagnosis relapse
Mutationi=open('./Home2/GAGRA1_GAGRA2_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./Home2/RESULTS_GAGRA1_GAGRA2_NoMismatchHla_MHCII_weak_binders_RelDiag').readlines()
Excel=open('./Home2/RESULTS_GAGRA1_GAGRA2_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapseOnly.xls').readlines()
out=open('./Home2/WeakBinders_Annotati_GAGRA_RelDiag_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()


Mutationi=open('./Home2/GAGRA1_GAGRA2_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./Home2/RESULTS_GAGRA1_GAGRA2_NoMismatchHla_MHCII_strong_binders_RelDiag').readlines()
Excel=open('./Home2/RESULTS_GAGRA1_GAGRA2_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapseOnly.xls').readlines()
out=open('./Home2/StrongBinders_Annotati_GAGRA_RelDiag_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()


#LUAN
#Common diagnosis relapse
Mutationi=open('./Home2/LUAN1_LUAN2_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./Home2/RESULTS_LUAN1_LUAN2_NoMismatchHla_MHCII_weak_binders_RelDiag').readlines()
Excel=open('./Home2/RESULTS_LUAN1_LUAN2_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapseOnly.xls').readlines()
out=open('./Home2/WeakBinders_Annotati_LUAN_RelOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()


#Common diagnosis relapse
Mutationi=open('./Home2/LUAN1_LUAN2_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./Home2/RESULTS_LUAN1_LUAN2_NoMismatchHla_MHCII_strong_binders_RelDiag').readlines()
Excel=open('./Home2/RESULTS_LUAN1_LUAN2_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapseOnly.xls').readlines()
out=open('./Home2/StrongBinders_Annotati_LUAN_RelOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()


#Relapse Only
Mutationi=open('./Home2/LUAN1_LUAN2_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./Home2/RESULTS_LUAN1_LUAN2_NoMismatchHla_MHCII_weak_binders_RelOnly').readlines()
Excel=open('./Home2/RESULTS_LUAN1_LUAN2_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./Home2/WeakBinders_Annotati_LUAN_RelOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()


Mutationi=open('./Home2/LUAN1_LUAN2_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./Home2/RESULTS_LUAN1_LUAN2_NoMismatchHla_MHCII_strong_binders_RelOnly').readlines()
Excel=open('./Home2/RESULTS_LUAN1_LUAN2_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./Home2/StrongBinders_Annotati_LUAN_RelOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()

#MABI
#Common diagnosis relapse
Mutationi=open('./Home2/MABI1_MABI2_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./Home2/RESULTS_MABI1_MABI2_NoMismatchHla_MHCII_weak_binders_RelDiag').readlines()
Excel=open('./Home2/RESULTS_MABI1_MABI2_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapseOnly.xls').readlines()
out=open('./Home2/WeakBinders_Annotati_MABI_RelOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\n'))



out.close()


#Common diagnosis relapse
Mutationi=open('./Home2/MABI1_MABI2_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./Home2/RESULTS_MABI1_MABI2_NoMismatchHla_MHCII_strong_binders_RelDiag').readlines()
Excel=open('./Home2/RESULTS_MABI1_MABI2_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapseOnly.xls').readlines()
out=open('./Home2/StrongBinders_Annotati_MABI_RelOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\n'))



out.close()





#MOGE
#Common diagnosis relapse
Mutationi=open('./Home2/MOGE1_MOGE2_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./Home2/RESULTS_MOGE1_MOGE2_NoMismatchHla_MHCII_weak_binders_RelDiag').readlines()
Excel=open('./Home2/RESULTS_MOGE1_MOGE2_NoMismatchHla_MHCII_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./Home2/WeakBinders_Annotati_MOGE_RelDiag_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()


Mutationi=open('./Home2/MOGE1_MOGE2_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./Home2/RESULTS_MOGE1_MOGE2_NoMismatchHla_MHCII_strong_binders_RelDiag').readlines()
Excel=open('./Home2/RESULTS_MOGE1_MOGE2_NoMismatchHla_MHCII_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./Home2/StrongBinders_Annotati_MOGE_RelDiag_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()



#Relapse Only
Mutationi=open('./Home2/MOGE1_MOGE2_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./Home2/RESULTS_MOGE1_MOGE2_NoMismatchHla_MHCII_weak_binders_RelOnly').readlines()
Excel=open('./Home2/RESULTS_MOGE1_MOGE2_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./Home2/WeakBinders_Annotati_MOGE_RelOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Relapse'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()


Mutationi=open('./Home2/MOGE1_MOGE2_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./Home2/RESULTS_MOGE1_MOGE2_NoMismatchHla_MHCII_strong_binders_RelOnly').readlines()
Excel=open('./Home2/RESULTS_MOGE1_MOGE2_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./Home2/StrongBinders_Annotati_MOGE_RelOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Relapse'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\t'+str(j.split('\t')[7])+'\n'))



out.close()



#PIAG
#Relapse Only
Mutationi=open('./Home2/PIAG1_PIAG2_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./Home2/RESULTS_PIAG1_PIAG2_NoMismatchHla_MHCII_weak_binders_RelOnly').readlines()
Excel=open('./Home2/RESULTS_PIAG1_PIAG2_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./Home2/WeakBinders_Annotati_PIAG_RelOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Relapse'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\n'))



out.close()


Mutationi=open('./Home2/PIAG1_PIAG2_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.txt').readlines()
Summary=open('./Home2/RESULTS_PIAG1_PIAG2_NoMismatchHla_MHCII_strong_binders_RelOnly').readlines()
Excel=open('./Home2/RESULTS_PIAG1_PIAG2_NoMismatchHla_MHCII_PeptidesPrediction_RelapseOnly.xls').readlines()
out=open('./Home2/StrongBinders_Annotati_PIAG_RelOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Relapse'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\n'))



out.close()


#Diagnosis Only
Mutationi=open('./Home2/PIAG1_PIAG2_NoMismatchHla_MHCII_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./Home2/RESULTS_PIAG1_PIAG2_NoMismatchHla_MHCII_weak_binders_DiagOnly').readlines()
Excel=open('./Home2/RESULTS_PIAG1_PIAG2_NoMismatchHla_MHCII_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./Home2/WeakBinders_Annotati_PIAG_DiagOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\n'))



out.close()


Mutationi=open('./Home2/PIAG1_PIAG2_NoMismatchHla_MHCII_PeptidesPrediction_DiagnosisOnly.txt').readlines()
Summary=open('./Home2/RESULTS_PIAG1_PIAG2_NoMismatchHla_MHCII_strong_binders_DiagOnly').readlines()
Excel=open('./Home2/RESULTS_PIAG1_PIAG2_NoMismatchHla_MHCII_PeptidesPrediction_DiagnosisOnly.xls').readlines()
out=open('./Home2/StrongBinders_Annotati_PIAG_DiagOnly_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Diagnosis'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\n'))



out.close()





#PRELU
#Common diagnosis relapse
Mutationi=open('./Home2/PRELU1_PRELU2_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./Home2/RESULTS_PRELU1_PRELU2_NoMismatchHla_MHCII_weak_binders_RelDiag').readlines()
Excel=open('./Home2/RESULTS_PRELU1_PRELU2_NoMismatchHla_MHCII_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./Home2/WeakBinders_Annotati_PRELU_RelDiag_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\n'))



out.close()


Mutationi=open('./Home2/PRELU1_PRELU2_NoMismatchHla_MHCII_PeptidesPrediction_CommonDiagnosisRelapse.txt').readlines()
Summary=open('./Home2/RESULTS_PRELU1_PRELU2_NoMismatchHla_MHCII_strong_binders_RelDiag').readlines()
Excel=open('./Home2/RESULTS_PRELU1_PRELU2_NoMismatchHla_MHCII_PeptidesPrediction_CommonRelapseDiagnosis.xls').readlines()
out=open('./Home2/StrongBinders_Annotati_PRELU_RelDiag_MHCII','w')
header=Summary[0:2]

out.write(''.join(header).replace(',DRB1','\tDRB1').replace('\t\t','\t').rstrip()+'\tMutation'+'\t'+'GeneSymbol'+'\t'+'Chrom'+'\t'+'Pos_NT'+'\t'+'WT_NT'+'\t'+'MUT_NT'+'\t'+'Expr_Relapse'+'\t'+'Expr_Diagnosis'+'\t'+'\t'.join(''.join(header).replace(' HLA','\tHLA').replace('\t\t','\t').split('\t')[2:]).replace('DRB1','IC50_DRB1').rstrip()+'\n')
binders=[]
for i in Summary[2:]:
	for j in Mutationi:
		if str(j.split('\t')[0].rstrip()) == str(i.split('\t')[0].split('_')[1].split(';')[0]).rstrip() and str(i.split('\t')[1]).rstrip() in str(j.split('\t')[15].rstrip()):
			binders.append(''.join(''.join(i).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip()))


for i in binders:
	for j in Excel[2:]:
		if str(i.split('\t')[1])==str(j.split('\t')[1]):
			out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\n'))



out.close()