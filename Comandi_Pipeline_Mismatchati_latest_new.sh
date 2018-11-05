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
			#Alfe_Binders_NoExprFilt_Expression=Alfe_Binders_NoExprFilt[,c(3,18,17)]
			#sel <- grepl("IC50_HLA",names(Alfe_Binders_NoExprFilt))
			#Alfe_Binders_NoExprFilt[sel] <- lapply(Alfe_Binders_NoExprFilt[sel], function(x) replace(x,x >500 , 0) )
			#Alfe_Binders_NoExprFilt=Alfe_Binders_NoExprFilt[,c(5,14:16)]
			#Alfe_Binders_NoExprFilt=melt(Alfe_Binders_NoExprFilt)
			#Alfe_Binders_NoExprFilt=Alfe_Binders_NoExprFilt[Alfe_Binders_NoExprFilt$value>0,]
			#
			#Alfe_Binders=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Alfe_Binders))
			#Alfe_Binders[sel] <- lapply(Alfe_Binders[sel], function(x) replace(x,x >500 , 0) )
			#Alfe_Binders=Alfe_Binders[,c(5,14:16)]
			#Alfe_Binders=melt(Alfe_Binders)
			#Alfe_Binders=Alfe_Binders[Alfe_Binders$value>0,]
			#
			#
			#Alfe_Binders=rbind(strong_filtered,weak_filtered)
			#
			#Alfe_Binders_Coverages=Alfe_Binders[,c(1,2,3,27,29)]
			#
			#Alfe_Binders_Genes=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Alfe_Binders_Genes))
			#Alfe_Binders_Genes[sel] <- lapply(Alfe_Binders_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Alfe_Binders_Genes=Alfe_Binders_Genes[,c(1,2,3,14:16)]
			#Alfe_Binders_Genes=data.table::as.data.table(Alfe_Binders_Genes)
			#Alfe_Binders_Genes=data.table::melt.data.table(data=Alfe_Binders_Genes,id=1:3)
			#Alfe_Binders_Genes=Alfe_Binders_Genes[Alfe_Binders_Genes$value>0,]
			#
			#
			#
			#Alfe_Binders_NotExpressed=rbind(strong,weak)
			#
			#Alfe_Binders_NotExpressed_Coverages=Alfe_Binders_NotExpressed[,c(1,2,3,27,29)]
			#
			#Alfe_Binders_NotExpressed_Genes=rbind(strong,weak)
			#sel <- grepl("IC50_HLA",names(Alfe_Binders_NotExpressed_Genes))
			#Alfe_Binders_NotExpressed_Genes[sel] <- lapply(Alfe_Binders_NotExpressed_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Alfe_Binders_NotExpressed_Genes=Alfe_Binders_NotExpressed_Genes[,c(1,2,3,14:16)]
			#Alfe_Binders_NotExpressed_Genes=data.table::as.data.table(Alfe_Binders_NotExpressed_Genes)
			#Alfe_Binders_NotExpressed_Genes=data.table::melt.data.table(data=Alfe_Binders_NotExpressed_Genes,id=1:3)
			#Alfe_Binders_NotExpressed_Genes=Alfe_Binders_NotExpressed_Genes[Alfe_Binders_NotExpressed_Genes$value>0,]
			#
			
			
			
			#Alfe_Binders_Genes_AlleleCoverage=unique(merge(Alfe_Binders_NotExpressed_Genes,Alfe_Binders_NotExpressed_Coverages,by=c('chr','pos','genes')))
			
			
			
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
			
			#sel <- grepl("IC50_HLA",names(Besu_Binders_NoExprFilt))
			#Besu_Binders_NoExprFilt[sel] <- lapply(Besu_Binders_NoExprFilt[sel], function(x) replace(x,x > 500 , 0) )
			#Besu_Binders_NoExprFilt=Besu_Binders_NoExprFilt[,c(5,17:22)]
			#Besu_Binders_NoExprFilt=melt(Besu_Binders_NoExprFilt)
			#Besu_Binders_NoExprFilt=Besu_Binders_NoExprFilt[Besu_Binders_NoExprFilt$value>0,]
			#
			#
			#Besu_Binders=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Besu_Binders))
			#Besu_Binders[sel] <- lapply(Besu_Binders[sel], function(x) replace(x,x > 500 , 0) )
			#Besu_Binders=Besu_Binders[,c(5,17:22)]
			#Besu_Binders=melt(Besu_Binders)
			#Besu_Binders=Besu_Binders[Besu_Binders$value>0,]
			#
			#Besu_Binders_Genes=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Besu_Binders_Genes))
			#Besu_Binders_Genes[sel] <- lapply(Besu_Binders_Genes[sel], function(x) replace(x,x > 500 , 0) )
			#Besu_Binders_Genes=Besu_Binders_Genes[,c(3,17:22)]
			#Besu_Binders_Genes=melt(Besu_Binders_Genes)
			#Besu_Binders_Genes=Besu_Binders_Genes[Besu_Binders_Genes$value>0,]
			#
			#
			#Besu_Binders_Coverages=Besu_Binders[,c(1,2,3,33,35)]
			#
			#Besu_Binders_Genes=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Besu_Binders_Genes))
			#Besu_Binders_Genes[sel] <- lapply(Besu_Binders_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Besu_Binders_Genes=Besu_Binders_Genes[,c(1,2,3,17:22)]
			#Besu_Binders_Genes=data.table::as.data.table(Besu_Binders_Genes)
			#Besu_Binders_Genes=data.table::melt.data.table(data=Besu_Binders_Genes,id=1:3)
			#Besu_Binders_Genes=Besu_Binders_Genes[Besu_Binders_Genes$value>0,]
			#
			#Besu_Binders_Genes_AlleleCoverage=unique(merge(Besu_Binders_Genes,Besu_Binders_Coverages,by=c('chr','pos','genes')))
			#
			#write.table(unique(besu_NeoAG_mhci_MinorAntigens_filtered[,c(1:24,33:36)]),"Summary_BESU_Annotato_MHCI_MinorAntigens_NoMismatches_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)
			#
			#
			#
			#Besu_Binders_NotExpressed=rbind(strong,weak)
			#Besu_Binders_NotExpressed_Coverages=Besu_Binders_NotExpressed[,c(1,2,3,33,35)]
			#
			#Besu_Binders_NotExpressed_Genes=rbind(strong,weak)
			#sel <- grepl("IC50_HLA",names(Besu_Binders_NotExpressed_Genes))
			#Besu_Binders_NotExpressed_Genes[sel] <- lapply(Besu_Binders_NotExpressed_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Besu_Binders_NotExpressed_Genes=Besu_Binders_NotExpressed_Genes[,c(1,2,3,17:22)]
			#Besu_Binders_NotExpressed_Genes=data.table::as.data.table(Besu_Binders_NotExpressed_Genes)
			#Besu_Binders_NotExpressed_Genes=data.table::melt.data.table(data=Besu_Binders_NotExpressed_Genes,id=1:3)
			#Besu_Binders_NotExpressed_Genes=Besu_Binders_NotExpressed_Genes[Besu_Binders_NotExpressed_Genes$value>0,]
			#
			#
			#Besu_Binders_NotExpressed_Genes_AlleleCoverage=unique(merge(Besu_Binders_NotExpressed_Genes,Besu_Binders_NotExpressed_Coverages,by=c('chr','pos','genes')))
			#
			
			
			
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
			
			
			
			
			
			
			
			
			
			
			#sel <- grepl("IC50_HLA",names(Calu_Binders_NoExprFilt))
			#Calu_Binders_NoExprFilt[sel] <- lapply(Calu_Binders_NoExprFilt[sel], function(x) replace(x,x > 500 , 0) )
			#Calu_Binders_NoExprFilt=Calu_Binders_NoExprFilt[,c(5,14:16)]
			#Calu_Binders_NoExprFilt=melt(Calu_Binders_NoExprFilt)
			#Calu_Binders_NoExprFilt=Calu_Binders_NoExprFilt[Calu_Binders_NoExprFilt$value>0,]
			#
			#
			#
			#Calu_Binders=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Calu_Binders))
			#Calu_Binders[sel] <- lapply(Calu_Binders[sel], function(x) replace(x,x > 500 , 0) )
			#Calu_Binders=Calu_Binders[,c(5,14:16)]
			#Calu_Binders=melt(Calu_Binders)
			#Calu_Binders=Calu_Binders[Calu_Binders$value>0,]
			#
			#Calu_Binders=rbind(strong_filtered,weak_filtered)
			#
			#Calu_Binders_Coverages=Calu_Binders[,c(1,2,3,27,29)]
			#
			#Calu_Binders_Genes=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Calu_Binders_Genes))
			#Calu_Binders_Genes[sel] <- lapply(Calu_Binders_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Calu_Binders_Genes=Calu_Binders_Genes[,c(1,2,3,14:16)]
			#Calu_Binders_Genes=data.table::as.data.table(Calu_Binders_Genes)
			#Calu_Binders_Genes=data.table::melt.data.table(data=Calu_Binders_Genes,id=1:3)
			#Calu_Binders_Genes=Calu_Binders_Genes[Calu_Binders_Genes$value>0,]
			#
			#Calu_Binders_Genes_AlleleCoverage=unique(merge(Calu_Binders_Genes,Calu_Binders_Coverages,by=c('chr','pos','genes')))
			#
			#
			#Calu_Binders_NotExpressed=rbind(strong,weak)
			#
			#Calu_Binders_NotExpressed_Coverages=Calu_Binders_NotExpressed[,c(1,2,3,27,29)]
			#
			#Calu_Binders_NotExpressed_Genes=rbind(strong,weak)
			#sel <- grepl("IC50_HLA",names(Calu_Binders_NotExpressed_Genes))
			#Calu_Binders_NotExpressed_Genes[sel] <- lapply(Calu_Binders_NotExpressed_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Calu_Binders_NotExpressed_Genes=Calu_Binders_NotExpressed_Genes[,c(1,2,3,14:16)]
			#Calu_Binders_NotExpressed_Genes=data.table::as.data.table(Calu_Binders_NotExpressed_Genes)
			#Calu_Binders_NotExpressed_Genes=data.table::melt.data.table(data=Calu_Binders_NotExpressed_Genes,id=1:3)
			#Calu_Binders_NotExpressed_Genes=Calu_Binders_NotExpressed_Genes[Calu_Binders_NotExpressed_Genes$value>0,]
			#
			#Calu_Binders_NotExpressed_Genes_AlleleCoverage=unique(merge(Calu_Binders_NotExpressed_Genes,Calu_Binders_NotExpressed_Coverages,by=c('chr','pos','genes')))
			#
			#
			#
			#
			#write.table(unique(calu_NeoAG_mhci_MinorAntigens_filtered[,c(1:18,27:30)]),"Summary_CALU_Annotato_MHCI_MinorAntigens_NoMismatches_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)
			
			
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
			
			
			
			
			#
			#Deiv_Binders_NoExprFilt=rbind(strong,weak)
			#sel <- grepl("IC50_HLA",names(Deiv_Binders_NoExprFilt))
			#Deiv_Binders_NoExprFilt[sel] <- lapply(Deiv_Binders_NoExprFilt[sel], function(x) replace(x,x > 500 , 0) )
			#Deiv_Binders_NoExprFilt=Deiv_Binders_NoExprFilt[,c(5,17:22)]
			#Deiv_Binders_NoExprFilt=melt(Deiv_Binders_NoExprFilt)
			#Deiv_Binders_NoExprFilt=Deiv_Binders_NoExprFilt[Deiv_Binders_NoExprFilt$value>0,]
			#
			#
			#Deiv_Binders=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Deiv_Binders))
			#Deiv_Binders[sel] <- lapply(Deiv_Binders[sel], function(x) replace(x,x > 500 , 0) )
			#Deiv_Binders=Deiv_Binders[,c(5,17:22)]
			#Deiv_Binders=melt(Deiv_Binders)
			#Deiv_Binders=Deiv_Binders[Deiv_Binders$value>0,]
			#
			#Deiv_Binders=rbind(strong_filtered,weak_filtered)
			#
			#Deiv_Binders_Coverages=Deiv_Binders[,c(1,2,3,33,35)]
			#
			#Deiv_Binders_Genes=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Deiv_Binders_Genes))
			#Deiv_Binders_Genes[sel] <- lapply(Deiv_Binders_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Deiv_Binders_Genes=Deiv_Binders_Genes[,c(1,2,3,17:22)]
			#Deiv_Binders_Genes=data.table::as.data.table(Deiv_Binders_Genes)
			#Deiv_Binders_Genes=data.table::melt.data.table(data=Deiv_Binders_Genes,id=1:3)
			#Deiv_Binders_Genes=Deiv_Binders_Genes[Deiv_Binders_Genes$value>0,]
			#
			#Deiv_Binders_Genes_AlleleCoverage=unique(merge(Deiv_Binders_Genes,Deiv_Binders_Coverages,by=c('chr','pos','genes')))
			#
			#
			#
			#Deiv_Binders_NotExpressed=rbind(strong,weak)
			#
			#Deiv_Binders_NotExpressed_Coverages=Deiv_Binders_NotExpressed[,c(1,2,3,33,35)]
			#
			#Deiv_Binders_NotExpressed_Genes=rbind(strong,weak)
			#sel <- grepl("IC50_HLA",names(Deiv_Binders_NotExpressed_Genes))
			#Deiv_Binders_NotExpressed_Genes[sel] <- lapply(Deiv_Binders_NotExpressed_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Deiv_Binders_NotExpressed_Genes=Deiv_Binders_NotExpressed_Genes[,c(1,2,3,17:22)]
			#Deiv_Binders_NotExpressed_Genes=data.table::as.data.table(Deiv_Binders_NotExpressed_Genes)
			#Deiv_Binders_NotExpressed_Genes=data.table::melt.data.table(data=Deiv_Binders_NotExpressed_Genes,id=1:3)
			#Deiv_Binders_NotExpressed_Genes=Deiv_Binders_NotExpressed_Genes[Deiv_Binders_NotExpressed_Genes$value>0,]
			#
			#Deiv_Binders_NotExpressed_Genes_AlleleCoverage=unique(merge(Deiv_Binders_NotExpressed_Genes,Deiv_Binders_NotExpressed_Coverages,by=c('chr','pos','genes'),allow.cartesian=T))
			#
			#Deiv_Binders_NotExpressed_Genes_AlleleCoverage=unique(Deiv_Binders_NotExpressed_Genes_AlleleCoverage)
			#
			#
			#
			#write.table(unique(dr1_NeoAG_mhci_MinorAntigens_filtered[,c(1:24,33:36)]),"Summary_DR1_Annotato_MHCI_MinorAntigens_NoMismatches_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)
			#
			
			
			
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
			
			
			
			
			#
			#Dest_Binders_NoExprFilt=rbind(strong,weak)
			#sel <- grepl("IC50_HLA",names(Dest_Binders_NoExprFilt))
			#Dest_Binders_NoExprFilt[sel] <- lapply(Dest_Binders_NoExprFilt[sel], function(x) replace(x,x > 500 , 0) )
			#Dest_Binders_NoExprFilt=Dest_Binders_NoExprFilt[,c(5,14:16)]
			#Dest_Binders_NoExprFilt=melt(Dest_Binders_NoExprFilt)
			#Dest_Binders_NoExprFilt=Dest_Binders_NoExprFilt[Dest_Binders_NoExprFilt$value>0,]
			#
			#
			#
			#Dest_Binders=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Dest_Binders))
			#Dest_Binders[sel] <- lapply(Dest_Binders[sel], function(x) replace(x,x > 500 , 0) )
			#Dest_Binders=Dest_Binders[,c(5,14:16)]
			#Dest_Binders=melt(Dest_Binders)
			#Dest_Binders=Dest_Binders[Dest_Binders$value>0,]
			#
			#
			#Dest_Binders=rbind(strong_filtered,weak_filtered)
			#
			#Dest_Binders_Coverages=Dest_Binders[,c(1,2,3,27,29)]
			#
			#Dest_Binders_Genes=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Dest_Binders_Genes))
			#Dest_Binders_Genes[sel] <- lapply(Dest_Binders_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Dest_Binders_Genes=Dest_Binders_Genes[,c(1,2,3,14:16)]
			#Dest_Binders_Genes=data.table::as.data.table(Dest_Binders_Genes)
			#Dest_Binders_Genes=data.table::melt.data.table(data=Dest_Binders_Genes,id=1:3)
			#Dest_Binders_Genes=Dest_Binders_Genes[Dest_Binders_Genes$value>0,]
			#
			#Dest_Binders_Genes_AlleleCoverage=unique(merge(Dest_Binders_Genes,Dest_Binders_Coverages,by=c('chr','pos','genes')))
			#
			#
			#Dest_Binders_NotExpressed=rbind(strong,weak)
			#
			#Dest_Binders_NotExpressed_Coverages=Dest_Binders_NotExpressed[,c(1,2,3,27,29)]
			#
			#Dest_Binders_NotExpressed_Genes=rbind(strong,weak)
			#sel <- grepl("IC50_HLA",names(Dest_Binders_NotExpressed_Genes))
			#Dest_Binders_NotExpressed_Genes[sel] <- lapply(Dest_Binders_NotExpressed_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Dest_Binders_NotExpressed_Genes=Dest_Binders_NotExpressed_Genes[,c(1,2,3,14:16)]
			#Dest_Binders_NotExpressed_Genes=data.table::as.data.table(Dest_Binders_NotExpressed_Genes)
			#Dest_Binders_NotExpressed_Genes=data.table::melt.data.table(data=Dest_Binders_NotExpressed_Genes,id=1:3)
			#Dest_Binders_NotExpressed_Genes=Dest_Binders_NotExpressed_Genes[Dest_Binders_NotExpressed_Genes$value>0,]
			#
			#Dest_Binders_NotExpressed_Genes_AlleleCoverage=unique(merge(Dest_Binders_NotExpressed_Genes,Dest_Binders_NotExpressed_Coverages,by=c('chr','pos','genes')))
			#
			#
			#
			#write.table(unique(dest_NeoAG_mhci_MinorAntigens_filtered[,c(1:18,27:30)]),"Summary_DEST_Annotato_MHCI_MinorAntigens_NoMismatches_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)
			#
			
			
			
			
			
			
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
			
			
			
			#sel <- grepl("IC50_HLA",names(Dr1_Binders_NoExprFilt))
			#Dr1_Binders_NoExprFilt[sel] <- lapply(Dr1_Binders_NoExprFilt[sel], function(x) replace(x,x >500 , 0) )
			#Dr1_Binders_NoExprFilt=Dr1_Binders_NoExprFilt[,c(5,17:22)]
			#Dr1_Binders_NoExprFilt=melt(Dr1_Binders_NoExprFilt)
			#Dr1_Binders_NoExprFilt=Dr1_Binders_NoExprFilt[Dr1_Binders_NoExprFilt$value>0,]
			#
			#
			#Dr1_Binders=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Dr1_Binders))
			#Dr1_Binders[sel] <- lapply(Dr1_Binders[sel], function(x) replace(x,x >500 , 0) )
			#Dr1_Binders=Dr1_Binders[,c(5,17:22)]
			#Dr1_Binders=melt(Dr1_Binders)
			#Dr1_Binders=Dr1_Binders[Dr1_Binders$value>0,]
			#
			#
			#Dr1_Binders=rbind(strong_filtered,weak_filtered)
			#
			#Dr1_Binders_Coverages=Dr1_Binders[,c(1,2,3,33,35)]
			#
			#Dr1_Binders_Genes=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Dr1_Binders_Genes))
			#Dr1_Binders_Genes[sel] <- lapply(Dr1_Binders_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Dr1_Binders_Genes=Dr1_Binders_Genes[,c(1,2,3,17:22)]
			#Dr1_Binders_Genes=data.table::as.data.table(Dr1_Binders_Genes)
			#Dr1_Binders_Genes=data.table::melt.data.table(data=Dr1_Binders_Genes,id=1:3)
			#Dr1_Binders_Genes=Dr1_Binders_Genes[Dr1_Binders_Genes$value>0,]
			#
			#Dr1_Binders_Genes=nrow(unique(Dr1_Binders_Genes)
			#Dr1_Binders_Genes_AlleleCoverage=unique(merge(Dr1_Binders_Genes,Dr1_Binders_Coverages,by=c('chr','pos','genes'),allow.cartesian=T))
			#
			#
			#
			#
			#Dr1_Binders_NotExpressed=rbind(strong,weak)
			#
			#Dr1_Binders_NotExpressed_Coverages=Dr1_Binders_NotExpressed[,c(1,2,3,33,35)]
			#
			#Dr1_Binders_NotExpressed_Genes=rbind(strong,weak)
			#sel <- grepl("IC50_HLA",names(Dr1_Binders_NotExpressed_Genes))
			#Dr1_Binders_NotExpressed_Genes[sel] <- lapply(Dr1_Binders_NotExpressed_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Dr1_Binders_NotExpressed_Genes=Dr1_Binders_NotExpressed_Genes[,c(1,2,3,17:22)]
			#Dr1_Binders_NotExpressed_Genes=data.table::as.data.table(Dr1_Binders_NotExpressed_Genes)
			#Dr1_Binders_NotExpressed_Genes=data.table::melt.data.table(data=Dr1_Binders_NotExpressed_Genes,id=1:3)
			#Dr1_Binders_NotExpressed_Genes=Dr1_Binders_NotExpressed_Genes[Dr1_Binders_NotExpressed_Genes$value>0,]
			#
			#Dr1_Binders_NotExpressed_Genes=nrow(unique(Dr1_Binders_NotExpressed_Genes))
			#Dr1_Binders_NotExpressed_Genes_AlleleCoverage=unique(merge(Dr1_Binders_NotExpressed_Genes,Dr1_Binders_NotExpressed_Coverages,by=c('chr','pos','genes'),allow.cartesian = T))
			#Dr1_Binders_NotExpressed_Genes_AlleleCoverage=unique(Dr1_Binders_NotExpressed_Genes_AlleleCoverage)
			#
			#
			#write.table(unique(dr1_NeoAG_mhci_MinorAntigens_filtered[,c(1:24,33:36)]),"Summary_DR1_Annotato_MHCI_MinorAntigens_NoMismatches_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)
			#
			#
			
			
			
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
			
			
			
			#selsel <- grepl("IC50_HLA",names(Dr4_Binders_NoExprFilt))
			#Dr4_Binders_NoExprFilt[sel] <- lapply(Dr4_Binders_NoExprFilt[sel], function(x) replace(x,x > 500 , 0) )
			#Dr4_Binders_NoExprFilt=Dr4_Binders_NoExprFilt[,c(5,16:20)]
			#Dr4_Binders_NoExprFilt=melt(Dr4_Binders_NoExprFilt)
			#Dr4_Binders_NoExprFilt=Dr4_Binders_NoExprFilt[Dr4_Binders_NoExprFilt$value>0,]
			#
			#
			#
			#Dr4_Binders=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Dr4_Binders))
			#Dr4_Binders[sel] <- lapply(Dr4_Binders[sel], function(x) replace(x,x > 500 , 0) )
			#Dr4_Binders=Dr4_Binders[,c(5,16:20)]
			#Dr4_Binders=melt(Dr4_Binders)
			#Dr4_Binders=Dr4_Binders[Dr4_Binders$value>0,]
			#
			#
			#Dr4_Binders=rbind(strong_filtered,weak_filtered)
			#
			#Dr4_Binders_Coverages=Dr4_Binders[,c(1,2,3,31,33)]
			#
			#Dr4_Binders_Genes=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Dr4_Binders_Genes))
			#Dr4_Binders_Genes[sel] <- lapply(Dr4_Binders_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Dr4_Binders_Genes=Dr4_Binders_Genes[,c(1,2,3,16:20)]
			#Dr4_Binders_Genes=data.table::as.data.table(Dr4_Binders_Genes)
			#Dr4_Binders_Genes=data.table::melt.data.table(data=Dr4_Binders_Genes,id=1:3)
			#Dr4_Binders_Genes=Dr4_Binders_Genes[Dr4_Binders_Genes$value>0,]
			#
			#Dr4_Binders_Genes=unique(Dr4_Binders_Genes))
			#Dr4_Binders_Genes_AlleleCoverage=merge(Dr4_Binders_Genes,Dr4_Binders_Coverages,by=c('chr','pos','genes'),)
			#
			#
			#
			#
			#Dr4_Binders_NotExpressed=rbind(strong,weak)
			#
			#Dr4_Binders_NotExpressed_Coverages=Dr4_Binders_NotExpressed[,c(1,2,3,31,33)]
			#
			#Dr4_Binders_NotExpressed_Genes=rbind(strong,weak)
			#sel <- grepl("IC50_HLA",names(Dr4_Binders_NotExpressed_Genes))
			#Dr4_Binders_NotExpressed_Genes[sel] <- lapply(Dr4_Binders_NotExpressed_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Dr4_Binders_NotExpressed_Genes=Dr4_Binders_NotExpressed_Genes[,c(1,2,3,16:20)]
			#Dr4_Binders_NotExpressed_Genes=data.table::as.data.table(Dr4_Binders_NotExpressed_Genes)
			#Dr4_Binders_NotExpressed_Genes=data.table::melt.data.table(data=Dr4_Binders_NotExpressed_Genes,id=1:3)
			#Dr4_Binders_NotExpressed_Genes=Dr4_Binders_NotExpressed_Genes[Dr4_Binders_NotExpressed_Genes$value>0,]
			#
			#Dr4_Binders_NotExpressed_Genes=unique(Dr4_Binders_NotExpressed_Genes)
			#Dr4_Binders_NotExpressed_Genes_AlleleCoverage=merge(Dr4_Binders_NotExpressed_Genes,Dr4_Binders_NotExpressed_Coverages,by=c('chr','pos','genes'),allow.cartesian=T)
			#
			#
			#
			#
			#write.table(unique(dr4_NeoAG_mhci_MinorAntigens_filtered[,c(1:22,31:34)]),"Summary_DR4_Annotato_MHCI_MinorAntigens_NoMismatches_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)
			
			
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
			
			
			
			#Dr5_Binders_NoExprFilt=rbind(strong,weak)
			#sel <- grepl("IC50_HLA",names(Dr5_Binders_NoExprFilt))
			#Dr5_Binders_NoExprFilt[sel] <- lapply(Dr5_Binders_NoExprFilt[sel], function(x) replace(x,x > 500 , 0) )
			#Dr5_Binders_NoExprFilt=Dr5_Binders_NoExprFilt[,c(5,17:22)]
			#Dr5_Binders_NoExprFilt=melt(Dr5_Binders_NoExprFilt)
			#Dr5_Binders_NoExprFilt=Dr5_Binders_NoExprFilt[Dr5_Binders_NoExprFilt$value>0,]
			#
			#
			#
			#
			#Dr5_Binders=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Dr5_Binders))
			#Dr5_Binders[sel] <- lapply(Dr5_Binders[sel], function(x) replace(x,x > 500 , 0) )
			#Dr5_Binders=Dr5_Binders[,c(5,17:22)]
			#Dr5_Binders=melt(Dr5_Binders)
			#Dr5_Binders=Dr5_Binders[Dr5_Binders$value>0,]
			#
			#Dr5_Binders_Genes=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Dr5_Binders_Genes))
			#Dr5_Binders_Genes[sel] <- lapply(Dr5_Binders_Genes[sel], function(x) replace(x,x > 500 , 0) )
			#Dr5_Binders_Genes=Dr5_Binders_Genes[,c(3,17:22)]
			#Dr5_Binders_Genes=melt(Dr5_Binders_Genes)
			#Dr5_Binders_Genes=Dr5_Binders_Genes[Dr5_Binders_Genes$value>0,]
			#
			#
			#
			#Dr5_Binders=rbind(strong_filtered,weak_filtered)
			#
			#Dr5_Binders_Coverages=Dr5_Binders[,c(1,2,3,33,35)]
			#
			#Dr5_Binders_Genes=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Dr5_Binders_Genes))
			#Dr5_Binders_Genes[sel] <- lapply(Dr5_Binders_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Dr5_Binders_Genes=Dr5_Binders_Genes[,c(1,2,3,17:22)]
			#Dr5_Binders_Genes=data.table::as.data.table(Dr5_Binders_Genes)
			#Dr5_Binders_Genes=data.table::melt.data.table(data=Dr5_Binders_Genes,id=1:3)
			#Dr5_Binders_Genes=Dr5_Binders_Genes[Dr5_Binders_Genes$value>0,]
			#
			#Dr5_Binders_Genes=unique(Dr5_Binders_Genes))
			#Dr5_Binders_Genes_AlleleCoverage=merge(Dr5_Binders_Genes,Dr5_Binders_Coverages,by=c('chr','pos','genes'),allow.cartesian = T)
			#
			#
			#
			#
			#Dr5_Binders_NotExpressed=rbind(strong,weak)
			#
			#Dr5_Binders_NotExpressed_Coverages=Dr5_Binders_NotExpressed[,c(1,2,3,33,35)]
			#
			#Dr5_Binders_NotExpressed_Genes=rbind(strong,weak)
			#sel <- grepl("IC50_HLA",names(Dr5_Binders_NotExpressed_Genes))
			#Dr5_Binders_NotExpressed_Genes[sel] <- lapply(Dr5_Binders_NotExpressed_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Dr5_Binders_NotExpressed_Genes=Dr5_Binders_NotExpressed_Genes[,c(1,2,3,17:22)]
			#Dr5_Binders_NotExpressed_Genes=data.table::as.data.table(Dr5_Binders_NotExpressed_Genes)
			#Dr5_Binders_NotExpressed_Genes=data.table::melt.data.table(data=Dr5_Binders_NotExpressed_Genes,id=1:3)
			#Dr5_Binders_NotExpressed_Genes=Dr5_Binders_NotExpressed_Genes[Dr5_Binders_NotExpressed_Genes$value>0,]
			#
			#Dr5_Binders_NotExpressed_Genes=unique(Dr5_Binders_NotExpressed_Genes)
			#Dr5_Binders_NotExpressed_Genes_AlleleCoverage=merge(Dr5_Binders_NotExpressed_Genes,Dr5_Binders_NotExpressed_Coverages,by=c('chr','pos','genes'),allow.cartesian = T)
			#
			#
			#
			#
			#write.table(unique(dr5_NeoAG_mhci_MinorAntigens_filtered[,c(1:24,33:36)]),"Summary_DR5_Annotato_MHCI_MinorAntigens_NoMismatches_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)
			#
			
			
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
			
			
			
			
			#sel <- grepl("IC50_HLA",names(Foca_Binders_NoExprFilt))
			#Foca_Binders_NoExprFilt[sel] <- lapply(Foca_Binders_NoExprFilt[sel], function(x) replace(x,x >500 , 0) )
			#Foca_Binders_NoExprFilt=Foca_Binders_NoExprFilt[,c(5,17:22)]
			#Foca_Binders_NoExprFilt=melt(Foca_Binders_NoExprFilt)
			#Foca_Binders_NoExprFilt=Foca_Binders_NoExprFilt[Foca_Binders_NoExprFilt$value>0,]
			#
			#
			#
			#
			#Foca_Binders=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Foca_Binders))
			#Foca_Binders[sel] <- lapply(Foca_Binders[sel], function(x) replace(x,x >500 , 0) )
			##Foca_Binders$IC50_Binder=rowSums(Foca_Binders[sel])
			##Foca_Binders=Foca_Binders[order(Foca_Binders$IC50_Binder),]
			##Foca_Binders=Foca_Binders[Foca_Binders$IC50_Binder>0,]
			#Foca_Binders=Foca_Binders[,c(5,17:22)]
			#Foca_Binders=melt(Foca_Binders)
			#Foca_Binders=Foca_Binders[Foca_Binders$value>0,]
			#
			#Foca_Binders_Genes=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Foca_Binders_Genes))
			#Foca_Binders_Genes[sel] <- lapply(Foca_Binders_Genes[sel], function(x) replace(x,x >500 , 0) )
			##Foca_Binders_Genes$IC50_Binder=rowSums(Foca_Binders_Genes[sel])
			##Foca_Binders_Genes=Foca_Binders_Genes[order(Foca_Binders_Genes$IC50_Binder),]
			##Foca_Binders_Genes=Foca_Binders_Genes[Foca_Binders_Genes$IC50_Binder>0,]
			#Foca_Binders_Genes=Foca_Binders_Genes[,c(3,17:22)]
			#Foca_Binders_Genes=melt(Foca_Binders_Genes)
			#Foca_Binders_Genes=Foca_Binders_Genes[Foca_Binders_Genes$value>0,]
			#
			#Foca_Binders=rbind(strong_filtered,weak_filtered)
			#
			#Foca_Binders_Coverages=Foca_Binders[,c(1,2,3,33,35)]
			#
			#Foca_Binders_Genes=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Foca_Binders_Genes))
			#Foca_Binders_Genes[sel] <- lapply(Foca_Binders_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Foca_Binders_Genes=Foca_Binders_Genes[,c(1,2,3,17:22)]
			#Foca_Binders_Genes=data.table::as.data.table(Foca_Binders_Genes)
			#Foca_Binders_Genes=data.table::melt.data.table(data=Foca_Binders_Genes,id=1:3)
			#Foca_Binders_Genes=Foca_Binders_Genes[Foca_Binders_Genes$value>0,]
			#
			#Foca_Binders_Genes=unique(Foca_Binders_Genes)
			#Foca_Binders_Genes_AlleleCoverage=merge(Foca_Binders_Genes,Foca_Binders_Coverages,by=c('chr','pos','genes'),allow.cartesian = T)
			#Foca_Binders_Genes_AlleleCoverage=unique(Foca_Binders_Genes_AlleleCoverage)
			#
			#
			#Foca_Binders_NotExpressed=rbind(strong,weak)
			#
			#Foca_Binders_NotExpressed_Coverages=Foca_Binders_NotExpressed[,c(1,2,3,33,35)]
			#
			#Foca_Binders_NotExpressed_Genes=rbind(strong,weak)
			#sel <- grepl("IC50_HLA",names(Foca_Binders_NotExpressed_Genes))
			#Foca_Binders_NotExpressed_Genes[sel] <- lapply(Foca_Binders_NotExpressed_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Foca_Binders_NotExpressed_Genes=Foca_Binders_NotExpressed_Genes[,c(1,2,3,17:22)]
			#Foca_Binders_NotExpressed_Genes=data.table::as.data.table(Foca_Binders_NotExpressed_Genes)
			#Foca_Binders_NotExpressed_Genes=data.table::melt.data.table(data=Foca_Binders_NotExpressed_Genes,id=1:3)
			#Foca_Binders_NotExpressed_Genes=Foca_Binders_NotExpressed_Genes[Foca_Binders_NotExpressed_Genes$value>0,]
			#
			#Foca_Binders_NotExpressed_Genes=unique(Foca_Binders_NotExpressed_Genes)
			#Foca_Binders_NotExpressed_Genes_AlleleCoverage=merge(Foca_Binders_NotExpressed_Genes,Foca_Binders_NotExpressed_Coverages,by=c('chr','pos','genes'),allow.cartesian = T)
			#Foca_Binders_NotExpressed_Genes_AlleleCoverage=unique(Foca_Binders_NotExpressed_Genes_AlleleCoverage)
			#
			#
			#
			#
			#
			#write.table(unique(foca_NeoAG_mhci_MinorAntigens_filtered[,c(1:26,35:38)]),"Summary_FOCA_Annotato_MHCI_MinorAntigens_NoMismatches_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)
			
			
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
			
			
			
			
			#Gagra_Binders_NoExprFilt=rbind(strong,weak)
			#sel <- grepl("IC50_HLA",names(Gagra_Binders_NoExprFilt))
			#Gagra_Binders_NoExprFilt[sel] <- lapply(Gagra_Binders_NoExprFilt[sel], function(x) replace(x,x >500 , 0) )
			#Gagra_Binders_NoExprFilt=Gagra_Binders_NoExprFilt[,c(5,17:22)]
			#Gagra_Binders_NoExprFilt=melt(Gagra_Binders_NoExprFilt)
			#Gagra_Binders_NoExprFilt=Gagra_Binders_NoExprFilt[Gagra_Binders_NoExprFilt$value>0,]
			#
			#
			#Gagra_Binders=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Gagra_Binders))
			#Gagra_Binders[sel] <- lapply(Gagra_Binders[sel], function(x) replace(x,x >500 , 0) )
			#Gagra_Binders=Gagra_Binders[,c(5,17:22)]
			#Gagra_Binders=melt(Gagra_Binders)
			#Gagra_Binders=Gagra_Binders[Gagra_Binders$value>0,]
			#
			#Gagra_Binders_Genes=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Gagra_Binders_Genes))
			#Gagra_Binders_Genes[sel] <- lapply(Gagra_Binders_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Gagra_Binders_Genes=Gagra_Binders_Genes[,c(3,17:22)]
			#Gagra_Binders_Genes=melt(Gagra_Binders_Genes)
			#Gagra_Binders_Genes=Gagra_Binders_Genes[Gagra_Binders_Genes$value>0,]
			#
			#
			#Gagra_Binders=rbind(strong_filtered,weak_filtered)
			#
			#Gagra_Binders_Coverages=Gagra_Binders[,c(1,2,3,33,35)]
			#
			#Gagra_Binders_Genes=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Gagra_Binders_Genes))
			#Gagra_Binders_Genes[sel] <- lapply(Gagra_Binders_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Gagra_Binders_Genes=Gagra_Binders_Genes[,c(1,2,3,17:22)]
			#Gagra_Binders_Genes=data.table::as.data.table(Gagra_Binders_Genes)
			#Gagra_Binders_Genes=data.table::melt.data.table(data=Gagra_Binders_Genes,id=1:3)
			#Gagra_Binders_Genes=Gagra_Binders_Genes[Gagra_Binders_Genes$value>0,]
			#
			#Gagra_Binders_Genes=unique(Gagra_Binders_Genes)
			#Gagra_Binders_Genes_AlleleCoverage=merge(Gagra_Binders_Genes,Gagra_Binders_Coverages,by=c('chr','pos','genes'),allow.cartesian = T)
			#Gagra_Binders_Genes_AlleleCoverage=unique(Gagra_Binders_Genes_AlleleCoverage)
			#
			#
			#Gagra_Binders_NotExpressed=rbind(strong,weak)
			#
			#Gagra_Binders_NotExpressed_Coverages=Gagra_Binders_NotExpressed[,c(1,2,3,33,35)]
			#
			#Gagra_Binders_NotExpressed_Genes=rbind(strong,weak)
			#sel <- grepl("IC50_HLA",names(Gagra_Binders_NotExpressed_Genes))
			#Gagra_Binders_NotExpressed_Genes[sel] <- lapply(Gagra_Binders_NotExpressed_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Gagra_Binders_NotExpressed_Genes=Gagra_Binders_NotExpressed_Genes[,c(1,2,3,17:22)]
			#Gagra_Binders_NotExpressed_Genes=data.table::as.data.table(Gagra_Binders_NotExpressed_Genes)
			#Gagra_Binders_NotExpressed_Genes=data.table::melt.data.table(data=Gagra_Binders_NotExpressed_Genes,id=1:3)
			#Gagra_Binders_NotExpressed_Genes=Gagra_Binders_NotExpressed_Genes[Gagra_Binders_NotExpressed_Genes$value>0,]
			#
			#Gagra_Binders_NotExpressed_Genes=unique(Gagra_Binders_NotExpressed_Genes)
			#Gagra_Binders_NotExpressed_Genes_AlleleCoverage=merge(Gagra_Binders_NotExpressed_Genes,Gagra_Binders_NotExpressed_Coverages,by=c('chr','pos','genes'),allow.cartesian = T)
			#Gagra_Binders_NotExpressed_Genes_AlleleCoverage=unique(Gagra_Binders_NotExpressed_Genes_AlleleCoverage)
			#
			#
			#
			#
			#
			#write.table(unique(gagra_NeoAG_mhci_MinorAntigens_filtered[,c(1:24,33:36)]),"Summary_GAGRA_Annotato_MHCI_MinorAntigens_NoMismatches_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)
			#
			
			
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
			
			
			#
			#
			#Luan_Binders_NoExprFilt=rbind(strong,weak)
			#sel <- grepl("IC50_HLA",names(Luan_Binders_NoExprFilt))
			#Luan_Binders_NoExprFilt[sel] <- lapply(Luan_Binders_NoExprFilt[sel], function(x) replace(x,x > 500 , 0) )
			#Luan_Binders_NoExprFilt=Luan_Binders_NoExprFilt[,c(5,16:20)]
			#Luan_Binders_NoExprFilt=melt(Luan_Binders_NoExprFilt)
			#Luan_Binders_NoExprFilt=Luan_Binders_NoExprFilt[Luan_Binders_NoExprFilt$value>0,]
			#
			#
			#
			#
			#Luan_Binders=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Luan_Binders))
			#Luan_Binders[sel] <- lapply(Luan_Binders[sel], function(x) replace(x,x > 500 , 0) )
			#Luan_Binders=Luan_Binders[,c(5,16:20)]
			#Luan_Binders=melt(Luan_Binders)
			#Luan_Binders=Luan_Binders[Luan_Binders$value>0,]
			#
			#
			#Luan_Binders_Genes=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Luan_Binders_Genes))
			#Luan_Binders_Genes[sel] <- lapply(Luan_Binders_Genes[sel], function(x) replace(x,x > 500 , 0) )
			#Luan_Binders_Genes=Luan_Binders_Genes[,c(3,16:20)]
			#Luan_Binders_Genes=melt(Luan_Binders_Genes)
			#Luan_Binders_Genes=Luan_Binders_Genes[Luan_Binders_Genes$value>0,]
			#
			#
			#
			#
			#Luan_Binders=rbind(strong_filtered,weak_filtered)
			#
			#Luan_Binders_Coverages=Luan_Binders[,c(1,2,3,31,33)]
			#
			#Luan_Binders_Genes=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Luan_Binders_Genes))
			#Luan_Binders_Genes[sel] <- lapply(Luan_Binders_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Luan_Binders_Genes=Luan_Binders_Genes[,c(1,2,3,16:20)]
			#Luan_Binders_Genes=data.table::as.data.table(Luan_Binders_Genes)
			#Luan_Binders_Genes=data.table::melt.data.table(data=Luan_Binders_Genes,id=1:3)
			#Luan_Binders_Genes=Luan_Binders_Genes[Luan_Binders_Genes$value>0,]
			#
			#Luan_Binders_Genes=unique(Luan_Binders_Genes)
			#Luan_Binders_Genes_AlleleCoverage=merge(Luan_Binders_Genes,Luan_Binders_Coverages,by=c('chr','pos','genes'),allow.cartesian = T)
			#Luan_Binders_Genes_AlleleCoverage=unique(Luan_Binders_Genes_AlleleCoverage)
			#
			#
			#Luan_Binders_NotExpressed=rbind(strong,weak)
			#
			#Luan_Binders_NotExpressed_Coverages=Luan_Binders_NotExpressed[,c(1,2,3,31,33)]
			#
			#Luan_Binders_NotExpressed_Genes=rbind(strong,weak)
			#sel <- grepl("IC50_HLA",names(Luan_Binders_NotExpressed_Genes))
			#Luan_Binders_NotExpressed_Genes[sel] <- lapply(Luan_Binders_NotExpressed_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Luan_Binders_NotExpressed_Genes=Luan_Binders_NotExpressed_Genes[,c(1,2,3,16:20)]
			#Luan_Binders_NotExpressed_Genes=data.table::as.data.table(Luan_Binders_NotExpressed_Genes)
			#Luan_Binders_NotExpressed_Genes=data.table::melt.data.table(data=Luan_Binders_NotExpressed_Genes,id=1:3)
			#Luan_Binders_NotExpressed_Genes=Luan_Binders_NotExpressed_Genes[Luan_Binders_NotExpressed_Genes$value>0,]
			#
			#Luan_Binders_NotExpressed_Genes=unique(Luan_Binders_NotExpressed_Genes)
			#Luan_Binders_NotExpressed_Genes_AlleleCoverage=merge(Luan_Binders_NotExpressed_Genes,Luan_Binders_NotExpressed_Coverages,by=c('chr','pos','genes'),allow.cartesian = T)
			#Luan_Binders_NotExpressed_Genes_AlleleCoverage=unique(Luan_Binders_NotExpressed_Genes_AlleleCoverage)
			#
			#
			#
			#
			#write.table(unique(luan_NeoAG_mhci_MinorAntigens_filtered[,c(1:22,31:34)]),"Summary_LUAN_Annotato_MHCI_MinorAntigens_NoMismatches_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)
			#
			#
			
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
			
			
			
			#
			#Mabi_Binders_NoExprFilt=rbind(strong,weak)
			#sel <- grepl("IC50_HLA",names(Mabi_Binders_NoExprFilt))
			#Mabi_Binders_NoExprFilt[sel] <- lapply(Mabi_Binders_NoExprFilt[sel], function(x) replace(x,x > 500 , 0) )
			#Mabi_Binders_NoExprFilt=Mabi_Binders_NoExprFilt[,c(5,14:16)]
			#Mabi_Binders_NoExprFilt=melt(Mabi_Binders_NoExprFilt)
			#Mabi_Binders_NoExprFilt=Mabi_Binders_NoExprFilt[Mabi_Binders_NoExprFilt$value>0,]
			#
			#Mabi_Binders=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Mabi_Binders))
			#Mabi_Binders[sel] <- lapply(Mabi_Binders[sel], function(x) replace(x,x > 500 , 0) )
			#Mabi_Binders=Mabi_Binders[,c(5,14:16)]
			#Mabi_Binders=melt(Mabi_Binders)
			#Mabi_Binders=Mabi_Binders[Mabi_Binders$value>0,]
			#
			#
			#Mabi_Binders_Genes=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Mabi_Binders_Genes))
			#Mabi_Binders_Genes[sel] <- lapply(Mabi_Binders_Genes[sel], function(x) replace(x,x > 500 , 0) )
			#Mabi_Binders_Genes=Mabi_Binders_Genes[,c(3,14:16)]
			#Mabi_Binders_Genes=melt(Mabi_Binders_Genes)
			#Mabi_Binders_Genes=Mabi_Binders_Genes[Mabi_Binders_Genes$value>0,]
			#
			#
			#
			#Mabi_Binders=rbind(strong_filtered,weak_filtered)
			#
			#Mabi_Binders_Coverages=Mabi_Binders[,c(1,2,3,27,29)]
			#
			#Mabi_Binders_Genes=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Mabi_Binders_Genes))
			#Mabi_Binders_Genes[sel] <- lapply(Mabi_Binders_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Mabi_Binders_Genes=Mabi_Binders_Genes[,c(1,2,3,14:16)]
			#Mabi_Binders_Genes=data.table::as.data.table(Mabi_Binders_Genes)
			#Mabi_Binders_Genes=data.table::melt.data.table(data=Mabi_Binders_Genes,id=1:3)
			#Mabi_Binders_Genes=Mabi_Binders_Genes[Mabi_Binders_Genes$value>0,]
			#
			#Mabi_Binders_Genes=unique(Mabi_Binders_Genes)
			#Mabi_Binders_Genes_AlleleCoverage=merge(Mabi_Binders_Genes,Mabi_Binders_Coverages,by=c('chr','pos','genes'),allow.cartesian = T)
			#Mabi_Binders_Genes_AlleleCoverage=unique(Mabi_Binders_Genes_AlleleCoverage)
			#
			#
			#
			#Mabi_Binders_NotExpressed=rbind(strong,weak)
			#
			#Mabi_Binders_NotExpressed_Coverages=Mabi_Binders_NotExpressed[,c(1,2,3,27,29)]
			#
			#Mabi_Binders_NotExpressed_Genes=rbind(strong,weak)
			#sel <- grepl("IC50_HLA",names(Mabi_Binders_NotExpressed_Genes))
			#Mabi_Binders_NotExpressed_Genes[sel] <- lapply(Mabi_Binders_NotExpressed_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Mabi_Binders_NotExpressed_Genes=Mabi_Binders_NotExpressed_Genes[,c(1,2,3,14:16)]
			#Mabi_Binders_NotExpressed_Genes=data.table::as.data.table(Mabi_Binders_NotExpressed_Genes)
			#Mabi_Binders_NotExpressed_Genes=data.table::melt.data.table(data=Mabi_Binders_NotExpressed_Genes,id=1:3)
			#Mabi_Binders_NotExpressed_Genes=Mabi_Binders_NotExpressed_Genes[Mabi_Binders_NotExpressed_Genes$value>0,]
			#
			#Mabi_Binders_NotExpressed_Genes=unique(Mabi_Binders_NotExpressed_Genes)
			#Mabi_Binders_NotExpressed_Genes_AlleleCoverage=merge(Mabi_Binders_NotExpressed_Genes,Mabi_Binders_NotExpressed_Coverages,by=c('chr','pos','genes'),allow.cartesian = T)
			#Mabi_Binders_NotExpressed_Genes_AlleleCoverage=unique(Mabi_Binders_NotExpressed_Genes_AlleleCoverage)
			#
			#
			#
			#write.table(unique(mabi_NeoAG_mhci_MinorAntigens_filtered[,c(1:18,27:30)]),"Summary_MABI_Annotato_MHCI_MinorAntigens_NoMismatches_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)
			#
			
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
			
			
			#Moge_Binders_NoExprFilt=rbind(strong,weak)
			#sel <- grepl("IC50_HLA",names(Moge_Binders_NoExprFilt))
			#Moge_Binders_NoExprFilt[sel] <- lapply(Moge_Binders_NoExprFilt[sel], function(x) replace(x,x >500 , 0) )
			#Moge_Binders_NoExprFilt=Moge_Binders_NoExprFilt[,c(5,16:20)]
			#Moge_Binders_NoExprFilt=melt(Moge_Binders_NoExprFilt)
			#Moge_Binders_NoExprFilt=Moge_Binders_NoExprFilt[Moge_Binders_NoExprFilt$value>0,]
			#
			#Moge_Binders=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Moge_Binders))
			#Moge_Binders[sel] <- lapply(Moge_Binders[sel], function(x) replace(x,x >500 , 0) )
			#Moge_Binders=Moge_Binders[,c(5,16:20)]
			#Moge_Binders=melt(Moge_Binders)
			#Moge_Binders=Moge_Binders[Moge_Binders$value>0,]
			#
			#
			#Moge_Binders_Genes=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Moge_Binders_Genes))
			#Moge_Binders_Genes[sel] <- lapply(Moge_Binders_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Moge_Binders_Genes=Moge_Binders_Genes[,c(3,16:20)]
			#Moge_Binders_Genes=melt(Moge_Binders_Genes)
			#Moge_Binders_Genes=Moge_Binders_Genes[Moge_Binders_Genes$value>0,]
			#
			#
			#
			#
			#Moge_Binders=rbind(strong_filtered,weak_filtered)
			#
			#Moge_Binders_Coverages=Moge_Binders[,c(1,2,3,31,33)]
			#
			#Moge_Binders_Genes=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Moge_Binders_Genes))
			#Moge_Binders_Genes[sel] <- lapply(Moge_Binders_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Moge_Binders_Genes=Moge_Binders_Genes[,c(1,2,3,16:20)]
			#Moge_Binders_Genes=data.table::as.data.table(Moge_Binders_Genes)
			#Moge_Binders_Genes=data.table::melt.data.table(data=Moge_Binders_Genes,id=1:3)
			#Moge_Binders_Genes=Moge_Binders_Genes[Moge_Binders_Genes$value>0,]
			#
			#Moge_Binders_Genes=unique(Moge_Binders_Genes)
			#Moge_Binders_Genes_AlleleCoverage=merge(Moge_Binders_Genes,Moge_Binders_Coverages,by=c('chr','pos','genes'),allow.cartesian = T)
			#Moge_Binders_Genes_AlleleCoverage=unique(Moge_Binders_Genes_AlleleCoverage)
			#
			#
			#
			#Moge_Binders_NotExpressed=rbind(strong,weak)
			#
			#Moge_Binders_NotExpressed_Coverages=Moge_Binders_NotExpressed[,c(1,2,3,31,33)]
			#
			#Moge_Binders_NotExpressed_Genes=rbind(strong,weak)
			#sel <- grepl("IC50_HLA",names(Moge_Binders_NotExpressed_Genes))
			#Moge_Binders_NotExpressed_Genes[sel] <- lapply(Moge_Binders_NotExpressed_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Moge_Binders_NotExpressed_Genes=Moge_Binders_NotExpressed_Genes[,c(1,2,3,16:20)]
			#Moge_Binders_NotExpressed_Genes=data.table::as.data.table(Moge_Binders_NotExpressed_Genes)
			#Moge_Binders_NotExpressed_Genes=data.table::melt.data.table(data=Moge_Binders_NotExpressed_Genes,id=1:3)
			#Moge_Binders_NotExpressed_Genes=Moge_Binders_NotExpressed_Genes[Moge_Binders_NotExpressed_Genes$value>0,]
			#
			#Moge_Binders_NotExpressed_Genes=unique(Moge_Binders_NotExpressed_Genes)
			#Moge_Binders_NotExpressed_Genes_AlleleCoverage=merge(Moge_Binders_NotExpressed_Genes,Moge_Binders_NotExpressed_Coverages,by=c('chr','pos','genes'),allow.cartesian = T)
			#Moge_Binders_NotExpressed_Genes_AlleleCoverage=unique(Moge_Binders_NotExpressed_Genes_AlleleCoverage)
			#
			#
			#
			#
			#write.table(unique(moge_NeoAG_mhci_MinorAntigens_filtered[,c(1:22,31:34)]),"Summary_MOGE_Annotato_MHCI_MinorAntigens_NoMismatches_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)
			
			
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
			
			
			
			
			
			
			
			
			#Piag_Binders_NoExprFilt=rbind(strong,weak)
			#sel <- grepl("IC50_HLA",names(Piag_Binders_NoExprFilt))
			#Piag_Binders_NoExprFilt[sel] <- lapply(Piag_Binders_NoExprFilt[sel], function(x) replace(x,x >500 , 0) )
			#Piag_Binders_NoExprFilt=Piag_Binders_NoExprFilt[,c(5,16:20)]
			#Piag_Binders_NoExprFilt=melt(Piag_Binders_NoExprFilt)
			#Piag_Binders_NoExprFilt=Piag_Binders_NoExprFilt[Piag_Binders_NoExprFilt$value>0,]
			#
			#Piag_Binders=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Piag_Binders))
			#Piag_Binders[sel] <- lapply(Piag_Binders[sel], function(x) replace(x,x >500 , 0) )
			#Piag_Binders=Piag_Binders[,c(5,16:20)]
			#Piag_Binders=melt(Piag_Binders)
			#Piag_Binders=Piag_Binders[Piag_Binders$value>0,]
			#
			#Piag_Binders_Genes=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Piag_Binders_Genes))
			#Piag_Binders_Genes[sel] <- lapply(Piag_Binders_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Piag_Binders_Genes=Piag_Binders_Genes[,c(3,16:20)]
			#Piag_Binders_Genes=melt(Piag_Binders_Genes)
			#Piag_Binders_Genes=Piag_Binders_Genes[Piag_Binders_Genes$value>0,]
			#
			#
			#Piag_Binders=rbind(strong_filtered,weak_filtered)
			#
			#Piag_Binders_Coverages=Piag_Binders[,c(1,2,3,31,33)]
			#
			#Piag_Binders_Genes=rbind(strong_filtered,weak_filtered)
			#sel <- grepl("IC50_HLA",names(Piag_Binders_Genes))
			#Piag_Binders_Genes[sel] <- lapply(Piag_Binders_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Piag_Binders_Genes=Piag_Binders_Genes[,c(1,2,3,16:20)]
			#Piag_Binders_Genes=data.table::as.data.table(Piag_Binders_Genes)
			#Piag_Binders_Genes=data.table::melt.data.table(data=Piag_Binders_Genes,id=1:3)
			#Piag_Binders_Genes=Piag_Binders_Genes[Piag_Binders_Genes$value>0,]
			#
			#Piag_Binders_Genes=unique(Piag_Binders_Genes)
			#Piag_Binders_Genes_AlleleCoverage=merge(Piag_Binders_Genes,Piag_Binders_Coverages,by=c('chr','pos','genes'),allow.cartesian = T)
			#Piag_Binders_Genes_AlleleCoverage=unique(Piag_Binders_Genes_AlleleCoverage)
			#
			#
			#Piag_Binders_NotExpressed=rbind(strong,weak)
			#
			#Piag_Binders_NotExpressed_Coverages=Piag_Binders_NotExpressed[,c(1,2,3,31,33)]
			#
			#Piag_Binders_NotExpressed_Genes=rbind(strong,weak)
			#sel <- grepl("IC50_HLA",names(Piag_Binders_NotExpressed_Genes))
			#Piag_Binders_NotExpressed_Genes[sel] <- lapply(Piag_Binders_NotExpressed_Genes[sel], function(x) replace(x,x >500 , 0) )
			#Piag_Binders_NotExpressed_Genes=Piag_Binders_NotExpressed_Genes[,c(1,2,3,16:20)]
			#Piag_Binders_NotExpressed_Genes=data.table::as.data.table(Piag_Binders_NotExpressed_Genes)
			#Piag_Binders_NotExpressed_Genes=data.table::melt.data.table(data=Piag_Binders_NotExpressed_Genes,id=1:3)
			#Piag_Binders_NotExpressed_Genes=Piag_Binders_NotExpressed_Genes[Piag_Binders_NotExpressed_Genes$value>0,]
			#
			#Piag_Binders_NotExpressed_Genes=unique(Piag_Binders_NotExpressed_Genes)
			#Piag_Binders_NotExpressed_Genes_AlleleCoverage=merge(Piag_Binders_NotExpressed_Genes,Piag_Binders_NotExpressed_Coverages,by=c('chr','pos','genes'),allow.cartesian = T)
			#Piag_Binders_NotExpressed_Genes_AlleleCoverage=unique(Piag_Binders_NotExpressed_Genes_AlleleCoverage)
			#
			#
			#
			#
			#
			#write.table(unique(piag_NeoAG_mhci_MinorAntigens_filtered[,c(1:22,31:34)]),"Summary_PIAG_Annotato_MHCI_MinorAntigens_NoMismatches_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)
			
			
			
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
			
			
			
			
			alfe_auc  =auc(1:length(Alfe_Binders_Ic50_NoExprFilt_ord$value), y=Alfe_Binders_Ic50_NoExprFilt_ord$value)
			besu_auc  =auc(1:length(Besu_Binders_Ic50_NoExprFilt_ord$value), y=Besu_Binders_Ic50_NoExprFilt_ord$value)
			calu_auc  =auc(1:length(Calu_Binders_Ic50_NoExprFilt_ord_homo$value), y=Calu_Binders_Ic50_NoExprFilt_ord_homo$value)
			deiv_auc  =auc(1:length(Deiv_Binders_Ic50_NoExprFilt_ord$value), y=Deiv_Binders_Ic50_NoExprFilt_ord$value)
			dest_auc  =auc(1:length(Dest_Binders_Ic50_NoExprFilt_ord$value), y=Dest_Binders_Ic50_NoExprFilt_ord$value)
			dr1_auc   =auc(1:length(Dr1_Binders_Ic50_NoExprFilt_ord$value), y=Dr1_Binders_Ic50_NoExprFilt_ord$value)
			dr4_auc   =auc(1:length(Dr4_Binders_Ic50_NoExprFilt_ord$value), y=Dr4_Binders_Ic50_NoExprFilt_ord$value)
			dr5_auc   =auc(1:length(Dr5_Binders_Ic50_NoExprFilt_ord$value), y=Dr5_Binders_Ic50_NoExprFilt_ord$value)
			foca_auc  =auc(1:length(Foca_Binders_Ic50_NoExprFilt_ord$value), y=Foca_Binders_Ic50_NoExprFilt_ord$value)
			gagra_auc =auc(1:length(Gagra_Binders_Ic50_NoExprFilt_ord$value), y=Gagra_Binders_Ic50_NoExprFilt_ord$value)
			luan_auc  =auc(1:length(Luan_Binders_Ic50_NoExprFilt_ord_homo$value), y=Luan_Binders_Ic50_NoExprFilt_ord_homo$value)
			mabi_auc  =auc(1:length(Mabi_Binders_Ic50_NoExprFilt_ord$value), y=Mabi_Binders_Ic50_NoExprFilt_ord$value)
			moge_auc  =auc(1:length(Moge_Binders_Ic50_NoExprFilt_ord_homo$value), y=Moge_Binders_Ic50_NoExprFilt_ord_homo$value)
			piag_auc  =auc(1:length(Piag_Binders_Ic50_NoExprFilt_ord_homo$value), y=Piag_Binders_Ic50_NoExprFilt_ord_homo$value)
			prelu_auc =auc(1:length(Prelu_Binders_Ic50_NoExprFilt_ord$value), y=Prelu_Binders_Ic50_NoExprFilt_ord$value)
			
			
			all_aucs=rbind(alfe_auc,besu_auc,calu_auc,deiv_auc,dest_auc,dr1_auc,dr4_auc,dr5_auc,foca_auc,gagra_auc,luan_auc,mabi_auc,moge_auc,piag_auc,prelu_auc)
			all_aucs=as.data.frame(all_aucs)
			all_aucs$genes=rownames(all_aucs)
			all_aucs$RelDon=c(1,0,1,0,1,0,1,1,0,0,1,1,0,0,1)
			#all_aucs$RelDon=c(0,1,0,1,0,1,0,0,1,1,0,0,1,1,0)
			all_aucs=all_aucs[order(all_aucs$V1),]
			colnames(all_aucs)=c('aucs','patients','RelDon')
			all_aucs$Donors=ifelse(all_aucs$RelDon==1,'Related','Unrelated')
			
			ggplot(all_aucs,aes(x=seq(1,length(all_aucs$aucs)),y=aucs,col=Donors)) +
			geom_point()  +
			geom_rug() +
			scale_color_manual(values=c( Related= "#E69F00", Unrelated="darkorchid2")) +
			geom_smooth( method = "lm",mapping = aes(fill=all_aucs$Donors)) +
			scale_fill_manual(values=c(Related= "#E69F00", Unrelated="darkorchid2")) +
			ggtitle('MRD vs MUD AUCs Binders\n linear regression') +
			xlab('Patients')
			
			
			plot(log(all_aucs$aucs),col=as.factor(all_aucs$Donor),lwd=3)
			text(log(all_aucs$aucs),labels=all_aucs$patients, cex= 0.7,adj=1,pos=3,col=color[1:15])
			
			summary(glm(factor(all_aucs$Donors) ~ all_aucs$aucs,family = binomial(link = "probit")))
			
			
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
			



#Prelu_Binders_NoExprFilt=rbind(strong,weak)
#sel <- grepl("IC50_HLA",names(Prelu_Binders_NoExprFilt))
#Prelu_Binders_NoExprFilt[sel] <- lapply(Prelu_Binders_NoExprFilt[sel], function(x) replace(x,x > 500 , 0) )
#Prelu_Binders_NoExprFilt=Prelu_Binders_NoExprFilt[,c(5,14:16)]
#Prelu_Binders_NoExprFilt=melt(Prelu_Binders_NoExprFilt)
#Prelu_Binders_NoExprFilt=Prelu_Binders_NoExprFilt[Prelu_Binders_NoExprFilt$value>0,]
#
#
#Prelu_Binders=rbind(strong_filtered,weak_filtered)
#sel <- grepl("IC50_HLA",names(Prelu_Binders))
#Prelu_Binders[sel] <- lapply(Prelu_Binders[sel], function(x) replace(x,x > 500 , 0) )
#Prelu_Binders=Prelu_Binders[,c(5,14:16)]
#Prelu_Binders=melt(Prelu_Binders)
#Prelu_Binders=Prelu_Binders[Prelu_Binders$value>0,]
#
#
#Prelu_Binders_Genes=rbind(strong_filtered,weak_filtered)
#sel <- grepl("IC50_HLA",names(Prelu_Binders_Genes))
#Prelu_Binders_Genes[sel] <- lapply(Prelu_Binders_Genes[sel], function(x) replace(x,x > 500 , 0) )
#Prelu_Binders_Genes=Prelu_Binders_Genes[,c(3,14:16)]
#Prelu_Binders_Genes=melt(Prelu_Binders_Genes)
#Prelu_Binders_Genes=Prelu_Binders_Genes[Prelu_Binders_Genes$value>0,]
#
#
#Prelu_Binders=rbind(strong_filtered,weak_filtered)
#
#Prelu_Binders_Coverages=Prelu_Binders[,c(1,2,3,27,29)]
#
#Prelu_Binders_Genes=rbind(strong_filtered,weak_filtered)
#sel <- grepl("IC50_HLA",names(Prelu_Binders_Genes))
#Prelu_Binders_Genes[sel] <- lapply(Prelu_Binders_Genes[sel], function(x) replace(x,x >500 , 0) )
#Prelu_Binders_Genes=Prelu_Binders_Genes[,c(1,2,3,14:16)]
#Prelu_Binders_Genes=data.table::as.data.table(Prelu_Binders_Genes)
#Prelu_Binders_Genes=data.table::melt.data.table(data=Prelu_Binders_Genes,id=1:3)
#Prelu_Binders_Genes=Prelu_Binders_Genes[Prelu_Binders_Genes$value>0,]
#
#Prelu_Binders_Genes=unique(Prelu_Binders_Genes)
#Prelu_Binders_Genes_AlleleCoverage=merge(Prelu_Binders_Genes,Prelu_Binders_Coverages,by=c('chr','pos','genes'),allow.cartesian = T)
#Prelu_Binders_Genes_AlleleCoverage=unique(Prelu_Binders_Genes_AlleleCoverage)
#
#
#Prelu_Binders_NotExpressed=rbind(strong,weak)
#
#Prelu_Binders_NotExpressed_Coverages=Prelu_Binders_NotExpressed[,c(1,2,3,27,29)]
#
#Prelu_Binders_NotExpressed_Genes=rbind(strong,weak)
#sel <- grepl("IC50_HLA",names(Prelu_Binders_NotExpressed_Genes))
#Prelu_Binders_NotExpressed_Genes[sel] <- lapply(Prelu_Binders_NotExpressed_Genes[sel], function(x) replace(x,x >500 , 0) )
#Prelu_Binders_NotExpressed_Genes=Prelu_Binders_NotExpressed_Genes[,c(1,2,3,14:16)]
#Prelu_Binders_NotExpressed_Genes=data.table::as.data.table(Prelu_Binders_NotExpressed_Genes)
#Prelu_Binders_NotExpressed_Genes=data.table::melt.data.table(data=Prelu_Binders_NotExpressed_Genes,id=1:3)
#Prelu_Binders_NotExpressed_Genes=Prelu_Binders_NotExpressed_Genes[Prelu_Binders_NotExpressed_Genes$value>0,]
#
#Prelu_Binders_NotExpressed_Genes=unique(Prelu_Binders_NotExpressed_Genes)
#Prelu_Binders_NotExpressed_Genes_AlleleCoverage=merge(Prelu_Binders_NotExpressed_Genes,Prelu_Binders_NotExpressed_Coverages,by=c('chr','pos','genes'),allow.cartesian = T)
#Prelu_Binders_NotExpressed_Genes_AlleleCoverage=unique(Prelu_Binders_NotExpressed_Genes_AlleleCoverage)
#
#





#Alfe_Binders_NoExprFilt_ord=Alfe_Binders_NoExprFilt[order(Alfe_Binders_NoExprFilt$value),]
#Besu_Binders_NoExprFilt_ord=Besu_Binders_NoExprFilt[order(Besu_Binders_NoExprFilt$value),]
#Calu_Binders_NoExprFilt_ord=Calu_Binders_NoExprFilt[order(Calu_Binders_NoExprFilt$value),]
#Deiv_Binders_NoExprFilt_ord=Deiv_Binders_NoExprFilt[order(Deiv_Binders_NoExprFilt$value),]
#Dest_Binders_NoExprFilt_ord=Dest_Binders_NoExprFilt[order(Dest_Binders_NoExprFilt$value),]
#Dr1_Binders_NoExprFilt_ord=Dr1_Binders_NoExprFilt[order(Dr1_Binders_NoExprFilt$value),]
#Dr4_Binders_NoExprFilt_ord=Dr4_Binders_NoExprFilt[order(Dr4_Binders_NoExprFilt$value),]
#Dr5_Binders_NoExprFilt_ord=Dr5_Binders_NoExprFilt[order(Dr5_Binders_NoExprFilt$value),]
#Foca_Binders_NoExprFilt_ord=Foca_Binders_NoExprFilt[order(Foca_Binders_NoExprFilt$value),]
#Gagra_Binders_NoExprFilt_ord=Gagra_Binders_NoExprFilt[order(Gagra_Binders_NoExprFilt$value),]
#Luan_Binders_NoExprFilt_ord=Luan_Binders_NoExprFilt[order(Luan_Binders_NoExprFilt$value),]
#Mabi_Binders_NoExprFilt_ord=Mabi_Binders_NoExprFilt[order(Mabi_Binders_NoExprFilt$value),]
#Moge_Binders_NoExprFilt_ord=Moge_Binders_NoExprFilt[order(Moge_Binders_NoExprFilt$value),]
#Piag_Binders_NoExprFilt_ord=Piag_Binders_NoExprFilt[order(Piag_Binders_NoExprFilt$value),]
#Prelu_Binders_NoExprFilt_ord=Prelu_Binders_NoExprFilt[order(Prelu_Binders_NoExprFilt$value),]
#




Calu_Binders_NoExprFilt_ord_homo=Calu_Binders_NoExprFilt_ord[Calu_Binders_NoExprFilt_ord$variable=='IC50_HLA-A02:01',]
Calu_Binders_NoExprFilt_ord_homo=rbind(Calu_Binders_NoExprFilt_ord_homo,Calu_Binders_NoExprFilt_ord)
Calu_Binders_NoExprFilt_ord_homo=Calu_Binders_NoExprFilt_ord_homo[order(Calu_Binders_NoExprFilt_ord_homo$value),]


Luan_Binders_NoExprFilt_ord_homo=Luan_Binders_NoExprFilt_ord[Luan_Binders_NoExprFilt_ord$variable=='IC50_HLA-C07:01',]
Luan_Binders_NoExprFilt_ord_homo=rbind(Luan_Binders_NoExprFilt_ord_homo,Luan_Binders_NoExprFilt_ord)
Luan_Binders_NoExprFilt_ord_homo=Luan_Binders_NoExprFilt_ord_homo[order(Luan_Binders_NoExprFilt_ord_homo$value),]

Moge_Binders_NoExprFilt_ord_homo=Moge_Binders_NoExprFilt_ord[Moge_Binders_NoExprFilt_ord$variable=='IC50_HLA-A24:02',]
Moge_Binders_NoExprFilt_ord_homo=rbind(Moge_Binders_NoExprFilt_ord_homo,Moge_Binders_NoExprFilt_ord)
Moge_Binders_NoExprFilt_ord_homo=Moge_Binders_NoExprFilt_ord_homo[order(Moge_Binders_NoExprFilt_ord_homo$value),]


Piag_Binders_NoExprFilt_ord_homo=Piag_Binders_NoExprFilt_ord[Piag_Binders_NoExprFilt_ord$variable=='IC50_HLA-C04:01',]
Piag_Binders_NoExprFilt_ord_homo=rbind(Piag_Binders_NoExprFilt_ord_homo,Piag_Binders_NoExprFilt_ord)
Piag_Binders_NoExprFilt_ord_homo=Piag_Binders_NoExprFilt_ord_homo[order(Piag_Binders_NoExprFilt_ord_homo$value),]




color=c25[c(1:11,13,18,22,24)]
plot(lines(loess.smooth(1:length(Alfe_Binders_NoExprFilt_ord$value), y=Alfe_Binders_NoExprFilt_ord$value)),xlim=c(0,3000),ylim=c(0,500),col=color[1],xlab='Number of Predicted mhAGs Binders',ylab='IC50 (nM)')
lines(loess.smooth(1:length(Alfe_Binders_NoExprFilt_ord$value), y=Alfe_Binders_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[1],xlab='',ylab='IC50 (nM)',lwd=3)
#polygon(1:length(Besu_Binders_ord$value), y=Besu_Binders_ord$value,col='skyblue')
par(new=TRUE)
lines(loess.smooth(1:length(Besu_Binders_NoExprFilt_ord$value), y=Besu_Binders_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[2],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Calu_Binders_NoExprFilt_ord_homo$value), y=Calu_Binders_NoExprFilt_ord_homo$value,degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[3],xlab='',ylab='',lwd=3)
par(new=TRUE)
lines(loess.smooth(1:length(Deiv_Binders_NoExprFilt_ord$value), y=Deiv_Binders_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[4],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Dest_Binders_NoExprFilt_ord$value), y=Dest_Binders_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[5],xlab='',ylab='',lwd=3)
par(new=TRUE)
lines(loess.smooth(1:length(Dr1_Binders_NoExprFilt_ord$value), y=Dr1_Binders_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[6],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Dr4_Binders_NoExprFilt_ord$value), y=Dr4_Binders_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[7],xlab='',ylab='',lwd=3,)
par(new=TRUE)
lines(loess.smooth(1:length(Dr5_Binders_NoExprFilt_ord$value), y=Dr5_Binders_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[8],xlab='',ylab='',lwd=3,)
par(new=TRUE)
lines(loess.smooth(1:length(Foca_Binders_NoExprFilt_ord$value), y=Foca_Binders_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[9],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Gagra_Binders_NoExprFilt_ord$value), y=Gagra_Binders_NoExprFilt_ord$value, degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[10],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Luan_Binders_NoExprFilt_ord_homo$value), y=Luan_Binders_NoExprFilt_ord_homo$value,  degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[11],xlab='',ylab='',lwd=3,)
par(new=TRUE)
lines(loess.smooth(1:length(Mabi_Binders_NoExprFilt_ord$value), y=Mabi_Binders_NoExprFilt_ord$value,  degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[12],xlab='',ylab='',lwd=3)
par(new=TRUE)
lines(loess.smooth(1:length(Moge_Binders_NoExprFilt_ord_homo$value), y=Moge_Binders_NoExprFilt_ord_homo$value,  degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[13],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Piag_Binders_NoExprFilt_ord_homo$value), y=Piag_Binders_NoExprFilt_ord_homo$value,  degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[14],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Prelu_Binders_NoExprFilt_ord$value), y=Prelu_Binders_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,3000),ylim=c(0,500),col=color[15],xlab='',ylab='',lwd=3)
legend("topright", inset=c(-0,0), c('alfe','besu','calu','deiv','dest','dr1','dr4','dr5','foca','gagra','luan','mabi','moge','piag','prelu'),col=color[1:15], title="Patients",cex=0.7,lwd=rep(2,15),lty=c(1,2,1,2,1,2,1,1,2,2,1,1,2,2,1))




alfe_auc  =auc(1:length(Alfe_Binders_NoExprFilt_ord$value), y=Alfe_Binders_NoExprFilt_ord$value)
besu_auc  =auc(1:length(Besu_Binders_NoExprFilt_ord$value), y=Besu_Binders_NoExprFilt_ord$value)
calu_auc  =auc(1:length(Calu_Binders_NoExprFilt_ord_homo$value), y=Calu_Binders_NoExprFilt_ord_homo$value)
deiv_auc  =auc(1:length(Deiv_Binders_NoExprFilt_ord$value), y=Deiv_Binders_NoExprFilt_ord$value)
dest_auc  =auc(1:length(Dest_Binders_NoExprFilt_ord$value), y=Dest_Binders_NoExprFilt_ord$value)
dr1_auc   =auc(1:length(Dr1_Binders_NoExprFilt_ord$value), y=Dr1_Binders_NoExprFilt_ord$value)
dr4_auc   =auc(1:length(Dr4_Binders_NoExprFilt_ord$value), y=Dr4_Binders_NoExprFilt_ord$value)
dr5_auc   =auc(1:length(Dr5_Binders_NoExprFilt_ord$value), y=Dr5_Binders_NoExprFilt_ord$value)
foca_auc  =auc(1:length(Foca_Binders_NoExprFilt_ord$value), y=Foca_Binders_NoExprFilt_ord$value)
gagra_auc =auc(1:length(Gagra_Binders_NoExprFilt_ord$value), y=Gagra_Binders_NoExprFilt_ord$value)
luan_auc  =auc(1:length(Luan_Binders_NoExprFilt_ord_homo$value), y=Luan_Binders_NoExprFilt_ord_homo$value)
mabi_auc  =auc(1:length(Mabi_Binders_NoExprFilt_ord$value), y=Mabi_Binders_NoExprFilt_ord$value)
moge_auc  =auc(1:length(Moge_Binders_NoExprFilt_ord_homo$value), y=Moge_Binders_NoExprFilt_ord_homo$value)
piag_auc  =auc(1:length(Piag_Binders_NoExprFilt_ord_homo$value), y=Piag_Binders_NoExprFilt_ord_homo$value)
prelu_auc =auc(1:length(Prelu_Binders_NoExprFilt_ord$value), y=Prelu_Binders_NoExprFilt_ord$value)

all_aucs=rbind(alfe_auc,besu_auc,calu_auc,deiv_auc,dest_auc,dr1_auc,dr4_auc,dr5_auc,foca_auc,gagra_auc,luan_auc,mabi_auc,moge_auc,piag_auc,prelu_auc)
all_aucs=as.data.frame(all_aucs)
all_aucs$genes=rownames(all_aucs)
all_aucs$RelDon=c(1,0,1,0,1,0,1,1,0,0,1,1,0,0,1)
#all_aucs$RelDon=c(0,1,0,1,0,1,0,0,1,1,0,0,1,1,0)
all_aucs=all_aucs[order(all_aucs$V1),]
colnames(all_aucs)=c('aucs','patients','RelDon')


plot(log(all_aucs$aucs),col=as.factor(all_aucs$Donor),lwd=3)
text(log(all_aucs$aucs),labels=all_aucs$patients, cex= 0.7,adj=1,pos=3,col=color[1:15])

summary(glm(all_aucs$Donors~ all_aucs$aucs,family = binomial(link = "probit")))

Call:
glm(formula = all_aucs$Donor ~ all_aucs$aucs, family = binomial(link = "probit"))

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-1.8418  -0.3673   0.2091   0.7152   1.4387  

Coefficients:
                Estimate Std. Error z value Pr(>|z|)  
(Intercept)    2.443e+00  1.144e+00   2.136   0.0327 *
all_aucs$aucs -1.178e-05  5.779e-06  -2.039   0.0415 *	
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 20.728  on 14  degrees of freedom
Residual deviance: 12.257  on 13  degrees of freedom
AIC: 16.257


#
#
#
#
#Alfe_Binders_Genes_AlleleCoverage_ord=Alfe_Binders_Genes_AlleleCoverage[order(Alfe_Binders_Genes_AlleleCoverage$value),]
#Besu_Binders_NotExpressed_Genes_AlleleCoverage_ord=Besu_Binders_NotExpressed_Genes_AlleleCoverage[order(Besu_Binders_NotExpressed_Genes_AlleleCoverage$value),]
#Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord=Calu_Binders_NotExpressed_Genes_AlleleCoverage[order(Calu_Binders_NotExpressed_Genes_AlleleCoverage$value),]
#Deiv_Binders_NotExpressed_Genes_AlleleCoverage_ord=Deiv_Binders_NotExpressed_Genes_AlleleCoverage[order(Deiv_Binders_NotExpressed_Genes_AlleleCoverage$value),]
#Dest_Binders_NotExpressed_Genes_AlleleCoverage_ord=Dest_Binders_NotExpressed_Genes_AlleleCoverage[order(Dest_Binders_NotExpressed_Genes_AlleleCoverage$value),]
#Dr1_Binders_NotExpressed_Genes_AlleleCoverage_ord=Dr1_Binders_NotExpressed_Genes_AlleleCoverage[order(Dr1_Binders_NotExpressed_Genes_AlleleCoverage$value),]
#Dr4_Binders_NotExpressed_Genes_AlleleCoverage_ord=Dr4_Binders_NotExpressed_Genes_AlleleCoverage[order(Dr4_Binders_NotExpressed_Genes_AlleleCoverage$value),]
#Dr5_Binders_NotExpressed_Genes_AlleleCoverage_ord=Dr5_Binders_NotExpressed_Genes_AlleleCoverage[order(Dr5_Binders_NotExpressed_Genes_AlleleCoverage$value),]
#Foca_Binders_NotExpressed_Genes_AlleleCoverage_ord=Foca_Binders_NotExpressed_Genes_AlleleCoverage[order(Foca_Binders_NotExpressed_Genes_AlleleCoverage$value),]
#Gagra_Binders_NotExpressed_Genes_AlleleCoverage_ord=Gagra_Binders_NotExpressed_Genes_AlleleCoverage[order(Gagra_Binders_NotExpressed_Genes_AlleleCoverage$value),]
#Luan_Binders_NotExpressed_Genes_AlleleCoverage_ord=Luan_Binders_NotExpressed_Genes_AlleleCoverage[order(Luan_Binders_NotExpressed_Genes_AlleleCoverage$value),]
#Mabi_Binders_NotExpressed_Genes_AlleleCoverage_ord=Mabi_Binders_NotExpressed_Genes_AlleleCoverage[order(Mabi_Binders_NotExpressed_Genes_AlleleCoverage$value),]
#Moge_Binders_NotExpressed_Genes_AlleleCoverage_ord=Moge_Binders_NotExpressed_Genes_AlleleCoverage[order(Moge_Binders_NotExpressed_Genes_AlleleCoverage$value),]
#Piag_Binders_NotExpressed_Genes_AlleleCoverage_ord=Piag_Binders_NotExpressed_Genes_AlleleCoverage[order(Piag_Binders_NotExpressed_Genes_AlleleCoverage$value),]
#Prelu_Binders_NotExpressed_Genes_AlleleCoverage_ord=Prelu_Binders_NotExpressed_Genes_AlleleCoverage[order(Prelu_Binders_NotExpressed_Genes_AlleleCoverage$value),]
#
#
#
#Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo=Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord[Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord$variable=='IC50_HLA-A02:01',]
#Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo=rbind(Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo,Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord)
#Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo=Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo[order(Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo$value),]
#
#
#Luan_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo=Luan_Binders_NotExpressed_Genes_AlleleCoverage_ord[Luan_Binders_NotExpressed_Genes_AlleleCoverage_ord$variable=='IC50_HLA-C07:01',]
#Luan_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo=rbind(Luan_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo,Luan_Binders_NotExpressed_Genes_AlleleCoverage_ord)
#Luan_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo=Luan_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo[order(Luan_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo$value),]
#
#Moge_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo=Moge_Binders_NotExpressed_Genes_AlleleCoverage_ord[Moge_Binders_NotExpressed_Genes_AlleleCoverage_ord$variable=='IC50_HLA-A24:02',]
#Moge_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo=rbind(Moge_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo,Moge_Binders_NotExpressed_Genes_AlleleCoverage_ord)
#Moge_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo=Moge_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo[order(Moge_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo$value),]
#
#
#Piag_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo=Piag_Binders_NotExpressed_Genes_AlleleCoverage_ord[Piag_Binders_NotExpressed_Genes_AlleleCoverage_ord$variable=='IC50_HLA-C04:01',]
#Piag_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo=rbind(Piag_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo,Piag_Binders_NotExpressed_Genes_AlleleCoverage_ord)
#Piag_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo=Piag_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo[order(Piag_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo$value),]
#
#
#pdf('Ic50_vs_NumberPeptides_MhAGs_NoExprFiltr.pdf')
#plot.new()
#color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
#par(mar=c(5.1,5.1,5.1,5.1))
##ets the bottom, left, top and right margins respectively of the plot region in number of lines of text.
#
#c25 <- c("dodgerblue2","#E31A1C", # red
#                "green4",
#                "#6A3D9A", # purple
#                "#FF7F00", # orange
#                "black","gold1",
#                "skyblue2","#FB9A99", # lt pink
#                "palegreen2",
#                "#CAB2D6", # lt purple
#                "#FDBF6F", # lt orange
#                "gray70", "khaki2",
#                "maroon","orchid1","deeppink1","blue1","steelblue4",
#                "darkturquoise","green1","yellow4","yellow3",
#                "darkorange4","brown")
#
#
#color=c25[c(1:11,13,18,22,24)]
#plot(lines(loess.smooth(1:length(Alfe_Binders_Genes_AlleleCoverage_ord$value), y=Alfe_Binders_Genes_AlleleCoverage_ord$value)),xlim=c(0,2500),ylim=c(0,500),col=color[1],xlab='Number of Predicted mhAGs Binders_NoExprFilt',ylab='IC50 (nM)')
#lines(loess.smooth(1:length(Alfe_Binders_Genes_AlleleCoverage_ord$value), y=Alfe_Binders_Genes_AlleleCoverage_ord$value,degree=1,span=0.8),xlim=c(0,2500),ylim=c(0,500),col=color[1],xlab='',ylab='IC50 (nM)',lwd=3)
##polygon(1:length(Besu_Binders_NotExpressed_Genes_AlleleCoverage_ord$value), y=Besu_Binders_NotExpressed_Genes_AlleleCoverage_ord$value,col='skyblue')
#par(new=TRUE)
#lines(loess.smooth(1:length(Besu_Binders_NotExpressed_Genes_AlleleCoverage_ord$value), y=Besu_Binders_NotExpressed_Genes_AlleleCoverage_ord$value,degree=1,span=0.8),xlim=c(0,2500),ylim=c(0,500),col=color[2],xlab='',ylab='',lwd=3,lty=2)
#par(new=TRUE)
#lines(loess.smooth(1:length(Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo$value), y=Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo$value,degree=1,span=0.8),xlim=c(0,2500),ylim=c(0,500),col=color[3],xlab='',ylab='',lwd=3)
#par(new=TRUE)
#lines(loess.smooth(1:length(Deiv_Binders_NotExpressed_Genes_AlleleCoverage_ord$value), y=Deiv_Binders_NotExpressed_Genes_AlleleCoverage_ord$value,degree=1,span=0.8),xlim=c(0,2500),ylim=c(0,500),col=color[4],xlab='',ylab='',lwd=3,lty=2)
#par(new=TRUE)
#lines(loess.smooth(1:length(Dest_Binders_NotExpressed_Genes_AlleleCoverage_ord$value), y=Dest_Binders_NotExpressed_Genes_AlleleCoverage_ord$value,degree=1,span=0.8),xlim=c(0,2500),ylim=c(0,500),col=color[5],xlab='',ylab='',lwd=3)
#par(new=TRUE)
#lines(loess.smooth(1:length(Dr1_Binders_NotExpressed_Genes_AlleleCoverage_ord$value), y=Dr1_Binders_NotExpressed_Genes_AlleleCoverage_ord$value,degree=1,span=0.8),xlim=c(0,2500),ylim=c(0,500),col=color[6],xlab='',ylab='',lwd=3,lty=2)
#par(new=TRUE)
#lines(loess.smooth(1:length(Dr4_Binders_NotExpressed_Genes_AlleleCoverage_ord$value), y=Dr4_Binders_NotExpressed_Genes_AlleleCoverage_ord$value,degree=1,span=0.8),xlim=c(0,2500),ylim=c(0,500),col=color[7],xlab='',ylab='',lwd=3,)
#par(new=TRUE)
#lines(loess.smooth(1:length(Dr5_Binders_NotExpressed_Genes_AlleleCoverage_ord$value), y=Dr5_Binders_NotExpressed_Genes_AlleleCoverage_ord$value,degree=1,span=0.8),xlim=c(0,2500),ylim=c(0,500),col=color[8],xlab='',ylab='',lwd=3,)
#par(new=TRUE)
#lines(loess.smooth(1:length(Foca_Binders_NotExpressed_Genes_AlleleCoverage_ord$value), y=Foca_Binders_NotExpressed_Genes_AlleleCoverage_ord$value,degree=1,span=0.8),xlim=c(0,2500),ylim=c(0,500),col=color[9],xlab='',ylab='',lwd=3,lty=2)
#par(new=TRUE)
#lines(loess.smooth(1:length(Gagra_Binders_NotExpressed_Genes_AlleleCoverage_ord$value), y=Gagra_Binders_NotExpressed_Genes_AlleleCoverage_ord$value, degree=1,span=0.8),xlim=c(0,2500),ylim=c(0,500),col=color[10],xlab='',ylab='',lwd=3,lty=2)
#par(new=TRUE)
#lines(loess.smooth(1:length(Luan_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo$value), y=Luan_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo$value,  degree=1,span=0.8),xlim=c(0,2500),ylim=c(0,500),col=color[11],xlab='',ylab='',lwd=3,)
#par(new=TRUE)
#lines(loess.smooth(1:length(Mabi_Binders_NotExpressed_Genes_AlleleCoverage_ord$value), y=Mabi_Binders_NotExpressed_Genes_AlleleCoverage_ord$value,  degree=1,span=0.8),xlim=c(0,2500),ylim=c(0,500),col=color[12],xlab='',ylab='',lwd=3)
#par(new=TRUE)
#lines(loess.smooth(1:length(Moge_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo$value), y=Moge_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo$value,  degree=1,span=0.8),xlim=c(0,2500),ylim=c(0,500),col=color[13],xlab='',ylab='',lwd=3,lty=2)
#par(new=TRUE)
#lines(loess.smooth(1:length(Piag_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo$value), y=Piag_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo$value,  degree=1,span=0.8),xlim=c(0,2500),ylim=c(0,500),col=color[14],xlab='',ylab='',lwd=3,lty=2)
#par(new=TRUE)
#lines(loess.smooth(1:length(Prelu_Binders_NotExpressed_Genes_AlleleCoverage_ord$value), y=Prelu_Binders_NotExpressed_Genes_AlleleCoverage_ord$value,degree=1,span=0.8),xlim=c(0,2500),ylim=c(0,500),col=color[15],xlab='',ylab='',lwd=3)
#legend("topright", inset=c(-0,0), c('alfe','besu','calu','deiv','dest','dr1','dr4','dr5','foca','gagra','luan','mabi','moge','piag','prelu'),col=color[1:15], title="Patients",cex=0.7,lwd=rep(2,15),lty=c(1,2,1,2,1,2,1,1,2,2,1,1,2,2,1))
#
#dev.off()
#
#
#
#
#
#Alfe_Binders_Genes_AlleleCoverage_ord,Besu_Binders_NotExpressed_Genes_AlleleCoverage_ord,Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo,Deiv_Binders_NotExpressed_Genes_AlleleCoverage_ord,Dest_Binders_NotExpressed_Genes_AlleleCoverage_ord,Dr1_Binders_NotExpressed_Genes_AlleleCoverage_ord,Dr4_Binders_NotExpressed_Genes_AlleleCoverage_ord,)
#known_Binders=c('HMHA1','MYO1G','AKAP13','KIAA0020','HMHB1','BCL2A1','CENPM','SP110','UGT2B17')
#
#
#
#
#
#
#
Alfe_Binders_Genes_AlleleCoverage_ord_Dx=Alfe_Binders_Genes_AlleleCoverage_ord[Alfe_Binders_Genes_AlleleCoverage_ord$Dx_Mut>=3,]
Besu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx=Besu_Binders_NotExpressed_Genes_AlleleCoverage_ord[Besu_Binders_NotExpressed_Genes_AlleleCoverage_ord$Dx_Mut>=3,]
Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx=Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo[Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo$Dx_Mut>=3,]
Deiv_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx=Deiv_Binders_NotExpressed_Genes_AlleleCoverage_ord[Deiv_Binders_NotExpressed_Genes_AlleleCoverage_ord$Dx_Mut>=3,]
Dest_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx=Dest_Binders_NotExpressed_Genes_AlleleCoverage_ord[Dest_Binders_NotExpressed_Genes_AlleleCoverage_ord$Dx_Mut>=3,]
Dr1_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx=Dr1_Binders_NotExpressed_Genes_AlleleCoverage_ord[Dr1_Binders_NotExpressed_Genes_AlleleCoverage_ord$Dx_Mut>=3,]
Dr4_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx=Dr4_Binders_NotExpressed_Genes_AlleleCoverage_ord[Dr4_Binders_NotExpressed_Genes_AlleleCoverage_ord$Dx_Mut>=3,]
Dr5_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx=Dr5_Binders_NotExpressed_Genes_AlleleCoverage_ord[Dr5_Binders_NotExpressed_Genes_AlleleCoverage_ord$Dx_Mut>=3,]
Foca_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx=Foca_Binders_NotExpressed_Genes_AlleleCoverage_ord[Foca_Binders_NotExpressed_Genes_AlleleCoverage_ord$Dx_Mut>=3,]
Gagra_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx=Gagra_Binders_NotExpressed_Genes_AlleleCoverage_ord[Gagra_Binders_NotExpressed_Genes_AlleleCoverage_ord$Dx_Mut>=3,]
Luan_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx=Luan_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo[Luan_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo$Dx_Mut>=3,]
Mabi_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx=Mabi_Binders_NotExpressed_Genes_AlleleCoverage_ord[Mabi_Binders_NotExpressed_Genes_AlleleCoverage_ord$Dx_Mut>=3,]
Moge_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx=Moge_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo[Moge_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo$Dx_Mut>=3,]
Piag_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx=Piag_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo[Piag_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo$Dx_Mut>=3,]
Prelu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx=Prelu_Binders_NotExpressed_Genes_AlleleCoverage_ord[Prelu_Binders_NotExpressed_Genes_AlleleCoverage_ord$Dx_Mut>=3,]



alfe_auc_Dx  =auc(1:length(Alfe_Binders_Genes_AlleleCoverage_ord_Dx$value), y=Alfe_Binders_Genes_AlleleCoverage_ord_Dx$value)
besu_auc_Dx  =auc(1:length(Besu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Besu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value)
calu_auc_Dx  =auc(1:length(Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value)
deiv_auc_Dx  =auc(1:length(Deiv_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Deiv_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value)
dest_auc_Dx  =auc(1:length(Dest_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Dest_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value)
dr1_auc_Dx   =auc(1:length(Dr1_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Dr1_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value)
dr4_auc_Dx   =auc(1:length(Dr4_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Dr4_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value)
dr5_auc_Dx   =auc(1:length(Dr5_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Dr5_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value)
foca_auc_Dx  =auc(1:length(Foca_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Foca_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value)
gagra_auc_Dx =auc(1:length(Gagra_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Gagra_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value)
luan_auc_Dx  =auc(1:length(Luan_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Luan_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value)
mabi_auc_Dx  =auc(1:length(Mabi_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Mabi_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value)
moge_auc_Dx  =auc(1:length(Moge_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Moge_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value)
piag_auc_Dx  =auc(1:length(Piag_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Piag_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value)
prelu_auc_Dx =auc(1:length(Prelu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Prelu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value)

all_aucs_Dx=rbind(alfe_auc_Dx,besu_auc_Dx,calu_auc_Dx,deiv_auc_Dx,dest_auc_Dx,dr1_auc_Dx,dr4_auc_Dx,dr5_auc_Dx,foca_auc_Dx,gagra_auc_Dx,luan_auc_Dx,mabi_auc_Dx,moge_auc_Dx,piag_auc_Dx,prelu_auc_Dx)
all_aucs_Dx=as.data.frame(all_aucs_Dx)
all_aucs_Dx$genes=rownames(all_aucs_Dx)
all_aucs_Dx$RelDon=c(1,0,1,0,1,0,1,1,0,0,1,1,0,0,1)

#all_aucs_Dx$RelDon=c(0,1,0,1,0,1,0,0,1,1,0,0,1,1,0)
all_aucs_Dx=all_aucs_Dx[order(all_aucs_Dx$V1),]
colnames(all_aucs_Dx)=c('aucs','patients','RelDon')
all_aucs_Dx$Donors=ifelse(all_aucs_Dx$RelDon==1,'Related','Unrelated')

plot(log(all_aucs_Dx$aucs),col=as.factor(all_aucs_Dx$Donor),lwd=3)
text(log(all_aucs_Dx$aucs),labels=all_aucs_Dx$patients, cex= 0.7,adj=1,pos=3,col=color[1:15])

summary(glm(all_aucs_Dx$Donors~ all_aucs_Dx$aucs,family = binomial(link = "probit")))

Call:
glm(formula = factor(all_aucs_Dx$Donors) ~ all_aucs_Dx$aucs, family = binomial(link = "probit"))

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-1.7572  -0.1077   0.1533   0.5689   1.3880  

Coefficients:
                   Estimate Std. Error z value Pr(>|z|)  
(Intercept)       2.876e+00  1.409e+00   2.041   0.0413 *
all_aucs_Dx$aucs -5.611e-05  2.970e-05  -1.890   0.0588 .
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 20.728  on 14  degrees of freedom
Residual deviance: 10.036  on 13  degrees of freedom
AIC: 14.036

Number of Fisher Scoring iterations: 7




ggplot(all_aucs_Dx,aes(x=seq(1,length(all_aucs_Dx$aucs)),y=aucs,col=Donors)) +
geom_point()  +
geom_rug() +
scale_color_manual(values=c( Related= "#E69F00", Unrelated="darkorchid2")) +
geom_smooth( method = "lm",mapping = aes(fill=all_aucs_Dx$Donors)) +
scale_fill_manual(values=c(Related= "#E69F00", Unrelated="darkorchid2")) +
ggtitle('MRD vs MUD AUCs Binders\n linear regression') +
xlab('Patients')




alfe_auc_Rel  =auc(1:length(Alfe_Binders_Genes_AlleleCoverage_ord_Rel$value), y=Alfe_Binders_Genes_AlleleCoverage_ord_Rel$value)
besu_auc_Rel  =auc(1:length(Besu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Besu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value)
calu_auc_Rel  =auc(1:length(Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value)
deiv_auc_Rel  =auc(1:length(Deiv_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Deiv_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value)
dest_auc_Rel  =auc(1:length(Dest_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Dest_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value)
dr1_auc_Rel   =auc(1:length(Dr1_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Dr1_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value)
dr4_auc_Rel   =auc(1:length(Dr4_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Dr4_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value)
dr5_auc_Rel   =auc(1:length(Dr5_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Dr5_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value)
foca_auc_Rel  =auc(1:length(Foca_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Foca_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value)
gagra_auc_Rel =auc(1:length(Gagra_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Gagra_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value)
luan_auc_Rel  =auc(1:length(Luan_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Luan_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value)
mabi_auc_Rel  =auc(1:length(Mabi_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Mabi_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value)
moge_auc_Rel  =auc(1:length(Moge_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Moge_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value)
piag_auc_Rel  =auc(1:length(Piag_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Piag_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value)
prelu_auc_Rel =auc(1:length(Prelu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Prelu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value)

all_aucs_Rel=rbind(alfe_auc_Rel,besu_auc_Rel,calu_auc_Rel,deiv_auc_Rel,dest_auc_Rel,dr1_auc_Rel,dr4_auc_Rel,dr5_auc_Rel,foca_auc_Rel,gagra_auc_Rel,luan_auc_Rel,mabi_auc_Rel,moge_auc_Rel,piag_auc_Rel,prelu_auc_Rel)
all_aucs_Rel=as.data.frame(all_aucs_Rel)
all_aucs_Rel$genes=rownames(all_aucs_Rel)
all_aucs_Rel$RelDon=c(1,0,1,0,1,0,1,1,0,0,1,1,0,0,1)

#all_aucs_Rel$RelDon=c(0,1,0,1,0,1,0,0,1,1,0,0,1,1,0)
all_aucs_Rel=all_aucs_Rel[order(all_aucs_Rel$V1),]
colnames(all_aucs_Rel)=c('aucs','patients','RelDon')
all_aucs_Rel$Donors=ifelse(all_aucs_Rel$RelDon==1,'Related','Unrelated')

plot(log(all_aucs_Rel$aucs),col=as.factor(all_aucs_Rel$Donor),lwd=3)
text(log(all_aucs_Rel$aucs),labels=all_aucs_Rel$patients, cex= 0.7,adj=1,pos=3,col=color[1:15])

summary(glm(factor(all_aucs_Rel$Donors)~ all_aucs_Rel$aucs,family = binomial(link = "probit")))

Call:
glm(formula = factor(all_aucs_Rel$Donors) ~ all_aucs_Rel$aucs, 
    family = binomial(link = "probit"))

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-1.6819  -0.7757  -0.3412   0.4579   1.7968  

Coefficients:
                    Estimate Std. Error z value Pr(>|z|)  
(Intercept)       -1.995e+00  9.005e-01  -2.216   0.0267 *
all_aucs_Rel$aucs  3.335e-05  1.459e-05   2.285   0.0223 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 20.728  on 14  degrees of freedom
Residual deviance: 13.978  on 13  degrees of freedom
AIC: 17.978

Number of Fisher Scoring iterations: 5


ggplot(all_aucs_Rel,aes(x=seq(1,length(all_aucs_Rel$aucs)),y=aucs,col=Donors)) +
geom_point()  +
geom_rug() +
scale_color_manual(values=c( Related= "#E69F00", Unrelated="darkorchid2")) +
geom_smooth( method = "lm",mapping = aes(fill=all_aucs_Rel$Donors)) +
scale_fill_manual(values=c(Related= "#E69F00", Unrelated="darkorchid2")) +
ggtitle('MRD vs MUD AUCs Binders\n linear regression') +
xlab('Patients')


#
#
#pdf('Ic50_vs_NumberPeptides_MhAGs_DiagnosisExprFiltr.pdf')
plot.new()
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
par(mar=c(5.1,5.1,5.1,5.1))
#ets the bottom, left, top and right margins respectively of the plot region in number of lines of text.

c25 <- c("dodgerblue2","#E31A1C", # red
                "green4",
                "#6A3D9A", # purple
                "#FF7F00", # orange
                "black","gold1",
                "skyblue2","#FB9A99", # lt pink
                "palegreen2",
                "#CAB2D6", # lt purple
                "#FDBF6F", # lt orange
                "gray70", "khaki2",
                "maroon","orchid1","deeppink1","blue1","steelblue4",
                "darkturquoise","green1","yellow4","yellow3",
                "darkorange4","brown")


alfe_auc  =auc(1:length(Alfe_Binders_NoExprFilt_ord$value), y=Alfe_Binders_NoExprFilt_ord$value)
besu_auc  =auc(1:length(Besu_Binders_NoExprFilt_ord$value), y=Besu_Binders_NoExprFilt_ord$value)
calu_auc  =auc(1:length(Calu_Binders_NoExprFilt_ord_homo$value), y=Calu_Binders_NoExprFilt_ord_homo$value)
deiv_auc  =auc(1:length(Deiv_Binders_NoExprFilt_ord$value), y=Deiv_Binders_NoExprFilt_ord$value)
dest_auc  =auc(1:length(Dest_Binders_NoExprFilt_ord$value), y=Dest_Binders_NoExprFilt_ord$value)
dr1_auc   =auc(1:length(Dr1_Binders_NoExprFilt_ord$value), y=Dr1_Binders_NoExprFilt_ord$value)
dr4_auc   =auc(1:length(Dr4_Binders_NoExprFilt_ord$value), y=Dr4_Binders_NoExprFilt_ord$value)
dr5_auc   =auc(1:length(Dr5_Binders_NoExprFilt_ord$value), y=Dr5_Binders_NoExprFilt_ord$value)
foca_auc  =auc(1:length(Foca_Binders_NoExprFilt_ord$value), y=Foca_Binders_NoExprFilt_ord$value)
gagra_auc =auc(1:length(Gagra_Binders_NoExprFilt_ord$value), y=Gagra_Binders_NoExprFilt_ord$value)
luan_auc  =auc(1:length(Luan_Binders_NoExprFilt_ord_homo$value), y=Luan_Binders_NoExprFilt_ord_homo$value)
mabi_auc  =auc(1:length(Mabi_Binders_NoExprFilt_ord$value), y=Mabi_Binders_NoExprFilt_ord$value)
moge_auc  =auc(1:length(Moge_Binders_NoExprFilt_ord_homo$value), y=Moge_Binders_NoExprFilt_ord_homo$value)
piag_auc  =auc(1:length(Piag_Binders_NoExprFilt_ord_homo$value), y=Piag_Binders_NoExprFilt_ord_homo$value)
prelu_auc =auc(1:length(Prelu_Binders_NoExprFilt_ord$value), y=Prelu_Binders_NoExprFilt_ord$value)


color=c25[c(1:11,13,18,22,24)]
plot(lines(loess.smooth(1:length(Alfe_Binders_Genes_AlleleCoverage_ord_Dx$value), y=Alfe_Binders_Genes_AlleleCoverage_ord_Dx$value)),xlim=c(0,1000),ylim=c(0,500),col=color[1],xlab='Number of Predicted mhAGs Binders Expressed at Diagnosis',ylab='IC50 (nM)')
lines(loess.smooth(1:length(Alfe_Binders_Genes_AlleleCoverage_ord_Dx$value), y=Alfe_Binders_Genes_AlleleCoverage_ord_Dx$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[1],xlab='',ylab='IC50 (nM)',lwd=3)
#polygon(1:length(Besu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Besu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value,col='skyblue')
par(new=TRUE)
lines(loess.smooth(1:length(Besu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Besu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[2],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[3],xlab='',ylab='',lwd=3)
par(new=TRUE)
lines(loess.smooth(1:length(Deiv_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Deiv_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[4],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Dest_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Dest_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[5],xlab='',ylab='',lwd=3)
par(new=TRUE)
lines(loess.smooth(1:length(Dr1_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Dr1_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[6],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Dr4_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Dr4_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[7],xlab='',ylab='',lwd=3,)
par(new=TRUE)
lines(loess.smooth(1:length(Dr5_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Dr5_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[8],xlab='',ylab='',lwd=3,)
par(new=TRUE)
lines(loess.smooth(1:length(Foca_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Foca_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[9],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Gagra_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Gagra_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value, degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[10],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Luan_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Luan_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value,  degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[11],xlab='',ylab='',lwd=3,)
par(new=TRUE)
lines(loess.smooth(1:length(Mabi_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Mabi_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value,  degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[12],xlab='',ylab='',lwd=3)
par(new=TRUE)
lines(loess.smooth(1:length(Moge_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Moge_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value,  degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[13],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Piag_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Piag_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value,  degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[14],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Prelu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value), y=Prelu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Dx$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[15],xlab='',ylab='',lwd=3)
legend("topright", inset=c(-0,0), c('alfe','besu','calu','deiv','dest','dr1','dr4','dr5','foca','gagra','luan','mabi','moge','piag','prelu'),col=color[1:15], title="Patients",cex=0.7,lwd=rep(2,15),lty=c(1,2,1,2,1,2,1,1,2,2,1,1,2,2,1))
#
#

Alfe_Binders_Genes_AlleleCoverage_ord_Rel=Alfe_Binders_Genes_AlleleCoverage_ord[Alfe_Binders_Genes_AlleleCoverage_ord$Rel_Mut>=3,]
Besu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel=Besu_Binders_NotExpressed_Genes_AlleleCoverage_ord[Besu_Binders_NotExpressed_Genes_AlleleCoverage_ord$Rel_Mut>=3,]
Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel=Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo[Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo$Rel_Mut>=3,]
Deiv_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel=Deiv_Binders_NotExpressed_Genes_AlleleCoverage_ord[Deiv_Binders_NotExpressed_Genes_AlleleCoverage_ord$Rel_Mut>=3,]
Dest_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel=Dest_Binders_NotExpressed_Genes_AlleleCoverage_ord[Dest_Binders_NotExpressed_Genes_AlleleCoverage_ord$Rel_Mut>=3,]
Dr1_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel=Dr1_Binders_NotExpressed_Genes_AlleleCoverage_ord[Dr1_Binders_NotExpressed_Genes_AlleleCoverage_ord$Rel_Mut>=3,]
Dr4_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel=Dr4_Binders_NotExpressed_Genes_AlleleCoverage_ord[Dr4_Binders_NotExpressed_Genes_AlleleCoverage_ord$Rel_Mut>=3,]
Dr5_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel=Dr5_Binders_NotExpressed_Genes_AlleleCoverage_ord[Dr5_Binders_NotExpressed_Genes_AlleleCoverage_ord$Rel_Mut>=3,]
Foca_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel=Foca_Binders_NotExpressed_Genes_AlleleCoverage_ord[Foca_Binders_NotExpressed_Genes_AlleleCoverage_ord$Rel_Mut>=3,]
Gagra_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel=Gagra_Binders_NotExpressed_Genes_AlleleCoverage_ord[Gagra_Binders_NotExpressed_Genes_AlleleCoverage_ord$Rel_Mut>=3,]
Luan_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel=Luan_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo[Luan_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo$Rel_Mut>=3,]
Mabi_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel=Mabi_Binders_NotExpressed_Genes_AlleleCoverage_ord[Mabi_Binders_NotExpressed_Genes_AlleleCoverage_ord$Rel_Mut>=3,]
Moge_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel=Moge_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo[Moge_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo$Rel_Mut>=3,]
Piag_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel=Piag_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo[Piag_Binders_NotExpressed_Genes_AlleleCoverage_ord_homo$Rel_Mut>=3,]
Prelu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel=Prelu_Binders_NotExpressed_Genes_AlleleCoverage_ord[Prelu_Binders_NotExpressed_Genes_AlleleCoverage_ord$Rel_Mut>=3,]


#pdf('Ic50_vs_NumberPeptides_MhAGs_RelapseExprFiltr.pdf')
plot.new()
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
par(mar=c(5.1,5.1,5.1,5.1))
#ets the bottom, left, top and right margins respectively of the plot region in number of lines of text.

c25 <- c("dodgerblue2","#E31A1C", # red
                "green4",
                "#6A3D9A", # purple
                "#FF7F00", # orange
                "black","gold1",
                "skyblue2","#FB9A99", # lt pink
                "palegreen2",
                "#CAB2D6", # lt purple
                "#FDBF6F", # lt orange
                "gray70", "khaki2",
                "maroon","orchid1","deeppink1","blue1","steelblue4",
                "darkturquoise","green1","yellow4","yellow3",
                "darkorange4","brown")


color=c25[c(1:11,13,18,22,24)]
plot(lines(loess.smooth(1:length(Alfe_Binders_Genes_AlleleCoverage_ord_Rel$value), y=Alfe_Binders_Genes_AlleleCoverage_ord_Rel$value)),xlim=c(0,1000),ylim=c(0,500),col=color[1],xlab='Number of Predicted mhAGs Binders Expressed at Relapse',ylab='IC50 (nM)')
lines(loess.smooth(1:length(Alfe_Binders_Genes_AlleleCoverage_ord_Rel$value), y=Alfe_Binders_Genes_AlleleCoverage_ord_Rel$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[1],xlab='',ylab='IC50 (nM)',lwd=3)
#polygon(1:length(Besu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Besu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value,col='skyblue')
par(new=TRUE)
lines(loess.smooth(1:length(Besu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Besu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[2],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[3],xlab='',ylab='',lwd=3)
par(new=TRUE)
lines(loess.smooth(1:length(Deiv_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Deiv_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[4],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Dest_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Dest_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[5],xlab='',ylab='',lwd=3)
par(new=TRUE)
lines(loess.smooth(1:length(Dr1_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Dr1_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[6],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Dr4_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Dr4_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[7],xlab='',ylab='',lwd=3,)
par(new=TRUE)
lines(loess.smooth(1:length(Dr5_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Dr5_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[8],xlab='',ylab='',lwd=3,)
par(new=TRUE)
lines(loess.smooth(1:length(Foca_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Foca_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[9],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Gagra_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Gagra_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value, degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[10],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Luan_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Luan_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value,  degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[11],xlab='',ylab='',lwd=3,)
par(new=TRUE)
lines(loess.smooth(1:length(Mabi_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Mabi_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value,  degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[12],xlab='',ylab='',lwd=3)
par(new=TRUE)
lines(loess.smooth(1:length(Moge_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Moge_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value,  degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[13],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Piag_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Piag_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value,  degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[14],xlab='',ylab='',lwd=3,lty=2)
par(new=TRUE)
lines(loess.smooth(1:length(Prelu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value), y=Prelu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[15],xlab='',ylab='',lwd=3)
legend("topright", inset=c(-0,0), c('alfe','besu','calu','deiv','dest','dr1','dr4','dr5','foca','gagra','luan','mabi','moge','piag','prelu'),col=color[1:15], title="Patients",cex=0.7,lwd=rep(2,15),lty=c(1,2,1,2,1,2,1,1,2,2,1,1,2,2,1))
#
#dev.off()
#
#
#
#all_binders=rbind(Alfe_Binders_Genes_AlleleCoverage_ord_Rel,Besu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel,Calu_Binders_NotExpressed_Genes_AlleleCoverage_ord_Rel,Deiv_Binders_Genes_ord,Dr1_Binders_Genes_ord,Dr1_Binders_Genes_ord,Dr4_Binders_Genes_ord,Dr5_Binders_Genes_ord,Foca_Binders_Genes_ord,Gagra_Binders_Genes_ord,Luan_Binders_Genes_ord,Mabi_Binders_Genes_ord,Moge_Binders_Genes_ord,Piag_Binders_Genes_ord,Moge_Binders_Genes_ord)
#known_Binders=c('HMHA1','MYO1G','AKAP13','KIAA0020','HMHB1','BCL2A1','CENPM','SP110','UGT2B17')
#
#
#
#
#
#
#
#
#
#Besu_Binders_NoExprFilt_ord=Besu_Binders_NoExprFilt[order(Besu_Binders_NoExprFilt$value),]
#Besu_Binders_NoExprFilt_ord=Besu_Binders_NoExprFilt[order(Besu_Binders_NoExprFilt$value),]
#Calu_Binders_NoExprFilt_ord=Calu_Binders_NoExprFilt[order(Calu_Binders_NoExprFilt$value),]
#Deiv_Binders_NoExprFilt_ord=Deiv_Binders_NoExprFilt[order(Deiv_Binders_NoExprFilt$value),]
#Dest_Binders_NoExprFilt_ord=Dest_Binders_NoExprFilt[order(Dest_Binders_NoExprFilt$value),]
#Dr1_Binders_NoExprFilt_ord=Dr1_Binders_NoExprFilt[order(Dr1_Binders_NoExprFilt$value),]
#Dr4_Binders_NoExprFilt_ord=Dr4_Binders_NoExprFilt[order(Dr4_Binders_NoExprFilt$value),]
#Dr5_Binders_NoExprFilt_ord=Dr5_Binders_NoExprFilt[order(Dr5_Binders_NoExprFilt$value),]
#Foca_Binders_NoExprFilt_ord=Foca_Binders_NoExprFilt[order(Foca_Binders_NoExprFilt$value),]
#Gagra_Binders_NoExprFilt_ord=Gagra_Binders_NoExprFilt[order(Gagra_Binders_NoExprFilt$value),]
#Luan_Binders_NoExprFilt_ord=Luan_Binders_NoExprFilt[order(Luan_Binders_NoExprFilt$value),]
#Mabi_Binders_NoExprFilt_ord=Mabi_Binders_NoExprFilt[order(Mabi_Binders_NoExprFilt$value),]
#Moge_Binders_NoExprFilt_ord=Moge_Binders_NoExprFilt[order(Moge_Binders_NoExprFilt$value),]
#Piag_Binders_NoExprFilt_ord=Piag_Binders_NoExprFilt[order(Piag_Binders_NoExprFilt$value),]
#Moge_Binders_NoExprFilt_ord=Moge_Binders_NoExprFilt[order(Moge_Binders_NoExprFilt$value),]
#
#
#
#Besu_Binders_ord=Besu_Binders[order(Besu_Binders$value),]
#Besu_Binders_ord=Besu_Binders[order(Besu_Binders$value),]
#Calu_Binders_ord=Calu_Binders[order(Calu_Binders$value),]
#Deiv_Binders_ord=Deiv_Binders[order(Deiv_Binders$value),]
#Dest_Binders_ord=Dest_Binders[order(Dest_Binders$value),]
#Dr1_Binders_ord=Dr1_Binders[order(Dr1_Binders$value),]
#Dr4_Binders_ord=Dr4_Binders[order(Dr4_Binders$value),]
#Dr5_Binders_ord=Dr5_Binders[order(Dr5_Binders$value),]
#Foca_Binders_ord=Foca_Binders[order(Foca_Binders$value),]
#Gagra_Binders_ord=Gagra_Binders[order(Gagra_Binders$value),]
#Luan_Binders_ord=Luan_Binders[order(Luan_Binders$value),]
#Mabi_Binders_ord=Mabi_Binders[order(Mabi_Binders$value),]
#Moge_Binders_ord=Moge_Binders[order(Moge_Binders$value),]
#Piag_Binders_ord=Piag_Binders[order(Piag_Binders$value),]
#Moge_Binders_ord=Moge_Binders[order(Moge_Binders$value),]
#
#
#
#
#Besu_Binders_Genes_ord=Besu_Binders_Genes[order(Besu_Binders_Genes$value),]
#Besu_Binders_Genes_ord=Besu_Binders_Genes[order(Besu_Binders_Genes$value),]
#Calu_Binders_Genes_ord=Calu_Binders_Genes[order(Calu_Binders_Genes$value),]
#Deiv_Binders_Genes_ord=Deiv_Binders_Genes[order(Deiv_Binders_Genes$value),]
#Dest_Binders_Genes_ord=Dest_Binders_Genes[order(Dest_Binders_Genes$value),]
#Dr1_Binders_Genes_ord=Dr1_Binders_Genes[order(Dr1_Binders_Genes$value),]
#Dr4_Binders_Genes_ord=Dr4_Binders_Genes[order(Dr4_Binders_Genes$value),]
#Dr5_Binders_Genes_ord=Dr5_Binders_Genes[order(Dr5_Binders_Genes$value),]
#Foca_Binders_Genes_ord=Foca_Binders_Genes[order(Foca_Binders_Genes$value),]
#Gagra_Binders_Genes_ord=Gagra_Binders_Genes[order(Gagra_Binders_Genes$value),]
#Luan_Binders_Genes_ord=Luan_Binders_Genes[order(Luan_Binders_Genes$value),]
#Mabi_Binders_Genes_ord=Mabi_Binders_Genes[order(Mabi_Binders_Genes$value),]
#Moge_Binders_Genes_ord=Moge_Binders_Genes[order(Moge_Binders_Genes$value),]
#Piag_Binders_Genes_ord=Piag_Binders_Genes[order(Piag_Binders_Genes$value),]
#Moge_Binders_Genes_ord=Moge_Binders_Genes[order(Moge_Binders_Genes$value),]
#
#
#
#
#Calu_Binders_NoExprFilt_ord_homo=Calu_Binders_NoExprFilt_ord[Calu_Binders_NoExprFilt_ord$variable=='IC50_HLA-A02:01',]
#Calu_Binders_NoExprFilt_ord_homo=rbind(Calu_Binders_NoExprFilt_ord_homo,Calu_Binders_NoExprFilt_ord)
#Calu_Binders_NoExprFilt_ord_homo=Calu_Binders_NoExprFilt_ord_homo[order(Calu_Binders_NoExprFilt_ord_homo$value),]
#
#Luan_Binders_NoExprFilt_ord_homo=Luan_Binders_NoExprFilt_ord[Luan_Binders_NoExprFilt_ord$variable=='IC50_HLA-C07:01',]
#Luan_Binders_NoExprFilt_ord_homo=rbind(Luan_Binders_NoExprFilt_ord_homo,Luan_Binders_NoExprFilt_ord)
#Luan_Binders_NoExprFilt_ord_homo=Luan_Binders_NoExprFilt_ord_homo[order(Luan_Binders_NoExprFilt_ord_homo$value),]
#
#Moge_Binders_NoExprFilt_ord_homo=Moge_Binders_NoExprFilt_ord[Moge_Binders_NoExprFilt_ord$variable=='IC50_HLA-A24:02',]
#Moge_Binders_NoExprFilt_ord_homo=rbind(Moge_Binders_NoExprFilt_ord_homo,Moge_Binders_NoExprFilt_ord)
#Moge_Binders_NoExprFilt_ord_homo=Moge_Binders_NoExprFilt_ord_homo[order(Moge_Binders_NoExprFilt_ord_homo$value),]
#
#
#Piag_Binders_NoExprFilt_ord_homo=Piag_Binders_NoExprFilt_ord[Piag_Binders_NoExprFilt_ord$variable=='IC50_HLA-C04:01',]
#Piag_Binders_NoExprFilt_ord_homo=rbind(Piag_Binders_NoExprFilt_ord_homo,Piag_Binders_NoExprFilt_ord)
#Piag_Binders_NoExprFilt_ord_homo=Piag_Binders_NoExprFilt_ord_homo[order(Piag_Binders_NoExprFilt_ord_homo$value),]
#
#
#
#
#Calu_Binders_ord_homo=Calu_Binders_ord[Calu_Binders_ord$variable=='IC50_HLA-A02:01',]
#Calu_Binders_ord_homo=rbind(Calu_Binders_ord_homo,Calu_Binders_ord)
#Calu_Binders_ord_homo=Calu_Binders_ord_homo[order(Calu_Binders_ord_homo$value),]
#
#Luan_Binders_ord_homo=Luan_Binders_ord[Luan_Binders_ord$variable=='IC50_HLA-C07:01',]
#Luan_Binders_ord_homo=rbind(Luan_Binders_ord_homo,Luan_Binders_ord)
#Luan_Binders_ord_homo=Luan_Binders_ord_homo[order(Luan_Binders_ord_homo$value),]
#
#Moge_Binders_ord_homo=Moge_Binders_ord[Moge_Binders_ord$variable=='IC50_HLA-A24:02',]
#Moge_Binders_ord_homo=rbind(Moge_Binders_ord_homo,Moge_Binders_ord)
#Moge_Binders_ord_homo=Moge_Binders_ord_homo[order(Moge_Binders_ord_homo$value),]
#
#
#Piag_Binders_ord_homo=Piag_Binders_ord[Piag_Binders_ord$variable=='IC50_HLA-C04:01',]
#Piag_Binders_ord_homo=rbind(Piag_Binders_ord_homo,Piag_Binders_ord)
#Piag_Binders_ord_homo=Piag_Binders_ord_homo[order(Piag_Binders_ord_homo$value),]
#
#
#
##write.table(unique(prelu_NeoAG_mhci_MinorAntigens_filtered[,c(1:18,27:30)]),"Summary_PRELU_Annotato_MHCI_MinorAntigens_NoMismatches_AllInfos.txt",sep="\t",col.names=T,row.names=F,quote=F)
#
#all_binders=rbind(Besu_Binders_Genes_ord,Besu_Binders_Genes_ord,Calu_Binders_Genes_ord,Deiv_Binders_Genes_ord,Dr1_Binders_Genes_ord,Dr1_Binders_Genes_ord,Dr4_Binders_Genes_ord,Dr5_Binders_Genes_ord,Foca_Binders_Genes_ord,Gagra_Binders_Genes_ord,Luan_Binders_Genes_ord,Mabi_Binders_Genes_ord,Moge_Binders_Genes_ord,Piag_Binders_Genes_ord,Moge_Binders_Genes_ord)
#known_Binders=c('HMHA1','MYO1G','AKAP13','KIAA0020','HMHB1','BCL2A1','CENPM','SP110','UGT2B17')
#
#
#color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
#par(mar=c(5.1,5.1,5.1,5.1))
##ets the bottom, left, top and right margins respectively of the plot region in number of lines of text.
#
#c25 <- c("dodgerblue2","#E31A1C", # red
#                "green4",
#                "#6A3D9A", # purple
#                "#FF7F00", # orange
#                "black","gold1",
#                "skyblue2","#FB9A99", # lt pink
#                "palegreen2",
#                "#CAB2D6", # lt purple
#                "#FDBF6F", # lt orange
#                "gray70", "khaki2",
#                "maroon","orchid1","deeppink1","blue1","steelblue4",
#                "darkturquoise","green1","yellow4","yellow3",
#                "darkorange4","brown")
#
#
#color=c25[c(1:11,13,18,22,24)]
#plot(lines(loess.smooth(1:length(Besu_Binders_NoExprFilt_ord$value), y=Besu_Binders_NoExprFilt_ord$value)),xlim=c(0,3500),ylim=c(0,500),col=color[1],xlab='Number of Predicted mhAGs Binders_NoExprFilt',ylab='IC50 (nM)')
#lines(loess.smooth(1:length(Besu_Binders_NoExprFilt_ord$value), y=Besu_Binders_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,3500),ylim=c(0,500),col=color[1],xlab='',ylab='IC50 (nM)',lwd=3)
##polygon(1:length(Besu_Binders_NoExprFilt_ord$value), y=Besu_Binders_NoExprFilt_ord$value,col='skyblue')
#par(new=TRUE)
#lines(loess.smooth(1:length(Besu_Binders_NoExprFilt_ord$value), y=Besu_Binders_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,3500),ylim=c(0,500),col=color[2],xlab='',ylab='',lwd=3,lty=2)
#par(new=TRUE)
#lines(loess.smooth(1:length(Calu_Binders_NoExprFilt_ord_homo$value), y=Calu_Binders_NoExprFilt_ord_homo$value,degree=1,span=0.8),xlim=c(0,3500),ylim=c(0,500),col=color[3],xlab='',ylab='',lwd=3)
#par(new=TRUE)
#lines(loess.smooth(1:length(Deiv_Binders_NoExprFilt_ord$value), y=Deiv_Binders_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,3500),ylim=c(0,500),col=color[4],xlab='',ylab='',lwd=3,lty=2)
#par(new=TRUE)
#lines(loess.smooth(1:length(Dest_Binders_NoExprFilt_ord$value), y=Dest_Binders_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,3500),ylim=c(0,500),col=color[5],xlab='',ylab='',lwd=3)
#par(new=TRUE)
#lines(loess.smooth(1:length(Dr1_Binders_NoExprFilt_ord$value), y=Dr1_Binders_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,3500),ylim=c(0,500),col=color[6],xlab='',ylab='',lwd=3,lty=2)
#par(new=TRUE)
#lines(loess.smooth(1:length(Dr4_Binders_NoExprFilt_ord$value), y=Dr4_Binders_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,3500),ylim=c(0,500),col=color[7],xlab='',ylab='',lwd=3,)
#par(new=TRUE)
#lines(loess.smooth(1:length(Dr5_Binders_NoExprFilt_ord$value), y=Dr5_Binders_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,3500),ylim=c(0,500),col=color[8],xlab='',ylab='',lwd=3,)
#par(new=TRUE)
#lines(loess.smooth(1:length(Foca_Binders_NoExprFilt_ord$value), y=Foca_Binders_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,3500),ylim=c(0,500),col=color[9],xlab='',ylab='',lwd=3,lty=2)
#par(new=TRUE)
#lines(loess.smooth(1:length(Gagra_Binders_NoExprFilt_ord$value), y=Gagra_Binders_NoExprFilt_ord$value, degree=1,span=0.8),xlim=c(0,3500),ylim=c(0,500),col=color[10],xlab='',ylab='',lwd=3,lty=2)
#par(new=TRUE)
#lines(loess.smooth(1:length(Luan_Binders_NoExprFilt_ord_homo$value), y=Luan_Binders_NoExprFilt_ord_homo$value,  degree=1,span=0.8),xlim=c(0,3500),ylim=c(0,500),col=color[11],xlab='',ylab='',lwd=3,)
#par(new=TRUE)
#lines(loess.smooth(1:length(Mabi_Binders_NoExprFilt_ord$value), y=Mabi_Binders_NoExprFilt_ord$value,  degree=1,span=0.8),xlim=c(0,3500),ylim=c(0,500),col=color[12],xlab='',ylab='',lwd=3)
#par(new=TRUE)
#lines(loess.smooth(1:length(Moge_Binders_NoExprFilt_ord_homo$value), y=Moge_Binders_NoExprFilt_ord_homo$value,  degree=1,span=0.8),xlim=c(0,3500),ylim=c(0,500),col=color[13],xlab='',ylab='',lwd=3,lty=2)
#par(new=TRUE)
#lines(loess.smooth(1:length(Piag_Binders_NoExprFilt_ord_homo$value), y=Piag_Binders_NoExprFilt_ord_homo$value,  degree=1,span=0.8),xlim=c(0,3500),ylim=c(0,500),col=color[14],xlab='',ylab='',lwd=3,lty=2)
#par(new=TRUE)
#lines(loess.smooth(1:length(Moge_Binders_NoExprFilt_ord$value), y=Moge_Binders_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,3500),ylim=c(0,500),col=color[15],xlab='',ylab='',lwd=3)
#legend("topright", inset=c(-0,0), c('alfe','besu','calu','deiv','dest','dr1','dr4','dr5','foca','gagra','luan','mabi','moge','piag','prelu'),col=color[1:15], title="Patients",cex=0.7,lwd=rep(2,15),lty=c(1,2,1,2,1,2,1,1,2,2,1,1,2,2,1))
#
#
#
#
#color=c25[c(1:11,13,18,22,24)]
#plot(lines(loess.smooth(1:length(Alfe_Binders_NoExprFilt_ord$value), y=Alfe_Binders_NoExprFilt_ord$value)),xlim=c(0,1000),ylim=c(0,500),col=color[1],xlab='Number of Predicted mhAGs Binders',ylab='IC50 (nM)')
#lines(loess.smooth(1:length(Alfe_Binders_NoExprFilt_ord$value), y=Alfe_Binders_NoExprFilt_ord$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[1],xlab='',ylab='IC50 (nM)',lwd=3)
##polygon(1:length(Besu_Binders_ord$value), y=Besu_Binders_ord$value,col='skyblue')
#par(new=TRUE)
#lines(loess.smooth(1:length(Besu_Binders_ord$value), y=Besu_Binders_ord$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[2],xlab='',ylab='',lwd=3,lty=2)
#par(new=TRUE)
#lines(loess.smooth(1:length(Calu_Binders_ord_homo$value), y=Calu_Binders_ord_homo$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[3],xlab='',ylab='',lwd=3)
#par(new=TRUE)
#lines(loess.smooth(1:length(Deiv_Binders_ord$value), y=Deiv_Binders_ord$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[4],xlab='',ylab='',lwd=3,lty=2)
#par(new=TRUE)
#lines(loess.smooth(1:length(Dest_Binders_ord$value), y=Dest_Binders_ord$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[5],xlab='',ylab='',lwd=3)
#par(new=TRUE)
#lines(loess.smooth(1:length(Dr1_Binders_ord$value), y=Dr1_Binders_ord$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[6],xlab='',ylab='',lwd=3,lty=2)
#par(new=TRUE)
#lines(loess.smooth(1:length(Dr4_Binders_ord$value), y=Dr4_Binders_ord$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[7],xlab='',ylab='',lwd=3,)
#par(new=TRUE)
#lines(loess.smooth(1:length(Dr5_Binders_ord$value), y=Dr5_Binders_ord$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[8],xlab='',ylab='',lwd=3,)
#par(new=TRUE)
#lines(loess.smooth(1:length(Foca_Binders_ord$value), y=Foca_Binders_ord$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[9],xlab='',ylab='',lwd=3,lty=2)
#par(new=TRUE)
#lines(loess.smooth(1:length(Gagra_Binders_ord$value), y=Gagra_Binders_ord$value, degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[10],xlab='',ylab='',lwd=3,lty=2)
#par(new=TRUE)
#lines(loess.smooth(1:length(Luan_Binders_ord_homo$value), y=Luan_Binders_ord_homo$value,  degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[11],xlab='',ylab='',lwd=3,)
#par(new=TRUE)
#lines(loess.smooth(1:length(Mabi_Binders_ord$value), y=Mabi_Binders_ord$value,  degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[12],xlab='',ylab='',lwd=3)
#par(new=TRUE)
#lines(loess.smooth(1:length(Moge_Binders_ord_homo$value), y=Moge_Binders_ord_homo$value,  degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[13],xlab='',ylab='',lwd=3,lty=2)
#par(new=TRUE)
#lines(loess.smooth(1:length(Piag_Binders_ord_homo$value), y=Piag_Binders_ord_homo$value,  degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[14],xlab='',ylab='',lwd=3,lty=2)
#par(new=TRUE)
#lines(loess.smooth(1:length(Moge_Binders_ord$value), y=Moge_Binders_ord$value,degree=1,span=0.8),xlim=c(0,1000),ylim=c(0,500),col=color[15],xlab='',ylab='',lwd=3)
#legend("topright", inset=c(-0,0), c('alfe','besu','calu','deiv','dest','dr1','dr4','dr5','foca','gagra','luan','mabi','moge','piag','prelu'),col=color[1:15], title="Patients",cex=0.7,lwd=rep(2,15),lty=c(1,2,1,2,1,2,1,1,2,2,1,1,2,2,1))
#
#
#
#alfe_auc  =auc(1:length(Besu_Binders_ord$value), y=Besu_Binders_ord$value)
#besu_auc  =auc(1:length(Besu_Binders_ord$value), y=Besu_Binders_ord$value)
#calu_auc  =auc(1:length(Calu_Binders_ord_homo$value), y=Calu_Binders_ord_homo$value)
#deiv_auc  =auc(1:length(Deiv_Binders_ord$value), y=Deiv_Binders_ord$value)
#dest_auc  =auc(1:length(Dest_Binders_ord$value), y=Dest_Binders_ord$value)
#dr1_auc   =auc(1:length(Dr1_Binders_ord$value), y=Dr1_Binders_ord$value)
#dr4_auc   =auc(1:length(Dr4_Binders_ord$value), y=Dr4_Binders_ord$value)
#dr5_auc   =auc(1:length(Dr5_Binders_ord$value), y=Dr5_Binders_ord$value)
#foca_auc  =auc(1:length(Foca_Binders_ord$value), y=Foca_Binders_ord$value)
#gagra_auc =auc(1:length(Gagra_Binders_ord$value), y=Gagra_Binders_ord$value)
#luan_auc  =auc(1:length(Luan_Binders_ord_homo$value), y=Luan_Binders_ord_homo$value)
#mabi_auc  =auc(1:length(Mabi_Binders_ord$value), y=Mabi_Binders_ord$value)
#moge_auc  =auc(1:length(Moge_Binders_ord_homo$value), y=Moge_Binders_ord_homo$value)
#piag_auc  =auc(1:length(Piag_Binders_ord_homo$value), y=Piag_Binders_ord_homo$value)
#prelu_auc =auc(1:length(Prelu_Binders_ord$value), y=Prelu_Binders_ord$value)
#
#
#
#legend("topright", inset=c(-0,0), c('alfe','besu','calu','deiv','dest','dr1','dr4','dr5','foca','gagra','luan','mabi','moge','piag','prelu'),col=color[1:15], title="Patients",cex=0.7,lwd=rep(2,15),lty=c(1,2,1,2,1,2,1,1,2,2,1,1,2,2,1)) 
# alfe_auc 
#[1] 18883.07
#> besu_auc 
#[1] 110785.9
#> calu_auc 
#[1] 31014.48
#> deiv_auc 
#[1] 119524.4
#> dest_auc 
#[1] 35478.98
#> dr1_auc  
#[1] 154515.3
#> dr4_auc  
#[1] 70518.68
#> dr5_auc  
#[1] 4815.255
#> foca_auc 
#[1] 75342.28
#> gagra_auc
#[1] 123475.1
#> luan_auc 
#[1] 35164.39
#> mabi_auc 
#[1] 54018.22
#> moge_auc 
#[1] 39343.16
#> piag_auc 
#[1] 38724.35
#> prelu_auc
#[1] 12145.08
#
#
#
#
#Besu_Binders$Fc=Besu_Binders$ALFE2_kallisto-Besu_Binders$ALFE1_kallisto
#Besu_Binders$Fc=Besu_Binders$BESU2_kallisto-Besu_Binders$BESU1_kallisto
#Calu_Binders$Fc=Calu_Binders$CALU2_kallisto-Calu_Binders$CALU1_kallisto
#Deiv_Binders$Fc=Deiv_Binders$DEIV2b_kallisto-Deiv_Binders$DEIV2_kallisto
#Dest_Binders$Fc=Dest_Binders$DEST2_kallisto-Dest_Binders$DEST1_kallisto
#Dr1_Binders$Fc=Dr1_Binders$DR1_2_kallisto-Dr1_Binders$DR1_1_kallisto
#Dr4_Binders$Fc=Dr4_Binders$DR42_kallisto-Dr4_Binders$DR41_kallisto
#Dr5_Binders$Fc=Dr5_Binders$DR52_kallisto-Dr5_Binders$DR51_kallisto
#Foca_Binders$Fc=Foca_Binders$FOCA2_kallisto-Foca_Binders$FOCA1_kallisto
#Gagra_Binders$Fc=Gagra_Binders$GAGRA2_kallisto-Gagra_Binders$GAGRA1_kallisto
#Luan_Binders$Fc=Luan_Binders$LUAN2_kallisto-Luan_Binders$LUAN1_kallisto
#Mabi_Binders$Fc=Mabi_Binders$MABI2_kallisto-Mabi_Binders$MABI1_kallisto
#Moge_Binders$Fc=Moge_Binders$MOGE2_kallisto-Moge_Binders$MOGE1_kallisto
#Piag_Binders$Fc=Piag_Binders$PIAG2_kallisto-Piag_Binders$PIAG1_kallisto
#Moge_Binders$Fc=Moge_Binders$PRELU2_kallisto-Moge_Binders$PRELU1_kallisto
#
#
##Trasformo tutto in TPM
#tpm <- function(counts, lengths) {
#rate <- counts / lengths
#rate / sum(rate) * 1e6
#}
#
#transcripts_lengths <- apply(txi_length$length, 1, max)
#transcripts_lengths=as.data.frame(transcripts_lengths)
#transcripts_lengths$genes=rownames(transcripts_lengths)
#
#countdata_filt <- countdata[rowSums(countdata[, 1:30]) > 0, ]
#countdata_filtDGE_Nothreshold=DGEList(countdata_filt[,1:30],genes=rownames(countdata_filt),group=condition)
#keep_nothreshold <- rowSums(cpm(countdata_filtDGE_Nothreshold)>1) >= 1 #PRENDO SOLO I GENI CON CPM>1 IN ALMENO 8 CAMPIONI
#countdata_filtDGE_keep_Nothreshold <- countdata_filtDGE_Nothreshold[keep_nothreshold, , keep.lib.sizes=FALSE]
#countdata_filtDGE_keep_Nothreshold=calcNormFactors(countdata_filtDGE_keep_Nothreshold,method="TMM")
#countdata_filtDGE_Nothreshold_CPM=as.data.frame(cpm(countdata_filtDGE_Nothreshold,log=T))
#countdata_filtDGE_Nothreshold_CPM$genes=rownames(countdata_filtDGE_Nothreshold_CPM)
#
#conte=as.data.frame(countdata_filtDGE_Nothreshold$counts)
#conte$genes=rownames(conte)
#transcripts_lengths=transcripts_lengths[transcripts_lengths$genes %in% conte$genes,]
#conte_tpm=tpm(countdata_filtDGE_Nothreshold$counts,transcripts_lengths$transcripts_lengths)
#
#conte_tpm=as.data.frame(conte_tpm)
##Ora su questi dati creo dei files RNK
##Inizio con i FCs
#conte_tpm$BesuFC  =  conte_tpm$ALFE2_kallisto   -   conte_tpm$ALFE1_kallisto 
#conte_tpm$BesuFC  =  conte_tpm$BESU2_kallisto   -   conte_tpm$BESU1_kallisto 
#conte_tpm$CaluFC  =  conte_tpm$CALU2_kallisto   -   conte_tpm$CALU1_kallisto 
#conte_tpm$DeivFC  =  conte_tpm$DEIV2b_kallisto  -   conte_tpm$DEIV2_kallisto
#conte_tpm$DestFC  =  conte_tpm$DEST2_kallisto   -   conte_tpm$DEST1_kallisto 
#conte_tpm$Dr1FC   =   conte_tpm$DR1_2_kallisto    -  conte_tpm$DR1_1_kallisto 
#conte_tpm$Dr4FC   =   conte_tpm$DR42_kallisto     -  conte_tpm$DR41_kallisto 
#conte_tpm$Dr5FC   =   conte_tpm$DR52_kallisto     -  conte_tpm$DR51_kallisto 
#conte_tpm$FocaFC  =  conte_tpm$FOCA2_kallisto   -   conte_tpm$FOCA1_kallisto 
#conte_tpm$GagraFC = conte_tpm$GAGRA2_kallisto -    conte_tpm$GAGRA1_kallisto
#conte_tpm$LuanC   =   conte_tpm$LUAN2_kallisto    -  conte_tpm$LUAN1_kallisto 
#conte_tpm$MabiFC  =  conte_tpm$MABI2_kallisto   -   conte_tpm$MABI1_kallisto 
#conte_tpm$MogeFC  =  conte_tpm$MOGE2_kallisto   -   conte_tpm$MOGE1_kallisto 
#conte_tpm$PiagFC  =  conte_tpm$PIAG2_kallisto   -   conte_tpm$PIAG1_kallisto 
#conte_tpm$MogeFC = conte_tpm$PRELU2_kallisto -    conte_tpm$PRELU1_kallisto
#
#
#conte_tpm$BesuRNK=   ifelse((conte_tpm$ALFE2_kallisto  >= 1.5  | conte_tpm$ALFE1_kallisto   >= 1.5 ) & ( conte_tpm$BesuFC   <= -1.5 | conte_tpm$BesuFC   >= 1.5),conte_tpm$BesuFC  ,NA) 
#conte_tpm$BesuRNK=   ifelse((conte_tpm$BESU2_kallisto  >= 1.5  | conte_tpm$BESU1_kallisto   >= 1.5 ) & ( conte_tpm$BesuFC   <= -1.5 | conte_tpm$BesuFC   >= 1.5),conte_tpm$BesuFC  ,NA) 
#conte_tpm$CaluRNK=   ifelse((conte_tpm$CALU2_kallisto  >= 1.5  | conte_tpm$CALU1_kallisto   >= 1.5 ) & ( conte_tpm$CaluFC   <= -1.5 | conte_tpm$CaluFC   >= 1.5),conte_tpm$CaluFC  ,NA) 
#conte_tpm$DeivRNK=   ifelse((conte_tpm$DEIV2b_kallisto >= 1.5  | conte_tpm$DEIV2_kallisto  >= 1.5 ) & (  conte_tpm$DeivFC   <= -1.5 | conte_tpm$DeivFC   >= 1.5),conte_tpm$DeivFC  ,NA) 
#conte_tpm$DestRNK=   ifelse((conte_tpm$DEST2_kallisto  >= 1.5  | conte_tpm$DEST1_kallisto   >= 1.5 ) & ( conte_tpm$DestFC   <= -1.5 | conte_tpm$DestFC   >= 1.5),conte_tpm$DestFC  ,NA) 
#conte_tpm$Dr1RNK=    ifelse((conte_tpm$DR1_2_kallisto  >= 1.5  | conte_tpm$DR1_1_kallisto   >= 1.5 ) & ( conte_tpm$Dr1FC    <= -1.5 | conte_tpm$Dr1FC    >= 1.5),conte_tpm$Dr1FC   ,NA) 
#conte_tpm$Dr4RNK=    ifelse((conte_tpm$DR42_kallisto   >= 1.5  | conte_tpm$DR41_kallisto   >= 1.5 ) & (  conte_tpm$Dr4FC    <= -1.5 | conte_tpm$Dr4FC    >= 1.5),conte_tpm$Dr4FC   ,NA) 
#conte_tpm$Dr5RNK=    ifelse((conte_tpm$DR52_kallisto   >= 1.5  | conte_tpm$DR51_kallisto   >= 1.5 ) & (  conte_tpm$Dr5FC    <= -1.5 | conte_tpm$Dr5FC    >= 1.5),conte_tpm$Dr5FC   ,NA) 
#conte_tpm$FocaRNK=   ifelse((conte_tpm$FOCA2_kallisto  >= 1.5  | conte_tpm$FOCA1_kallisto   >= 1.5 ) & ( conte_tpm$FocaFC   <= -1.5 | conte_tpm$FocaFC   >= 1.5),conte_tpm$FocaFC  ,NA) 
#conte_tpm$GagraRNK = ifelse((conte_tpm$GAGRA2_kallisto >= 1.5  | conte_tpm$GAGRA1_kallisto  >= 1.5 ) & ( conte_tpm$GagraFC  <= -1.5 | conte_tpm$GagraFC  >= 1.5),conte_tpm$GagraFC ,NA) 
#conte_tpm$LuanRNK =  ifelse((conte_tpm$LUAN2_kallisto  >= 1.5  | conte_tpm$LUAN1_kallisto   >= 1.5 ) & ( conte_tpm$LuanC    <= -1.5 | conte_tpm$LuanC    >= 1.5),conte_tpm$LuanC   ,NA) 
#conte_tpm$MabiRNK =  ifelse((conte_tpm$MABI2_kallisto  >= 1.5  | conte_tpm$MABI1_kallisto   >= 1.5 ) & ( conte_tpm$MabiFC   <= -1.5 | conte_tpm$MabiFC   >= 1.5),conte_tpm$MabiFC  ,NA) 
#conte_tpm$MogeRNK =  ifelse((conte_tpm$MOGE2_kallisto  >= 1.5  | conte_tpm$MOGE1_kallisto   >= 1.5 ) & ( conte_tpm$MogeFC   <= -1.5 | conte_tpm$MogeFC   >= 1.5),conte_tpm$MogeFC  ,NA) 
#conte_tpm$PiagRNK =  ifelse((conte_tpm$PIAG2_kallisto  >= 1.5  | conte_tpm$PIAG1_kallisto   >= 1.5 ) & ( conte_tpm$PiagFC   <= -1.5 | conte_tpm$PiagFC   >= 1.5),conte_tpm$PiagFC  ,NA) 
#conte_tpm$MogeRNK = ifelse((conte_tpm$PRELU2_kallisto >= 1.5  | conte_tpm$PRELU1_kallisto  >= 1.5 ) & ( conte_tpm$MogeFC  <= -1.5 | conte_tpm$MogeFC  >= 1.5),conte_tpm$MogeFC ,NA) 
#
#
#conte_tpm$genes=rownames(conte_tpm)
##Analisi GSEA
#
#
#
#ALFE_RNK_tpm=conte_tpm[conte_tpm$ALFE2_kallisto>= 1.5 | conte_tpm$ALFE1_kallisto>=1.5,c(61,2,1)]
#BESU_RNK_tpm=conte_tpm[conte_tpm$BESU2_kallisto>= 1.5 | conte_tpm$BESU1_kallisto>=1.5,c(61,4,3)]
#CALU_RNK_tpm=conte_tpm[conte_tpm$CALU2_kallisto>= 1.5 | conte_tpm$CALU1_kallisto>=1.5,c(61,6,5)]
#DEIV_RNK_tpm=conte_tpm[conte_tpm$DEIV2b_kallisto>= 1.5 | conte_tpm$DEIV2_kallisto>=1.5,c(61,8,7)]
#DEST_RNK_tpm=conte_tpm[conte_tpm$DEST2_kallisto>= 1.5 | conte_tpm$DEST1_kallisto>=1.5,c(61,10,9)]
#DR1_RNK_tpm=conte_tpm[conte_tpm$DR1_2_kallisto>= 1.5 | conte_tpm$DR1_1_kallisto>=1.5,c(61,12,11)]
#DR4_RNK_tpm=conte_tpm[conte_tpm$DR42_kallisto>= 1.5 | conte_tpm$DR41_kallisto>=1.5,c(61,14,13)]
#DR5_RNK_tpm=conte_tpm[conte_tpm$DR52_kallisto>= 1.5 | conte_tpm$DR51_kallisto>=1.5,c(61,16,15)]
#FOCA_RNK_tpm=conte_tpm[conte_tpm$FOCA2_kallisto>= 1.5 | conte_tpm$FOCA1_kallisto>=1.5,c(61,18,17)]
#GAGRA_RNK_tpm=conte_tpm[conte_tpm$GAGRA2_kallisto>= 1.5 | conte_tpm$GAGRA1_kallisto>=1.5,c(61,20,19)]
#LUAN_RNK_tpm=conte_tpm[conte_tpm$LUAN2_kallisto>= 1.5 | conte_tpm$LUAN1_kallisto>=1.5,c(61,22,21)]
#MABI_RNK_tpm=conte_tpm[conte_tpm$MABI2_kallisto>= 1.5 | conte_tpm$MABI1_kallisto>=1.5,c(61,24,23)]
#MOGE_RNK_tpm=conte_tpm[conte_tpm$MOGE2_kallisto>= 1.5 | conte_tpm$MOGE1_kallisto>=1.5,c(61,26,25)]
#PIAG_RNK_tpm=conte_tpm[conte_tpm$PIAG2_kallisto>= 1.5 | conte_tpm$PIAG1_kallisto>=1.5,c(61,28,27)]
#PRELU_RNK_tpm=conte_tpm[conte_tpm$PRELU2_kallisto>= 1.5 | conte_tpm$PRELU1_kallisto>=1.5,c(61,30,29)]
#
#
#
#pdf('HumanProteinAtlas_Expression_mhAGs_15PAtients.pdf')
#
#alfe_tissue=merge(Besu_Binders,rna_tissue,by='genes')
#alfe_tissue=alfe_tissue[(alfe_tissue$Sample=='bone marrow'| alfe_tissue$Sample=='skin'| alfe_tissue$Sample =='testis'| alfe_tissue$Sample=='liver'| alfe_tissue$Sample=='small intestine') & alfe_tissue$Value>=1, ]
#ggplot(alfe_tissue, aes(x=genes,y=Value, fill=Sample)) +geom_bar(stat="identity") + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) + ggtitle("Expression in tissues: ALFE") + scale_fill_brewer(palette="Dark2")
#
#
#besu_tissue=merge(Besu_Binders,rna_tissue,by='genes')
#besu_tissue=besu_tissue[(besu_tissue$Sample=='bone marrow'| besu_tissue$Sample=='skin'| besu_tissue$Sample =='testis'| besu_tissue$Sample=='liver'| besu_tissue$Sample=='small intestine') & besu_tissue$Value>=1, ]
#ggplot(besu_tissue, aes(x=genes,y=Value, fill=Sample)) +geom_bar(stat="identity") + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) + ggtitle("Expression in tissues: BESU") + scale_fill_brewer(palette="Dark2")
#
#
#calu_tissue=merge(Calu_Binders,rna_tissue,by='genes')
#calu_tissue=calu_tissue[(calu_tissue$Sample=='bone marrow'| calu_tissue$Sample=='skin'| calu_tissue$Sample =='testis'| calu_tissue$Sample=='liver'| calu_tissue$Sample=='small intestine') & calu_tissue$Value>=1, ]
#ggplot(calu_tissue, aes(x=genes,y=Value, fill=Sample)) +geom_bar(stat="identity") + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) + ggtitle("Expression in tissues: CALU") + scale_fill_brewer(palette="Dark2")
#
#deiv_tissue=merge(Deiv_Binders,rna_tissue,by='genes')
#deiv_tissue=deiv_tissue[(deiv_tissue$Sample=='bone marrow'| deiv_tissue$Sample=='skin'| deiv_tissue$Sample =='testis'| deiv_tissue$Sample=='liver'| deiv_tissue$Sample=='small intestine') & deiv_tissue$Value>=1, ]
#ggplot(deiv_tissue, aes(x=genes,y=Value, fill=Sample)) +geom_bar(stat="identity") + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) + ggtitle("Expression in tissues: DEIV") + scale_fill_brewer(palette="Dark2")
#
#
#dest_tissue=merge(Dest_Binders,rna_tissue,by='genes')
#dest_tissue=dest_tissue[(dest_tissue$Sample=='bone marrow'| dest_tissue$Sample=='skin'| dest_tissue$Sample =='testis'| dest_tissue$Sample=='liver'| dest_tissue$Sample=='small intestine') & dest_tissue$Value>=1, ]
#ggplot(dest_tissue, aes(x=genes,y=Value, fill=Sample)) +geom_bar(stat="identity") + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) + ggtitle("Expression in tissues: DEST") + scale_fill_brewer(palette="Dark2")
#
#dr1_tissue=merge(Dr1_Binders,rna_tissue,by='genes')
#dr1_tissue=dr1_tissue[(dr1_tissue$Sample=='bone marrow'| dr1_tissue$Sample=='skin'| dr1_tissue$Sample =='testis'| dr1_tissue$Sample=='liver'| dr1_tissue$Sample=='small intestine') & dr1_tissue$Value>=1, ]
#ggplot(dr1_tissue, aes(x=genes,y=Value, fill=Sample)) +geom_bar(stat="identity") + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) + ggtitle("Expression in tissues: DR1") + scale_fill_brewer(palette="Dark2")
#
#dr4_tissue=merge(Dr4_Binders,rna_tissue,by='genes')
#dr4_tissue=dr4_tissue[(dr4_tissue$Sample=='bone marrow'| dr4_tissue$Sample=='skin'| dr4_tissue$Sample =='testis'| dr4_tissue$Sample=='liver'| dr4_tissue$Sample=='small intestine') & dr4_tissue$Value>=1, ]
#ggplot(dr4_tissue, aes(x=genes,y=Value, fill=Sample)) +geom_bar(stat="identity") + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) + ggtitle("Expression in tissues: DR4") + scale_fill_brewer(palette="Dark2")
#
#dr5_tissue=merge(Dr5_Binders,rna_tissue,by='genes')
#dr5_tissue=dr5_tissue[(dr5_tissue$Sample=='bone marrow'| dr5_tissue$Sample=='skin'| dr5_tissue$Sample =='testis'| dr5_tissue$Sample=='liver'| dr5_tissue$Sample=='small intestine') & dr5_tissue$Value>=1, ]
#ggplot(dr5_tissue, aes(x=genes,y=Value, fill=Sample)) +geom_bar(stat="identity") + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) + ggtitle("Expression in tissues: DR5") + scale_fill_brewer(palette="Dark2")
#
#foca_tissue=merge(Foca_Binders,rna_tissue,by='genes')
#foca_tissue=foca_tissue[(foca_tissue$Sample=='bone marrow'| foca_tissue$Sample=='skin'| foca_tissue$Sample =='testis'| foca_tissue$Sample=='liver'| foca_tissue$Sample=='small intestine') & foca_tissue$Value>=1, ]
#ggplot(foca_tissue, aes(x=genes,y=Value, fill=Sample)) +geom_bar(stat="identity") + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) + ggtitle("Expression in tissues: FOCA") + scale_fill_brewer(palette="Dark2")
#
#gagra_tissue=merge(Gagra_Binders,rna_tissue,by='genes')
#gagra_tissue=gagra_tissue[(gagra_tissue$Sample=='bone marrow'| gagra_tissue$Sample=='skin'| gagra_tissue$Sample =='testis'| gagra_tissue$Sample=='liver'| gagra_tissue$Sample=='small intestine') & gagra_tissue$Value>=1, ]
#ggplot(gagra_tissue, aes(x=genes,y=Value, fill=Sample)) +geom_bar(stat="identity") + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) + ggtitle("Expression in tissues: GAGRA") + scale_fill_brewer(palette="Dark2")
#
#luan_tissue=merge(Luan_Binders,rna_tissue,by='genes')
#luan_tissue=luan_tissue[(luan_tissue$Sample=='bone marrow'| luan_tissue$Sample=='skin'| luan_tissue$Sample =='testis'| luan_tissue$Sample=='liver'| luan_tissue$Sample=='small intestine') & luan_tissue$Value>=1, ]
#ggplot(luan_tissue, aes(x=genes,y=Value, fill=Sample)) +geom_bar(stat="identity") + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) + ggtitle("Expression in tissues: LUAN") + scale_fill_brewer(palette="Dark2")
#
#mabi_tissue=merge(Mabi_Binders,rna_tissue,by='genes')
#mabi_tissue=mabi_tissue[(mabi_tissue$Sample=='bone marrow'| mabi_tissue$Sample=='skin'| mabi_tissue$Sample =='testis'| mabi_tissue$Sample=='liver'| mabi_tissue$Sample=='small intestine') & mabi_tissue$Value>=1, ]
#ggplot(mabi_tissue, aes(x=genes,y=Value, fill=Sample)) +geom_bar(stat="identity") + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) + ggtitle("Expression in tissues: MABI") + scale_fill_brewer(palette="Dark2")
#
#moge_tissue=merge(Moge_Binders,rna_tissue,by='genes')
#moge_tissue=moge_tissue[(moge_tissue$Sample=='bone marrow'| moge_tissue$Sample=='skin'| moge_tissue$Sample =='testis'| moge_tissue$Sample=='liver'| moge_tissue$Sample=='small intestine') & moge_tissue$Value>=1, ]
#ggplot(moge_tissue, aes(x=genes,y=Value, fill=Sample)) +geom_bar(stat="identity") + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) + ggtitle("Expression in tissues: MOGE") + scale_fill_brewer(palette="Dark2")
#
#piag_tissue=merge(Piag_Binders,rna_tissue,by='genes')
#piag_tissue=piag_tissue[(piag_tissue$Sample=='bone marrow'| piag_tissue$Sample=='skin'| piag_tissue$Sample =='testis'| piag_tissue$Sample=='liver'| piag_tissue$Sample=='small intestine') & piag_tissue$Value>=1, ]
#ggplot(piag_tissue, aes(x=genes,y=Value, fill=Sample)) +geom_bar(stat="identity") + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) + ggtitle("Expression in tissues: PIAG") + scale_fill_brewer(palette="Dark2")
#
#prelu_tissue=merge(Moge_Binders,rna_tissue,by='genes')
#prelu_tissue=prelu_tissue[(prelu_tissue$Sample=='bone marrow'| prelu_tissue$Sample=='skin'| prelu_tissue$Sample =='testis'| prelu_tissue$Sample=='liver'| prelu_tissue$Sample=='small intestine') & prelu_tissue$Value>=1, ]
#ggplot(prelu_tissue, aes(x=genes,y=Value, fill=Sample)) +geom_bar(stat="identity") + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) + ggtitle("Expression in tissues: PRELU") + scale_fill_brewer(palette="Dark2")
#
#dev.off()