#!/bin/sh
#GLI ARGOMENTI SONO:
#INPUT_FILE=$1
#OUTPUT_FILE=$2
#LISTA_ALLELI=$3
#

name_input=$(printf $1 | awk -F/ '{print $NF}' | sed 's/\..*//g') #prende solo la radice del nome del file  del sample (tutto quello che viene prima del punto ), che poi mi servira per rinominare i nuovi files
name_output=$(printf $2 | awk -F/ '{print $NF}' | sed 's/\..*//g')
input_dir=$(dirname $1)
output_dir=$(dirname $2)

echo $name_input
echo $name_output

echo $input_dir
echo $output_dir

echo $input_dir"/"$name_input
echo $output_dir"/"$name_output

haplotype=$(echo  $3 | tr ',' '\t' ) 

echo $haplotype

sed -e 's/\"//g' $1 | awk '{OFS="\t"; if ($15!="." && length($3)==1 && length($4)==1 && $29 != "0/0" && $30 =="0/0") print $15,$13,$11,$28,$25,$26,$1,$2,$3,$4}' | sed 's/p\.//g' | sed -e 's/^[[:alpha:]]\{3\}/&\t/g'  |  sed  's/\s[[:digit:]]*/&\t/g'  - > $input_dir"/"$name_input"_ForEpitopes_DiagnosisOnly.txt"
sed -e 's/\"//g' $1 | awk '{OFS="\t"; if ($15!="." && length($3)==1 && length($4)==1 && $29 == "0/0" && $30 !="0/0") print $15,$13,$11,$28,$25,$26,$1,$2,$3,$4}'| sed 's/p\.//g' | sed -e 's/^[[:alpha:]]\{3\}/&\t/g'  |  sed  's/\s[[:digit:]]*/&\t/g'  - > $input_dir"/"$name_input"_ForEpitopes_RelapseOnly.txt"
sed -e 's/\"//g' $1 | awk '{OFS="\t"; if ($15!="." && length($3)==1 && length($4)==1 && $29 != "0/0" && $30 !="0/0") print $15,$13,$11,$28,$25,$26,$1,$2,$3,$4}' | sed 's/p\.//g' | sed -e 's/^[[:alpha:]]\{3\}/&\t/g'  |  sed  's/\s[[:digit:]]*/&\t/g'  - > $input_dir"/"$name_input"_ForEpitopes_CommonDiagnosisRelapse.txt"

sed -i '/stop_gain/d'  $input_dir"/"$name_input"_ForEpitopes_DiagnosisOnly.txt" 
sed -i '/stop_gain/d'  $input_dir"/"$name_input"_ForEpitopes_RelapseOnly.txt"
sed -i '/stop_gain/d'  $input_dir"/"$name_input"_ForEpitopes_CommonDiagnosisRelapse.txt"

Rscript /home/fsantaniello/snp2epi_new_MHCII.R $input_dir"/"$name_input"_ForEpitopes_DiagnosisOnly.txt" $4 $output_dir"/"$name_output"_PeptidesPrediction_DiagnosisOnly.txt"
Rscript /home/fsantaniello/snp2epi_new_MHCII.R $input_dir"/"$name_input"_ForEpitopes_RelapseOnly.txt" $4 $output_dir"/"$name_output"_PeptidesPrediction_RelapseOnly.txt"
Rscript /home/fsantaniello/snp2epi_new_MHCII.R $input_dir"/"$name_input"_ForEpitopes_CommonDiagnosisRelapse.txt" $4 $output_dir"/"$name_output"_PeptidesPrediction_CommonDiagnosisRelapse.txt"

awk '{if (($8~"MODERATE"||$8~"HIGH")&&$17~"True") print ">WT_"$1"; Transcript: "$5"; Effect: "$7"\n"$15"\n>MUT_"$1"; Transcript: "$5"; Effect: "$7"\n"$16"\n"}' $output_dir"/"$name_output"_PeptidesPrediction_DiagnosisOnly.txt" > $output_dir"/FastaSequences_Peptides"$name_output"_DiagnosisOnly.fa"
awk '{if (($8~"MODERATE"||$8~"HIGH")&&$17~"True") print ">WT_"$1"; Transcript: "$5"; Effect: "$7"\n"$15"\n>MUT_"$1"; Transcript: "$5"; Effect: "$7"\n"$16"\n"}' $output_dir"/"$name_output"_PeptidesPrediction_RelapseOnly.txt" > $output_dir"/FastaSequences_Peptides"$name_output"_RelapseOnly.fa"
awk '{if (($8~"MODERATE"||$8~"HIGH")&&$17~"True") print ">WT_"$1"; Transcript: "$5"; Effect: "$7"\n"$15"\n>MUT_"$1"; Transcript: "$5"; Effect: "$7"\n"$16"\n"}' $output_dir"/"$name_output"_PeptidesPrediction_CommonDiagnosisRelapse.txt" > $output_dir"/FastaSequences_Peptides"$name_output"_CommonRelapseDiagnosis.fa"

grep -v "\*" $output_dir"/FastaSequences_Peptides"$name_output"_DiagnosisOnly.fa" > tt && mv tt $output_dir"/FastaSequences_Peptides"$name_output"_DiagnosisOnly.fa"	
grep -v "\*"  $output_dir"/FastaSequences_Peptides"$name_output"_RelapseOnly.fa" > tt && mv tt  $output_dir"/FastaSequences_Peptides"$name_output"_RelapseOnly.fa"	
grep -v "\*" $output_dir"/FastaSequences_Peptides"$name_output"_CommonRelapseDiagnosis.fa" > tt && mv tt $output_dir"/FastaSequences_Peptides"$name_output"_CommonRelapseDiagnosis.fa"

echo -e "*************************************************************************************************\nCALCULATING BINDING AFFINITY FOR"$name_input" DIAGNOSIS ONLY\n*************************************************************************************************"

/home/fsantaniello/netMHC2pan/netMHCIIpan-3.1/netMHCIIpan -length 15 -filter 1 -xls -xlsfile $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_DiagnosisOnly.xls" -f $output_dir"/FastaSequences_Peptides"$name_output"_DiagnosisOnly.fa" -a $3  


echo -e "*************************************************************************************************\nCALCULATING BINDING AFFINITY FOR"$name_input" RELAPSE ONLY\n*************************************************************************************************"

/home/fsantaniello/netMHC2pan/netMHCIIpan-3.1/netMHCIIpan -length 15 -filter 1 -xls -xlsfile $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_RelapseOnly.xls" -f $output_dir"/FastaSequences_Peptides"$name_output"_RelapseOnly.fa" -a $3  


echo -e "*************************************************************************************************\nCALCULATING BINDING AFFINITY FOR"$name_input" RELAPSE DIAGNOSIS COMMON\n**********************************************************************************************"

/home/fsantaniello/netMHC2pan/netMHCIIpan-3.1/netMHCIIpan -length 15 -filter 1 -xls -xlsfile $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_CommonRelapseDiagnosis.xls" -f $output_dir"/FastaSequences_Peptides"$name_output"_CommonRelapseDiagnosis.fa" -a $3 

condeactivate
#
#
##Faccio anche un summary dei binders
haplotype=$(echo  $3 | tr ',' '\t' ) 
length=$(echo $haplotype | wc -w)

 if [ "$length" = 1 ]
     then 
	echo "*************************************************************************************************Haplotype is "$haplotype" *************************************************************************************************"
        #Summary Binders Diagnosis Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\t" $haplotype)
        more +2 $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_DiagnosisOnly.xls"| awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1}' - | (echo -e $weak_binders; cat -) > $output_dir"/RESULTS_"$name_output"_weak_binders_DiagOnly"
         #/home/fsantaniello/Summary_Epitopi.py $input_dir"/"$name_input"_ForEpitopes_DiagnosisOnly.txt" /home/fsantaniello/weak_binders_DiagOnly $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_DiagnosisOnly.xls" $output_dir"/RESULTS_"$name_output"_weak_binders_annotated_DiagOnly"
        strong_binders=$(echo "Strong binders\nGene\tPeptide\t" $haplotype)
        more +2 $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_DiagnosisOnly.xls" | awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1}' - | (echo -e $strong_binders; cat -) > $output_dir"/RESULTS_"$name_output"_strong_binders_DiagOnly"
         #/home/fsantaniello/Summary_Epitopi.py $input_dir"/"$name_input"_ForEpitopes_DiagnosisOnly.txt" /home/fsantaniello/strong_binders_DiagOnly $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_DiagnosisOnly.xls" /home/fsantaniello/strong_binders_annotated_DiagOnly

        #cat /home/fsantaniello/strong_binders_annotated_DiagOnly /home/fsantaniello/weak_binders_annotated_DiagOnly > $output_dir"/Summary_Binders_RESULTS_"$name_output"__PeptidesPrediction_DiagnosisOnly.txt"



        #Summary Binders Relapse Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\t" $haplotype)
        more +2 $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_RelapseOnly.xls"| awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1}' - | (echo -e $weak_binders; cat -) > $output_dir"/RESULTS_"$name_output"_weak_binders_RelOnly"
         #/home/fsantaniello/Summary_Epitopi.py $input_dir"/"$name_input"_ForEpitopes_RelapseOnly.txt" /home/fsantaniello/weak_binders_RelOnly $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_RelapseOnly.xls" /home/fsantaniello/weak_binders_annotated_RelOnly

        strong_binders=$(echo "Strong binders\nGene\tPeptide\t" $haplotype)
        more +2 $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_RelapseOnly.xls" | awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1}' - | (echo -e $strong_binders; cat -) > $output_dir"/RESULTS_"$name_output"_strong_binders_RelOnly"
         #/home/fsantaniello/Summary_Epitopi.py $input_dir"/"$name_input"_ForEpitopes_RelapseOnly.txt" /home/fsantaniello/strong_binders_RelOnly $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_RelapseOnly.xls" /home/fsantaniello/strong_binders_annotated_RelOnly
        
        cat /home/fsantaniello/strong_binders_annotated_RelOnly /home/fsantaniello/weak_binders_annotated_RelOnly > $output_dir"/Summary_Binders_RESULTS_"$name_output"__PeptidesPrediction_RelapseOnly.txt"
        
        
        #Summary Binders Common Diagnosis Relapse
        weak_binders=$(echo "Weak binders\nGene\tPeptide\t" $haplotype)
        more +2 $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_CommonRelapseDiagnosis.xls"| awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1}' - | (echo -e $weak_binders; cat -) > $output_dir"/RESULTS_"$name_output"_weak_binders_RelDiag"
         #/home/fsantaniello/Summary_Epitopi.py $input_dir"/"$name_input"_ForEpitopes_CommonDiagnosisRelapse.txt" /home/fsantaniello/weak_binders_RelDiag $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_CommonRelapseDiagnosis.xls" /home/fsantaniello/weak_binders_annotated_RelDiag
        strong_binders=$(echo "Strong binders\nGene\tPeptide\t" $haplotype)
        more +2 $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_CommonRelapseDiagnosis.xls" | awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1}' - | (echo -e $strong_binders; cat -) > $output_dir"/RESULTS_"$name_output"_strong_binders_RelDiag"
         #/home/fsantaniello/Summary_Epitopi.py $input_dir"/"$name_input"_ForEpitopes_CommonDiagnosisRelapse.txt" /home/fsantaniello/strong_binders_RelDiag $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_CommonRelapseDiagnosis.xls" /home/fsantaniello/strong_binders_annotated_RelDiag
       
        cat /home/fsantaniello/strong_binders_annotated_RelDiag /home/fsantaniello/weak_binders_annotated_RelDiag > $output_dir"/Summary_Binders_RESULTS_"$name_output"__PeptidesPrediction_CommonRelapseDiagnosisOnly.txt"

else

	echo "*************************************************************************************************Haplotype is "$haplotype" *************************************************************************************************"
        #Summary Binders Diagnosis Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\t" $haplotype)
        more +2 $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_DiagnosisOnly.xls"| awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0;if ($8>50&&$8<=500) print $3,$2,0,1}' - | (echo -e $weak_binders; cat -) > $output_dir"/RESULTS_"$name_output"_weak_binders_DiagOnly"
         #/home/fsantaniello/Summary_Epitopi.py $input_dir"/"$name_input"_ForEpitopes_DiagnosisOnly.txt" /home/fsantaniello/weak_binders_DiagOnly $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_DiagnosisOnly.xls" $output_dir"/RESULTS_"$name_output"_weak_binders_annotated_DiagOnly"
        strong_binders=$(echo "Strong binders\nGene\tPeptide\t" $haplotype)
        more +2 $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_DiagnosisOnly.xls" | awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0;if ($8>0&&$8<=50) print $3,$2,0,1}' - | (echo -e $strong_binders; cat -) > $output_dir"/RESULTS_"$name_output"_strong_binders_DiagOnly"
         #/home/fsantaniello/Summary_Epitopi.py $input_dir"/"$name_input"_ForEpitopes_DiagnosisOnly.txt" /home/fsantaniello/strong_binders_DiagOnly $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_DiagnosisOnly.xls" /home/fsantaniello/strong_binders_annotated_DiagOnly

        #cat /home/fsantaniello/strong_binders_annotated_DiagOnly /home/fsantaniello/weak_binders_annotated_DiagOnly > $output_dir"/Summary_Binders_RESULTS_"$name_output"__PeptidesPrediction_DiagnosisOnly.txt"



        #Summary Binders Relapse Only
        weak_binders=$(echo "Weak binders\nGene\tPeptide\t" $haplotype)
        more +2 $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_RelapseOnly.xls"| awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0;if ($8>50&&$8<=500) print $3,$2,0,1}' - | (echo -e $weak_binders; cat -) > $output_dir"/RESULTS_"$name_output"_weak_binders_RelOnly"
         #/home/fsantaniello/Summary_Epitopi.py $input_dir"/"$name_input"_ForEpitopes_RelapseOnly.txt" /home/fsantaniello/weak_binders_RelOnly $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_RelapseOnly.xls" /home/fsantaniello/weak_binders_annotated_RelOnly

        strong_binders=$(echo "Strong binders\nGene\tPeptide\t" $haplotype)
        more +2 $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_RelapseOnly.xls" | awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0;if ($8>0&&$8<=50) print $3,$2,0,1}' - | (echo -e $strong_binders; cat -) > $output_dir"/RESULTS_"$name_output"_strong_binders_RelOnly"
         #/home/fsantaniello/Summary_Epitopi.py $input_dir"/"$name_input"_ForEpitopes_RelapseOnly.txt" /home/fsantaniello/strong_binders_RelOnly $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_RelapseOnly.xls" /home/fsantaniello/strong_binders_annotated_RelOnly
        
        cat /home/fsantaniello/strong_binders_annotated_RelOnly /home/fsantaniello/weak_binders_annotated_RelOnly > $output_dir"/Summary_Binders_RESULTS_"$name_output"__PeptidesPrediction_RelapseOnly.txt"
        
        
        #Summary Binders Common Diagnosis Relapse
        weak_binders=$(echo "Weak binders\nGene\tPeptide\t" $haplotype)
        more +2 $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_CommonRelapseDiagnosis.xls"| awk '{OFS="\t"; if ($5>50&&$5<=500) print $3,$2,1,0;if ($8>50&&$8<=500) print $3,$2,0,1}' - | (echo -e $weak_binders; cat -) > $output_dir"/RESULTS_"$name_output"_weak_binders_RelDiag"
         #/home/fsantaniello/Summary_Epitopi.py $input_dir"/"$name_input"_ForEpitopes_CommonDiagnosisRelapse.txt" /home/fsantaniello/weak_binders_RelDiag $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_CommonRelapseDiagnosis.xls" /home/fsantaniello/weak_binders_annotated_RelDiag
        strong_binders=$(echo "Strong binders\nGene\tPeptide\t" $haplotype)
        more +2 $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_CommonRelapseDiagnosis.xls" | awk '{OFS="\t"; if ($5>0&&$5<=50) print $3,$2,1,0;if ($8>0&&$8<=50) print $3,$2,0,1}' - | (echo -e $strong_binders; cat -) > $output_dir"/RESULTS_"$name_output"_strong_binders_RelDiag"
         #/home/fsantaniello/Summary_Epitopi.py $input_dir"/"$name_input"_ForEpitopes_CommonDiagnosisRelapse.txt" /home/fsantaniello/strong_binders_RelDiag $output_dir"/RESULTS_"$name_output"_PeptidesPrediction_CommonRelapseDiagnosis.xls" /home/fsantaniello/strong_binders_annotated_RelDiag
       
        cat /home/fsantaniello/strong_binders_annotated_RelDiag /home/fsantaniello/weak_binders_annotated_RelDiag > $output_dir"/Summary_Binders_RESULTS_"$name_output"__PeptidesPrediction_CommonRelapseDiagnosisOnly.txt"

 fi
#
#
#
