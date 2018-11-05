#!/home/fsantaniello/.ctgb-sistemistica/bin/python
import sys
Mutazioni=open(str(sys.argv[1])).readlines()
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
				binders.append(str(j.split('\t')[8]).rstrip()+'\t'+str(j.split('\t')[9]).rstrip()+'\t'+str(i.split('\t')[0])+'\t'+str(i.split('\t')[1])+'\t'+str(i.split('\t')[2]).rstrip()+'\t'+str(j.split('\t')[6]).rstrip()+'\t'+str(j.split('\t')[0]).rstrip()+'\t'+str(j.split('\t')[11]).rstrip()+'\t'+str(j.split('\t')[10]).rstrip()+'\t'+str(j.split('\t')[17]).rstrip()+'\t'+str(j.split('\t')[18]).rstrip())
	for i in binders:
		for j in Excel_mut[2:]:
			if str(i.split('\t')[3])==str(j.split('\t')[1]):
				out.write(str(''.join(i)+'\t'+str(j.split('\t')[4])+'\n'))

	out.close()


