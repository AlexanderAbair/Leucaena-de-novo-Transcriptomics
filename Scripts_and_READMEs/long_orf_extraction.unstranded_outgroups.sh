#6/18/2018
#This for loop extracts long orfs from all transcriptome files (including outgroups)

for i in Albizia_julibrissin Entada_abyssinica Microlobius_foetidus
	do
		/Programs/TransDecoder-3.0.1/TransDecoder.LongOrfs -t $i*fasta
	done

