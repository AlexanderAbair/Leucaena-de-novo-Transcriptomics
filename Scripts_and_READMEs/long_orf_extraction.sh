#6/18/2018
#This for loop extracts long orfs from all transcriptome files (including outgroups)

for i in Albizia_julibrissin Entada_abyssinica L_collinsii L_confertiflora-adeno L_cruziana L_cuspidata L_diversifolia L_esculenta L_greggii L_involucrata L_lanceolata L_lempirana L_leuco-glabrata L_leuco-leuco L_macrophylla_ist L_macrophylla_mac L_magni L_matudae L_multicapitula L_pallida L_pulverulenta L_retusa L_salvadorensis L_shannonii L_trichandra L_trichodes L_zacapana Microlobius_foetidus
	do
		/Programs/TransDecoder-3.0.1/TransDecoder.LongOrfs -t $i*fasta -S
	done

