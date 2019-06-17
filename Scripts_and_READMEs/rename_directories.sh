#6/18/2018
#Cleaning up directory names

	for i in L_collinsii L_confertiflora-adeno L_cruziana L_cuspidata L_diversifolia L_esculenta L_greggii L_involucrata L_lanceolata L_lempirana L_leuco-glabrata L_leuco-leuco L_macrophylla_ist L_macrophylla_mac L_magni L_matudae L_multicapitula L_pallida L_pulverulenta L_retusa L_salvadorensis L_shannonii L_trichandra L_trichodes L_zacapana Albizia_julibrissin Entada_abyssinica Microlobius_foetidus

		do mv ./$i*.namedContigs.fasta.transdecoder_dir ./$i.transdecoder_dir

	done
