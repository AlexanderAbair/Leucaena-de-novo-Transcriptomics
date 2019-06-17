#6/21/2018
#This script retains only the orfs from each taxon with blast hits.

	for i in L_collinsii L_confertiflora-adeno L_cruziana L_cuspidata L_diversifolia L_esculenta L_greggii L_involucrata L_lanceolata L_lempirana L_leuco-glabrata L_leuco-leuco L_macrophylla_ist L_macrophylla_mac L_magni L_matudae L_multicapitula L_pallida L_pulverulenta L_retusa L_salvadorensis L_shannonii L_trichandra L_trichodes L_zacapana Albizia_julibrissin Entada_abyssinica Microlobius_foetidus 

		do /Programs/TransDecoder-3.0.1/TransDecoder.Predict -t $i.Trinity.fasta --retain_blastp_hits $i.blastp.outfmt6 --cpu 60

	done
