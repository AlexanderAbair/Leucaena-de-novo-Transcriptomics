#6/25/2018
#This script cleans sequence ID tags in the cds files

	for i in L_collinsii L_confertiflora-adeno L_cruziana L_cuspidata L_diversifolia L_esculenta L_greggii L_involucrata L_lanceolata L_lempirana L_leuco-glabrata L_leuco-leuco L_macrophylla-ist L_macrophylla-mac L_magni L_matudae L_multicapitula L_pallida L_pulverulenta L_retusa L_salvadorensis L_shannonii L_trichandra L_trichodes L_zacapana Albizia_julibrissin Entada_abyssinica Microlobius_foetidus

		do cat /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/2_sequence_processing/$i.Trinity.fasta.transdecoder.cds | sed -e 's/>[^)]*) />/' > $i.temp.cds
	           cat $i.temp.cds | sed 's/(-)/_neg/' > $i.temp.negs.cds
		   cat $i.temp.negs.cds | sed 's/:/_/' > $i.character_swap_1.cds
		   cat $i.character_swap_1.cds | sed 's/-/_/' > $i.character_swap_2.cds
		   cat $i.character_swap_2.cds | sed -e 's/\_/\@/2' > $i.arroba.cds
		   cat $i.arroba.cds | sed 's/(+)/_pos/' > /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/2_sequence_processing/Cleaned_sequence_tags/$i.Trinity.fasta.transdecoder.cleaned.cds
		   rm $i.arroba.cds
		   rm $i.character_swap_1.cds
		   rm $i.character_swap_2.cds
		   rm $i.temp.cds
		   rm $i.temp.negs.cds
		done

