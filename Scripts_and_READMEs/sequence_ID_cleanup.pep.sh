#6/25/2018
#This script cleans sequence ID tags in the pep files

	for i in L_collinsii L_confertiflora-adeno L_cruziana L_cuspidata L_diversifolia L_esculenta L_greggii L_involucrata L_lanceolata L_lempirana L_leuco-glabrata L_leuco-leuco L_macrophylla-ist L_macrophylla-mac L_magni L_matudae L_multicapitula L_pallida L_pulverulenta L_retusa L_salvadorensis L_shannonii L_trichandra L_trichodes L_zacapana Albizia_julibrissin Entada_abyssinica Microlobius_foetidus

		do cat /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/2_sequence_processing/$i.Trinity.fasta.transdecoder.pep | sed -e 's/>[^)]*) />/' > $i.temp.pep
	           cat $i.temp.pep | sed 's/(-)/_neg/' > $i.temp.negs.pep
		   cat $i.temp.negs.pep | sed 's/:/_/' > $i.character_swap_1.pep
		   cat $i.character_swap_1.pep | sed 's/-/_/' > $i.character_swap_2.pep
		   cat $i.character_swap_2.pep | sed -e 's/\_/\@/2' > $i.arroba.pep
		   cat $i.arroba.pep | sed 's/(+)/_pos/' > /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/2_sequence_processing/Cleaned_sequence_tags/$i.Trinity.fasta.transdecoder.cleaned.pep
		   rm $i.arroba.pep
		   rm $i.character_swap_1.pep
		   rm $i.character_swap_2.pep
		   rm $i.temp.pep
		   rm $i.temp.negs.pep
		done

