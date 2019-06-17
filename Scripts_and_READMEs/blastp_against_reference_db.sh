#6/18/2018
#This for loop blasts each taxon against the reference database created from Aribidopsis thaliana and Glycine max proteomes.

	for i in L_collinsii L_confertiflora-adeno L_cruziana L_cuspidata L_diversifolia L_esculenta L_greggii L_involucrata L_lanceolata L_lempirana L_leuco-glabrata L_leuco-leuco L_macrophylla_ist L_macrophylla_mac L_magni L_matudae L_multicapitula L_pallida L_pulverulenta L_retusa L_salvadorensis L_shannonii L_trichandra L_trichodes L_zacapana Albizia_julibrissin Entada_abyssinica Microlobius_foetidus
		do blastp -query /media/Scratch8TB/Aabair/Projects/de_novo_Leucaeana_phylogenomic_inference/2_sequence_processing/$i.transdecoder_dir/longest_orfs.pep -db /media/Scratch8TB/Aabair/Projects/de_novo_Leucaeana_phylogenomic_inference/1_reference_proteomes/db -max_target_seqs 1 -outfmt 6 -evalue 10 -num_threads 60 > $i.blastp.outfmt6
		done

