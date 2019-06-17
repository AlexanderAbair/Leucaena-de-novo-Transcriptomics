#6/25/2018
#Reducing redundancy in the Leucaena cds files with cleaned sequence ID tags. Outgroups are in a different script.

	for i in L_collinsii L_confertiflora-adeno L_cruziana L_cuspidata L_diversifolia L_esculenta L_greggii L_involucrata L_lanceolata L_lempirana L_leuco-glabrata L_leuco-leuco L_macrophylla-ist L_macrophylla-mac L_magni L_matudae L_multicapitula L_pallida L_pulverulenta L_retusa L_salvadorensis L_shannonii L_trichandra L_trichodes L_zacapana

		do /Programs/TransDecoder-3.0.1/util/bin/cd-hit-est -i ../2_sequence_processing/Cleaned_sequence_tags/$i.Trinity.fasta.transdecoder.cleaned.cds -o ./Data/$i.fa.cds.cdhitest -c 0.99 -n 10 -r 0 -T 2

		done
