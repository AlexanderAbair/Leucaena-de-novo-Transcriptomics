#   Preparing pueblana for incorporation into transcriptome analysis
#
#   2/13/2019
#
#   I transfered the assembled pueblana DNA to the current working directory

	cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/pueblana_sequence_data/assembled_pueblana_DNA
	rsync -arvz -e 'ssh -p 222' abair@128.123.95.198:/data/projects/Lpueblana_genome01_10_16/graph_prefix.scafSeq .

#   Extracting protein sequences from assembled pueblana genomic DNA
#
#   Using the emboss package, I used the getorf function to find orfs in the assembled scaffolds of pueblana and translate them to protein sequences.
#   I specified a minimum orf size of 100 amino acids (300 nucleotides) for the second run.

	getorf -minsize 300

	Input sequence(s): graph_prefix.scafSeq
	Output sequence [scaffold1.orf]: pueblana_minsize_300.pep

#   I also extracted long orfs and output them as cds

	 getorf -minsize 300 -find 2

        Input sequence(s): graph_prefix.scafSeq
        Output sequence [scaffold1.orf]: pueblana_minsize_300.cds


#   2/16/2019
#   Sequence ID tags were cleaned.

	sed 's/ /_/g' pueblana_minsize_300.pep > test.pep
	sed 's/-/DASH/g' test.pep > test.2.pep
	sed 's/[[]//g' test.2.pep > test.3.pep
	sed 's/[]]//g' test.3.pep > test.4.pep
	sed 's/[.]//g' test.4.pep > test.5.pep
	sed 's/>s/>L_pueblana@s/g' test.5.pep > test.6.pep
	sed 's/[(]//g' test.6.pep > test.7.pep
        sed 's/[)]//g' test.7.pep > test.8.pep
        sed 's/>C/>L_pueblana@s/g' test.8.pep > test.9.pep
	mv test.9.pep pueblana_minsize_300.cleaned_seq_ID_tags.pep
	rm test.*

#   The same cleaning process was done with the cds sequences

        sed 's/ /_/g' pueblana_minsize_300.cds > test.cds
        sed 's/-/DASH/g' test.cds > test.2.cds
        sed 's/[[]//g' test.2.cds > test.3.cds
        sed 's/[]]//g' test.3.cds > test.4.cds
        sed 's/[.]//g' test.4.cds > test.5.cds
        sed 's/>s/>L_pueblana@s/g' test.5.cds > test.6.cds
        sed 's/[(]//g' test.6.cds > test.7.cds
        sed 's/[)]//g' test.7.cds > test.8.cds
        sed 's/>C/>L_pueblana@s/g' test.8.cds > test.9.cds
        mv test.9.cds pueblana_minsize_300.cleaned_seq_ID_tags.cds
	rm test.*

#   3/12/19
#   Blasting against the reference proteome database

	blastp -query /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/1_pueblana_sequence_data/assembled_pueblana_DNA/pueblana_minsize_300.cleaned_seq_ID_tags.pep -db /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/1_reference_proteomes/db -max_target_seqs 1 -outfmt 6 -evalue 10 -num_threads 60 > L_pueblana.blastp.outfmt6

#   cds, peptide, and blastp files were coppied to their respective folders for further processing

	mkdir /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/2_sequence_processing
	mkdir /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/2_sequence_processing/cds_and_pep_sequences
	mkdir /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/2_sequence_processing/blastp_results

	cp pueblana_minsize_300.cleaned_seq_ID_tags.cds ../../2_sequence_processing/cds_and_pep_sequences/L_pueblana.long_orfs.cds
        cp pueblana_minsize_300.cleaned_seq_ID_tags.pep ../../2_sequence_processing/cds_and_pep_sequences/L_pueblana.long_orfs.pep
	cp L_pueblana.blastp.outfmt6 ../../2_sequence_processing/blastp_results/

#   Blastp output, cds, and pep files were all brought over to the 2_sequence_processing subfolders from the failed pueblana exclusion experiment.
#   See the README files in their respective folders for information about how sequence IDs were cleaned up.


#   I reduced redundancy in the cds files using th following for loop

        for i in Albizia_julibrissin Entada_abyssinica L_collinsii L_confertiflora_adenotheloidea L_cruziana L_cuspidata L_diversifolia L_esculenta L_greggii L_involucrata L_lanceolata L_lempirana L_leucocephala_glabrata L_leucocephala_leucocephala L_macrophylla_istmensis L_macrophylla_macrophylla L_magnifica L_matudae L_multicapitula L_pallida L_pueblana L_pulverulenta L_retusa L_salvadorensis L_shannonii L_trichandra L_trichodes L_zacapana Microlobius_foetidus

                do /Programs/TransDecoder-3.0.1/util/bin/cd-hit-est -i /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/2_sequence_processing/cds_and_pep_sequences/$i.long_orfs.cds -o /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/2_sequence_processing/cds_and_pep_sequences/$i.cds.cdhitest -c 0.99 -n 10 -r 0 -T 2
                done

#   Concatenating cdhitest files

        mkdir /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/
        mkdir /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/1_clustering

        touch 1_clustering/diploids_and_outgroups.fa
        touch 1_clustering/diploids_and_outgroups.without_pueblana.fa

        for i in Albizia_julibrissin Entada_abyssinica L_collinsii L_cruziana L_cuspidata L_esculenta L_greggii L_lanceolata L_lempirana L_macrophylla_istmensis L_macrophylla_macrophylla L_magnifica L_matudae L_multicapitula L_pueblana L_pulverulenta L_retusa L_salvadorensis L_shannonii L_trichandra L_trichodes L_zacapana Microlobius_foetidus

                do cat $i.*.cdhitest >> ../../3_homology/1_clustering/diploids_and_outgroups.fa
                done

        for i in Albizia_julibrissin Entada_abyssinica L_collinsii L_cruziana L_cuspidata L_esculenta L_greggii L_lanceolata L_lempirana L_macrophylla_istmensis L_macrophylla_macrophylla L_magnifica L_matudae L_multicapitula L_pulverulenta L_retusa L_salvadorensis L_shannonii L_trichandra L_trichodes L_zacapana Microlobius_foetidus

                do cat $i.*.cdhitest >> ../../3_homology/1_clustering/diploids_and_outgroups.without_pueblana.fa
                done

        cat *.cdhitest >> ../../3_homology/1_clustering/all_taxa_and_outgroups.fa

#   Making a blast database for all diploids and the outgroups

        makeblastdb -in diploids_and_outgroups.fa -parse_seqids -dbtype nucl -out diploids_and_outgroups.fa

#   The completed program displayed the following:

	[...]added 1501607 sequences in 88.4408 seconds.

#   Making a blast database for all diploids (except pueblana) and the outgroups

	makeblastdb -in diploids_and_outgroups.without_pueblana.fa -parse_seqids -dbtype nucl -out diploids_and_outgroups.without_pueblana.fa

#   The completed program displayed the following:

	[...]added 1175709 sequences in 70.5755 seconds.

#   An all-by-all blast of the diploids and outgroups was done using the following command:

        blastn -db diploids_and_outgroups.fa -query diploids_and_outgroups.fa -evalue 10 -num_threads 60 -max_target_seqs 30 -out diploids_and_outgroups.rawblast -outfmt '6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore'

#   The blastn results from the diploids_and_outgroups.without_pueblana.fa file were copied over from the previous experiment, since they will not be effected by the previous errors due to mislabeled contigs in the pueblana file

	cp /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion_FAILED/3_homology/1_clustering/diploids_and_outgroups.without_pueblana.rawblast* /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/1_clustering/
	cp /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion_FAILED/3_homology/1_clustering/diploids_and_outgroups.without_pueblana.hit-frac0.4_I1.4_e5 /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/1_clustering/

#   Filtering raw blast output from all diploids and outgroups by hit fraction and prepare input file for mcl

        python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/blast_to_mcl.py /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/1_clustering/diploids_and_outgroups.rawblast 0.4

#   Running mcl on all diploids and outgroups:

        mcl diploids_and_outgroups.rawblast.hit-frac0.4.minusLogEvalue --abc -te 32 -tf 'gq(5)' -I 1.4 -o hit-frac0.4_I1.4_e5

#   The caption at the end of the mcl command prompt:

	[mcl] jury pruning marks: <99,94,95>, out of 100
	[mcl] jury pruning synopsis: <97.2 or superb> (cf -scheme, -do log)
	[mcl] output is in hit-frac0.4_I1.4_e5
	[mcl] 80701 clusters found
	[mcl] output is in hit-frac0.4_I1.4_e5

#   Renaming the hit fraction file

	mv hit-frac0.4_I1.4_e5 diploids_and_outgroups.hit-frac0.4_I1.4_e5

#   Writing fasta files for each cluster from mcl output that have all diploid and outgroup taxa.

        mkdir /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/2_clusters
        mkdir /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/2_clusters/diploids_and_outgroups

        python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/write_fasta_files_from_mcl.py /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/1_clustering/diploids_and_outgroups.fa /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/1_clustering/diploids_and_outgroups.hit-frac0.4_I1.4_e5 23 /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/2_clusters/diploids_and_outgroups/

#   The program printed the following message upon completion

	11555 clusters with at least 23 taxa read

#   Splitting up clusters for parallelizing RaxML

        cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/2_clusters/diploids_and_outgroups

        for i in {00..12} ; do mkdir $i.clusters ; done

        mv cluster?.fa ./00.clusters/
        mv cluster??.fa ./00.clusters/
        mv cluster???.fa ./00.clusters/
 
        for i in {1..9} ; do mv cluster$i???.fa 0$i.clusters ; done
        for i in 10 11 ; do mv cluster$i???.fa $i.clusters ; done

#   Align each cluster, trim alignment, and infer a tree.
        
	for i in 00 01 
		do sudo python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/fasta_to_tree.py /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/2_clusters/diploids_and_outgroups/$i.clusters 12 dna n
                done

	for i in 02 03 
		do sudo python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/fasta_to_tree.py /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/2_clusters/diploids_and_outgroups/$i.clusters 12 dna n
                done

	for i in 04 05 
		do sudo python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/fasta_to_tree.py /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/2_clusters/diploids_and_outgroups/$i.clusters 10 dna n
                done

	for i in 06 07 
		do sudo python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/fasta_to_tree.py /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/2_clusters/diploids_and_outgroups/$i.clusters 10 dna n
                done

	for i in 08 09 
		do sudo python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/fasta_to_tree.py /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/2_clusters/diploids_and_outgroups/$i.clusters 10 dna n
                done

	for i in 10 11 
		do sudo python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/fasta_to_tree.py /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/2_clusters/diploids_and_outgroups/$i.clusters 10 dna n
                done

#   3/14/2019
#   Trimming tips that are longer than 0.2 and more than 10 times longer than its sister. Also trim tips that are longer than 0.4

	for i in {00..11}
                do python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/trim_tips.py /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/2_clusters/diploids_and_outgroups/$i.clusters .tre 0.2 0.4
                done

#   Setting up directories for building homolog groups

	mkdir /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/3_build_homolog_trees 
        cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/3_build_homolog_trees

        mkdir diploids_and_outgroups
        cd diploids_and_outgroups

        for i in {00..11}
                do mkdir $i.clusters
                done

#   Copying the necessary files

        for i in {00..11} 
                do cp /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/2_clusters/diploids_and_outgroups/$i.clusters/* /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/3_build_homolog_trees/diploids_and_outgroups/$i.clusters
		done

#   Masking both mono- and (optional) paraphyletic tips that belong to the same taxon. Keep the tip that has the most un-ambiguous charactors in the trimmed alignment.

        for i in {00..11}
                do python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/mask_tips_by_taxonID_transcripts.py /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/3_build_homolog_trees/diploids_and_outgroups/$i.clusters /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/3_build_homolog_trees/diploids_and_outgroups/$i.clusters y
                done

#   Cutting long internal branches

        mkdir /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/4_homologs
	mkdir /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/4_homologs/diploids_and_outgroups

        for i in {00..11}
                do sudo python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/cut_long_internal_branches.py /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/3_build_homolog_trees/diploids_and_outgroups/$i.clusters .mm 0.3 20 /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/4_homologs/diploids_and_outgroups
                done

#   Writing fasta files from .subtree files

        python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/write_fasta_files_from_trees.py /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/1_clustering/diploids_and_outgroups.fa /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/4_homologs/diploids_and_outgroups .subtree /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/4_homologs/diploids_and_outgroups

#   Setting up directories and moving files for parallelized computing

        cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/4_homologs/diploids_and_outgroups

        for i in {00..12}
                do mkdir $i.clusters
                done

        mv cluster?_* 00.clusters
        mv cluster??_* 00.clusters
        mv cluster???_* 00.clusters

        for i in {1..9}
                do mv cluster$i???_* ./0$i.clusters
                done

        for i in {0..2}
                do mv cluster1$i???_* ./1$i.clusters
                done

#   Calculating the final homolog trees and bootstrap from these new fasta files
#   This script was run simultaneously in separate screen sessions with varying core numbers (the "$i" represents the different number prefixes on the cluster directories).  The lower number folders contain larger clusters, so were given more cores.

        sudo python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/fasta_to_tree.py /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/4_homologs/diploids_and_outgroups/$i.clusters [core#] dna y

____________________________________________________________________________________________________________________________________________

4/4/2019

#   Setting up directories for orthology assessment and tree estimation.

	cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/
	mkdir 5_pruning
	mkdir 5_pruning/complete_diploid_dataset
	cd 5_pruning/complete_diploid_dataset
	mkdir 1_tre_files 2_filtered_trees 3_fasta_from_filtered_trees 4_new_alignments 5_RAxML_species_tree 6_astral

#   Copying homolog trees for orthology assessment

	cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/4_homologs/diploids_and_outgroups
	for i in {00..11} ; do cp $i.clusters/*tre ../../5_pruning_trees/complete_diploid_dataset/1_tre_files/; done

#   Filtering for strictly one-to-one orthologs 

	cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/
        python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/filter_1to1_orthologs.py 1_tre_files/ .tre 23 2_filtered_trees

#   Writing fasta files from the ortholog trees.

        python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/write_ortholog_fasta_files.py /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/1_clustering/diploids_and_outgroups.fa /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/2_filtered_trees /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/3_fasta_from_filtered_trees 23

#   Generating new alignments.

        python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/prank_wrapper.py /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/3_fasta_from_filtered_trees /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/4_new_alignments .fa dna

#   Trimming alignments using 0.3 for MIN_COLUMN_OCCUPANCY

        python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/phyutility_wrapper.py /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/4_new_alignments 0.3 dna

#   Choosing the minimal cleaned alignment length (300 nucleotides) and minimal number of taxa filters for whether to include an ortholog in the supermatrix. Concatenating selected cleaned matrices.

        python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/concatenate_matrices.py /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/4_new_alignments 300 23 dna one-to-one_filter300-23

#   Estimating a species tree with RAxML

        raxml -T 2 -p 12345 -m GTRCAT -q one_to_one_filter300-23.model -s one_to_one_filter300-23.phy -n _filter300-23

#   Running 1000 jackknife replicates resampling 63% genes each time

        cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/5_RAxML_species_tree
	python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/jackknife_by_percent_genes.py 10 0.63 dna
        cat *result* >../1000_reps_at_63_percent

#   Mapping jackknife results to the best tree.
        
        raxml -f b -t /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/5_RAxML_species_tree/RAxML_bestTree.one-to-one_filter300-23 -z 1000_reps_at_63_percent -T 2 -m GTRCAT Leucaena_diploids.1000-63.tre

_____________________________________________________________________________________________________________________________________

#   4/5/2019
#   Astral Species Tree Estimation

#   Setting up directories for generating gene trees

	mkdir /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/6_astral/gene_trees_for_astral

	mkdir /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/6_astral/gene_trees_for_astral/alignments

#   Copying alignments for generating gene trees

	cp /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/4_new_alignments/*aln-cln /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/6_astral/gene_trees_for_astral/alignments

#   Running RAxML to generate gene trees

	python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/raxml_wrapper.py /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/6_astral/gene_trees_for_astral/alignments 15 dna

#   Preparing gene trees for Astral

	mv /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/6_astral/gene_trees_for_astral/alignments/*tre /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/6_astral/gene_trees_for_astral/

	cat /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/6_astral/gene_trees_for_astral/*tre >> /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/6_astral/gene_trees_for_astral/pueblana_inclusion.190_gene_trees.tre

#   Running Astral

	java -jar /Programs/Astral/astral.5.6.2.jar -i /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/6_astral/gene_trees_for_astral/pueblana_inclusion.190_gene_trees.tre -o /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/6_astral/astral.Leucaena_diploid_species.tre

__________________________________________________________________________________________________________________________________________________________

#   Jackknife support for RAxML species tree

#   Preparing for jackknifing on the RAxML species tree

	cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/5_RAxML_species_tree

	mkdir jackknife_1000-63
	for i in 001-200 201-400 401-600 601-800 801-1000 ; do mkdir jackknife_1000-63/reps.$i ; done

	cp ./*jk* jackknife_1000-63/reps.001-200
	cp jacknife* jackknife_1000-63/reps.001-200

	for i in 001-200 201-400 401-600 601-800 801-1000 ; do cp one-to-one_filter300-23.* jackknife_1000-63/reps.$i ; done

#   Setting up separate screens for each new set of 200 reps

        screen -S jk2
                cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/5_RAxML_species_tree/jackknife_1000-63/reps.201-400
                python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/jackknife_by_percent_genes.py 10 0.63 dna

        screen -S jk3
                cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/5_RAxML_species_tree/jackknife_1000-63/reps.401-600
                python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/jackknife_by_percent_genes.py 10 0.63 dna

        screen -S jk4
                cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/5_RAxML_species_tree/jackknife_1000-63/reps.601-800
                python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/jackknife_by_percent_genes.py 10 0.63 dna

        screen -S jk5
                cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/5_RAxML_species_tree/jackknife_1000-63/reps.801-1000
                python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/jackknife_by_percent_genes.py 10 0.63 dna

#   Concatenating the results from jacknife resampling

        cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/5_RAxML_species_tree/jackknife_1000-63/
        cat ./*/*result* >> JK1000-63_trees

#   Mapping jackknife results to the best tree

        raxml -f b -t /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/5_RAxML_species_tree/RAxML_bestTree.one-to-one_filter300-23 -z JK1000-63_trees -T 2 -m GTRCAT -n Leucaena_diploid.jackknife_1000-63.tre

_____________________________________________________________________________________________________________________________________________________________________

#   Jackknife support for Astral tree

#   Setting up directories

	cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/pueblana_inclusion/3_homology/5_pruning_trees/complete_diploid_dataset/6_astral/
	mkdir resampling
	cd resampling 
	for i in {0001..1000}
	        do mkdir $i.rep
	        done

#   Sampling

	for i in {0001..1000}
	        do ls ../gene_trees_for_astral | grep cluster | sort -R | head -120 > $i.rep/rep.$i.list_of_resampled_trees.txt
	           samples=$( cat $i.rep/rep.$i.list_of_resampled_trees.txt )
                   for a in $samples
                        do cp ../gene_trees_for_astral/$a $i.rep
                        done
	        done

#   Resample concatenation

        for i in {0001..1000}
                do cat $i.rep/*tre >> $i.rep/concatenated_trees_for_astral.rep_$i.tre
                done

#   Running Astral

        for i in {0001..1000}
                do java -jar /Programs/Astral/astral.5.6.2.jar -i $i.rep/concatenated_trees_for_astral.rep_$i.tre -o $i.rep/jacknife_pseudoreplicate_$i.tre
                done

#   Setting up directories, and copying rerooted trees

        cd sumtrees/
        mkdir 1000_astral_tree_reps
        for i in {0001..1000};  do nw_reroot ../$i.rep/jackknife* Entada_abyssinica > 1000_astral_tree_reps/jackknife_pseudoreplicate_$i.rerooted.tre ; done

#   Mapping results to the species and calculating support values

	/Programs/DendroPy-4.4.0/applications/sumtrees/sumtrees.py --target=../../astral.Leucaena_diploid_species.tre --decimals=0 --percentages --output=Leucaena_diploids_without_pueblana.jackknife1000-63.tre 1000_astral_tree_reps/*
