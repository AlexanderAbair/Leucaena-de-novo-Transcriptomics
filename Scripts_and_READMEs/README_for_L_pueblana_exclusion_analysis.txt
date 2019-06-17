#6/16/2018
#De novo Phylogenomic Inference Using Transcriptome Data from 24 Leucaena Taxa

#Homologies will first be infered from proteomes from Arabidopsis thaliana and Glycine max.
#Proteomes were downloaded from the internet.
#The Glycine max proteome came from http://www.uniprot.org/proteomes/UP000008827
#The proteome was published by Schmutz et al. (2010) "Genome sequence of the palaeopolyploid soybean."
#The Arabidopsis thaliana proteome came from https://www.arabidopsis.org/download_files/Proteins/Araport11_protein_lists/Araport11_genes.201606.pep.fasta.gz
#The proteome was published on TAIR (https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FProteins)[date accessed 6/18/18] as "Araport11 protein lists" 

#The proteomes are stored in the directory /media/Scratch8TB/Aabair/Projects/de_novo_Leucaeana_phylogenomic_inference/1_reference_proteomes
#The transcriptome files for each Leucaena taxon (and 3 outgroups) previously generated using Trinity are stored in /media/Scratch8TB/Aabair/Projects/de_novo_Leucaeana_phylogenomic_inference/2_sequence_processing

#The proteome files were concatenated into one database file

	cat *fa > db

#A blastp database was made from this db file

	makeblastdb -in ./db -parse_seqids -dbtype prot -out db

#Long orfs were pulled from the L_collinsia trinity file

	/Programs/TransDecoder-3.0.1/TransDecoder.LongOrfs -t L_collinsii.Trinity.3seedling.Summer2015.namedContigs.fasta -S 

#Long orfs were blasted against the references

	blastp -query /media/Scratch8TB/Aabair/Projects/de_novo_Leucaeana_phylogenomic_inference/2_sequence_processing/L_collinsii.Trinity.3seedling.Summer2015.namedContigs.fasta.transdecoder_dir/longest_orfs.pep -db /media/Scratch8TB/Aabair/Projects/de_novo_Leucaeana_phylogenomic_inference/1_reference_proteomes/db -max_target_seqs 1 -outfmt 6 -evalue 10 -num_threads 60 > taxonID.blastp.outfmt6

#This process took about two hours with 60 threads.
______________________________________________________________________________________________________________________________________________________________________________________

#6/18/2018

#Applying this process to all taxa
#Since this process worked for L_collinsii, a for loop will now be made to run on all taxa (including outgroups)
#The for loop is called long_orf_extraction.sh
#The output directories were cleaned using the for loop called renaming_directories.sh
#Long orfs from all samples were blasted against the reference db using a for loop called blastp_against_reference_db.sh*
#This blast script started running at 3:00 and is expected to finish in about 24-30 hours.
______________________________________________________________________________________________________________________________________________________________________________________

#6/21/2018
#For the next step, I renamed the trinity files and the transdecoder directories using scripts in the /2_sequence_processing/Scripts directory.
#The orfs with blast hits were retained using the script filtering_for_blast_hits.sh
______________________________________________________________________________________________________________________________________________________________________________________

#6/25/2018

#Renaming sequence ID tags in cds and peptide files
#I started to experiment with the sed command to clean up redundant and unnecesaary bits of information from the sequence tag IDs.
#I want to have the tags in the following format: taxonID@sequenceID
#The first sed command reduced the sequence IDs to: Genus_species_sequence_ID

        cat Taxon.Trinity.fasta.transdecoder.cds | sed -e 's/::L_[^)]*)[^)]*)//' > test.cds ##This did not work.  See the sequence_renaming scripts for details.

#I want to replace the second of three underscores with an @ symbol
#I used the following sed command to replace every second underscore in a new line with an @ symbol:

        cat test.cds | sed -e 's/\_/\@/2' > test.with_@.cds

#I made for loop to clean all sequence ID tags for cds and peptide files called seqeunce_ID_cleanup.cds.sh and seqeunce_ID_cleanup.pep.sh
#The second underscore in the L_macrophylla subspecies had to be changed a hyphen throughout the two files, and the file names had to be changed as well
#I used individual sed commands on the outgroups

#I reduced redundancy in the cds files using the basic command provided by the Ya Yang tutorial in a for loop called reducing_redundancy_in_cds_files.sh
#The basic command is as follows:

	/Programs/TransDecoder-3.0.1/util/bin/cd-hit-est -i ../2_sequence_processing/Cleaned_sequence_tags/$i.Trinity.fasta.transdecoder.cds -o ./Data/$i.fa.cds.cdhitest -c 0.99 -n 10 -r 0 -T 2

________________________________________________________________________________________________________________________________________________________________________________________

#6/28/2018

#Concatenating cdhitest files
#One file containing all Leucaena taxa plus outgroups was made:

	cat *cdhitest > all_taxa_and_outgroups.fa

#Another file just containing the diploid Leucaena and outgroups was mde with the following for loop:

	for i in L_collinsii L_cruziana L_cuspidata L_esculenta L_greggii L_lanceolata L_lempirana L_macrophylla-ist L_macrophylla-mac L_magni L_matudae L_multicapitula L_pulverulenta L_retusa L_salvadorensis L_shannonii L_trichandra L_trichodes L_zacapana Albizia_julibrissin Entada_abyssinica Microlobius_foetidus  

		do cat Data/$i.*.cdhitest >> 1_clustering/diploids_and_outgroups.fa

	done

#A blast database was made for the diploids and outgroups using:

	makeblastdb -in diploids_and_outgroups.fa -parse_seqids -dbtype nucl -out diploids_and_outgroups.fa

#"added 1173453 sequences"
_________________________________________________________________________________________________________________________________________________________________________________________

7/6/2018

#An all-by-all blast of the diploids and outgroups was done using the following command:

	blastn -db diploids_and_outgroups.fa -query diploids_and_outgroups.fa -evalue 10 -num_threads 32 -max_target_seqs 30 -out diploids_and_outgroups.rawblast -outfmt '6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore'

#Filtering raw blast output by hit fraction and prepare input file for mcl

	python ../../../phylogenomic_dataset_construction/blast_to_mcl.py diploids_and_outgroups.rawblast 0.4

#"Output written to diploids_and_outgroups.rawblast.hit-frac0.4.minusLogEvalue
#Highly similar interspecific hits written to diploids_and_outgroups.rawblast.ident
#Check for possible contamination"

#Running mcl:

	mcl diploids_and_outgroups.rawblast.hit-frac0.4.minusLogEvalue --abc -te 32 -tf 'gq(5)' -I 1.4 -o hit-frac0.4_I1.4_e5

#The caption at the end of the mcl command prompt:
#[mcl] jury pruning marks: <99,98,98>, out of 100
#[mcl] jury pruning synopsis: <98.6 or marvelous> (cf -scheme, -do log)
#[mcl] output is in hit-frac0.4_I1.4_e5
#[mcl] 48957 clusters found

#Writing fasta files for each cluster from mcl output that have all 5 taxa.

	mkdir ../2_clusters
	python ../../../phylogenomic_dataset_construction/write_fasta_files_from_mcl.py diploids_and_outgroups.fa hit-frac0.4_I1.4_e5 22 ../2_clusters/

#14381 clusters with at least 22 taxa read were found
_________________________________________________________________________________________________________________________________________________________________________________________

#7/7/2018
#Align each cluster, trim alignment, and infer a tree.

	cd /Programs/Pasta/pasta   #necessary in order for the following command to run properly
	sudo python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/fasta_to_tree.py /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/2_clusters 32 dna n
________________________________________________________________________________________________________________________________________________________________________________________

#7/17/18
#Some of the files from the fasta_to_tree.py process ended up in the /Programs/Pasta/pasta directory
#I moved all the RAxML_bestTree*, RAxML_parsimonyTree*, RAxML_result*, RAxML_log*, and RAxML_info* files into the /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology/2_clusters directory
#The *.tre files do not seem to have been generated.
#I need to look through the log, info, and result files to see were an error may have occured (if it did occur).
#I suspect the problem stems from having to run the fasta_to_tree.py from the /Programs/Pasta/pasta directory.

#UPDATE: naming convention with RAxML made RAxML_* rather than *.tre

________________________________________________________________________________________________________________________________________________________________________________________

#7/30/18

#Moving RAxML_bestTree* files to a new folder and renaming them with a .tre extension
	
	mv /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/2_clusters/RAxML_bestTree* /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/test_trees
	cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/test_trees
	for i in *; do mv "$i" "$i.tre"; done
	for i in *; rename s/RAxML_bestTree.// $i; done
	for i in *; rename s/fa.mafft.aln-cln/raxml/ $i ; done

#Trimming tips that are longer than 0.2 and more than 10 times longer than its sister. Also trim tips that are longer than 0.4.

	python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/trim_tips.py /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/test_trees .tre 0.2 0.4

#Renaming test_trees directory to 3_build_homolog_trees

        mv test_trees 3_build_homolog_trees

_______________________________________________________________________________________________________________________________________________________________________________________

#7/31/18
#Copying necessary files to the 3_build_homolog_trees directory

	for i in * ; cp 2_clusters/cluster*.fa 3_build_homolog_trees ; done
        for i in * ; cp 2_clusters/cluster*.fa.mafft.aln 3_build_homolog_trees ; done
        for i in * ; cp 2_clusters/cluster*.fa.mafft.aln-cln 3_build_homolog_trees ; done

#Masking both mono- and (optional) paraphyletic tips that belong to the same taxon. Keep the tip that has the most un-ambiguous charactors in the trimmed alignment.
	
	python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/mask_tips_by_taxonID_transcripts.py /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/3_build_homolog_trees /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/3_build_homolog_trees y

#Cutting long internal branches 

	mkdir /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/4_homolog
	python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/cut_long_internal_branches.py /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/3_build_homolog_trees .mm 0.3 5 /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/4_homolog

#Writing fasta files from .subtree files

	python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/write_fasta_files_from_trees.py /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/1_clustering/diploids_and_outgroups.fa /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/4_homolog/ .subtree /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/4_homolog

#Splitting up data for faster downstream processing

	cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/4_homolog
	for i in 1 2 3 4 5 6 7 8 9 ; do mkdir $i ; done
	for i in 1 2 3 4 5 6 7 8 9 ; do mv cluster$i* ./$i ; done

#Calculating the final homolog trees and bootstrap from these new fasta files

	#note the $i variable which was replaced with 1-9 in 9 separate screens
	python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/fasta_to_tree.py /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/4_homolog/$i 2 dna y

________________________________________________________________________________________________________________________________________________________________________________________

#8/1/2018
#TESTING
#At 12:45, I copied all of the comlpleted trees to a directory called test_prune

	for i in 2 3 4 5 6 7 8 9 1/0 1/1 1/2 1/3 1/4 1/5 1/6 1/7 1/8 1/9 ; do cp $i/*.tre test_prune/trees ; done

#I made a subdirectory in test_prune to prune for only 1-to-1 homologs.

	mkdir ortho_121
	mkdir ortho_121/tre

#I ran the python script from the test_prune directory to do the actual pruning 

	python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/filter_1to1_orthologs.py trees .tre 22 ortho_121/tre

#Working from the ortho_121/tre directory, I wrote fasta files from the ortholog trees and estimated new alignments.

	python ../../../../../../phylogenomic_dataset_construction/write_ortholog_fasta_files.py ../../../../1_clustering/diploids_and_outgroups.fa . ../fasta_from_ortho/ 22
	cd ../fasta_from_ortho
	mkdir ../new_estimated_alignments
	python ../../../../../../phylogenomic_dataset_construction/write_alignments_from_orthologs.py . ../new_estimated_alignments .fa dna

#Trimming alignments working from the new_estimated_alignments folder

        python ../../../../../../phylogenomic_dataset_construction/phyutility_wrapper.py . 0.3 dna	
________________________________________________________________________________________________________________________________________________________________________________________

#8/2/18
#Fixing species names in sequence ID tags

        mkdir ../fixed_species_names
        for i in cluster*; do sed 's/>c.*/>L_magni/g' $i | sed 's/L_macrophylla/L_macrophylla_ist/g' | sed 's/L_macrophlla/L_macrophylla_mac/g' > ../fixed_species_names/$i; done

#Choose the minimal cleaned alignment length (300 nucleotides) and minimal number of taxa filters for whether to include an ortholog in the supermatrix. Concatenate selected cleaned matrices:

	python ../../../../../../phylogenomic_dataset_construction/concatenate_matrices.py . 300 22 dna 121_filter300-22

#Run raxml with each ortholog as a separate partition.

	raxml -T 2 -p 12345 -m GTRCAT -q 121_filter300-22.model -s 121_filter300-22.phy -n 121_filter300-22

#View the tree using Geneious
#The tree structure exactly matches the reference-guided phylogeny!
________________________________________________________________________________________________________________________________________________________________________________________

#8/3/18
#Continuing to test a preliminary phylogeny using ~175 ortholog groups.

#Running jackknife replicates
#200 jackknife replicates resampling 10% genes each time

	#working from /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homolog_inference/test_prune/ortho_121/
	mkdir jackknife_200_reps_10_percent
	cp new_estimated_alignments/121_filter300-22.model jackknife_200_reps_10_percent/
	cp new_estimated_alignments/121_filter300-22.phy jackknife_200_reps_10_percent/
	cd jackknife_200_reps_10_percent/
	python ../../../../../../phylogenomic_dataset_construction/jackknife_by_percent_genes.py 2 0.1 dna
	cat *result* >../JK10_trees

#Mapping jackknife results to the best tree.
	
	raxml -f b -t new_estimated_alignments/RAxML_bestTree.121_filter300-22 -z JK10_trees -T 2 -m GTRCAT -n 121_filter300-5_JK10
___________________________________________________________________________________________________________________________________

#8/6/18
#Exploring homolog tree discordance in preliminary phylogeny

	cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homolog_inference/test_prune/ortho_121/
	mkdir inclades
	mkdir inclades_phyparts
	python ../../../../../phylogenomic_dataset_construction/extract_clades.py /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/4_homolog/test_prune/trees .tre inclades 19 in_out inclades_phyparts 

________________________________________________________________________________________________________________________________________________________________________________________

#Fixing species names in fasta sequence tags
#For some reason, the sequence tags for L_magni Trinity files were not properly edited early on in the sequence processing stage of the project.  They did not contaain the species name before the @ symbol.  Also the subspecies for L_macrophylla were not edited properly due to a second underscore and a misspelling of "macrophylla" in L_macrophylla_mac.  The following command fixes these errors:
#scripts for sequence ID tag editing can be found in /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/4_homolog 
#Scripts and seq_id.txt files were copied to each of the subdirectories where blocks of homolog trees were produced, and sequence ID tags were all edited simultaneously using separate screen sessions.
________________________________________________________________________________________________________________________________________________________________________________________

#8/14/18
#Moving ahead with full dataset now that the test set was successfully completed
#Since the trees were generated in blocks to save time, all tre files need to be consolidated into one location for further analysis.

	mkdir /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/5_pruning_trees/unpruned_trees
	for i in 2 3 4 5 6 7 8 9 1/0 1/1 1/2 1/3 1/4 1/5_to_9 ; do cp /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/4_homolog/$i/*updated* /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/5_pruning_trees/unpruned_trees ; done

#PARALOGY PRUNING TO INFER ORTHOLOGS

#1-to-1: Only look at homologs that are strictly one-to-one

	cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/5_pruning_trees/
	mkdir ortho_one-to-one
	mkdir ortho_one-to-one/tre
	python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/filter_1to1_orthologs.py /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/5_pruning_trees/unpruned_trees .tre 22 /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/5_pruning_trees/ortho_one-to-one/tre

#Cleaning up concatenated fasta file to be used for writing fasta files from sequence IDs

	sed -i 's/>c/>L_magnifica@c/g' diploids_and_outgroups.fa
	sed -i 's/@/_/2' diploids_and_outgroups.fa
	sed -i 's/>L_macrophlla@mac_/>L_macrophylla_mac@/g' diploids_and_outgroups.fa
	sed -i 's/>L_macrophylla@ist_/>L_macrophylla_ist@/g' diploids_and_outgroups.fa

________________________________________________________________________________________________________________________________________________________________________________________

#8/15/18
#Writing ortholog clusters to fasta files

	python ../../../../../phylogenomic_dataset_construction/write_ortholog_fasta_files.py diploids_and_outgroups.fa . . 22

#Realigning ortholog clusters using prank

	python ../../../../../phylogenomic_dataset_construction/prank_wrapper.py . ../new_alignments/ fa dna

#Trimming alignments in the new_alignments directory using 0.3 for MIN_COLUMN_OCCUPANCY

	python ../../../../../phylogenomic_dataset_construction/phyutility_wrapper.py . 0.3 dna

#Choosing the minimal cleaned alignment length (300 nucleotides) and minimal number of taxa filters for whether to include an ortholog in the supermatrix. Concatenate selected cleaned matices.

	python ../../../../../phylogenomic_dataset_construction/concatenate_matrices.py . 300 22 dna 121_filter300-22

#Estimating a species tree with RAxML

	raxml -T 2 -p 12345 -m GTRCAT -q 121_filter300-22.model -s 121_filter300-22.phy -n 121_filter300-22

#Running 200 jackknife replicates resampling 10% of genes each time

	cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/5_pruning_trees/ortho_one-to-one/
	mkdir jackknife_resampling
	cp new_alignments/121_filter300-22.model jackknife_resampling/
	cp new_alignments/121_filter300-22.phy jackknife_resampling/
	cd jackknife_resampling
	python ../../../../../phylogenomic_dataset_construction/jackknife_by_percent_genes.py 2 0.1 dna
___________________________________________________________________________________________________________________________________________________________________________________________

#8/16/18
#Concatenating the results from jacknife resampling

	cat *result* >../JK10_trees

#Mapping jackknife results to the best tree

	raxml -f b -t new_alignments/RAxML_bestTree.121_filter300-22 -z JK10_trees -T 2 -m GTRCAT -n 121_filter300-22_JK10

#Exploring homolog tree discordance using phyparts
#Prepare directorie and input files

	cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/5_pruning_trees/ortho_one-to-one/
	mkdir inclades inclades_phyparts
	cp ../../4_homolog/test_prune/ortho_121/in_out .

#Extract clades

	python ../../../../phylogenomic_dataset_construction/extract_clades.py ../unpruned_trees .tre inclades 19 in_out inclades_phyparts/

#Determine conflict and concordance

	phyparts -a 1 -d inclades_phyparts -m new_alignments/RAxML_bestTree.121_filter300-22 -o concon -s 50 -v

#Infer duplications

	phyparts -a 2 -d inclades_phyparts -m new_alignments/RAxML_bestTree.121_filter300-22 -o dup -s 50 -v

________________________________________________________________________________________________________________________________

#Experimenting with Astral

#Estimating a tree
	
	#example command
	java -jar astral.5.6.2.jar -i in.tree -o out.tre

#Preparing data for astral

	cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/5_pruning_trees/ortho_one-to-one/tre
	mkdir ../astral
	cat cluster*tre > ../astral/all_ortholog_clusters.tree
	

#Bootstrapped tree

	#example command
	java -jar astral.5.6.2.jar -i best_ml -b bs_paths -r 100

________________________________________________________________________________________________________________________________________________________________________________________

#8/20/18
#Running 200 jackknife replicates resampling 75% of genes each time

        cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/5_pruning_trees/ortho_one-to-one/
        mkdir jackknife_resampling_75
        cp new_alignments/121_filter300-22.model jackknife_resampling_75/
        cp new_alignments/121_filter300-22.phy jackknife_resampling_75/
        cd jackknife_resampling_75
        python ../../../../../phylogenomic_dataset_construction/jackknife_by_percent_genes.py 2 0.75 dna

#Concatenating the results from jacknife resampling

        cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/5_pruning_trees/ortho_one-to-one/jackknife_sampling_75
	cat *result* >../JK75_trees

#Mapping jackknife results to the best tree

        raxml -f b -t new_alignments/RAxML_bestTree.121_filter300-22 -z JK75_trees -T 2 -m GTRCAT -n 121_filter300-22_JK75

_______________________________________________________________________________________________________________________________________________________________________________________

#8/21/18
#Testing RAxML-ng 
#Copying necessary files

	cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/5_pruning_trees/ortho_one-to-one/
	mkdir RAxMLng
	cp new_alignments/121_filter300-22.model RAxMLng/
        cp new_alignments/121_filter300-22.phy RAxMLng/
	cd RAxMLng
	/Programs/RAXML-NG/bin/raxml-ng --msa 121_filter300-22.phy --model 121_filter300-22.model --threads 25

_______________________________________________________________________________________________________________________________________________________________________________________

#8/25/18
#Calculating the best tree and 1000 bootstraps

	raxmlHPC-PTHREADS -f a -s 121_filter300-22.phy.jk99_rep1 -x 12345 -p 12345 -T 32 -# 1000 -o Entada_abyssinica -m GTRCAT -n 121_filter300-22.phy.jk99_rep1.out

_______________________________________________________________________________________________________________________________________________________________________________________

#4/26/2019
#Redoing jackknife support for the RAxML diploid species tree.  1000 pseudorplicates at 63%

	cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/5_pruning_trees/ortho_one-to-one/
        mkdir jackknife_resampling_1000-63
        cp new_alignments/121_filter300-22.model jackknife_resampling_1000-63
        cp new_alignments/121_filter300-22.phy jackknife_resampling_1000-63
        cd jackknife_resampling_1000-63
        python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/jackknife_by_percent_genes.py 10 0.63 dna 

___________________________________________________________________________________________________________________________________________

#4/28/2019
#For some reason, only 200 replicates were made, so I need to make more directories.

	cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/5_pruning_trees/ortho_one-to-one/
	mv jackknife_resampling_1000-63 jackknife_resampling_1000-63.001-200_reps
	mkdir jackknife_resampling_1000-63
	mv jackknife_resampling_1000-63.001-200_reps jackknife_resampling_1000-63
	cd jackknife_resampling_1000-63
	mkdir jackknife_resampling_1000-63.201-400_reps/
	mkdir jackknife_resampling_1000-63.401-600_reps/
	mkdir jackknife_resampling_1000-63.601-800_reps/
	mkdir jackknife_resampling_1000-63.801-1000_reps/

#Copying starting files

	for i in 2 4 6 8
		do cp jackknife_resampling_1000-63.001-200_reps/121_filter300-22.model jackknife_resampling_1000-63.$i*
		done

	for i in 2 4 6 8
		do cp jackknife_resampling_1000-63.001-200_reps/121_filter300-22.phy jackknife_resampling_1000-63.$i*
		done

#Setting up separate screens for each new set of 200 reps

	screen -S jk2
		cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/5_pruning_trees/ortho_one-to-one/jackknife_resampling_1000-63/jackknife_resampling_1000-63.201-400_reps/
		python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/jackknife_by_percent_genes.py 10 0.63 dna

        screen -S jk3
                cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/5_pruning_trees/ortho_one-to-one/jackknife_resampling_1000-63/jackknife_resampling_1000-63.401-600_reps/
                python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/jackknife_by_percent_genes.py 10 0.63 dna

        screen -S jk4
                cd 
		/media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/5_pruning_trees/ortho_one-to-one/jackknife_resampling_1000-63/jackknife_resampling_1000-63.601-800_reps/
                python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/jackknife_by_percent_genes.py 10 0.63 dna

        screen -S jk5
                cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/5_pruning_trees/ortho_one-to-one/jackknife_resampling_1000-63/jackknife_resampling_1000-63.801-1000_reps/
                python /media/Scratch8TB/Aabair/Projects/phylogenomic_dataset_construction/jackknife_by_percent_genes.py 10 0.63 dna	 

#Concatenating the results from jacknife resampling

	cd /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/5_pruning_trees/ortho_one-to-one/jackknife_resampling_1000-63
        cat ./*/*result* >> JK1000-63_trees

#Mapping jackknife results to the best tree

        raxml -f b -t /media/Scratch8TB/Aabair/Projects/de_novo_Leucaena_phylogenomic_inference/3_homology_inference/5_pruning_trees/ortho_one-to-one/new_alignments/Phylip_format/RAxML/RAxML_bestTree.121_filter300-22 -z JK1000-63_trees -T 2 -m GTRCAT -n Leucaena_diploid_excluding_publana_species_tree.jackknife_1000-63.tre
