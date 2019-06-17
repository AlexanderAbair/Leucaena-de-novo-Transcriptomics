#6/25/2018
#Reducing redundancy in the cds files with cleaned sequence ID tags only on outgroups from unstranded data.

	for i in Albizia_julibrissin Entada_abyssinica Microlobius_foetidus

		do /Programs/TransDecoder-3.0.1/util/bin/cd-hit-est -i ../2_sequence_processing/Cleaned_sequence_tags/$i.Trinity.fasta.transdecoder.cleaned.cds -o ./Data/$i.fa.cds.cdhitest -c 0.99 -n 10 -r 1 -T 2

		done
