seq=$( cat seq_ids.txt )
cluster=$( cat cluster_numbers.txt )

for f in $cluster
	do cat $f > $f.updated_taxa.tre
	done

for i in $seq
	do sed -i "s/$i/L_magnifica@$i#/g" *.updated_taxa.tre
	done

sed -i 's/@#/_/g' *.updated_taxa.tre

sed -i 's/L_macrophlla/L_macrophylla_mac/g' *.updated_taxa.tre
sed -i 's/L_macrophylla@/L_macrophylla_ist@/g' *.updated_taxa.tre
sed -i 's/mac_c/c/g' *.updated_taxa.tre
sed -i 's/ist_c/c/g' *.updated_taxa.tre

