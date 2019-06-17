for i in {0001..1000}
	do mkdir $i.rep
	done

for i in {0001..1000}
	do ls ../gene_trees_for_astral | grep cluster | sort -R | head -120 > $i.rep/rep.$i.list_of_resampled_trees.txt
	   samples=$( cat $i.rep/rep.$i.list_of_resampled_trees.txt )
		for a in $samples
			do cp ../gene_trees_for_astral/$a $i.rep
			done
	done
