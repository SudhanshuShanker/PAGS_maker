
for i in $1/*.pdb
do
	cat $i | grep -E "ATOM|HETATM" | awk '{print $12,"  ",$7,"  ",$8,"  ",$9}' > ${i::-3}xyz

done

