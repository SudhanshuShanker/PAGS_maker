for i in $1/*.xyz
do
	cat $2 $i > ${i::-3}inp
	echo "*" >> ${i::-3}inp

done
