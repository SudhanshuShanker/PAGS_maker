for i in $1/*.job
do
	condor_submit $i
done
