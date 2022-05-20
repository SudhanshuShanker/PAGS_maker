counter=0
for i in $1/*.inp
do
        counter=$(($counter+1))
        echo $counter
        outdir=$(echo ${i}_output | sed 's_/_\\/_g')
	mkdir ${i}_output 
        cp ./tools/ex_job $1/ex_job_${counter}.job
        #echo "-s $i" >> flag_gly${counter}
        #echo "-out:path:pdb ./output_glycomutagenesis_${counter}" >> flag_gly_${counter} 
        #echo "-out:path:score ./output_glycomutagenesis_${counter}" >> flag_gly_${counter}
        #echo " " >> flag_gly_${counter}
	ii=$(echo ${i} | sed 's_/_\\/_g')
	sed -i "s/REPLACE/$ii/g" $1/ex_job_${counter}.job
        sed -i "s/output/$outdir/g" $1/ex_job_${counter}.job
        #sed -i "s/flag_gly/flag_gly_${counter}/g" job_gly_${counter}


done

