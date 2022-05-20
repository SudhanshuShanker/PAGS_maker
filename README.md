Steps To Run PAGS-Maker Pipeline
===============================

**Requirements**: AMBERTOOLS, ORCA, Python3

*After each step, analyze output files carefully*  

Edit the *settings.py* file for ORCA and AMBERHOME and other IO settings.  
$~$

**Step1: Make a simplified PDB file for a linkage**
*  Make a list file with linkage names (like ./bin/model_lists/list_example) 
	   [follow the extra/linkage_grammer_file.dat]  
*  Modify "for molecular maker" block in settings.py for the list path (model_name_list) and  pdb_out_dir if rquired.  
*  Run molecule_maker.py (./molecule_maker.py)  
``` sh
>> ./molecule_maker.py
```
This step will generate a simplified pdb molecule in global_settings.simplified_pdb_output_dir (settings.py)  
	
$~$

**Step2: Generate a rotational ensemble and QM input files**  
*  Modify settings.py in "for ensemble making and input file generation" block.  
*  run make_ensemble_and_qm_inputs.py  
``` sh
>> ./make_ensemble_and_qm_inputs.py
```
Use input files generated in *global_settings.output_ensemble_dir/project_****/QM_inps* for QM calculation (ORCA.)
I used condor for the calculation (see bin/tools/QM_tools file for all raw files).  

You can use your schedular and make directories like *global_settings.QM_out_dir_ends_with* (settings.py), 
and output_files named as *global_settings.QM_out_dir_file_name* (settings.py). It looks like:
	
	ensemble_0_0.inp                       <-- Input INP file for ORCA
	ensemble_0_0.inp_output/               <-- Output Directory
	ensemble_0_0.inp_output/condor.out_0   <-- Score file in Output Directory
	ensemble_0_45.inp                       <-- Input INP file for ORCA
	ensemble_0_45.inp_output/               <-- Output Directory
	ensemble_0_45.inp_output/condor.out_0   <-- Score file in Output Directory
	.
	.
	
$~$

**Step3: Extract QM energy for each structure**
*  Run "QM_energy extractor.py". This program will read QM energy values from score files and write final values in
"*.runconfig" files in *global_settings.output_ensemble_dir* (settings.py).  
The format is simple, any other QM/DFT approach can be used to generate similar runconfig files.  

Format:  
	
	completed
	project_8000
	180 180      -538.421270588952
	0 120      -538.430477861295
	300 300      -538.435997315140
	60 120      -538.426056330737
	300 0      -538.426373919027
	240 180      -538.432283903238
	.
	.
	
**1st line**: "completed" or "submit" shows calculation status.  
**2nd line**: Project directory name.  
**Next lines**:  "$\phi$($^o$)  $\psi$($^o$)   ENEGRY(in Hartree)"  

$~$

**Step4: Generate_energy histograms for $\phi$, $\psi$ and $\omega$.**
*  Modify settings.py "#For making Histogram Files in kcal/mol" block, if needed.  
*  Run "runconfig_to_phi_psi_histogram.py"  
``` sh
>> ./runconfig_to_phi_psi_histogram.py
```

*  It will generate histogram files in global_settings.QM_histogram_file_out_dir
CHEKC fFOR OMEGA

$~$

**Step5: Make independent PAGS files**  
Use *gaussian_fitter_independent.py*   
*This program uses sklearn to fit gaussian expansion terms with the histogram data*  
For very rouge energy surface, you may need to modify _bound parameters_ in gaussian_fitter_independent.py.

*  Go to *global_settings.QM_histogram_file_out_dir* (output histogram files)and run:  
``` sh
>> ../../gaussian_fitter_independent.py <exphist file>  <number of gaussian terms> <1:for plot (default), 0: for silent>
``` 
eg:  
``` sh
>> ../../gaussian_fitter_independent.py FNA2toFNA3W0_phi.exphist 4 
``` 

This code will try to fit the gaussian expansion equations with the histogram parameters. It will output multiple (default 6) GNU-plot windows with the calculation number in the title, and will ask for a selection:
``` sh
>> Input selection number: (0 [default] for no selection) 
```

Identify the best fitting and enter the best calculation figure number to generate parameters in ./params directory.

*  In the case of bad fitting try different higher number for gaussian expansion equations (like, 5,6,7..).  

_Note:  It may take minutes for gaussian terms more than 6.  

To save time and try many parameters, run this code in silent mode for using diffent gaussian expansion terms, like:
```sh
>> for i in 3 4 5 6 7; do ../../gaussian_fitter_independent FNA2toFNA3W0_phi.exphist $i 0; done
```

This step will perform all calculations for given set of gaussian expansion numbers and generate required parameters in _.data_parameters directory_ for next faster access. So you do not need to spend time for running same calculations twice.  

Now running the same command will take a few seconds:
```sh
>> ../../gaussian_fitter_independent.py FNA2toFNA3W0_phi.exphist 4 
```
***For the given example, best fitting can be obtained for g=7 and calculation number 1,3,4 or 5.*  

Generated PAGS file looks like: (for g=4)
	
	FUNCTION  GAUSSIAN
	LINKID s?5nA3_s?5nA2_0_PHI
	a     6.600016      6.979097      3.821354      9.828471      # magnitudes of the distributions
	b    -1.882325e2    3.191396e1   -1.280898e2    1.315602e2    # midpoints of the distributions
	c     2.168577e3    3.600000e3    9.139045e2    3.600000e3    # twice the squares of the widths of the distributions
	d    -0.4806548   # intercept (coefficient of zeroth order term)

$~$

**Step7: Cluster similar potential functions**


	


