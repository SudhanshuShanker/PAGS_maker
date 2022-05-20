STEPS TO RUN PAGSMAKER PIPELINE
===============================

**Requirements**: AMBERTOOLS, ORCA, Python3

*After each step, analyze output files carefully*  

**Step1**: edit the settings.py file for ORCA and AMBERHOME

**Step2**: To make a simplified pdb file for a linkage:  
*  Make a list file with linkage names (like ./bin/model_lists/list_example) 
	   [follow the extra/linkage_grammer_file.dat]  
*  Modify "for molecular maker" block in settings.py for the list path (model_name_list) and  pdb_out_dir if rquired.  
*  Run molecule_maker.py (./molecule_maker.py)  
``` sh
>> ./molecule_maker.py
```
This step will generate a simplified pdb molecule in global_settings.simplified_pdb_output_dir (settings.py)  
	
	
**Step3**: To generate rotational ensemble and QM input files:
	a: modify settings.py in "#for ensemble making and input file generation" block.
	b: run make_ensemble_and_qm_inputs.py
	use input files generated in global_settings.output_ensemble_dir/project_****/QM_inps for QM calculation using orca.
	I used condor for the calculation (see bin/tools/QM_tools file for all raw files)
	You can use your schedular and make directories like global_settings.QM_out_dir_ends_with (settings.py)
	and output_files with same names as global_settings.QM_out_dir_file_name (settings.py)
	looks like:
	
	ensemble_0_0.inp                       <-- Input INP file for ORCA
	ensemble_0_0.inp_output/               <-- Output Directory
	ensemble_0_0.inp_output/condor.out_0   <-- Score file in Output Directory
	ensemble_0_45.inp                       <-- Input INP file for ORCA
	ensemble_0_45.inp_output/               <-- Output Directory
	ensemble_0_45.inp_output/condor.out_0   <-- Score file in Output Directory
	.
	.
	
Step4: QM energy extraction for each structure.
	a: run "QM_energy extractor.py"
	this program will read QM energy values from score files and write final values in
	"*.runconfig" files in global_settings.output_ensemble_dir (settings.py)
	The format is simple so you can use any other QM/DFT approach to generate similar runconfig files.
	format:
	
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
	
first line : completed or submit shows calculation status
2nd line : project directory names
Next lines:  <PHI>  <PSI>   <ENEGRY IN Hartree>

	
Step5: Generate_energy histograms for phi, psi and omega.
	a: modify settings.py "#For making Histogram Files in kcal/mol" block, if needed.
	b: run "runconfig_to_phi_psi_histogram.py"	
	It will generate histogram files in global_settings.QM_histogram_file_out_dir
	CHEKC fFOR OMEGA


Step6: Making independent potential functions
	Application: gaussian_fitter_independent.py
	Note: This is done using sklearn fit for multiple gaussian expansion terms.
	For very rouge energy surface, you may need to modify bound parameters in gaussian_fitter_independent.py.
	To RUN:
	a): go to global_settings.QM_histogram_file_out_dir and run
	>> ../../gaussian_fitter_independent.py <exphist file>  <number of gaussian terms> <1: (default) for plot 0: for silent>
	eg:
	>> ../../gaussian_fitter_independent.py FNA2toFNA3W0_phi.exphist 4 
	
	This code will try to fit the gaussian expansion equation with the histogram parameters. It will show multiple GNU-plot windows with calculation number in the title, and will ask for a selection:
	>> Input selection number: (0 [default] for no selection) 

	Identify the best fitting and enter the best calculation figure number to generate parameters in ./params directory.
	
	b: in case of bad fitting try different number for gaussian expansion equations (like, 5,6,7..)
	
	Note:  It may take minutes for gaussian terms more than 6.
	
	To save time and try all parameters, run this code in silent mode for different gaussian equation numbers, like:
	>> for i in 3 4 5 6 7; do ../../gaussian_fitter_independent FNA2toFNA3W0_phi.exphist $i 0; done
	
	This step will perform all calculations for given series of gaussian equation numbers and generate useful parameters in 
	.data_parameters directory for next faster access.
	
	now running 
	>> ../../gaussian_fitter_independent.py FNA2toFNA3W0_phi.exphist 4 
	command will take few seconds to show all calculations.
	
	for the given example, best fitting is obtained for g=7 and 1,3,4 or 5.
	
	generated parameter looks like:
	
	FUNCTION  GAUSSIAN
	LINKID s?5nA3_s?5nA2_0_PHI
	a     6.600016      6.979097      3.821354      9.828471      # magnitudes of the distributions
	b    -1.882325e2    3.191396e1   -1.280898e2    1.315602e2    # midpoints of the distributions
	c     2.168577e3    3.600000e3    9.139045e2    3.600000e3    # twice the squares of the widths of the distributions
	d    -0.4806548   # intercept (coefficient of zeroth order term)


Step7: Clustering potential functions


	


