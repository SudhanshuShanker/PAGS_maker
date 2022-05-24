#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 17:58:13 2022

@author: sudhanshu
"""

class Data:
    pass

global_settings = Data()

# for ORCA dependency for QM calculation
global_settings.orca_path = '/home/sudhanshu/bin/orca_4_2_0_linux_x86-64_openmpi314/orca'

#for AmberTools
global_settings.AMBERSETTINGS = Data()
global_settings.AMBERSETTINGS.AMBERHOME = '/home/sudhanshu/HDD3/bin_rare/amber/amber16/'
global_settings.AMBERSETTINGS.param_dir = "./bin/tools"
global_settings.AMBERSETTINGS.temp_dir = "./temp/"



#for molecule maker
#global_settings.list_of_linkage_names = "./bin/model_lists/list_example" #for linkage without omega
global_settings.list_of_linkage_names = "./bin/model_lists/list_example_with_omega" #for linkage with omega
global_settings.simplified_pdb_output_dir = "./output/1_molecule_maker_pdbs/"

#for ensemble making and input file generation
global_settings.output_ensemble_dir = "./output/2_ensemble/"
global_settings.dihedral_rotation_step_size = 15  # for serious calculations make it 10 or 15



#for reading QM energies
#I used Condor and output dir name ends with ".inp_output"
#make your own rules and change in file QM_energy_extractor.py
global_settings.QM_out_dir_ends_with = ".inp_output"
global_settings.QM_out_dir_file_name = "condor.out_0"
        

#For making Histogram Files in kcal/mol
global_settings.QM_histogram_file_out_dir = "./output/3_histograms/"


