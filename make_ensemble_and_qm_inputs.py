#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 12:07:14 2021

@author: sudhanshu
"""
from src.one_param_converter_v1_rotation_type import run_pipeline
from settings import global_settings
import importlib



def read_simplified_pdbs_and_generate_ensemble(simplified_pdb_dir, step_size):

    py_path = (simplified_pdb_dir+".").replace("/",'.').replace("..","")+".connection_list_check"
    file_data = importlib.import_module(py_path) 
    
    files_and_torsions = file_data.files_and_torsions # information of pdb_file and phi_psi variables
    omega_tors = file_data.omega_tors # information of omega variable
 
    extra_tors_angs = "[ 0 ]" # depricated
    for index, i in enumerate(files_and_torsions):
        flnm= path_fl + i[0]
        tors_v = i[1]    
        extra_tors = omega_tors[index]
        print((flnm, tors_v, step_size, extra_tors, extra_tors_angs))
        run_pipeline(flnm, tors_v, step_size, extra_tors, extra_tors_angs,8000)



step_size = global_settings.dihedral_rotation_step_size
path_fl = global_settings.simplified_pdb_output_dir   # path of the previously generated simplified pdb files

read_simplified_pdbs_and_generate_ensemble(path_fl, step_size)