#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 21:56:06 2021

@author: sudhanshu
"""
from src.read_and_fit_from_data_file import (
    read_runconfig_file_and_make_phi_psi_histfiles, 
    read_runconfig_file_and_make_phi_psi_histfiles_no_norm)

import os
from settings import global_settings



def make_histogram_files():
    runconfig_files_dir = global_settings.output_ensemble_dir
    hist_file_out_dir = global_settings.QM_histogram_file_out_dir
    runconfig_files = []
    files_all = os.listdir(runconfig_files_dir)
    
    for i in files_all:
        if i.endswith("runconfig"):
            if i.startswith("p_"):  #p_9 for 2-->3  p_3 for 2-->4 # to control from many files 
                runconfig_files.append(i)
            
    counter = 1        
            
    for rcf in runconfig_files:
        runconfig_file = runconfig_files_dir + rcf
        fid1 = open(runconfig_file,"r")
        data_rc_file = fid1.readline()
        fid1.close()
        print("Treating file#", counter, ":" , rcf)
        if data_rc_file.count("completed") == 1:    
            counter += 1
            name_linakge = read_runconfig_file_and_make_phi_psi_histfiles(
                runconfig_file, hist_file_out_dir )
            if name_linakge.count("W1_")>0:
                if not os.path.exists(hist_file_out_dir + 'for_omega/'):
                    os.mkdir( hist_file_out_dir +'for_omega/')
                read_runconfig_file_and_make_phi_psi_histfiles_no_norm(
                    runconfig_file, hist_file_out_dir+'for_omega/')
        print("Completed!")

make_histogram_files()
