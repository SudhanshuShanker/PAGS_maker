#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 21:56:06 2021

@author: sudhanshu
"""
from src.read_and_fit_from_data_file import (
    read_runconfig_file_and_make_phi_psi_histfiles, 
    read_runconfig_file_and_make_phi_psi_histfiles_no_norm,
    identify_linkage_name, 
    data_from_multiple_exphist_files,
    make_exphist_data_from_x_y_data)

import os
from settings import global_settings
import numpy as np


def make_histogram_files():
    runconfig_files_dir = global_settings.output_ensemble_dir
    hist_file_out_dir = global_settings.QM_histogram_file_out_dir
    
    if not os.path.exists(hist_file_out_dir):
        os.mkdir(hist_file_out_dir)
    
    runconfig_files = []
    files_all = os.listdir(runconfig_files_dir)
    
    for i in files_all:
        if i.endswith("runconfig"):
            if i.startswith("p_"):  #p_9 for 2-->3  p_3 for 2-->4 # to control from many files 
                runconfig_files.append(i)
            
    counter = 1        
    w1_file_list = []     
    for rcf in runconfig_files:
        runconfig_file = runconfig_files_dir + rcf
        fid1 = open(runconfig_file,"r")
        data_rc_file = fid1.readline()
        fid1.close()
        print("Treating file#", counter, ":" , rcf)
        
        if data_rc_file.count("completed") == 1:    
            counter += 1
            
            identified_linkage_name = identify_linkage_name(runconfig_file)
            
            if identified_linkage_name.count("W0_")>0:
                read_runconfig_file_and_make_phi_psi_histfiles(
                    runconfig_file, hist_file_out_dir )
            elif identified_linkage_name.count("W1_")>0:   
                w1_file_list.append(identified_linkage_name.replace("60","").replace(
                    "180","").replace("300",""))
                if not os.path.exists(hist_file_out_dir + 'for_omega/'):
                    os.mkdir( hist_file_out_dir +'for_omega/')
                read_runconfig_file_and_make_phi_psi_histfiles_no_norm(
                    runconfig_file, hist_file_out_dir+'for_omega/')
                
                
                
                
                
                
        print("Completed! Omega part")
        
    # multiple files to one
    print("Combining for Phi Psi...")
    if len(w1_file_list) > 0:
        
        #w1_file_list = np.unique(w1_file_list)
        
        start_points =[]
        second_res_point=[]
        for i in w1_file_list:
            start_points.append(i.split("to")[0])
            second_res_point.append(i.split("to")[1][:-6] + i.split("to")[1][-5:-4])
            
        start_points = np.unique(start_points)
       
        
        all_exphist = os.listdir(hist_file_out_dir+'for_omega/')
        
        
        ends = ["_phi.exphist","_psi.exphist"]
        a_or_e = np.unique(second_res_point)
        
        #print(a_or_e)
          
        for i in start_points:
            
            for ae in a_or_e:
            
                for end_ in ends:
                    same_treat_files =[]
                    for fl in all_exphist:                
                        
                        if ( fl.startswith(i) and fl.endswith(end_) and 
                            ( (fl.count("to"+ae[:2]+"A"+ae[-1]) > 0) or (fl.count("to"+ae[:2]+"E"+ae[-1]) > 0)) ):
                            same_treat_files.append(hist_file_out_dir+'for_omega/'+fl)
                        
                    
                    #print(same_treat_files)
                    
                    if len(same_treat_files) <1:
                        continue
                    if len(ae) == 3:
                        ae = ae[:2] + "?" +ae[-1]
                    xh,yh=data_from_multiple_exphist_files(same_treat_files)
                    make_exphist_data_from_x_y_data(xh, yh, hist_file_out_dir+
                                                    i+ "to"+ae +"W1_?" +end_)
                                

    print("Done!")   
            
            
        #     print(i)
    
    
    
    

make_histogram_files()
