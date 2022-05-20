#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 12:24:16 2020

@author: sudhanshu
"""

import os
import time
import subprocess
from settings import global_settings

def ss_importdata(*inputs):
    data_only_mode = 0;
    if len (inputs) > 1:
        data_only_mode = inputs[1]
    flnm = inputs[0]  
    fid = open(flnm,'r')
    
    data = fid.readlines()
    fid.close()
    return data


def ss_write_in_file(fl, msg):
    fid2 = open(fl,"w+")
    for i in msg:
        fid2.write(i+"\n")
    fid2.close()


def read_energies(dirnm):
    all_fl = os.listdir(dirnm)
    dir_interest = []
    for i in all_fl:
        if i.endswith(global_settings.QM_out_dir_ends_with):
            dir_interest.append(i)
    all_data = []
    for i in dir_interest:
        d=subprocess.Popen(["grep","FINAL" , dirnm +"/" + i+  "/" + 
                            global_settings.QM_out_dir_file_name, ]
                           ,stdout=subprocess.PIPE)
        dd = d.communicate()[0]
        if len(dd) < 1:
            continue
        score_d = str(dd.splitlines()[-1])
        #final_data = i[9:-11] + score_d[27:-1]
        final_data = i.replace('ensemble_','').replace('.inp_output','').replace('_',' ') + score_d[27:-1]
        all_data.append(final_data)
        
    return all_data    
        

def main():
   
    base_dir = global_settings.output_ensemble_dir
    run_inp_files=os.listdir(base_dir)
   
    for r_file in run_inp_files:
        if not r_file.endswith("runconfig"):
            continue
        
        current_r_file = base_dir +  r_file
        r_file_data = ss_importdata(current_r_file)
        input_dir = r_file_data[1]
        input_dir = input_dir.rstrip()
        if r_file_data[0].count("submit") == 1:   
            data = read_energies(global_settings.output_ensemble_dir + 
                                 input_dir + "/QM_inps/")        
            ss_write_in_file( current_r_file, ["completed",input_dir]+data)
        
            print("QM energy data is written in:", current_r_file )
        print("All QM energy extraction written in runconfig files in", base_dir )
        time.sleep(1)
            
            

main()       
    #break
