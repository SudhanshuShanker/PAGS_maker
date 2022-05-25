#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 16:52:27 2021

@author: sudhanshu
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from src.read_and_fit_from_data_file import exphist_file_read
from settings import global_settings
from src.clustering_parameter_functions import link_name_converter


def link_name_converter_for_omega(link_name):
    lnk = link_name_converter(link_name+"_phi.params")    
    lnk =lnk.replace("_1E_","_E1_").replace("_1A_","_A1_")
    # as current support is for one omega only, we can add OM1 manualy for now
    lnk =lnk +"OM1"
    return lnk


def estimate_omega_minima_for_all_linkages(omega_dir, param_out_dir):

    all_files = os.listdir(omega_dir)

    #identify linkage names
    names_linkages = []
    exp_hist =[]
    for i in all_files:
        if not i.endswith('exphist'):
            continue
        exp_hist.append(i)
        if i.count('60_NoNorm_psi.exphist') > 0:
            names_linkages.append( i.replace('60_NoNorm_psi.exphist','') )
       
    #identify omeag values for all identified linkages
    all_out = dict()       
    for  i in names_linkages:        
        for j in exp_hist:
            
            if j.count(i+"60_NoNorm_psi") > 0:
                x,y = exphist_file_read(omega_dir+j)
      
                min_y_60 = min(y)
            elif j.count(i+"180_NoNorm_psi") > 0:
                x,y = exphist_file_read(omega_dir+j)
                min_y_180 = min(y)
                
            elif j.count(i+"300_NoNorm_psi") > 0:
                x,y = exphist_file_read(omega_dir+j)
                min_y_300 = min(y)
        
        d= np.array([min_y_60, min_y_180, min_y_300]) 
        d = d - min(d)
        
        all_out[i]= d
    
    #writing parameters
    if not os.path.exists(param_out_dir):
        os.mkdir(param_out_dir)
    

    for i in names_linkages:
        out_param_file_name = i + "_om1.params"
        fid = open(param_out_dir + out_param_file_name,"w+")

        fid.write("# PAGS Parameters for linkID "+ i + '_om1' +"\n")
        fid.write( "LINKID  "+ link_name_converter_for_omega(i) +"  # Linkage Name\n")
        fid.write( "FUNCTION    SQUARE     # the type of equation GAUSSIAN: gaussian, and SQUARE: quadratic\n")
        fid.write( "a     "+ "%4.2f" % (all_out[i][0])  +"    # Minima for 60 Degree\n")
        fid.write( "b     "+ "%4.2f" % (all_out[i][1])  +"    # Minima for 180 Degree\n")
        fid.write( "c     "+ "%4.2f" % (all_out[i][2])  +"    # Minima for 300 Degree\n")
        fid.write( "d     "+ "%4.2f" % 0  +"    # Unused value\n")
        fid.close()
          
        print(">------------------------------------------" )
        print("Omega params file for", i, "linkage save in", param_out_dir)
        print("params file name:",out_param_file_name )
        print("------------------------------------------!" )

            
if __name__ == "__main__":
    omega_dir = global_settings.QM_histogram_file_out_dir + "/for_omega/"
    out_dir = global_settings.QM_histogram_file_out_dir + "params/"    
    estimate_omega_minima_for_all_linkages(omega_dir,out_dir)
    
            
            
            
            
            