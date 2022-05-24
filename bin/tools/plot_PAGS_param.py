#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 15:25:52 2021

@author: sudhanshu
"""
import os
import sys
sys.path.insert(1, 'use path for this pipeline directory') ## MODIFY THIS

from src.ss_utils import CHI_energy_from_gaussian_file, link_id_handler

import matplotlib.pyplot as plt
import numpy as np

max_filter = 6.5


def plot_param_file(param_file_nm):


    if param_file_nm.upper().find('_OM') > 0:
        print("$\omega$ not analyzed here")
        return
    
    elif param_file_nm.upper().find('_PHI') > 0:
        phi_treatment = 1
        title = "$\phi$"
        
    elif param_file_nm.upper().find('_PSI') > 0:
        phi_treatment = 0
        title = "$\psi$"
    
    else:
        print("can not decide the phi or psi type. please modify the code.")
        return
    
    handle = CHI_energy_from_gaussian_file(param_file_nm) 
    
    handle.label_on=1
    handle.step_angle = 1
    handle.remap_360=0
    handle.phi_treat = phi_treatment
    handle.print_me = 0
    handle.run_me()
    # plt.title(title+"(" + res_V[num]+")", fontsize=12)
    # plt.plot([values[num],values[num]],[0,9])
    #print(handle.x_s)
    plt.xlabel("angle[$^o$]", fontsize=16)

    plt.ylabel("kcal/mol", fontsize=16)
    # plt.show(block=False)
    plt.title(param_file_nm)
    
    plt.xlim([handle.x_s[0], handle.x_s[-1]-handle.step_angle])
    plt.xticks(range(handle.x_s[0], handle.x_s[-1]+60, 60))
    plt.plot([handle.x_s[0], handle.x_s[-1]+60],[0,0],'b--',label = '_nolegend_')
    # plt.show()
    return handle.x_s, handle.y_s
    
mean_xx = []
mean_arr = []
legends = []
for i in range(1,len(sys.argv)): 
    
    input_file = sys.argv[i]
    xx, yy = plot_param_file(input_file)
    
    if len(sys.argv)>2:
        if max(yy) > max_filter:
            mean_arr.append(yy)
            mean_xx.append(xx)


mean_xx = np.array(mean_xx)
mean_arr = np.array(mean_arr)
plt.legend(sys.argv[1:])

if len(sys.argv) > 2:
    #print("avarage exphist parameter:")
    if np.sum(np.std(mean_xx,0)) == 0.0:
        yy = np.mean(mean_arr,0)
        yy = yy - np.min(yy)
        plt.plot(np.mean(mean_xx,0), yy , 'k--')
        for i,j in zip(xx,yy/627.509 ):
            if i%10 == 0:
                print(i,j)
#print("\n")   
plt.show()
