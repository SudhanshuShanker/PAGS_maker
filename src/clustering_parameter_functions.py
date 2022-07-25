#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 15:43:34 2021

@author: sudhanshu
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.insert(1, '/home/sudhanshu/bin/my_scripts/my_python_tools')

def make_exphist_data(use_x, use_y, outfile=""):
    
    if outfile == "":        
        for i,j in zip(use_x, use_y):
            print(("%5.1f " % i) + ("%21.19f" % ( j/627.509) ) )
    else:
        
        fid_out = open(outfile,"w+")
        
        for i,j in zip(use_x, use_y):
            fid_out.write(("%5.1f " % i) + ("%21.19f" % ( j/627.509) ) +"\n")
        fid_out.close()

def my_clustering(all_profiles, cut_off = 6):
    local_profiles = np.copy(all_profiles)
    y_s_mse_mat = np.zeros([40,40])
    for i in range(len(local_profiles)-1):
        i_ys = local_profiles[i]
        i_ys[i_ys>cut_off] = cut_off
        for j in range(i+1,len(local_profiles)):
            j_ys = local_profiles[j]
            j_ys[j_ys>cut_off] = cut_off
            val = np.std((i_ys - j_ys))
            y_s_mse_mat[i][j] = val
            #y_s_mse_mat[j][i] = val
    
    rep_mat = np.copy(y_s_mse_mat)

    rep_mat[rep_mat==0] = 300
    
    clusters = [ ]
    
    for i in range(80000000):
        min_val = np.min(rep_mat)
        
        if min_val == 300:
            break   
        idx = np.where(rep_mat==min_val)
        rep_mat[idx]=300
        if min_val < 0.2:
            
            
            same_pair = [idx[0][0],idx[1][0]]
            same_pair.sort()
            #print(same_pair)
            
        else:
            continue
        
        if len(clusters) > 0:    
            for i in range(len(clusters)):    
                match_found = 0
                print("::",clusters)
                for j in same_pair:
                    
                    if clusters[i].count(j) > 0:
                        match_found = 1
                        print("match found")
                        same_pair.remove(j)
                        
                        
                if match_found == 1:
                    break
                
            if match_found == 1:
                if len(same_pair)>0:
                    if clusters[i].count(same_pair[0]) == 0:
                        clusters[i].append(same_pair[0])
                continue
            else:
                if len(same_pair)>0:
                    clusters.append(same_pair)
                
            
            if len(clusters) > 20:
                break
                
                    
                    
                
                    
                
                
        else:
            
            clusters.append(same_pair)
                        
                        
            
        
        
        
       
        
        #print(clusters)



def my_clusers2(all_profiles,  cluster_base , cut_off = 2.):
    local_profiles = np.copy(all_profiles)
    
    # cluster_base = [0,1,2,3,6,7,8,11,24,25]
    # #cluster_base = [0,1,2,3,4,5,6,7,8,9,10,11,12]
    # cluster_base = [ 0,2,3,5,9 ]
    # cluster_base = [ 28,23,1,17,8 ]
    cluster_dict = dict()
    
    for i, prof in enumerate( cluster_base ):
        cluster_dict[prof] = [prof]
    
    
    for i, prf in enumerate( local_profiles ):
        #cut_off = 5
        i_ys = prf
        i_ys[i_ys>cut_off] = cut_off
        score = []
        for j in cluster_base:
            if i ==j:
                score.append(0)
                continue

            j_ys = local_profiles[j]
            j_ys[j_ys>cut_off] = cut_off
            val = np.std((i_ys - j_ys))
            score.append(val)
            
        score = np.array(score)
        
        position = np.where(score == np.min(score))
        #print(position[0][0])
        
        if cluster_dict[cluster_base[position[0][0]]].count(i) == 0:        
            cluster_dict[cluster_base[position[0][0]]].append(i)
            
   
    return cluster_dict
    

def compare_with_profiles(input_prof,  base_profiles , cut_off = 2.):
  
    all_vals = []
  
    base_profiles = np.array(base_profiles)
    for i, prf in enumerate( base_profiles ):
        #cut_off = 5
        i_ys = prf
        i_ys[i_ys>cut_off] = cut_off    
        
        j_ys = np.array(input_prof)
        j_ys[j_ys>cut_off] = cut_off
        val = np.std((i_ys - j_ys))
        all_vals.append(val)
    return np.array(all_vals)


def compare_with_profiles2(input_prof,  base_profiles , cut_off = 4.):
  
    all_vals = []
  
    base_profiles = np.array(base_profiles)
    
    for i, prf in enumerate( base_profiles ):
        #cut_off = 5
        i_ys = prf       
        j_ys = np.array(input_prof)
        
        
                
        diffs1 = np.std(j_ys[j_ys>cut_off] - i_ys[i_ys>cut_off])
        diffs2 = np.std(j_ys[j_ys>cut_off] - i_ys[i_ys>cut_off])
        
        val = np.std(diffs)
        all_vals.append(val)
    return np.array(all_vals)


    
    
def link_name_converter(link_name):
    if link_name.endswith(".params"):
        link_name = link_name.replace("?","$")
        vals=link_name.replace(".params","")[:-3].split("_")        
        phi_psi = vals[1].upper()
        A_or_E = ""
        if ( ( phi_psi == "E" ) or ( phi_psi == "A") or ( phi_psi == "$")):
            A_or_E = phi_psi
            phi_psi = vals[2].upper()
            
        omega = vals[0][-2:]+A_or_E
        vals = vals[0][:-2]
        vals = vals.replace("P", "6c").replace("FN", "5n").replace("FS", "5s")
        #vals = vals.replace("F", "5$")
        
        vals= vals.split('to')
        v1 = vals[1]
        v2 = vals[0]
        v1 = v1.replace("P","6c")
        vals = "s?"+v1+"_s?" + v2 + "_" +omega.replace("W","") + "_" +phi_psi
        
        return vals.replace("$","?")
    else:        
        vals = link_name.split("_")
        
        
        
        
      
    
    