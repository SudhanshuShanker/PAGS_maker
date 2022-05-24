#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 14:02:58 2021

@author: sudhanshu
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 14:19:04 2021

@author: sudhanshu
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os


from src.ss_utils import ( PAGS_energy_from_gaussian_file, 
                          ss_list_all_files_with_ext, 
                          ss_exist )
from src.clustering_parameter_functions import ( compare_with_profiles,
                                                link_name_converter,
                                                make_exphist_data)

from sklearn.cluster import DBSCAN




def list_phi_or_psi_files(in_file_dir, phi_or_psi_out):
    all_files = ss_list_all_files_with_ext(in_file_dir, '.params')       
    # PSI files
    psi_files = []
    phi_files = []
    for i in all_files:
        if i.count('_psi_g') > 0 : # identifier in the output param files
            psi_files.append(i)
        elif i.count('_phi_g') > 0 :
            phi_files.append(i)                
                
    psi_files.sort()
    phi_files.sort()
    
    if phi_or_psi_out == "PHI":
        use_files = phi_files
    elif phi_or_psi_out == "PSI":
        use_files = psi_files

    return use_files

def join_clusters(cluster_list,i,j):
    main = j
    merge = i
    if i < j:
        main = i
        merge = j
    
    return_v =[]
    for v in cluster_list:
        if v== merge:
            return_v.append(main)
        else:
            shifter=0
            if v > merge:
                shifter = 1
                
            return_v.append(v-shifter)
            
    
    
    
    return np.array(return_v)




def main(phi_or_psi_out, in_file_dir, out_file_dir):
     
    ## Setup the PAGS parameter loader 
    handle = PAGS_energy_from_gaussian_file(" ")
    handle.step_angle = 1
    handle.remap_360=0
    handle.print_me = 0
    handle.plot = False
    
    xticks = range(0,361,120)
    yticks = range(0,15,2)
    columns = 10
    rows = 4
    if phi_or_psi_out == "PHI":  #setting is loading PHI files
        handle.phi_treat = 1  # as phi ranges from -180 to 180
        cut_off = 8 # ignore high energy values more than
        eps_v = 2.   # for DBSCAN
        columns = 5
        rows = 1
        xticks = (range(-180,181,120)) # as phi ranges from -180 to 180
        yticks = (range(0,15,2))
        fig_size = (12,3)
    else:                       #setting is loading PSI files
        handle.phi_treat = 0
        cut_off = 2.  # ignore high energy values more than
        eps_v = .4 # for DBSCAN
        rows = 3
        columns = 5
        yticks = (range(0,9,2))
        fig_size = (12,12)

  
    # reading param files from directory 
    use_param_files = list_phi_or_psi_files(in_file_dir, phi_or_psi_out)
    

    # generating y values for each param file
    all_y_s =[]
    for i in use_param_files:      
        handle.filen = in_file_dir + i 
        handle.run_me() 
        all_y_s.append(handle.y_s) 
    
    profs = np.array(all_y_s)    
    
    
    # similarity score matrix for each profile using standard deviation 
    # between y points. High value for dissimal and low value for similar
    #profiles.   
    
    scores = np.zeros([len(profs),len(profs)])
    for i,prf1 in enumerate( profs):
        for j,prf2 in enumerate(profs):
            scores[i,j] = compare_with_profiles(prf1,[prf2],cut_off)[0]

    

    # clustering profiles using DBSCAN 
    clustering = DBSCAN(eps=eps_v, min_samples=1).fit(scores)
    
    
    manual_clusters = clustering.labels_[:]
    if phi_or_psi_out == 'PSI':
        manual_clusters = join_clusters(manual_clusters,17,26) # visually identified
        manual_clusters = join_clusters(manual_clusters,15,21) 
        manual_clusters = join_clusters(manual_clusters,21,26) 
        manual_clusters = join_clusters(manual_clusters,20,25) 
        manual_clusters = join_clusters(manual_clusters,10,16)
        manual_clusters = join_clusters(manual_clusters,2,21)
        manual_clusters = join_clusters(manual_clusters,15,18)
        manual_clusters = join_clusters(manual_clusters,10,13)
        manual_clusters = join_clusters(manual_clusters,17,16)
        manual_clusters = join_clusters(manual_clusters,3,19)
        manual_clusters = join_clusters(manual_clusters,2,19)
        manual_clusters = join_clusters(manual_clusters,1,6)
        manual_clusters = join_clusters(manual_clusters,1,17)
        manual_clusters = join_clusters(manual_clusters,10,15)
        manual_clusters = join_clusters(manual_clusters,0,15)
       
    # calculating mean profiles and names of linkages in each cluster
    cluster_files =dict()
    cluster_links = dict()
    mean_profs = dict()
    number_of_clusters = np.max(manual_clusters)+1  
    for i in range(number_of_clusters):              
        
        curr_files = []  # for files of same clusters        
        curr_links = "LINKID " + "PFW" + str(i+1) + phi_or_psi_out[1] +"  " 
        # for link names of same cluster
        idx = np.where( manual_clusters == i ) # members of same cluster
        for j in idx[0]:
            curr_files.append(use_param_files[j])
            curr_links = curr_links + link_name_converter(use_param_files[j]) +"  "
            
        cluster_files[i] = curr_files
        cluster_links[i] = curr_links

        # mean profiles of each cluster
        mean_prof = np.mean(profs[np.where(manual_clusters==i)[0]],0) 
        mean_profs[i]= mean_prof


    # Plotting similar profiles in same subplots   
    fig=plt.figure(figsize=fig_size)
    #plt.subplots(1,columns)
    for i in range(np.max(manual_clusters)+1):
        plt.subplot(rows,columns,i+1)
        
        
        if i==0:
            plt.ylabel("kcal/mol", fontsize = 12)
        if phi_or_psi_out == "PHI":
            if i==2:
                plt.xlabel("Angle [Degree]", fontsize = 12)
        else:
            if [5,10,15].count(i) == 1:
                plt.ylabel("kcal/mol", fontsize = 12)
            if i >14:
                plt.xlabel("Angle [Degree]", fontsize = 12)
                


        for j in np.where(manual_clusters==i)[0]:
            plt.plot(handle.x_s,profs[j])
        plt.title(i)
        plt.plot(handle.x_s,mean_profs[i],'k--')
        plt.yticks(yticks)
        plt.xticks(xticks)
        plt.grid(1)

    
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.3, wspace = 0.3)   
    plt.show()
        
            
    # writing avarage files
    if not ss_exist(out_file_dir,"dir"):
        os.mkdir(out_file_dir)
    
    for i in mean_profs.keys():
        outfile = out_file_dir+ phi_or_psi_out + "_" +str(i+1) +".exphist"
        link_detail_file = out_file_dir+ phi_or_psi_out + "_" +str(i+1) +".links"
        param_links = cluster_links[i]    
        fid2 = open(link_detail_file,"w+")
        fid2.write(param_links+"\n")
        fid2.close()        
        make_exphist_data(handle.x_s, mean_profs[i], outfile)
        
        
    #saving figure
    fig.savefig(out_file_dir+ phi_or_psi_out + ".png", dpi=300)
   

if __name__ == '__main__':
    phi_or_psi_out = "PSI"
    in_file_dir = "./bin/parameter_files/Pyranose_Furanose_W0/"    
    out_file_dir = "./bin/parameter_files/Pyranose_Furanose_W0/clustered_pots/"    
    
    print ("Making mean exphistogram files for", phi_or_psi_out,".")
    print (".params files from directory", in_file_dir, "are used." )
    print ("Mean exphist output directory is", out_file_dir)
    main(phi_or_psi_out, in_file_dir, out_file_dir)
    
    