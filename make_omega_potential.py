#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 16:52:27 2021

@author: sudhanshu
"""
import os
import numpy as np
import matplotlib.pyplot as plt

def exphist_file_read( data_p,two_col=0):
    data= open(data_p,"r").readlines() 
    all_data=[]
    for i in data:
        d = i.split()
        all_data.append(d)     
    all_data=np.array(all_data).astype(float)
    if all_data.shape[1]==1:
        x = np.array(range(len(all_data[:,0])))*15
    else:
        x = all_data[:,0]
        
    data_np = all_data[:,-1]*627.509    
    #data_np = data_np - np.min(data_np)   
    y = data_np 

    return x,y     




out_dir = "/home/sudhanshu/HDD3/pipeline_output_and_mount_jazz_local/histogram_and_potential_files/W1_with_side_OH/for_omega/"

all_files = os.listdir(out_dir)

exp_hist = []

for i in all_files:
    if i.endswith('exphist'):
        exp_hist.append(i)
        

names_linakages = []

for i in exp_hist:
    if i.count('60_NoNorm_psi.exphist') > 0:
        # if i.startswith("P"):
        #     continue
        names_linakages.append( i.replace('60_NoNorm_psi.exphist','') )
        
all_out = []        
        
for  i in names_linakages:
    
    for j in exp_hist:
        
        if j.count(i+"60_NoNorm_psi") > 0:
            x,y = exphist_file_read(out_dir+j)
  
            min_y_60 = min(y)
        elif j.count(i+"180_NoNorm_psi") > 0:
            x,y = exphist_file_read(out_dir+j)
            min_y_180 = min(y)
            
        elif j.count(i+"300_NoNorm_psi") > 0:
            x,y = exphist_file_read(out_dir+j)
            min_y_300 = min(y)
    
    d= np.array([min_y_60, min_y_180, min_y_300]) 
    d = d - min(d)
    
    all_out.append([i, d])
            
    
options_for_omega = [
    [0,1,2],
    [0,2,1],
    
    [1,0,2],
    [2,0,1],
    
    [1,2,0],
    [2,1,0],   
    
    ] 



# types_list =   [ 'toFNA1W1_a' , 'toFNE1W1_a', 'toFNA6W1_a' , 'toFNE6W1_a',
#                  'toFSA1W1_a' , 'toFSE1W1_a', 'toFSA6W1_a' , 'toFSE6W1_a',
#                  'toFNA1W1_e' , 'toFNE1W1_e', 'toFNA6W1_e' , 'toFNE6W1_e',
#                  'toFSA1W1_e' , 'toFSE1W1_e', 'toFSA6W1_e' , 'toFSE6W1_e']   

types_list =   [ '1W1_a' , '6W1_a' , 
                  '1W1_e' , '6W1_e' , ]      







# types_list =   [ 'A1W1_a' , 'A6W1_a' , 
#                   'A1W1_e' , 'A6W1_e' , 
#                   'E1W1_a' , 'E6W1_a' , 
#                   'E1W1_e' , 'E6W1_e' , ]    


# types_list = ['toFNE', 'toFNA']
 
pots = dict()
for i in types_list:
    pots[i]=[]
    

    
for t in types_list:
    plt.figure()
    for i in all_out:
        col = 'y'
        types= t
        if i[0].startswith(''):
            if i[0].count(types) > 0:
                col = 'k'
                if i[0].startswith('P'):
                    col = 'm'
                #print(i)
                pots[t].append(i[1])
                plt.plot(i[1], col)
        plt.title(types)
            
            
for i in types_list:
    pots[i] = np.mean(pots[i],0)
    plt.figure()
    plt.plot(pots[i], col)
    plt.ylim([0,8])
      
            
            
dets = [ 'FOA1', 'FOA6', 'FOE1' ,'FOE6' ]         
            
counter=-1
for i in types_list:
    counter +=1
    print (i,"\n\n")
    print( "LINKID  " + dets[counter]+"  "+ i +"  # Linkage Name")
    print( "FUNCTION    SQUARE     # the type of equation GAUSSIAN :gaussian and VONMISES :Von-Mises")
    print( "a     "+ "%4.2f" % (pots[i][0])  +"    # Minima for 60 Degree")
    print( "b     "+ "%4.2f" % (pots[i][1])  +"    # Minima for 180 Degree")
    print( "c     "+ "%4.2f" % (pots[i][2])  +"    # Minima for 300 Degree")
    print( "d     "+ "%4.2f" % 0  +"    # Unused value")

            
# {'1W1_a': array([2.68127559, 0.25625587, 2.67291683]),
#  '6W1_a': array([1.08540334, 4.90912885, 1.00665386]),
#  '1W1_e': array([4.5392214 , 0.33532454, 4.53653213]),
#  '6W1_e': array([1.04698457, 5.06644359, 0.94707239])}
            
            
            
            
            
            
            
            
            
            
            