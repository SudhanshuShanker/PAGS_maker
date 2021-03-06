#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 12:55:05 2021

@author: sudhanshu
"""
import numpy as np
from src.ss_utils import ss_exist


if not __name__ == "__main__":
    from settings import global_settings
    step_size = global_settings.dihedral_rotation_step_size
else:
    step_size = 15
    
class extract_multidimentional_data():
    '''hello'''
    def __init__(self, file_name, dim_, change_180 = 0, norm = 1):
        self.file_name = file_name
        self.dimention = dim_
        self.norm = norm
        self.change_180 = change_180        
        self.no_modify = 0
        if not file_name.strip() == "":
            self.data_file_read_server()
        
        
        
        
    def data_file_read_server(self ):
        data_p = self.file_name 
        data= open(data_p,"r").readlines() 
        all_data=[]
        counter = 0
        for i in data:
            counter +=1
            if counter <3:
                continue
            d = i.split()
            if (len(d) > 3):
                d= d[0:3]
            
            all_data.append(d)     
        
        #print(all_data)
        all_data=np.array(all_data).astype(float)
        
        # cleaning bad values
        
        indexes_to_del = list(np.where((all_data[:,2]+2 - np.median(all_data[:,2])) < 0))
        indexes_to_del.sort()
        indexes_to_del.reverse()
        
        for del_row in indexes_to_del:
            all_data = np.delete(all_data, del_row,0)
        
        
        
        
        #self.all_data = all_data
        
        sort_dim = self.dimention - 1
        
        
        all_data = all_data[all_data[:,sort_dim].argsort()]
        self.all_data = all_data
        if all_data.shape[1]==1:
            temp_x = np.array(range(len(all_data[:,0])))*15
        else:
            temp_x = all_data[:,sort_dim]
            
        data_np = all_data[:,-1]*627.509    
        if self.norm == 1:
            data_np = data_np - np.min(data_np)   
        temp_y = data_np    
        
        self.x = np.unique(temp_x)
        y=[]
        for i in self.x:
            y.append(np.min(temp_y[temp_x==i]))
            
            # if i == 240.:
            #     print(temp_y)
            
        self.y = np.array(y)
        
        if self.x[-1] != 360.0:
            self.x = np.append(self.x, 360.0)
            self.y = np.append(self.y, self.y[0])
            
        #return(x,data_np)
        return ("X and Y initiated")
    
    def change_to_180(self, x,y):
        x = np.array(x)
        # print(x)
        # print(y)
        id1 = np.where(x>=180)[0]
        id2 = np.where(x<180)[0]
        ids = np.concatenate((id1[:-1],id2,[id1[0]]))
        x[x>=180] =  x[x>=180] - 360
        x[-1] = -x[-1]
        # print(x)
        # print(x[ids])
        # print(y[ids])
        x =x[ids]
        x[-1]=180
        return [ x , y[ids]]
    
    
    def make_exphist_data(self,outfile=""):
        use_x = self.x
        use_y = self.y
        
        if self.no_modify == 0:
            if self.change_180 == 1:
                [use_x, use_y] = self.change_to_180(use_x, use_y)
            
        
        if outfile == "":        
            for i,j in zip(use_x, use_y):
                print(("%5.1f " % i) + ("%21.19f" % ( j/627.509) ) )
        else:
            fid_out = open(outfile,"w+")
            for i,j in zip(use_x, use_y):
                fid_out.write(("%5.1f " % i) + ("%21.19f" % ( j/627.509) ) +"\n")
            fid_out.close()
            



def identify_linkage_name(file):
    dir_path = file[:[i for i in range(len(file)) if file[i]=="/"][-1]+1]  
    dir_name = file.split("/")[-1].replace(".runconfig","").replace("p_","project_")
    
    if ss_exist(dir_path+dir_name+"/README", "file"):
        read_me_data = []
        fid=open(dir_path+dir_name+"/README","r")
        read_me_data = fid.readlines()
        fid.close()
        
        for i in read_me_data:
            if i.startswith("Input pdb file"):
                pdb_input_name = i.split("/")[-1].strip()[:-4]
                
    
    else:
        pdb_input_name = file.split("/")[-1].replace(".runconfig","")
    
    return pdb_input_name



def read_runconfig_file_and_make_phi_psi_histfiles(file, out_dir="./"):   
    
    pdb_input_name = identify_linkage_name(file)        

    gg = extract_multidimentional_data(file, 1,1)
    gg.make_exphist_data(out_dir+pdb_input_name+"_phi.exphist")    
    
    gg = extract_multidimentional_data(file, 2,0)
    gg.make_exphist_data(out_dir+pdb_input_name+"_psi.exphist")
    return 0


def read_runconfig_file_and_make_phi_psi_histfiles_no_norm(file, out_dir="./"):    
    #file = "/home/sudhanshu/HDD3/mount_jazz/p_7000.runconfig"

    
    pdb_input_name = identify_linkage_name(file) 
        
    nomalize = 0
    gg = extract_multidimentional_data(file, 1,1, nomalize)
    
    gg.make_exphist_data(out_dir+pdb_input_name+"_NoNorm_phi.exphist")
    
    
    gg = extract_multidimentional_data(file, 2,0,nomalize)
   
    gg.make_exphist_data(out_dir+pdb_input_name+"_NoNorm_psi.exphist")


def exphist_file_read( data_p,two_col=0):
    data= open(data_p,"r").readlines() 
    all_data=[]
    for i in data:
        d = i.split()
        all_data.append(d)     
    all_data=np.array(all_data).astype(float)
    if all_data.shape[1]==1:
        x = np.array(range(len(all_data[:,0])))*step_size
    else:
        x = all_data[:,0]
        
    data_np = all_data[:,-1]*627.509    
    #data_np = data_np - np.min(data_np)   
    y = data_np 

    return x,y     

def exphist_file_read_hartree( data_p):
    data= open(data_p,"r").readlines() 
    all_data=[]
    for i in data:
        d = i.split()
        all_data.append(d)     
    all_data=np.array(all_data).astype(float)
    if all_data.shape[1]==1:
        x = np.array(range(len(all_data[:,0])))*step_size
    else:
        x = all_data[:,0]
        
    data_np = all_data[:,-1]  
    #data_np = data_np - np.min(data_np)   
    y = data_np 

    return x,y     


def data_from_multiple_exphist_files(files_):
    x_vals =np.array([])
    for fl in files_:
        x,_ = exphist_file_read_hartree(fl)       
        x_vals = np.concatenate((x_vals,x))
    
    x_vals = np.unique(x_vals)
    xx= np.array([])
    yy= np.array([])
    
    for fl in files_:
        x,y = exphist_file_read_hartree(fl)   
        xx=np.concatenate((xx,x))
        yy=np.concatenate((yy,y))
        
    final_out = []    
    xy = np.array([xx,yy]).transpose()
    #print(xy)
    for i in x_vals:
        final_out.append( np.min(xy[xy[:,0]==i,1])*627.509 )
    
    final_out = np.array(final_out) -min(final_out)
    
    return x_vals, final_out
        
    
def make_exphist_data_from_x_y_data(xx, yy, outfile):
    gg = extract_multidimentional_data("", 1,0, norm=1)
    gg.no_modify = 1
    gg.x = xx
    gg.y = yy
    gg.make_exphist_data(outfile)
    