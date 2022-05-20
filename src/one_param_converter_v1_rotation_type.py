#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 15:46:23 2020
V 2.0 10_March_2020 :: ERRORR 
@author: sudhanshu
"""
import os
import shutil
from string import *
import numpy as np
import pandas as pd
import subprocess
from src.rotation_functions_v1 import *
from settings import global_settings


class my_pipeline():
    def __init__(self):
        self.tors = [["C4'1","C5'1","O5'1","C1'2"],]#["C5'1","O5'1","C1'2","C2'2"],]
        #tors = [["C3'1","C4'1","C5'1","O5'1"]]
        self.tors = [[22,1,2,3]]
        self.step_size = 45
        self.input_pdb = './ad.pdb'
        self.pdb_for_ensemble = self.input_pdb
        self.orca = global_settings.orca_path
        self.ensemble_dir = global_settings.output_ensemble_dir
        self.ring_torsions = []
        self.extra_torsions= []
        self.extra_torsion_angles = []
        self.project_count = 1000 ## change this for specific number
            
    def ss_exist(self,nm,typ):
        out=False
        if os.path.exists(nm):
            if typ=='dir':
                if os.path.isdir(nm):
                    out=True
            elif typ == 'file':
                if os.path.isfile(nm):
                    out=True
                                    
        return out
    
    def make_project_dir(self):
        working_dir = self.ensemble_dir
        if not self.ss_exist(working_dir,"dir"):
            os.mkdir(working_dir)
        
        project_count = self.project_count
        project_init = "project_"
        
        while 1:
            project_name = working_dir + "/" + project_init + ("%4d" % project_count).replace(" ","0")
            if self.ss_exist(project_name,"dir"):
                project_count +=1
            else:
                break
            
        os.mkdir(project_name)
        os.mkdir(project_name + "/pdbs")
        os.mkdir(project_name + "/QM_inps")
        return project_name
    
    
    def torsion_from_string2(self,tors_string):
        out_v = tors_string.replace(" ","")
        out_v = out_v.split("],[")
        out_v2 = [ ]
        print(out_v)
        for i in out_v:
            i = i.replace("],]","]]")
            i = i.replace("[","")
            i = i.replace("]","")        
            i = i.split(",")
            out_arr = []
            for ii in i:                
                # check is int each value
                if ii == '':
                    print("No torsion given")
                
                elif not ii=='NA':
                    out_arr.append(int(ii))
                else:
                    out_arr.append((ii))
        
            out_v2.append(out_arr)        
        
        return(out_v2)
    
    
    def torsion_from_string(self, tors_string):
        tors_arr = self.torsion_from_string2(tors_string) 
        self.tors = tors_arr
        print("Torsions",self.tors)
        
    def extra_torsion_from_string(self, extra_tors_string):
        tors_arr = self.torsion_from_string2(extra_tors_string)
        if not tors_arr == [[0,0,0,0]]:
            self.extra_torsions = tors_arr
            
    def extra_torsion_angles_from_string(self, extra_tors_angles):
        ang_arr = self.torsion_from_string2(extra_tors_angles)
        if not ang_arr ==[[0]]:
            self.extra_torsion_angles = ang_arr
        
    
    def pdb_splitter(self,pdbf):       
        '''
        
        # 
        #  1 -  6        Record name     "ATOM  "  
        #  7 - 11        Integer         Atom serial number.  
        # 13 - 16        Atom            Atom name.   
        # 17             Character       Alternate location indicator.
        # 18 - 20        Residue name    Residue name.  
        # 22             Character       Chain identifier.   
        # 23 - 26        Integer         Residue sequence number. 
        # 27             AChar           Code for insertion of residues.  
        # 31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms. 
        # 39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms. 
        # 47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.  	
        
        '''                   
        fl=open(pdbf,'r')
        data=fl.readlines()	
        fl.close()
        out=[];
        atm_name=[ ]
        xyz_all = [ ]
        for l in data:
            if len(l)>10:
                if l[0:4]=='ATOM' or l[0:6]=='HETATM':
                    srn=int(l[6:11].strip())
                    atmi=l[12:17].strip()
                    resi=l[17:20].strip()
                    chain=l[21].strip()
                    resno=int(l[22:26].strip())
                    xi=float(l[30:38].strip())
                    yi=float(l[38:46].strip())
                    zi=float(l[46:54].strip())
                    out.append([srn,atmi,resi,chain,resno,xi,yi,zi])
                    atm_name.append(atmi)
                    xyz_all.append([xi,yi,zi])
        return out,atm_name,xyz_all
    
    def euc_distance(self,a1,a2):
        return np.sqrt(np.sum((a1-a2)**2,1))
    
    
    
    def renamed_atom_indexer(self,pdbfl1, pdbfl2):
        data_fl1 = self.pdb_splitter(pdbfl1)
        data_fl2 = self.pdb_splitter(pdbfl2)
        xyz1 = np.array(data_fl1[2])
        xyz2 = np.array(data_fl2[2])
        return_dict = dict()
        
        euc_dict_array=[]
        for i in range(len(xyz1)):
            xyz_r1 = xyz1[i]
            euc_dists = self.euc_distance(xyz2, xyz_r1)
            indexes = np.where(euc_dists<=0.01)
            if len(indexes[0]) > 1:
                print ("clash present!")           
            euc_dict_array.append(indexes[0][0])
            return_dict[data_fl1[1][i]+str(data_fl1[0][i][4])]=data_fl2[1][indexes[0][0]]
            #print(data_fl1[1][indexes[0][0]])
        return return_dict
            
        
        #print(euc_dict_array)
            
    
    def renamed_atom_indexer2(self, pdbfl1, pose):
        data_fl1 = self.pdb_splitter(pdbfl1)
        xyz1 = np.array(data_fl1[2])
        
        atm_nm_from_pose = []
        xyz2_from_pose = []
        
        for i in range(pose.residue(1).natoms()):
            xyz2_from_pose.append( pose.conformation().residue(1).xyz(i+1) ) 
            atm_nm_from_pose.append( pose.conformation().residue(1).atom_name(i+1).strip()  )
        
        xyz2 = np.array(xyz2_from_pose)
        #print(atm_nm_from_pose)
        return_dict = dict()
        
        euc_dict_array=[]
        for i in range(len(xyz1)):
            xyz_r1 = xyz1[i]
            euc_dists = self.euc_distance(xyz2, xyz_r1)
            indexes = np.where(euc_dists<=0.1)
            if len(indexes[0]) > 1:
                print ("clash present!")           
            euc_dict_array.append(indexes[0][0])
            return_dict[data_fl1[1][i]+str(data_fl1[0][i][4])] = atm_nm_from_pose[indexes[0][0]]
        self.atm_nm_from_pose = atm_nm_from_pose
        return return_dict

    
    def one_working_dir_make(self):
        if not self.ss_exist(self.input_pdb, 'file'):
            print("PDB FILE DOES NOT EXIST")
            
        self.dir_name = self.make_project_dir()
        print(self.dir_name)
        print("CTPPM: Working Directory made.")
        return "Working Directory made."
    
    def adjust_extra_torsions(self):
        self.pdb_for_ensemble = self.input_pdb
        
        if (len(self.extra_torsions) > 0) and (len(self.extra_torsion_angles) > 0): 
            tmp_pdb_name = self.input_pdb[:-4]+"_temp.pdb"
            
            pdb_handler1 = pdb_functions()
            pdb_handler1.read_pdb(self.input_pdb)
            pdb_handler1.write_pdb(tmp_pdb_name)
            
            pdb_handler2 = pdb_functions()
           # gg2 = geometry_adjustment()
            # print("CURRENT TORS:", self.extra_torsions)
            # print("CURRENT ANG:", self.extra_torsion_angles)
            
            for i in range(len(self.extra_torsions)):
                curr_tors = self.extra_torsions[i]
                curr_ang = self.extra_torsion_angles[i]
                # print("CURR TORS:", curr_tors)
                # print("CURR ANG:", curr_ang)
                #print("CURRENT ANG ", curr_ang)
                if curr_ang == ['NA'] :
                    continue

                pdb_handler2 = pdb_functions()
                pdb_handler2.read_pdb(tmp_pdb_name)                
                gg2 = geometry_adjustment()
                gg2.input_xyz(pdb_handler1.xyz_data)
                gg2.input_torsion(curr_tors)
                gg2.detect_rotation_atoms()
                gg2.input_xyz(pdb_handler2.xyz_data)
                gg2.set_torsion(curr_ang[0])
                pdb_handler1.pdb_write_from_xyz(gg2.xyz,tmp_pdb_name)
            
            del gg2
            del pdb_handler1
            del pdb_handler2
            self.pdb_for_ensemble = tmp_pdb_name
                
                
    def two_ensemble_maker_old(self):   
        self.adjust_extra_torsions()
        input_pdb = self.pdb_for_ensemble
        pdb_for_quality_check = self.input_pdb
        tors = self.tors
        step_size = self.step_size
        dir_name = self.dir_name  + "/pdbs"
        
        num_torsions = len(self.tors)
        
        pdb_file = os.path.join(os.getcwd(), input_pdb)
        pdb_quality_check = os.path.join(os.getcwd(), pdb_for_quality_check)
        pdb_handler = pdb_functions()
        pdb_handler.read_pdb(pdb_file)
        
        ga = geometry_adjustment_2()
        ga.input_xyz(pdb_handler.xyz_data)

        qc = quality_check()
        qc.load_first_pdb(pdb_quality_check)
        
        
        self.score_file = dir_name + "/clash_score.dat"
        fid = open(self.score_file,"w+")
        
        if num_torsions == 0:
            out_pdb = dir_name+"/ensemble_0" +".pdb"
            
        
        elif num_torsions == 1:
            ga.input_torsion(tors[0])
            ga.detect_rotation_atoms(0)
            

            for i in range(0,360,step_size):
                ga.set_torsion(i)
                out_pdb = dir_name+"/ensemble_"+str(i) +".pdb"
                pdb_handler.pdb_write_from_xyz(ga.xyz,out_pdb)
                score = 1
                if qc.is_bond_changed_from_first_pdb(out_pdb):
                    score = 0

                fid.write( ("%2d" % score) + "  " + "ensemble_"+str(i) +".pdb\n")
            fid.close()
        elif num_torsions == 2:
            ga.input_torsion(tors[0])
            ga.detect_rotation_atoms(0) 
            ga.input_torsion(tors[1])
            ga.detect_rotation_atoms(1) 
            
            for i in range(0,360,step_size):
                ga.input_torsion(tors[0])
                ga.set_torsion(i,0)
                
                for j in range(0,360,step_size):                    
                    ga.input_torsion(tors[1])
                    #ga.detect_rotation_atoms()   
                    ga.set_torsion(j,0)                   
                
                    out_pdb = dir_name+"/ensemble_"+str(i)+ "_"+str(j) +".pdb"
                    pdb_handler.pdb_write_from_xyz(ga.xyz,out_pdb)
                    score = 1
                    if qc.is_bond_changed_from_first_pdb(out_pdb):
                        score = 0
        
            
                    fid.write( ("%2d" % score) + "  " + "ensemble_"+str(i)+"_"+str(j) +".pdb\n")
            fid.close()
            
            
        h = amber_tools()     
        h.run_analysis(pdb_file)      
        frt = find_ring_torsions(pdb_handler.atm_array, h._atom_type_dict_, h.get_bond_information()) 
        self.ring_torsions = frt.identify_ring_torsions(1)

        
        file_data = pd.read_fwf(self.score_file,header=None)
       # self.tors_idx_array = tors_idx_array
        
        self.write_README()
        
        print("CTPPM: Ensemble for given torsion is generated")
        return ["Ensemble for given torsion is generated", file_data]
 
    
    def three_filter_good_ensemble_points(self):
        dir_name = self.dir_name         
        df = pd.read_csv(self.score_file,header=None, delimiter=r"\s+")

        df2=df.loc[df.index[df[0]>0]][1].tolist()
     
        os.mkdir(self.dir_name+"/pdbs/use_ensemble")     
        number_of_torsions = len(self.tors)
        
        print("Number Torsions: ", number_of_torsions)
        
        print("CTPPM: Filtered Structures for Quantum Calculation")
        if number_of_torsions == 0:
            shutil.copyfile(dir_name+"/pdbs/ensemble_0.pdb",dir_name+"/pdbs/use_ensemble/ensemble_"+str(i)+".pdb")
            
 
        elif number_of_torsions == 1:        
            for i in df2:
                shutil.copyfile( dir_name+"/pdbs/"+i, dir_name+"/pdbs/use_ensemble/" + i)
                
 
        elif number_of_torsions == 2:
            
            for i in df2:
                shutil.copyfile( dir_name+"/pdbs/"+i, dir_name+"/pdbs/use_ensemble/" + i)
        
            # for i in range(0, 360, self.step_size): 
            #     for j in range(0, 360, self.step_size): 
            #         shutil.copyfile( dir_name+"/ensemble_"+str(i)+"_"+str(j)+".pdb",
            #                         dir_name+"/use_ensemble/ensemble_"+str(i)+"_"+str(j)+".pdb")
        else:
            print("Torsion count:", number_of_torsions, "is not supported")
                    
    
    def clean_atom_type(self, name):
        repl1 = "0123456789.'*"
        repl2 = '"'
        for i in repl1:
            name = name.replace(i,'')
     
        name = name.replace('"','')
        if len(name) == 1:  
            return name
        else:
            return name[0]
        
    
    def inp_file_header(self, torsion_angle):
         
        input_v = ""
        
        for torses in torsion_angle:    
            input_v_tmp = ""
            for i in torses:
                # nw_name = self.atm_dict[i]
                # position = self.atm_nm_from_pose.index(nw_name)
                input_v_tmp = input_v_tmp+ str(i-1) + "  "
                
            input_v = input_v + "{ D "+ input_v_tmp + " C } # D for Dihedral angle\n"
            
            
        
        if len(self.ring_torsions) > 0:
            for rts in self.ring_torsions:
                input_v_tmp = ""
                for i in rts:
                    # nw_name = self.atm_dict[i]
                    # position = self.atm_nm_from_pose.index(nw_name)
                    input_v_tmp = input_v_tmp+ str(i-1) + "  "
                    
                input_v = input_v +"{ D "+ input_v_tmp + " C } # D for Dihedral angle\n"

                    #xyz_fid.writelines(self.inp_file_header(rts))
        
        if len(self.extra_torsions) > 0:
            for ets in self.extra_torsions:
                input_v_tmp = ""
                for i in ets:
                    # nw_name = self.atm_dict[i]
                    # position = self.atm_nm_from_pose.index(nw_name)
                    input_v_tmp = input_v_tmp+ str(i-1) + "  "
                    
                input_v = input_v +"{ D "+ input_v_tmp + " C } # D for Dihedral angle\n"
            
        
        
        return '''# Test redundant internal optimization
#
! RKS B3LYP  Opt TightSCF SmallPrint

%basis
        Basis "6-31G"
end

%geom
Constraints
''' + input_v + '''end
end


#%pal 
#nprocs 8
#end

* xyz 0 1
'''
    def run_subprocess(self, command):
        out = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout,stderr = out.communicate()
        return stdout,stderr
    
    def four_create_quantum_inputs(self):
        dir_name = self.dir_name + "/pdbs/use_ensemble/"
        all_files = os.listdir(dir_name)
        
        # Making inp file.
        for i in all_files:
            if i.startswith("ensemble_") and i.endswith("pdb"):
                angles_array = np.array(i.replace("ensemble_","").replace(".pdb","").split("_")).astype(float)
                
                
                #angle_v = float(i.replace("ensemble_","").replace(".pdb","")) + 1 # NEW 
                
                pdb_file = dir_name + i
                pdb_data = self.pdb_splitter(pdb_file)
                xyz_fid = open( pdb_file[:-3]+"inp",'w+')
                
                input_tors=[]
                for num, torsionz in enumerate(self.tors):
                    #print (torsionz)
                    input_tors_temp= list(torsionz)
                    input_tors_temp.append(angles_array[num]+1)
                    input_tors.append(input_tors_temp)
                
                xyz_fid.writelines(self.inp_file_header(input_tors))
                

                    
                
                
                for j in range(len(pdb_data[1])):
                    print_line = (self.clean_atom_type(pdb_data[1][j]) +
                                  " " + ("%8.3f" % pdb_data[2][j][0]) + 
                                  " " + ("%8.3f" % pdb_data[2][j][1]) + 
                                  " " + ("%8.3f" % pdb_data[2][j][2]) + "\n")
                    xyz_fid.write(print_line)
                xyz_fid.write("*\n")
                xyz_fid.close()
                
        # submitting to orca
        inp_files = [ ]
        all_files = os.listdir(dir_name)
        for i in all_files:
            if (i.startswith("ensemble_") and i.endswith("inp")):
                inp_files.append(i)
                
         
        # for i in inp_files[:2]:   
        #     inpfile = dir_name + i
        #     os.system("bin/tools/run_orca.sh" + " " + inpfile+"&" )
        #     #subprocess.Popen([self.orca, inpfile, '>', inpfile+".out", "&"]) 
        #     # print (str(msg))
        
        dir_end = [i for i in range(len(self.dir_name)) if self.dir_name[i] =="_"][-1]+1
        quantum_run_dir = self.dir_name + "/QM_inps"
        #os.mkdir(quantum_run_dir)
        for fl in inp_files:
            shutil.copyfile(dir_name + fl, quantum_run_dir +"/" + fl)
            
        ## Entry in runfile for qunatum calculation on jazz
        
        fid = open(self.ensemble_dir + "p_" + self.dir_name[dir_end:] +".runconfig", 'w+' )
        fid.write("submit\n")
        fid.write("project_" + self.dir_name[dir_end:]+"\n")
        fid.close()
        shutil.copyfile(self.dir_name +"/README", quantum_run_dir +"/README")
            
        
        
        
        return "submitted"
        
    def write_README(self):
        fid2 = open(self.dir_name+"/README","w+")
        fid2.write("Used torsion:   "+ str(self.tors) + "\n")
        fid2.write("Used step size: "+ str(self.step_size) + "\n")
        fid2.write("Input pdb file: "+self.input_pdb + "\n")
        fid2.write("Orca exe:       "+self.orca + "\n")
        fid2.write("Ensemble output:"+self.ensemble_dir + "\n")
        fid2.close()


        
        
        

def run_pipeline(pdb_fl, torsion_v, step_sz, extra_tors, extra_tors_angs, project_count = 1000):
    gg = my_pipeline()
    gg.step_size = step_sz
    gg.project_count = project_count
    gg.input_pdb = pdb_fl
    gg.torsion_from_string(torsion_v)
    gg.extra_torsion_from_string(extra_tors) 
    gg.extra_torsion_angles_from_string(extra_tors_angs) 
    #gg.tors = torsion_v
    gg.one_working_dir_make()
    gg.two_ensemble_maker_old()
    gg.three_filter_good_ensemble_points()
    gg.four_create_quantum_inputs()
    

