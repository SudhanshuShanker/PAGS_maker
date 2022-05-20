#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 11:35:08 2020
V 1.0 27_Jan_2020
V 1.1 10_March_2020
@author: sudhanshu
"""

import numpy as np
#from string import *
import math
import os
import shutil
from settings import global_settings

### PDB Functions


    
class pdb_functions:
    def __init__(self):
        self.pdbfl = ""
        self.pdb_data = []
        self.xyz_data = []
        self.atm_array =[]
        self.atm_type = []
        self.info = ''' some function to read and write pdbfiles and
        converstion to different formats.'''
        
    def pdb_file(self,pdbf):
        self.pdbfl = pdbf
        
    def pdb_splitter(self):       
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
        fl=open(self.pdbfl,'r')
        data=fl.readlines()	
        fl.close()
        out=[];
        for l in data:
            if len(l)>5:
                if l[0:4]=='ATOM' or l[0:6]=='HETATM' :
                    srn = int(l[6:11].strip())
                    atmi = l[12:17].strip()
                    resi = l[17:20].strip()
                    chain = l[21].strip()
                    resno = int(l[22:26].strip())
                    xi = float(l[30:38].strip())
                    yi = float(l[38:46].strip())
                    zi = float(l[46:54].strip())
                    atm_id = l[77:].strip()
                    out.append([srn,atmi,resi,chain,resno,xi,yi,zi,atm_id])
        
        return out
    
    def read_pdb(self,pdbf):
        self.pdbfl = pdbf
        self.pdb_data = self.pdb_splitter() 
        self._xyz_data_()
        self._atm_array_()

    def _xyz_data_(self):         
        xyz_d = []
        for i in self.pdb_data:
            xyz_d.append(i[5:8])        
        self.xyz_data = xyz_d
        
    def _atm_array_(self):
        self.atm_array = [i[1] for i in self.pdb_data]
        self.atm_type = [i[1][0] for i in self.pdb_data]
    
    
   
    def import_pdb_models(pdbfl):
        fidmp=open(pdbfl,'r')
        data_mp=fidmp.readlines()
        fidmp.close()
        coord_mp=[]
        tmp_frame=[]
        for i in data_mp:
            if i.find('ENDMDL')>-1:
                # print tmp_frame
                coord_mp.append(tmp_frame)
                tmp_frame=[]
            if len(i)>5:
                if i[0:4]=='ATOM':
                    # srn=int(strip(l[6:11]))
                    # atmi=strip(l[12:17])
                    # resi=strip(l[17:20])
                    # chain=strip(l[21])
                    # resno=int(strip(l[22:26]))
                    xi=float(i[30:38].strip())
                    yi=float(i[38:46].strip())
                    zi=float(i[46:54].strip())
                    tmp_frame.append([xi,yi,zi])
        return coord_mp
         
   

    def pdb2atmdata(self,atms,pdbfile):
        fl=open(pdbfile,'r')
        data=fl.readlines()	
        fl.close()
        out=[];
        for line in data:
            spl_line=line.split()
            if len(spl_line)>6:
                atmind=2
            else:
                atmind=1	
            if spl_line[0]== 'ATOM':
                if spl_line[atmind]==atms:
                    if len(spl_line)>6:
                        out.append([float(spl_line[5]),float(spl_line[6]),float(spl_line[7])])
                    else:
                        out.append([float(spl_line[4]),float(spl_line[5]),float(spl_line[6])])
        return np.array(out)

    
    def p_data_string(self,en):
        strn=('ATOM  '+repr(en[0]).rjust(5)+" "+en[1].ljust(4)+" "+en[2].ljust(3)+" "+repr(en[4]).rjust(5) + 
              "    "+repr(en[5]).rjust(8)+repr(en[6]).rjust(8)+repr(en[7]).rjust(8)+"                     "+en[8].rjust(2)+"\n")
        return strn
    
    def pdb_write_from_xyz(self, xyz_cord, pdb_file_name):
        en = self.pdb_data
        counter = -1
        fid=open(pdb_file_name,'w+')
        for i in xyz_cord:
            counter=counter+1
            en_nm=en[counter][1]+" "
            if len(en_nm)>4:
                justi=5
                atmn_nm=" "+en_nm.ljust(justi)
            else:
                justi=4
                atmn_nm="  "+en_nm.ljust(justi)
            strn='ATOM  '+repr(en[counter][0]).rjust(5)+atmn_nm+""+en[counter][2].ljust(3)+" "+repr(en[counter][4]).rjust(5)+"    "+('%.3f' % i[0]).rjust(8)+('%.3f' % i[1]).rjust(8)+('%.3f' % i[2]).rjust(8)+" "*22+en[counter][8].rjust(2)+"\n"
        # print strn
            fid.write(strn)
        fid.write("END\n")
        fid.close()
        
    def write_pdb (self,pname):        
        xyz_cord = self.xyz_data
        self.pdb_write_from_xyz(xyz_cord,pname)
            
class geometry_adjustment_2:
    def __init__(self):
        self.xyz = []
        self.torsion_select = [] #pdb numbering
        self.rotation_step = 10
        self.rotation_atom_index = dict()
        self.info = ''' This class takes a pdb file as input,
        with a torsion angle [atom1 , atom2, atom3, atom4 ] as 
        integer and starts from 1.
        
        It rotates given pdb by given angle at defined torsion.  
        '''
        

    def input_torsion(self,t_array):
        if len(t_array) == 0:
            return
        if isinstance(t_array[0],int):       
            self.torsion_select = t_array
        
            
        
    def input_xyz(self,xyz_in):
        self.xyz = np.array(xyz_in)
        
    def detect_rotation_atoms(self, torsion_number = 0):
        trial_input = self.torsion_select #[22,1,2,3] #PDB numbering
        xyz_data = self.xyz  
        
        curr_atm = trial_input[2] - 1
        prev_atm = trial_input[1] - 1
        
        check_done = [ ]    
        up_side = []
        
        start_flag = 0
        atm_counter = -1
        
        while (1):
            if start_flag == 0 :
                curr_xyz = xyz_data[curr_atm]
                prev_xyz = xyz_data[prev_atm]
                start_flag = 1
            else:
                atm_counter += 1
                if atm_counter >= len(up_side):
                    break
                curr_atm = up_side[atm_counter]
                curr_xyz = xyz_data[curr_atm]
                    
                if check_done.count(curr_atm)>0:
                    continue
            
            
            distances_from_current = self.distn(xyz_data, curr_xyz)
            connections = np.where( (distances_from_current>0) &  (distances_from_current<1.7))[0]
            connections = connections [connections != prev_atm]  #####
            #print("connections", np.array(connections)+1)
            
            for i in connections:
                if up_side.count(i) == 0:
                    up_side.append(i)
            
            check_done.append(curr_atm) 
        self.rotation_atom_index[torsion_number] = up_side
        #return up_side  # 0 index numbering
            
    def distn(self, F,G):
        if ((len(np.shape(F)) > 1) and (len(np.shape(G)) > 1) ):
            return 0        
        elif ((len(np.shape(F)) > 1) or (len(np.shape(G)) > 1) ):
            
            dim = 1
        else:
            dim = 0
        distns = np.sqrt(np.sum(np.multiply(F-G,F-G),dim))        
        return distns 
    
    def localNewAngle(self, x,y):
        ang = 0; # This is the default value. In case y == 0 and x ~= 0, ang = 0.
        if (x != 0) & (y!= 0):
            c = x/np.sqrt(x*x + y*y)
            ang = np.sign(y)*np.arccos(c)
        elif (x == 0):
            if (y > 0):
                ang = np.pi/2
            elif (y < 0):
                ang = -np.pi/2
        return ang
    
    
    def calc_current_torsion(self):
        A = self.xyz[self.torsion_select[0] - 1 ]
        B = self.xyz[self.torsion_select[1] - 1]
        C = self.xyz[self.torsion_select[2] - 1]
        D = self.xyz[self.torsion_select[3] - 1]
        a=(A-B)/self.distn(A,B)
        b=(C-B)/self.distn(C,B)
        c=(D-B)/self.distn(D,B)
        b_cross_c = np.cross(b,c)
        x = -np.sum(np.conj(a)*c) + np.multiply([np.sum(np.conj(a)*b)],[np.sum(np.conj(b)*c)])
        y = np.sum(np.conj(a)*b_cross_c)
        out=180-np.rad2deg(self.localNewAngle(x,y))
        return 360 - float(out)
    

    
    def req_tors_ang(self,t_ang):
        orig_t=self.calc_current_torsion()
        req_ang=orig_t-t_ang
        return req_ang
    


    def rotation_mat(self,theta,axis):
        pi=math.pi
        theta=theta*pi/180
    #     line_cords=np.array(line_cords)
    #     origin=line_cords[1,:]
    #     vector=vector-origin
    #     axis=line_cords[0,:]-origin
    #     axis = axis/np.sqrt(np.dot(axis,axis))
        a = np.cos(theta/2)
        b,c,d = -axis*np.sin(theta/2)
        rot_mat=np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                         [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                         [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])
                         
        #new_vector=np.dot(rot_mat,vector)   
        return (rot_mat)#+origin)   

   
    
    def rotate_molecule_by(self,alfa,torsion_id=0):
        #data in xyz format, index of axis atoms [2,3],rotatioan starting atom index [4],angles[30]
        
        if len(self.rotation_atom_index[torsion_id]) == 0:
            self.detect_rotation_atoms()   
        rotation_xyz = np.copy(self.xyz[self.rotation_atom_index[torsion_id]][:])   
        
        #data=np.copy(data2)
        #for i in range(np.shape(rot_init)[0]):
        
        #line_cords= data[axisi[i]]
        origin = self.xyz[self.torsion_select[1] - 1 ]
        axist= self.xyz[self.torsion_select[2] - 1] - origin
    
        # origin=line_cords[0,:]  #O5'
        # axist=line_cords[1,:]-origin #P
        axist = axist/np.sqrt(np.dot(axist,axist))
        rot_m = self.rotation_mat(alfa,axist)
        t_data = rotation_xyz-origin
        for j in range (np.shape(t_data)[0]):
            c_data=t_data[j]
            t_data[j]=np.dot(rot_m,c_data)           
        self.xyz[self.rotation_atom_index[torsion_id]]=t_data +origin 
        
    def rotate_molecule_to(self,angle, torsion_id=0 ):
        alpha = self.req_tors_ang(angle)
        self.rotate_molecule_by(alpha, torsion_id)
        return 0
        
    def get_torsion(self):
        return self.calc_current_torsion()
    
    def set_torsion(self,angle,torsion_id=0):
        self.rotate_molecule_to(angle,torsion_id)  
        
class ss_functions:    
    def __init__(self):
        self.__name__ = "sudhanshu's functions"
        self._name_len_ = 6
        
    def ss_randname(self, sz):
        s=list(range(48,58))
        s.extend([95])
        s.extend(list(range(65,91)))
        s.extend(list(range(97,123)))	
        rv=np.random.rand(sz)
        rv=rv*len(s)
        rv =list(rv)
        rname=""
        for i in range(0,len(rv)):
            rname=rname+chr(s[int(rv[i])])
        return(rname)
    
    def ss_char(self,n_string):
        lens=len(n_string)
        cstring=""
        for i in range(0,lens):
            cstring=cstring+chr(int(n_string[i]))
        return cstring
        
    def ss_exist(self,nm,typ):
        out=0
        if os.path.exists(nm):
            if typ=='dir':
                if os.path.isdir(nm):
                    out = True
            elif typ == 'file':
                if os.path.isfile(nm):
                    out = False
                                    
        return out
    
    def ss_make_random_dir(self,dir_path):
        while(1):
            trial_name = "ss_" + self.ss_randname(self._name_len_)
            trial_dir_name = dir_path + "/" + trial_name
            if not self.ss_exist(trial_dir_name ,"dir" ):
                break
            
        os.mkdir(trial_dir_name)
        return trial_dir_name
    
    def ss_distance(self, F,G):
        if ((len(np.shape(F)) > 1) and (len(np.shape(G)) > 1) ):
            return 0        
        elif ((len(np.shape(F)) > 1) or (len(np.shape(G)) > 1) ):
            
            dim = 1
        else:
            dim = 0
        distns = np.sqrt(np.sum(np.multiply(F-G,F-G),dim))        
        return distns 
        
class amber_tools:
    def __init__(self):
        self._atom_type_dict_ = dict()
        self.gaff_bond_type_dict = dict()
        self.AMBERHOME = global_settings.AMBERSETTINGS.AMBERHOME
        self.param_dir = global_settings.AMBERSETTINGS.param_dir
        self.temp_dir = global_settings.AMBERSETTINGS.temp_dir
        self._temp_working_dir_ = ""
        self._pdb_file_name_ = ""
        self._bond_info_ = []
        
        self._load_gaff_params_()
        
        self.ssf = ss_functions()
        self._temp_working_dir_ = self.ssf.ss_make_random_dir(self.temp_dir)
        os.environ['AMBERHOME'] = self.AMBERHOME
        self.info = ''' Reads pdb  and finds amber atomtypes and bonds using AMBERTOOLS '''
        
        
        
    def _load_gaff_params_(self):
        file = self.param_dir +"/gaff.dat"
        #fid = open(file,"r")
        
        with open(file) as fp:
            for line in fp:                
                tmp_l = line.split("  ")
                if len(tmp_l) > 4:
                    if tmp_l[0].count("-") == 1:
                        if tmp_l[3].count(".") == 1:
                            self.gaff_bond_type_dict[tmp_l[0]] = float(tmp_l[3])
                            
        self.gaff_bond_type_dict['o -h '] = 0.973 #ho-oh
        self.gaff_bond_type_dict['h -h '] = 0.740 #ho-oh        
                            
    def clear_all(self):
        if len(self._temp_working_dir_) > len(self.temp_dir)+1:
            if self._temp_working_dir_.count(self.temp_dir) == 1: 
                shutil.rmtree(self._temp_working_dir_)
      
    def input_pdb_file(self,pdbfl):
        self._pdb_file_name_ = pdbfl         # if not same than previous clean atm and bond info        
        
    def run_antechamber(self):
        os.system(self.AMBERHOME + "bin/antechamber -pf y -fi pdb -fo mol2 -i " + self._pdb_file_name_ + " -o " +
                  self._temp_working_dir_ + "/temp.mol2")
        
        
    def _read_mol2_(self, is_first=1):
        fl = self._temp_working_dir_ + "/temp.mol2"
        
        atmnm_dict = dict()
        bond_info = []
        
        
        if os.path.exists(fl):
            atom_flag = 0
            bond_flag = 0
            
            with open(fl) as fp:
                for line in fp:
                    tmp_l = line.split()
                    if len(tmp_l) > 0:
                        if tmp_l[0] == "@<TRIPOS>ATOM":
                            atom_flag = 1
                            continue
                            
                        elif tmp_l[0] == "@<TRIPOS>BOND":
                            atom_flag = 0
                            bond_flag = 1
                            continue
                        elif tmp_l[0] == "@<TRIPOS>SUBSTRUCTURE":
                            bond_flag = 0
                            continue
                            
                        if atom_flag == 1:
                            atmnm_dict[tmp_l[1]] = tmp_l[5]
                            
                        elif bond_flag == 1:
                            bond_info.append([int(tmp_l[1]),int(tmp_l[2])])
                            
            #print(atmnm_dict)
            bond_info = np.array(bond_info)
            bond_info = bond_info[bond_info[:,1].argsort()]
            bond_info = bond_info[bond_info[:,0].argsort()]
        if is_first == 1:
            self._atom_type_dict_ = atmnm_dict   
            self._bond_info_ = bond_info
        else:
            return (atmnm_dict , bond_info)
                    
   
    def run_analysis(self, pdbfl, is_first=1):
        if len(self._pdb_file_name_) == 0:
            self.input_pdb_file(pdbfl)
        self.run_antechamber()
        if is_first == 1:
            self._read_mol2_(is_first)
            self.clear_all()
        else:
            a,b = self._read_mol2_(is_first)
            self.clear_all()
            return a,b
        
            
    def get_bond_information(self):
        return self._bond_info_
        
    def get_atm_types(self,atm_name):
        return self._atom_type_dict_[atm_name]
    
    def get_gaff_bond_length(self, atm1,atm2):
        atm1_opts = [atm1, atm1[0], atm1[0]+"x"]
        atm2_opts = [atm2, atm2[0], atm2[0]+"x"]
        for i in atm1_opts:
            for j in atm2_opts:
                key1 = i.ljust(2) + "-" + j.ljust(2)
                key2 = j.ljust(2) + "-" + i.ljust(2)
                if key1 in self.gaff_bond_type_dict.keys():
                    return self.gaff_bond_type_dict[key1]
                elif key2 in self.gaff_bond_type_dict.keys():
                    return self.gaff_bond_type_dict[key2]
        return -1  
   
    
   
class quality_check:
    def __init__(self):
        self.amber_tools = amber_tools()
        self._first_pdb_bonds_ = []
        self._first_atom_type_ = dict()
        self._pdb_atms_ = []
        self._if_first_pdb_done_ = 0   
        self.ssf = ss_functions()
    
    def load_first_pdb(self,pdbfl):
        if self._if_first_pdb_done_ == 0:
            tmp_amb = amber_tools()
            tmp_amb.run_analysis(pdbfl)
            self._first_pdb_bonds_ = tmp_amb.get_bond_information()
            self._first_atom_type_ = tmp_amb._atom_type_dict_
            self._if_first_pdb_done_ = 1
            del tmp_amb
            
            tmp_pdb_hand = pdb_functions()
            tmp_pdb_hand.read_pdb(pdbfl)
            self._pdb_atms_ = tmp_pdb_hand.atm_array
            del tmp_pdb_hand
        else:
            print("already calculated values for parent pdb")


    def is_bond_changed_from_first_pdb(self,pdbfl2):
        tmp_amb2 = amber_tools()
        atm_dict, bond_info = tmp_amb2.run_analysis(pdbfl2,0)
        del tmp_amb2
        if not len(self._first_pdb_bonds_) ==  len(bond_info):
            return 1
        
        if not self._first_atom_type_ == atm_dict:
            return 1
        
        for i in bond_info:
            if not sum(np.sum(self._first_pdb_bonds_==i,1) ==2 ) == 1:
                return 1

        return 0
   
        
    def check_if_this_bond_exist(self,bond):
        if sum(np.sum(self._first_pdb_bonds_==bond,1)==2):
            return True
        bond.reverse()
        if sum(np.sum(self._first_pdb_bonds_==bond,1)==2):
            return True
        return False
        
    
    
    
    
class find_ring_torsions:
    def __init__(self,atms, gaff_type_atm_dict, bonds):
        self.atm_arr = atms
        self.bond_info = bonds
        self.gaff_atm_dict = gaff_type_atm_dict
        self.gaff_type_atm_arr = self.conv_atms_to_gaff_array()
        self.rings = []
     
    def conv_atms_to_gaff_array(self):
        out = []
        for i in self.atm_arr:
            out.append(self.gaff_atm_dict[i])
        return out
            
    def identify_rings(self):
        
        
        usable_bonds = []
        for i in self.bond_info:
            if (self.gaff_type_atm_arr[i[0]-1][0] != 'h' ) & (self.gaff_type_atm_arr[i[1]-1][0] != 'h' ):
                usable_bonds.append(i)
                
        usable_bonds = np.array(usable_bonds)
        
        connect_paths = []
        start_v = usable_bonds[0]
        #usable_bonds = np.delete(usable_bonds,0,0)        
        #print(usable_bonds)

        connect = np.array(usable_bonds)
        connect2 = np.array(usable_bonds)
        one_finder = []

        num_nodes =len( set(np.reshape(connect,-1)))
        
        # Identifying possible rings        
        num_rings = 0        
        if num_nodes-1 ==len(connect):
            num_rings = 0
        else:
            num_rings = len(connect) - num_nodes +1
        
            
        # removing heavy atom terminal bonds : single bonds
        for i in range(np.max(connect)):
            count = sum(sum(connect==i+1))
            #print(i+1,count)
            if count == 1:
                index_v = np.where(np.sum(connect == i+1,1))[0][0]
                one_finder.append(i+1)
        
        one_finder.sort()
        print(connect2)
        one_finder.reverse()
        #print(one_finder)
        for i in one_finder:
            indx_del = np.sum(connect2 == i,1)
            indx_del = np.where(indx_del)[0][0]
            connect2 = np.delete(connect2,indx_del,0)
         
       
        
        #print("coneft2",connect2)
            
        #return connect2
        # Identifying ring start atoms.
        count_keeper = []
        ring_start_atms = []
        for i in range(np.max(connect2)):
            count = sum(sum(connect2==i+1))
            count_keeper.append(count)
            if count > 2:
                ring_start_atms.append(i+1)
        
        
        max_carb_ring_size = 9
        #print(ring_start_atms)
        
        if len(ring_start_atms) == 0:
            if num_rings == 1:
                if count_keeper.count(2) == len(count_keeper):
                    ring_start_atms.append(connect2[0][0])
        
        
        all_rings = [ ]
        for i in ring_start_atms:
            tmp_v = self.find_ring_from_start_atm(connect2,i)
            for j in tmp_v:
                all_rings.append(j)
         
        self.rings = all_rings
            
            #for j in all_connects:
      
    def identify_ring_torsions(self,format_ = 1):
        if len(self.rings)==0:
            self.identify_rings()
        return_v = [] 
        for i in self.rings:
            temp_v = i[:-1]*2
            
            if format_ == 0:
                temp_v = [i-1 for i in temp_v]
            
            
            for j in range(len(i)-1):
                return_v.append(temp_v[j:j+4])
        
        return (return_v)
            
    
            
    def all_connectivity_atm(self,atmnum,bonds):
        out = []
        for i in bonds:
            if np.sum(i==atmnum):
                out.append(i[i!=atmnum][0])
        
        return out
                
            
    def find_ring_from_start_atm(self, bonds, start_atm = -1 ):
        #BETA version: may contain some errors
        # Convert it to full tree
        
        not_all = 1
        if start_atm == -1:
            start_atm = 1
            not_all = 0
        
        paths = [[start_atm]]
        curr_atm = start_atm
        done_atm = []
        trial_atms=[]
        rings = []
        
        
        
        while_check = 1
        while (while_check == 1) :    
            
            all_conct_curr_atm = self.all_connectivity_atm(curr_atm, bonds)  
            
            counts_done = 0
            #print (done_atm )
    
            for i in paths:            
                if i[-1] == curr_atm:     #end path atom           
                    for j in all_conct_curr_atm: 
                        # initialize error handle
                        if len(i) > 1:    
                            # ignore back propogation
                            if i[-2] == j:                            
                                continue
                        new_path = i + [j]
                        
                        # kill if all paths are found
                        
                        #print(new_path, paths.count(new_path))
                        if paths.count(new_path)>0:
                            #return rings #paths for more path
                            while_check = 0
                            continue
                        
                        # check repeats 
                        # crpts = 0
                        # for rv2 in new_path:
                        #     if not_all == 1:
                        #         if rv2 == start_atm:
                        #             continue
                            
                        #     if new_path.count(rv2)>1:
                        #         crpts = 1
                        #         break
                            
                        #if crpts == 0:
                        paths.append(new_path)
                        
                        # Ring identification by atom repeat 
                        for rv in new_path:
                            
                            if not_all == 1:
                                if rv != start_atm:
                                    continue
                            
                            if new_path.count(rv)>1:
                                rings.append(  new_path[new_path.index(rv):] )
                                break
                        
                        # atoms for next trial
                        trial_atms.append(j)
            # ignoring treated atoms
            if done_atm.count(curr_atm) > 0:
                    while_check = 0
            done_atm.append(curr_atm)
            #print(paths)
            # Selection of next current atom        
            for ii in trial_atms:
                if done_atm.count(ii)>0:
                    continue
                curr_atm = ii
                
        #print(paths)             
        return rings
            
        
            
            

        
# h= amber_tools()     
# h.run_analysis("af.pdb")        
 
# pdb_handler = pdb_functions()
# pdb_handler.read_pdb('bin/ad.pdb')


# gg = geometry_adjustment()
# gg.input_xyz(pdb_handler.xyz_data)
# gg.input_torsion([22,1,2,3])
# gg.detect_rotation_atoms()


# qc = quality_check()
# qc.load_first_pdb('bin/ad.pdb')

# frt = find_ring_torsions(pdb_handler.atm_array, h._atom_type_dict_, h.get_bond_information()) 


# trial_dir = "/home/sudhanshu/Desktop/projects/pipeline/my_web/bin/working_dir/project_9"
# all_file = os.listdir(trial_dir)

# all_pdb = []
# for i in all_file:
#     if i.startswith("ensemble_") & i.endswith(".pdb"):
#         all_pdb.append(i)

# for i in all_pdb:
#     fl = trial_dir + "/" +i
#     if qc.find_new_possible_bonds(fl):
#         print("clash present", fl)
    



# for i in range(0,360,45):
#     gg.set_torsion(i)
#     out_pdb = "check_"+str(i)+".pdb"
#     pdb_handler.pdb_write_from_xyz(gg.xyz,out_pdb)
#     if qc.find_new_possible_bonds(out_pdb):
#         print("clash present", out_pdb)

    
        




        
        
# def torsion_connectivity(torsion_input, xyz_data):
#     trial_input = torsion_input #[22,1,2,3] #PDB numbering
#     xyz_data = np.array(xyz_data)    
    
#     curr_atm = trial_input[2] - 1
#     prev_atm = trial_input[1] - 1
    
#     check_done = [ ]    
#     up_side = []
    
#     start_flag = 0
#     atm_counter = -1
    
#     while (1):
#         if start_flag == 0 :
#             curr_xyz = xyz_data[curr_atm]
#             prev_xyz = xyz_data[prev_atm]
#             start_flag = 1
#         else:
#             atm_counter += 1
#             if atm_counter >= len(up_side):
#                 break
#             curr_atm = up_side[atm_counter]
#             curr_xyz = xyz_data[curr_atm]
                
#             # print (np.array(check_done)+1)  
#             # print(curr_atm+1)
#             if check_done.count(curr_atm)>0:
#                 continue
        
        
#         distances_from_current = distn(xyz_data, curr_xyz)
#         connections = np.where( (distances_from_current>0) &  (distances_from_current<1.7))[0]
#         connections = connections [connections != prev_atm]  #####
#         #print("connections", np.array(connections)+1)
        
#         for i in connections:
#             if up_side.count(i) == 0:
#                 up_side.append(i)
        
#         check_done.append(curr_atm) 
        
#         #print("upside", np.array(up_side)+1)
        
#         # out_atm =""
#         # for j in up_side:
#         #     out_atm = out_atm + "," + np.array(atms)[j]
        
#         #print(out_atm)
        
        
#         # print (check_done)
#         # print(up_side)
#     return up_side  # 0 index numbering
        



