import os

'''
to import this

import sys
sys.path.insert(1, '/home/sudhanshu/bin/my_scripts/my_python_tools')

import ss_utils

'''
from Bio.PDB.Polypeptide import index_to_one
import math
import pandas as pd

from sklearn.linear_model import LinearRegression


def rosetta_score_file_read(scorefile):
    fid = open(scorefile,'r')
    data1 = fid.readline()
    fid.close()
    skip_rows = 1
    if data1.startswith("SCORE: "):
        skip_rows = 0
    
      
    data = pd.read_csv(scorefile,skiprows=skip_rows,delim_whitespace=True)
    return data


def aa_1_letter_code():
    seq=''
    
    for i in range(20):
        seq = seq+ index_to_one(i)
        
    return seq
    

 

def ss_randname( sz ):
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

def ss_char( n_string):
    #Converts integer array to characters
    lens=len(n_string)
    cstring=""
    for i in range(0,lens):
        #cstring=cstring+chr(int(n_string[i]))
        cstring=cstring+str(int(n_string[i]))
    return cstring
    
def ss_exist( nm,typ):
    out=0
    if os.path.exists(nm):
        if typ=='dir':
            if os.path.isdir(nm):
                out=1
        elif typ == 'file':
            if os.path.isfile(nm):
                out=1
                                
    return out
    
def ss_print( data, yon):
    if (yon==1):
        print(data)

def ss_importdata( *inputs):
    data_only_mode = 0;
    if len (inputs) > 1:
        data_only_mode = inputs[1]
    flnm = inputs[0]  
    fid = open(flnm,'r')
    data = fid.readlines()
    
    if data_only_mode > 0:
        if ss_filetype(flnm) == "pdb":
            return extract_pdb_data(data)
        
        elif ss_filetype(flnm) == "xyz":
            return extract_xyz_data(data, data_only_mode)
         
    
    
    
    return data

def ss_isnumber( val):
    val = val.replace('E','e')
    if ((val.count('.')>1) or (val.count('e')>1) or ((val[-1]).isnumeric() == False)):
        return False
    if ((val[0]=='+') | (val[0]=='-')):
        val = val[1:]
    val = val.replace('.','')
    val = val.replace('e-','')
    val = val.replace('e+','')
    val = val.replace('e','')    
    return val.isnumeric()
    
    

def ss_filetype( filenm):
    if filenm.count('.') != 0:
        return filenm[len(filenm)-filenm[:0:-1].index('.'):]
    return "NAN"


def ss_list_all_files_with_ext(dir_path, extn):
    files = os.listdir(dir_path)
    out_list = []
    for  i in files:
        if i.endswith(extn):
            out_list.append(i)
            
    return out_list
            


def extract_pdb_data( all_data):
    return_block =[]
    for i in all_data:
        if len(i) > 6:
            if (i[:4] == 'ATOM') or (i[:6] == 'HETATM'):
                return_block.append(i.strip())
    
    return return_block

def mode8point3(val):
    return "%8.3f" % float(val)

def extract_xyz_data( all_data,types):
    return_block =[]
    for i in all_data:
        i_spl = i.split()
        if len (i_spl) > 3:
            if ((i_spl[0].isalpha()) and (ss_isnumber(i_spl[1])) and 
                (ss_isnumber(i_spl[2])) and (ss_isnumber(i_spl[3]))):
                if types == 2:
                    return_block.append([i_spl[0],mode8point3(i_spl[1])+mode8point3(i_spl[2])+mode8point3(i_spl[3])])
                else:
                    return_block.append(i.strip())
                    

    return return_block
                
    
def ss_distance(F,G):
    if ((len(np.shape(F)) > 1) and (len(np.shape(G)) > 1) ):
        return 0        
    elif ((len(np.shape(F)) > 1) or (len(np.shape(G)) > 1) ):
        
        dim = 1
    else:
        dim = 0
    distns = np.sqrt(np.sum(np.multiply(F-G,F-G),dim))        
    return distns 

import numpy as np
import matplotlib.pyplot as plt

class PAGS_energy_from_gaussian_file():
    def __init__(self,filen):
        self.filen = filen
        self.abc=[]
        self.d=[]
        self.phi_treat = 0
        self.step_angle = 10
        self.remap_360 = 1
        self.print_me = 1
        self.label_on = 0
        self.plot = True
        self.link_list = []
        self.function_type = ""
        self.link_identifier=""
        
    
    def run_me(self, colors=""):
        self.read_chi_file(self.filen)        
        self.plot_me(colors)
        
    
    def gaus1(self,x,a,b,c):
        return a*np.exp((-1*(x-b)**2)/c)    
    
    #def gaus1(self, x,a,b,c):
    #    #return a*np.cos(np.deg2rad(x-b))
    #    return a*np.exp(c*np.cos(np.deg2rad(x-b)))/(np.pi*2*c) 
    
    def gausx(self,x,abc):
        count = len(abc)//3
        #print(abc)
        val = 0
        for i in range(count): 
            #print(abc[i*3],abc[i*3+1],abc[i*3+2])
            val = val+ self.gaus1(x,abc[i*3],abc[i*3+1],abc[i*3+2])
        return val
    
    
    def read_chi_file(self,flnm=""):
        if flnm == "":
            flnm = self.filen
        
        fid = open(flnm,"r")
        data = fid.readlines();
        fid.close()
        req=[]
        
        for i in data:
            if len(i) > 5:
                if (i.startswith('a') or i.startswith('b') or  i.startswith('c') or i.startswith('d')):
                    curr_line=[]
                    d2 = i.split()[1:]
                    
                    for j in d2:
                        if j[0]=="#":
                            break
                        curr_line.append(float(j))
                    
                    req.append(curr_line)
                    
                    #print(req)
                elif i.startswith('LINKID'):
                    link_list = []
                    for j in i.split():
                        if j.startswith('s'):
                            link_list.append(j)
                    
                    self.link_identifier = i.split()[1]
                    
                    self.link_list = link_list
                elif i.startswith('FUNCTION'):
                    self.function_type = i.split()[1]
                            
                
        
        self.abc = np.array(req[:-1])
        self.d=req[-1]
        
        #print(data)
        
        
        
    def plot_me(self,colors):
        x_s = range(0,360+self.step_angle,self.step_angle)
        if self.phi_treat == 1:
            x_s = range(-180,181 ,self.step_angle)
        y_s=[]
        
        for i in x_s:
            y_v = self.gausx(i, (self.abc.transpose().reshape(1,-1)[0]))
            #print(i,y_v)
            y_s.append(y_v+self.d)
        
        
        y_s = np.concatenate((y_s))
        y_s = list(y_s)
        
        # print(y_s)
        # print(x_s)
        
        if self.remap_360 == 1:
            x_s2 = range(0,360+self.step_angle,self.step_angle)
            x_s = np.array(x_s)
            y_s = np.array(y_s)
            ids= np.where(x_s<0)
            y_s = np.concatenate( (  y_s[np.where(x_s>=0)],y_s[ids[0][1:]]))
            y_s = np.concatenate( (y_s,[y_s[0]]))
            #print(y_s)
            x_s = list(x_s2)
            y_s = list(y_s)            
        
        labels = ''
        if self.label_on == 0:
            labels='_nolegend_'
            
        if self.plot == True:    
            if colors == "":
                plt.plot(x_s,y_s, label=labels)
            else:
                plt.plot(x_s,y_s, label=labels,color=colors)
                
        #plt.plot(x_s,y_s)
        self.x_s = x_s
        self.y_s = y_s
        
        if self.print_me == 1:
            for i,j in zip(x_s,y_s): 
                print(i,j/627.509)
        #self.write_dist(x_s,y_s)

    def value_for_angle(self, angle):
        y_v = self.gausx(angle, (self.abc.transpose().reshape(1,-1)[0]))
            #print(i,y_v)
        
        return y_v+self.d
        



class pdb_functions:
    def __init__(self):
        self.pdbfl = ""
        self.pdb_data = []
        self.xyz_data = []
        self.atm_array = []
        self.info = ''' some function to read and write pdbfiles and
        converstion to different formats.'''
        self.pdb_data_continuous = []
    def _pdb_file_(self,pdbf):
        self.pdbfl = pdbf
        
    def _pdb_splitter_(self):       
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
                    
                    if atmi.count(" ") > 0:
                        if atmi.endswith("A"):
                            new_atom_i = atmi.split(" ")[0]
                            atmi = new_atom_i
                    
                    resi = l[17:20].strip()
                    chain = l[21].strip()
                    resno = int(l[22:26].strip())
                    xi = float(l[30:38].strip())
                    yi = float(l[38:46].strip())
                    zi = float(l[46:54].strip())
                    atm_id = l[77:].strip()
                    out.append([srn,atmi,resi,chain,resno,xi,yi,zi,atm_id])
        
        return out
    
    def dummy_pdb_data_CB(self,size):
        data =[]
        for i in range(size):
            data.append([i+1,'CB','ALA','A',i+1,0,0,0,'C'])
        
        self.dummy_pdb_data = data
    
    
    
    def _pdb_continuous_residues_(self):
        """Removes chain change so residues number are different"""
        current_res = self.pdb_data[0][4]
        counter = 1
        pdb_data_continuous=[]
        for i in self.pdb_data:
            if not i[4] == current_res:
                counter +=1
                current_res=i[4]
        
            pdb_data_continuous.append(i[0:4]+[counter]+i[5:])
            
        self.pdb_data_continuous = pdb_data_continuous
                
    def renumber_atoms(self):
        for i, line in enumerate(self.pdb_data):
            self.pdb_data[i][0] = i+1
            
            
    def renumber_residues(self):
        prev_res = self.pdb_data[0][4]
        curr_res = 1
        
        for i in range(len(self.pdb_data)):
            if self.pdb_data[i][4] == prev_res:
                self.pdb_data[i][4] = curr_res
            else:
                prev_res = self.pdb_data[i][4]
                curr_res +=1
                self.pdb_data[i][4] = curr_res
                
           
    
    def read_pdb(self,pdbf):
        self.pdbfl = pdbf
        self.pdb_data = self._pdb_splitter_() 
        self._pdb_continuous_residues_()
        self._xyz_data_()
        self._atm_array_()

    def _xyz_data_(self):         
        xyz_d = []
        for i in self.pdb_data:
            xyz_d.append(i[5:8])        
        self.xyz_data = xyz_d
        
    def _atm_array_(self):
        self.atm_array = [i[1] for i in self.pdb_data]
    
    
   
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
            strn='ATOM  '+repr(en[counter][0]).rjust(5)+atmn_nm+""+en[counter][2].ljust(3)+en[counter][3].rjust(2)+repr(en[counter][4]).rjust(4)+"    "+('%.3f' % i[0]).rjust(8)+('%.3f' % i[1]).rjust(8)+('%.3f' % i[2]).rjust(8)+" "*22+en[counter][8].rjust(2)+"\n"
        # print strn
            fid.write(strn)
        fid.write("END\n")
        fid.close()
        
    def write_pdb (self,pname):        
        xyz_cord = self.xyz_data
        self.pdb_write_from_xyz(xyz_cord,pname)
        
    def save_pdb(self,pname):
        self.write_pdb(pname)
        
    def dump_pdb(self,pname):
        self.write_pdb(pname)
        
    def pdb_write(self,pname):
        self.write_pdb(pname)
        
    def coord_of_atom_from_xyz_data(self, atom_name):
        counter = 0
        for i in self.pdb_data:
            if i[1]==atom_name:
                return self.xyz_data[counter]
            counter +=1
        return 0
    
    def index_number_of_atom(self, atom_name):
        counter = 0
        for i in self.pdb_data:
            if i[1]==atom_name:
                return counter
            counter +=1
        return -1
    
    def remove_line_x(self, line_num):
        print("line number starts from 0")
        if line_num == -1:
            return
        vac_pdb = []
        for i in range(len(self.pdb_data)):
            if i == line_num :
                continue
            vac_pdb.append(self.pdb_data[i])
            
        self.pdb_data = vac_pdb
        del vac_pdb
        self._xyz_data_()
        
    def refresh(self):
        self._xyz_data_()
        self._atm_array_()
        #renumbder_atom had to add
        
    def use_continuous_data(self):
        self.pdb_data = self.pdb_data_continuous
        
    def refresh_from_pdb_data(self):
        self.refresh()
        
    
    def refresh_from_xyz(self):
        for i in range(len(self.xyz_data)):
            self.pdb_data[i][5] = self.xyz_data[i][0]
            self.pdb_data[i][6] = self.xyz_data[i][1]
            self.pdb_data[i][7] = self.xyz_data[i][2]
            
        self._atm_array_()
            
    def rename_atom(self, old_name, new_name):
        for i in range(len(self.pdb_data)):
            if self.pdb_data[i][1] == old_name:
                self.pdb_data[i][1] = new_name
                return
            
    def rename_atom_type (self, index_id, new_type):
        self.pdb_data[index_id][8] = new_type
        
    def add_new_atom(self, xyz_cord, name , type_atom):
        self.xyz_data = np.append(self.xyz_data, [xyz_cord], axis=0)
        #print(self.xyz_data)
        len_pdb_data = len(self.pdb_data)
        prev_pdb_line = self.pdb_data[-1]
        self.pdb_data.append([ prev_pdb_line[0]+1, name , prev_pdb_line[2],
                              prev_pdb_line[3],  prev_pdb_line[4],
                              xyz_cord[0], xyz_cord[1], xyz_cord[2], type_atom] )
        
        self._atm_array_()
            
            
        
    def residue_data(self,residue_num):
        return_arr = []
        for i in self.pdb_data:
            if i[4] == residue_num:
                return_arr.append(i)
        return return_arr
        
    def coord_of_atom_of_residue(self, atom_name,residue_num):
        data = self.residue_data(residue_num)
        for i in data:
           # print(i)
            if i[1] == atom_name:
                return np.array(i[5:8])
            







def extract_CB(pdb_file, outfile=""):
    
    #pdb_file = "/home/sudhanshu/HDD2/projects2/utils/PDB_from_Morgan/MLN_pyranose_benchmark_set/only_protein/3OEB_protein.pdb"
    if outfile=="":
        outfile = pdb_file[:-4]+"_CB.pdb"
    
    pdb_f = pdb_functions()
    pdb_f.read_pdb(pdb_file)  
    pdb_f.use_continuous_data()
    
    cb_pdb_data = []
    
    for i in range(50000):
        d = pdb_f.residue_data(i+1)
        CB = 'CB'
        if len(d) == 0:
            break
        
        #print (CB)
        for j in d:
            if j[2] == 'GLY':
                CB = 'CA'
            
            if j[1] == CB:
                cb_pdb_data.append(j)
                break
            
    
    pdb_f.pdb_data = cb_pdb_data
    pdb_f.refresh_from_pdb_data()
    pdb_f.write_pdb(outfile)
    






class geometry_adjustment:
    def __init__(self):
        self.xyz = []
        self.atm_index = []
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
            
            #print(curr_atm)
            
            inter_atomic_dist = 1.68
            if len(self.atm_index ) > 0:
                if self.atm_index[curr_atm].startswith('H'):
                    inter_atomic_dist = 1.2
                    
                #print(self.atm_index[curr_atm])
                
                
            
            
            connections = np.where( (distances_from_current>0) &  (distances_from_current<inter_atomic_dist))[0]
            connections = connections [connections != prev_atm]  #####
            #print("connections", np.array(connections)+1)
            
            for i in connections:
                if up_side.count(i) == 0:
                    up_side.append(i)
            
            check_done.append(curr_atm) 
        self.rotation_atom_index[torsion_number] = up_side
        print (up_side)
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
        
        if len(self.rotation_atom_index) == 0:
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






class link_id_handler():
    def __init__(self):
        self.links_path = "/home/sudhanshu/HDD2/rosetta/Rosetta2/Rosetta/main/database/scoring/score_functions/carbohydrates/"
        self.param_files = []
        self.link_details=[]
        self._load_all_link_details()
          
    def _load_all_link_details(self):
        files = os.listdir(self.links_path)
        param_files = []
        for i in files:
            if i.endswith('.params'):
                param_files.append(i)   
                
        link_details=[]
        for i in param_files:
            handle = CHI_energy_from_gaussian_file(self.links_path+i)  
            handle.read_chi_file()
            link_details.append([i,handle.link_list, handle.function_type, handle.link_identifier])
        self.param_files = param_files
        self.link_details = link_details
        
        
    def if_link_ids_match (self,link_id1, link_id2):
       
        if (len(link_id1) != len(link_id2) ):
            return False
        
        qm_at_1=[i for  i,j in enumerate(link_id1) if j=="?"]
        qm_at_2=[i for  i,j in enumerate(link_id2) if j=="?"]
        all_qm = qm_at_1 +qm_at_2
        
        link_mod1 = ""
        link_mod2 = ""
        
        for i in range(len(link_id1)):
            if all_qm.count(i) == 0:
                link_mod1 = link_mod1 + link_id1[i]
                link_mod2 = link_mod2 + link_id2[i]
            else:
                link_mod1 = link_mod1 + "?"
                link_mod2 = link_mod2 + "?"
        return link_mod1 == link_mod2
    
    
    def find_link_id_matching_files(self, link_id):
        out_file_names = []
        for i in self.link_details:
            for j in i[1]:
                if self.if_link_ids_match(link_id, j):
                    out_file_names.append([i[0],i[3]])
                    
        return out_file_names




def rsquare_data(dx,dy): #d11 : experminetal data X data
    d11 = np.copy(dx)
    d22 = np.copy(dy)
    nan_vals = d11+d22    
    d1 = d11[~np.isnan(nan_vals)]
    d2 = d22[~np.isnan(nan_vals)]    
    d1 = d1.reshape((-1,1))    
    model = LinearRegression()
    model.fit(d1, d2)
    r_sq = model.score(d1, d2)
    slope = model.coef_
    y_pred = model.predict(d1)  
    
    d1_min_pos = np.where(d1==min(d1))[0][0]
    d1_max_pos = np.where(d1==max(d1))[0][0]
    # print(d1_max_pos)
    
    x_points = np.array([np.min(d1), np.max(d1)])
    y_points = np.array([y_pred[d1_min_pos], y_pred[d1_max_pos]])
    
    xy_out = [x_points, y_points]
    
    #print(r_sq)
    
    return r_sq,slope,xy_out,~np.isnan(nan_vals)




def normal_vector_of_the_plane(xx,yy,zz):
    points = np.array([xx,yy,zz])
    # now find the best-fitting plane for the test points
    
    # subtract out the centroid and take the SVD
    svd = np.linalg.svd(points - np.mean(points, axis=1, keepdims=True))
    
    # Extract the left singular vectors
    left = svd[0]
    return left[:, -1]



def angle_to_plane(xyz_atom, acceptor_C_cord, xyz_ring):
    xx= xyz_ring[:,0]
    yy= xyz_ring[:,1]
    zz= xyz_ring[:,2]
    
    nv_plane = normal_vector_of_the_plane(xx,yy,zz)
    
    atom_cord = xyz_atom
    
    unit_vector = np.array(acceptor_C_cord) - np.array(atom_cord)
    unit_vector = unit_vector / np.sqrt(np.sum(unit_vector**2))
        
    dot_product = np.dot(unit_vector, nv_plane)
    angle =  np.arccos(dot_product)*180/np.pi
    angle = min(angle , 180 -angle)
    
    return angle
        


def aa_parameters_all():

    aa_seq = aa_1_letter_code()
    
    aa_parameters = {}
    
    
    in_file = "/home/sudhanshu/HDD3/Data/CAPSIF_data/aa_properties/aa1.csv"
    
    data =pd.read_csv(in_file,delimiter=",")
    
    correct_aa_seq = [np.where(data.aa1 == i)[0][0] for i in aa_seq ]
    
    data = pd.DataFrame(data, index= correct_aa_seq)
    data = data.reset_index(drop=True)
    
    # normalization in 0 to 1 range
    for i in data.columns[2:]:
        data[i] = (data[i] - min(data[i]))/max(data[i])
    
    
    # for i in data.columns
    aa_parameters['info'] = ['Hydropathy', 'radius', 'Aromaphilicity', 'Hbond_D', 'Hbond_A']  
    
    
    
    for i in range(20):
        vector =[]
        for j in ["Hydropathy" ,"Volume(A3)" , "Aromaphilicity","H_bond_Doner","H_bond_Acceptor"]:
            val = data[j][i]
            if j == "Volume(A3)":
                val = val**(1/3)
                    
            vector.append(val)
        
        aa_parameters[data.aa3[i]] = vector
     
    aa_parameters['MSE'] = aa_parameters['MET']
    return aa_parameters
     
class scores:
    pass

def TM_score(pdb1,pdb2):
    tm_align_exe = "/home/sudhanshu/bin/my_scripts/TMalign"
    data = os.popen(tm_align_exe + " "+pdb1 + " "+pdb2).read()
    data = data.splitlines()
    out = scores()
    counter = 0
    for i in data:
        if i.startswith('Aligned'):
            isplit = i.split()
            out.rmsd= float(isplit[4][:-1])
            
        elif i.startswith('TM-score'):
            counter += 1
            isplit = i.split()
            if counter == 1:
                out.TM_score1 = float(isplit[1])
            elif counter == 2:
                out.TM_score2 = float(isplit[1])
            

        if counter == 2:
            break
            
    return out
    