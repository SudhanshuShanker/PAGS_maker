
import numpy as np
import math

from src.ss_utils import  pdb_functions, geometry_adjustment, ss_distance
import copy
import os


class atomic_data:
    def __init__(self):
        self.distance = 0
        self.angle = 0
 
# Omega 0  
COC= atomic_data()
COC.distance = 1.437
COC.angle = 114.35

# for omega 1
COCC = atomic_data()
COCC.distance_CO = 1.437
COCC.distance_OC = 1.444
COCC.distance_CC = 1.553
COCC.angle_COC = 115.163
COCC.angle_OCC = 111.364
COCC.angle_COH = 109.65


class Data:
    pass


def settings():
    handle = Data()
    # location of component pdbs to make simple models
    handle.component_dir = "./bin/pdb_components/"       
    return handle



def pdb_component_data(connection_data):
    pdb_components = dict()
    
    component_data = Data()
    component_data.pdb = "PYR_b3lyp_1C4.pdb"    
    pdb_components['PYRANOSE_1C4'] = component_data 
    
    component_data = Data()
    component_data.pdb = "PYR_b3lyp_4C1.pdb"    
    pdb_components['PYRANOSE_4C1'] = component_data 
    
    component_data = Data()
    component_data.pdb = 'furanose_north_b3lyp.pdb' 
    pdb_components['FURANOSE_NORTH'] = component_data
    
    component_data = Data()
    component_data.pdb = 'furanose_south_b3lyp.pdb'
    pdb_components['FURANOSE_SOUTH'] = component_data
    
    component_data = Data()
    component_data.pdb = 'COC.pdb'
    pdb_components['LINK_OMEGA0'] = component_data
    
    component_data = Data()
    component_data.pdb = 'OCC.pdb'
    pdb_components['LINK_OMEGA1'] = component_data
    
    
    
    # first res
    res1 = connection_data.first_res    
    static_head_atom = connection_data.first_hyd
    
    if static_head_atom[-1] == "1":
        static_head_atom = static_head_atom[:-1] + "2"
    else:
        static_head_atom = static_head_atom[:-1] + "1"
      
    pdb_components[res1].static_head_atom = static_head_atom
    pdb_components[res1].static_carbon = connection_data.first_carbon
    
    
    #last residue
    res2 = connection_data.second_res    
    movable_tail_atom = connection_data.second_hyd
    
    if movable_tail_atom[-1] == "1":
        movable_tail_atom = movable_tail_atom[:-1] + "2"
    else:
        movable_tail_atom = movable_tail_atom[:-1] + "1"
      
    pdb_components[res2].movable_tail_atom = movable_tail_atom
    pdb_components[res2].movable_carbon = connection_data.second_carbon    
    inbetween_seg = connection_data.omega
    
    if inbetween_seg[-1] == "1":
        pdb_components[inbetween_seg].movable_tail_atom = "CY"
        pdb_components[inbetween_seg].movable_carbon = "CY"
        pdb_components[inbetween_seg].static_carbon = "CY"
        pdb_components[inbetween_seg].head_atom = "CX"
        pdb_components[inbetween_seg].static_head_atom = "CY"
        pdb_components[inbetween_seg].head_type = "C"
        
    return pdb_components




def connection_decoder(connect_str):
 
    #path fpr P => P4
    if connect_str.count('P') == 1:
        if connect_str.count('P1') +  connect_str.count('P4') == 0:
            connect_str = connect_str.replace('P','P4')
            
    elif connect_str.count('P') == 2:
        if connect_str.count('P1') +  connect_str.count('P4') == 0:
            connect_str = connect_str.replace('P','P4')
        elif connect_str.count('P1') +  connect_str.count('P4') == 1:
            print("problem with the name P4,P1,and P combination is wrong!")
            return 0
            
      
    print(connect_str)
    
    data = connect_str.replace("to",'')
    out_d = [data[i] for i in range(len(data)) ]


    expand_dict = dict()
    #expand_dict['P'] = 'PYRANOSE_4C1'
    expand_dict['P1'] = 'PYRANOSE_1C4'
    expand_dict['P4'] = 'PYRANOSE_4C1'
    expand_dict['FN'] = 'FURANOSE_NORTH'
    expand_dict['FS'] = 'FURANOSE_SOUTH'

    out_x = []
    prev_i = ""
    
    output = Data()
    output.status = True  # if there is any error it will become False
    output._stage_ = -1  # for error handing, it shows the block
    output._errormsg_ = "" # for error handing, it shows the problem
    
    for i in out_d:
        if i =='P':
            prev_i = i            
            continue
        if i == 'F':
            prev_i = i
            continue
        if prev_i == "P":  # stage 1
            if ((i == '1') or (i == '4')):
                out_x.append(expand_dict[prev_i+i])
                prev_i = ""
                continue
            else:
                out_x.append("NAME_ERROR")
                prev_i = ""
                output.status = False
                output._stage_ = 1
                output._errormsg_ = "pyranose can be either P1: 1P4 or P4: 41P"
                return output
                continue
        
        if prev_i == "F":  # stage 2
            if ((i == 'N') or (i == 'S')):
                out_x.append(expand_dict[prev_i+i])
                prev_i = ""
                continue
            else:
                out_x.append("NAME_ERROR")
                prev_i = ""
                output.status = False
                output._stage_ = 2
                output._errormsg_ = "furanose can be either FN: F-north or FS: F-south"
                return output
                continue
        out_x.append(i)
        
    print (out_x)
            
    if ( out_x.count('NAME_ERROR') > 0 ): # stage 3
        output.status = False
        output._stage_ = 3
        print("NAME_ERROR")
        return output
       
    first_carbon = "C"+out_x[2]
    second_carbon = "C"+out_x[5]  
    new_oxygen = "O" + out_x[5]  
    first_hyd = ""
    second_hyd = ""
    

    if out_x[1] == 'A':
        first_hyd = "H" + out_x[2] + "1"
        
    elif  out_x[1] == 'E':
        first_hyd = "H" + out_x[2] + "2"
        
    # print(out_d)
    if out_x[4] == 'A':
        second_hyd = "H" + out_x[5] + "1"
        
    elif  out_x[4] == 'E':
        second_hyd = "H" + out_x[5] + "2"
    
    #print(new_oxygen,first_hyd, second_hyd)
    
    
    output.first_carbon = first_carbon
    output.second_carbon = second_carbon
    output.first_hyd = first_hyd
    output.second_hyd = second_hyd
    output.new_oxygen = new_oxygen
    output.first_res = out_x[0]
    output.second_res = out_x[3]
    output.omega = "LINK_OMEGA"+out_x[7]
    output.distance_interaction = +1
    output.head_atom_type = "O"
    output.original_second_carbon = second_carbon
    
    
    
    #omega related special treatment
    
    if out_x[7] == "1":
        if output.second_res.startswith("FUR"):
            if output.second_carbon == "C6":
                output.second_carbon = "C5"
                output.second_hyd = "H5" + output.second_hyd[-1]
                output.gt_carbon = "C4"
            elif output.second_carbon == "C1":
                output.second_carbon = "C2"
                output.second_hyd = "H2" + output.second_hyd[-1]
                output.gt_carbon = "C3"
                
        if output.second_res.startswith("PYR"):
            if output.second_carbon == "C6":
                output.second_carbon = "C5"
                output.second_hyd = "H5" + output.second_hyd[-1]
                output.gt_carbon = "C4"
            
       
    
    return output
    #return [first_carbon, first_hyd, second_carbon, second_hyd, new_oxygen, ]
    


class rotate_translate_molecule:
    def __init__(self, mol_handle, axis_rotation, angle_of_rotation):
        self.coord = mol_handle.xyz_data
        self.axis_rotation = np.array(axis_rotation)
        self.angle_of_rotation = angle_of_rotation
        self.mol_handle = mol_handle
        #self.rotation_mat(self.angle_of_rotation, self.axis_rotation)
        self.rotate_molecule_by()
        
        
    def rotation_mat(self,theta,axis):
        pi=math.pi
        theta=theta*pi/180
        a = np.cos(theta/2)
        bcd = -axis*np.sin(theta/2)
        #print(bcd)
        b= bcd[0]
        c = bcd[1]
        d = bcd[2]
        rot_mat=np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                         [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                         [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])
                         
        #new_vector=np.dot(rot_mat,vector)   
        return (rot_mat)#+origin)   

   
    
    def rotate_molecule_by(self):              
        origin = self.axis_rotation[0]              
        axist= self.axis_rotation[1] - origin
        rotation_xyz = np.copy(self.coord[:])   
        axist = axist/np.sqrt(np.dot(axist,axist))
        #print("acc", axist)
        
        rot_m = self.rotation_mat( self.angle_of_rotation, axist)
        t_data = rotation_xyz-origin
        for j in range (np.shape(t_data)[0]):
            c_data=t_data[j]
            t_data[j]=np.dot(rot_m,c_data)  
            
        self.xyz=t_data +origin 
        
        self.mol_handle.xyz_data_rotated = self.xyz
        

class rotate_translate_atom:
    def __init__(self, cord1, cord2, cord3, axis_rotation, angle_of_rotation):
        self.cord1 = cord1
        self.cord2 = cord2
        self.cord3 = cord3
        self.axis_rotation = np.array(axis_rotation)
        self.angle_of_rotation = angle_of_rotation        
        #self.rotation_mat(self.angle_of_rotation, self.axis_rotation)
        self.rotate_atom_by()
        
        
    def rotation_mat(self,theta,axis):
        pi=math.pi
        theta=theta*pi/180
        a = np.cos(theta/2)
        bcd = -axis*np.sin(theta/2)
        #print(bcd)
        b= bcd[0]
        c = bcd[1]
        d = bcd[2]
        rot_mat=np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                         [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                         [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])
                         
        #new_vector=np.dot(rot_mat,vector)   
        return (rot_mat)#+origin)   

   
    
    def rotate_atom_by(self):      
        
        origin = self.axis_rotation[0]              
        axist= self.axis_rotation[1] - origin
        rotation_xyz = np.copy(self.cord3[:])           

        axist = axist/np.sqrt(np.dot(axist,axist))
        
        #print("acc", axist)
        
        rot_m = self.rotation_mat( self.angle_of_rotation, axist)
        t_data = np.copy(rotation_xyz-origin)
        #for j in range (np.shape(t_data)[0]):
        c_data=t_data
        t_data=np.dot(rot_m,c_data)  
            
        self.xyz=t_data +origin 
        
        #self.mol_handle.xyz_data_rotated = self.xyz
        




def simple_sd_minimize_distance(target_atom_cord, data_handler, iteration_atom_id, axis_for_rot, distance_interaction = 1 ):
    direction = +1
    angle_val = 1
    
    prev_distance = ss_distance( target_atom_cord, data_handler.xyz_data_rotated[iteration_atom_id] )
    
    print("optimizing distance")
    
    direction_change_count = 0
    for i in range(1000):
        rotate_translate_molecule(data_handler, axis_for_rot, direction*angle_val)
        
        new_distance = ss_distance( target_atom_cord, data_handler.xyz_data_rotated[iteration_atom_id] )
        
        
        
        if (distance_interaction*new_distance) < (distance_interaction*prev_distance):
            data_handler.xyz_data = np.copy(data_handler.xyz_data_rotated[:])
            prev_distance = new_distance
            
            #print(new_distance, direction)
            angle_val = angle_val*2
        else:
            data_handler.xyz_data_rotated = np.copy(data_handler.xyz_data)
            direction = direction*-1
            direction_change_count +=1
            angle_val = angle_val*0.2
        if direction_change_count > 20:
            break
    
    print("optimization done!")
        
        
def angle_calculate(cord1,cord2,cord3):
    cord1 = np.array(cord1)
    cord2 = np.array(cord2)
    cord3 = np.array(cord3)
    vector_1 = cord1 - cord2
    vector_1 = vector_1 /np.sqrt(sum((vector_1)**2))  
    
    vector_2 = cord3 - cord2
    vector_2 = vector_2 /np.sqrt(sum((vector_2)**2)) 
    
    angle = np.arccos(np.dot(vector_1,vector_2))
    return np.rad2deg(angle)

    

def simple_sd_minimize_angle(target_atom_cord, target_atom_coord_midle, final_target_angle, data_handler, iteration_atom_id, axis_for_rot ):
   
    direction = +1
    angle_val = 1
    
    prev_angle_diff = np.abs(angle_calculate(target_atom_cord, target_atom_coord_midle, data_handler.xyz_data_rotated[iteration_atom_id]) - final_target_angle)
    
    print("optimizing angle")   
    direction_change_count = 0
    for i in range(1000):
        rotate_translate_molecule(data_handler, axis_for_rot, direction*angle_val)
        
        new_angle_diff = np.abs(angle_calculate(target_atom_cord, target_atom_coord_midle, data_handler.xyz_data_rotated[iteration_atom_id]) - final_target_angle)
        

        if new_angle_diff < prev_angle_diff :
            data_handler.xyz_data = np.copy(data_handler.xyz_data_rotated[:])
            prev_angle_diff  = new_angle_diff
            #print(prev_angle_diff, direction)
            angle_val = angle_val*2
        else:
            data_handler.xyz_data_rotated = np.copy(data_handler.xyz_data)
            direction = direction*-1
            direction_change_count +=1
            angle_val = angle_val*0.1
        if direction_change_count > 20:
            break        
        #print(i, new_angle_diff, prev_angle_diff, direction)
   
    print("optimization done!")
    


def text_file_reader(flnm):
    fid = open(flnm,"r")
    data = fid.readlines()
    fid.close()
    
    connections = []
    for i in data:
        if len(i.strip()) < 1:
            continue
        if i.strip().startswith("#"):
            continue
        else:
            connections.append(i.strip())
    return connections
    
    

def phi_psi_atom_decoder_from_connection_named_pdb(pdbfl):
    
    connection_name = pdbfl.split("/")[-1][:-4] # pdb file without extension
    connection_name = connection_name.split("_")[0] # pdb file remove extra detail after underscore    
    
    connects = connection_decoder(connection_name)   
    F_ring_arrange = ["O5", "C2", "C3", "C4", "C5", "O5"]
    P_ring_arrange = ["O5", "C1", "C2", "C3", "C4", "C5", "O5"]
    
    first_carbon = connects.first_carbon
    second_carbon = connects.original_second_carbon
    
    phi_psi_atoms = []
    
    if connects.omega[-1] == "0":        
        if connects.first_res.startswith("PYR"):
            index_c = P_ring_arrange.index(first_carbon)
            phi_psi_atoms.append(P_ring_arrange[index_c -1])
            phi_psi_atoms.append(P_ring_arrange[index_c ])
        elif connects.first_res.startswith("FUR"):
            index_c = F_ring_arrange.index(first_carbon)
            phi_psi_atoms.append(F_ring_arrange[index_c -1])
            phi_psi_atoms.append(F_ring_arrange[index_c ])
            
        else:
            phi_psi_atoms.append(0)
            
        phi_psi_atoms.append(second_carbon.replace("C","O"))
            
        if connects.second_res.startswith("PYR"):
            index_c = P_ring_arrange.index(second_carbon)
            phi_psi_atoms.append(P_ring_arrange[index_c ])
            phi_psi_atoms.append(P_ring_arrange[index_c - 1])
        elif connects.second_res.startswith("FUR"):
            index_c = F_ring_arrange.index(second_carbon)
            phi_psi_atoms.append(F_ring_arrange[index_c ])
            phi_psi_atoms.append(F_ring_arrange[index_c - 1 ])
            
        else:
            phi_psi_atoms.append(0)
            
        print(phi_psi_atoms)
            
        
            
        pdb_handle = pdb_functions()
        pdb_handle.read_pdb(pdbfl)
        res_relation = [1,1,2,2,2]
        
        index_arr = []
        counter = 0
        for step,i in enumerate(phi_psi_atoms):
            counter = 0
            temp_count = []            
            for j in pdb_handle.pdb_data:
                counter += 1
                if i == j[1]:
                    temp_count.append(counter)
            
            if len(temp_count) > 1:
                temp_count = [temp_count[res_relation[step]-1]]
            
            index_arr.append([i,temp_count])
                    
        
        #print(index_arr)
        
        
        return [i[1][0] for i in index_arr]
    
    if connects.omega[-1] == "1":   
        F_ring_arrange = ["O5", "C2", "C3", "C4", "C5", "C6", "O6", "O5","C2","C1","O1"]
        P_ring_arrange = [ "O5", "C1", "C2", "C3", "C4", "C5", "O5",  "C5", "C6", "O6",]
        if connects.first_res.startswith("PYR"):
            index_c = P_ring_arrange.index(first_carbon)
            phi_psi_atoms.append(P_ring_arrange[index_c -1])
            phi_psi_atoms.append(P_ring_arrange[index_c ])
        elif connects.first_res.startswith("FUR"):
            index_c = F_ring_arrange.index(first_carbon)
            phi_psi_atoms.append(F_ring_arrange[index_c -1])
            phi_psi_atoms.append(F_ring_arrange[index_c ])
            
        else:
            phi_psi_atoms.append(0)
            
        phi_psi_atoms.append(second_carbon.replace("C","O"))
            
        if connects.second_res.startswith("PYR"):
            index_c = P_ring_arrange.index(second_carbon)
            phi_psi_atoms.append(P_ring_arrange[index_c ])
            phi_psi_atoms.append(P_ring_arrange[index_c - 1])
            phi_psi_atoms.append(P_ring_arrange[index_c - 2])
        elif connects.second_res.startswith("FUR"):
            index_c = F_ring_arrange.index(second_carbon)
            phi_psi_atoms.append(F_ring_arrange[index_c ])
            phi_psi_atoms.append(F_ring_arrange[index_c - 1 ])
            phi_psi_atoms.append(F_ring_arrange[index_c - 2])
            
        else:
            phi_psi_atoms.append(0)
            
        print(phi_psi_atoms)
            
        
            
        pdb_handle = pdb_functions()
        pdb_handle.read_pdb(pdbfl)
        res_relation = [1,1,2,2,2,2]
        
        index_arr = []
        counter = 0
        for step,i in enumerate(phi_psi_atoms):
            counter = 0
            temp_count = []    
            
            for j in pdb_handle.pdb_data:
                
                counter += 1
                if i == j[1]:
                    temp_count.append(counter)
            
            if len(temp_count) > 1:
                temp_count = [temp_count[res_relation[step]-1]]
            
            #print(index_arr)
            index_arr.append([i,temp_count])
                    
        
        #print(index_arr)
        
        
        return [i[1][0] for i in index_arr]
            
            
    return [0]
            
        
        
                    
def generate_molecule_for_connection(connection = "P1A1toFSE4W0", model_dir = "/home/sudhanshu/HDD2/projects2/CHI_potential_data/bin/models_omega1/"):

    #connection = "PA1toFNA6W0"
    setting_ = settings()
    gt_arr = [""]    
    if connection.endswith("W1"):
        gt_arr = ["a",'e']  
        
    
    for gt_val in gt_arr:  
    
        connect_data = connection_decoder(connection)  # identify components required to make molecule
        
        pdb_components = pdb_component_data(connect_data)
        component_seq = [connect_data.first_res, connect_data.omega, connect_data.second_res]
        
        
        component_dir = setting_.component_dir  
        #model_dir = "/home/sudhanshu/HDD2/projects2/CHI_potential_data/bin/models_omega1/"
        
        
        component_handles = []
        
        for  i in component_seq:
            temp_hand = pdb_functions()
            temp_hand.read_pdb(component_dir + pdb_components[i].pdb)
            component_handles.append(temp_hand)
            del temp_hand
            
            
  
        
    # first component no move.  # first residue
        first_handle = component_handles[0]    
        for i in first_handle.pdb_data:
            if i[1] == connect_data.first_carbon:
                c_data = i
                c_cord = np.array(i[5:8])          #coordinates   
            
            if i[1] == connect_data.first_hyd:
                h_data = i
                h_cord = np.array(i[5:8])
                
        
     
    # Converts Hydrogen to Oxygen for the selected position for linkage making
        direction_vector = (h_cord - c_cord)/np.sqrt(np.sum((h_cord -c_cord)**2))    
        O_pos = direction_vector* COC.distance + c_cord
        
        for i in range(len(first_handle.pdb_data)):
            if first_handle.pdb_data[i][1] == connect_data.first_hyd:  # identify removable H
                
                first_handle.xyz_data[i] = list( O_pos)  # replace removable H with O coords
                first_handle.pdb_data[i][8] = 'O'       # replace atom type for O
                
                first_handle.pdb_data[i][1] = connect_data.new_oxygen  # giving atom name of O
                pdb_components[component_seq[0]].head_atom = connect_data.new_oxygen
        
        
        
        #first_handle.pdb_write_from_xyz(first_handle.xyz_data, "/home/sudhanshu/HDD2/projects2/CHI_potential_data/bin/x3.pdb")
        prev_handle = copy.deepcopy(first_handle)  # to keep it secure and work on next fragment
        combined_handle = copy.deepcopy(first_handle)  # to keep pdb data from all components in sequencial order
        
    
    #working on second residue and linker
        prev_res = component_seq[0]
        prev_static_hyd = connect_data.first_hyd    
        combined_handle.remove_line_x(prev_handle.index_number_of_atom(pdb_components[prev_res].head_atom))
        
   
        for ii, handle in enumerate(component_handles):    
            curr_res = component_seq[ii]
            
            if ii ==0:
                continue   # first residue is already done
          
            # attached tail to previous head  
            
            if ii == 1:
                if connect_data.omega[-1] == "0":
                    continue   # no special treatment for linker as no omega
                else:
                    
                    connect_backup = copy.copy(connect_data)
                    connect_data.second_hyd = "OY"
                    connect_data.second_carbon = "CY"
                    connect_data.movable_carbon = "CY"
                    connect_data.distance_interaction = -1
                    connect_data.movable_tail_atom = "CY"
                    # connect_data.head_atom_type = "C"
                    # connect_data.new_oxygen = "C5"
                    
    
      
            
            iteration_atom_index = -1
            iteration_atom_index2 = -1
            counter = -1
            
            for i in handle.pdb_data:
                counter += 1
                if i[1] == connect_data.second_carbon:
                    c_data2 = i
                    c_cord2 = np.array(i[5:8])
                    iteration_atom_index = counter  # used for 2nd rotation bond adjestment
                    
                
                if i[1] == connect_data.second_hyd:
                    h_data2 = i
                    h_cord2 = np.array(i[5:8])
                    
                if i[1] == pdb_components[component_seq[ii]].movable_tail_atom:
                    target_atom_coord2 = np.array(i[5:8])
                    iteration_atom_index2 = counter
            
            
            direction_vector2 = (h_cord2 - c_cord2)/np.sqrt(np.sum((h_cord2 -c_cord2)**2))            
            O_pos2 = direction_vector2* COC.distance + c_cord2
            
            for i in range(len(handle.pdb_data)):
                if handle.pdb_data[i][1] == connect_data.second_hyd:  # identify removable H                
                    handle.xyz_data[i] = list( O_pos2 )  # replace removable H with O coords
                    handle.pdb_data[i][8] = connect_data.head_atom_type       # replace atom type for O                
                    handle.pdb_data[i][1] = connect_data.new_oxygen   # giving atom name of O
        
        
            
        
            O_position_diff = O_pos2 - O_pos
            handle.xyz_data = handle.xyz_data - O_position_diff
            
            
            # angle between binding bonds
            vector_1 = c_cord - O_pos
            vector_1 = vector_1 /np.sqrt(sum((vector_1)**2))         
            
            vector_2 = c_cord2 - O_position_diff -O_pos
            vector_2 = vector_2 /np.sqrt(sum((vector_2)**2))         
            
            normal_vector = np.cross(vector_1, vector_2)        
            normal_vector = normal_vector/np.sqrt(sum((np.array(normal_vector))**2))
            
                    
            # to adjust angle between components
            
            axis_of_first_rotation = [O_pos, normal_vector + O_pos]
            target_atom_coord_midle = O_pos
            handle.xyz_data_rotated = np.copy(handle.xyz_data[:])
            simple_sd_minimize_angle( prev_handle.coord_of_atom_from_xyz_data(pdb_components[prev_res].static_carbon),
                                      prev_handle.coord_of_atom_from_xyz_data(pdb_components[prev_res].head_atom), 
                                      COC.angle, 
                                      handle, 
                                      handle.index_number_of_atom(pdb_components[curr_res].movable_carbon),
                                      axis_of_first_rotation)
            
            
          
            
            axis_of_second_rotation = [c_cord, O_pos]           
            simple_sd_minimize_distance( prev_handle.coord_of_atom_from_xyz_data(pdb_components[prev_res].static_head_atom), 
                                        handle, 
                                        handle.index_number_of_atom(pdb_components[curr_res].movable_tail_atom), 
                                        axis_of_second_rotation )
            
            
            if ii == 1:
                if connect_data.omega[-1] == "1":
                    pdb_components[curr_res].movable_tail_atom = "CX"
                
            handle.xyz_data = np.copy(handle.xyz_data_rotated[:])
            
            
            axis_of_third_rotation = [O_pos, handle.coord_of_atom_from_xyz_data(pdb_components[curr_res].movable_carbon)]
            print("starting 3rd min")
            simple_sd_minimize_distance( prev_handle.coord_of_atom_from_xyz_data(pdb_components[prev_res].static_head_atom), 
                                        handle, 
                                        handle.index_number_of_atom(pdb_components[curr_res].movable_tail_atom), 
                                        axis_of_third_rotation, connect_data.distance_interaction)
            
            #handle.refresh_from_xyz()
                    
            # rorate for gt cases
            
            if ii ==2:
                if gt_val == "a":
                    gt_C = connect_data.gt_carbon
                    gt_H = gt_C.replace('C','H') + "1"
                    gt_H_not = gt_C.replace('C','H') + "2"
                    
                    
                elif gt_val == "e":
                    gt_C = connect_data.gt_carbon
                    gt_H = gt_C.replace('C','H') + "2"
                    gt_H_not = gt_C.replace('C','H') + "1"
                    
                if connect_data.omega[-1] == "1":
                    gt_C_cord = handle.coord_of_atom_from_xyz_data(gt_C)
                    gt_H_cord = handle.coord_of_atom_from_xyz_data(gt_H)
                    index_gt_h = handle.index_number_of_atom(gt_H)
                    
                    direction_vector = (gt_H_cord - gt_C_cord)/np.sqrt(np.sum((gt_H_cord -gt_C_cord)**2))    
                    new_O_pos = direction_vector* COC.distance + gt_C_cord
                    new_HO_position = direction_vector* (COC.distance + 1.096 )+ gt_C_cord
                    new_O_name = gt_C.replace('C','O')
                    new_HO_name = "H"+new_O_name
                    
                    print('>>', gt_val, gt_C, gt_H, handle.pdb_data[index_gt_h],gt_C_cord,gt_H_cord )
                    handle.pdb_data[index_gt_h][5] = new_O_pos[0]
                    handle.pdb_data[index_gt_h][6] = new_O_pos[1]
                    handle.pdb_data[index_gt_h][7] = new_O_pos[2]
                    handle.xyz_data[index_gt_h] = new_O_pos
                    handle.rename_atom(gt_H, new_O_name)
                    handle.rename_atom_type(index_gt_h, "O")
                    handle.add_new_atom(new_HO_position, new_HO_name, 'H')
                    
                    
                    third_cord = handle.coord_of_atom_from_xyz_data(gt_H_not)
                    
                    vector_1a = new_HO_position - handle.coord_of_atom_from_xyz_data(gt_C)
                    vector_1 = vector_1a /np.sqrt(sum((vector_1a)**2))                     
                    vector_2a = third_cord - handle.coord_of_atom_from_xyz_data(gt_C)
                    vector_2a = vector_2a /np.sqrt(sum((vector_2a)**2))   
                    normal_vectora = np.cross(vector_1a, vector_2a)        
                    normal_vectora = normal_vectora/np.sqrt(sum((np.array(normal_vectora))**2))
                    index_new_ho = handle.index_number_of_atom(new_HO_name)
                    
                    # to adjust angle between components
                    
                    axis_of_h_rotation = [new_O_pos, normal_vectora + new_O_pos]
                            
                    
                    new_cordd = rotate_translate_atom(handle.coord_of_atom_from_xyz_data(gt_C), new_O_pos, 
                                              new_HO_position, axis_of_h_rotation, 180 - COCC.angle_COH)
                    
                    
                    
                    
                    #print(d.xyz)
                    new_cordd = np.copy(new_cordd.xyz)
                    
                    handle.pdb_data[index_new_ho][5] = new_cordd[0]
                    handle.pdb_data[index_new_ho][6] = new_cordd[1]
                    handle.pdb_data[index_new_ho][7] = new_cordd[2]
                    
                    print(new_cordd)
                    
                    handle.xyz_data[index_new_ho] = new_cordd
                    #handle.refresh()
                    handle.refresh_from_xyz()
                    
                    print(">>", connection + gt_val )
                    
                    combined_handle.write_pdb("/home/sudhanshu/HDD2/projects2/CHI_potential_data/bin/models_omega1/chk/x6.pdb")  
                    
                    #handle.refresh()
            
            handle.refresh_from_xyz()           
            
            
            for hdl in handle.pdb_data:
                combined_handle.pdb_data.append(hdl)
                combined_handle.refresh()
            
            
            if ii ==1:
                if connect_data.omega[-1] == "1":        
                    #combined_handle.remove_line_x(combined_handle.index_number_of_atom("CX"))
                    prev_handle = copy.copy(handle)
                    connect_data = copy.copy(connect_backup)
                    connect_data.first_carbon = "CY"
                    O_pos = handle.coord_of_atom_from_xyz_data("CY")
                    c_cord = handle.coord_of_atom_from_xyz_data("CX")
                    prev_res = component_seq[ii]
                    connect_data.head_atom_type = "C"
                    pdb_components[curr_res].static_head_atom = connect_data.new_oxygen
                    connect_data.new_oxygen = "CX"
                    
                                        
                    COC.angle = 115.163
                    COC.distance = 1.553
                    
                
        if connect_data.omega[-1] == "1": # extra vertual atom removal and renaming
            combined_handle.remove_line_x(combined_handle.index_number_of_atom("CX"))   
            combined_handle.remove_line_x(combined_handle.index_number_of_atom("CX"))  
            combined_handle.refresh()    
            combined_handle.rename_atom("CY", connect_data.original_second_carbon)
            combined_handle.rename_atom("HY1", "H"+connect_data.original_second_carbon[1:]+"1")
            combined_handle.rename_atom("HY2", "H"+connect_data.original_second_carbon[1:]+"2")
            combined_handle.refresh()  
        
        #
        
          
        #combined_handle.pdb_write_from_xyz(combined_handle.xyz_data, "/home/sudhanshu/HDD2/projects2/CHI_potential_data/bin/x5.pdb
        
        
        if connect_data.omega[-1] == "0":
            combined_handle.pdb_write_from_xyz(combined_handle.xyz_data, model_dir +connection +".pdb")   
    
        if connect_data.omega[-1] == "1":   
            
            
            
            combined_handle.pdb_write_from_xyz(combined_handle.xyz_data, model_dir + connection + gt_val +"_TEMP.pdb")
            atoms = phi_psi_atom_decoder_from_connection_named_pdb(model_dir + connection + gt_val +"_TEMP.pdb") 
            print(atoms)
            omega_tors = atoms[2:6]
            geom_adj = geometry_adjustment()
            geom_adj.input_xyz(combined_handle.xyz_data)
            geom_adj.input_torsion(omega_tors)
            geom_adj.atm_index = combined_handle.atm_array
            
            os.remove(model_dir + connection +gt_val+"_TEMP.pdb")
            for j in [60,180,300]:
                geom_adj.rotate_molecule_to(j)
                combined_handle.pdb_write_from_xyz(geom_adj.xyz, model_dir + connection +"_"+gt_val+str(j)+".pdb" )
                
            
                
            
        gt_val = ""
            
            
        
                
                
            
    
    
        
        
        
    
    
    
    
    