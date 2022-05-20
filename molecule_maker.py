#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 15:48:18 2021
#UPDATE MAY 19 2021
# This program creates simplified pdb file for given linkage name.
#Linkage name(s) can be given is a list file (see bin/model_lists )
@author: sudhanshu
"""
from src.molecule_maker_functions import (text_file_reader, 
                                          generate_molecule_for_connection,
                                          phi_psi_atom_decoder_from_connection_named_pdb,
                                          )

from settings import global_settings


def list_to_molecule_make(list_name_file, output_path):    
    connection_names = text_file_reader(list_name_file)
    fid = open( output_path+"connection_list_check.py","w+") # for linkage atoms details. 
    fid.write("files_and_torsions = [\n")
    omega_tors = []
    for i in connection_names:
        generate_molecule_for_connection(i, output_path) 
    
        if i[-1] == "0":        
            atoms = phi_psi_atom_decoder_from_connection_named_pdb(output_path+i+".pdb")   
            fid.write("[ '"+i+".pdb', '[" + str(atoms[:4]).replace(" ","") +"," +str(atoms[1:5]).replace(" ","") +"]'],\n")
            omega_tors.append("'[[0,0,0,0]]',")
        if i[-1] == "1":        # when omega present it creat tg,gt,gg forms
            for j in [60,180,300]:                
                for gt_type in ["a","e"]: # for axial and equitorial connections
                
                    atoms = phi_psi_atom_decoder_from_connection_named_pdb(output_path+i+"_"+gt_type+str(j)+".pdb")   
                    fid.write("[ '"+i+'_'+gt_type+str(j)+".pdb', '[" + str(atoms[:4]).replace(" ","") +"," +str(atoms[1:5]).replace(" ","") +"]'],\n")
                    omega_tors.append("'["+str(atoms[2:6]).replace(" ","") +"]',")   
                    

    fid.write("]\n")       
    # if i[-1] == "1":
    fid.write("omega_tors = [\n")
    for om in omega_tors:
        fid.write(om+"\n")
    fid.write("]\n")
    
    fid.close()
    


def main():
    model_name_list = global_settings.list_of_linkage_names
    pdb_out_dir = global_settings.simplified_pdb_output_dir    
    list_to_molecule_make( model_name_list, pdb_out_dir )
    

if __name__ == "__main__":
    main()