#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 18 11:19:11 2022

@author: sudhanshu
"""
from src.one_param_converter_v1_rotation_type import *

q = quality_check()

dir_path = "/home/sudhanshu/Desktop/projects/pipeline/2022/output/2_ensemble/project_8000/"

best_pdb = dir_path + "ensemble_180_135.pdb"
q.load_first_pdb(best_pdb)

array_ =[]
for i in os.listdir(dir_path):
    if not i.endswith('.pdb'):
        continue
    
    val = q.is_bonds_changed_from_first_pdb(dir_path + i)
    array_.append([i,val])

print(array_)