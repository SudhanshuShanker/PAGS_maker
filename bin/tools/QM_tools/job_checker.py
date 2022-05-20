#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 12:24:16 2020

@author: sudhanshu
"""

import os
import time
import subprocess
import datetime

def ss_importdata(*inputs):
    data_only_mode = 0;
    if len (inputs) > 1:
        data_only_mode = inputs[1]
    flnm = inputs[0]  
    fid = open(flnm,'r')
    
    data = fid.readlines()
    fid.close()
    return data

def ss_write_in_file(fl, msg):
    fid2 = open(fl,"w+")
    for i in msg:
        fid2.write(i+"\n")
    fid2.close()


def read_energies(dirnm):
    all_fl = os.listdir(dirnm)
    dir_interest = []
    for i in all_fl:
        if i.endswith("inp_output"):
            dir_interest.append(i)
    all_data = []
    for i in dir_interest:
        d=subprocess.Popen(["grep","FINAL" , dirnm +"/" + i+  "/condor.out_0", ],stdout=subprocess.PIPE)
        dd = d.communicate()[0]
        if len(dd) < 1:
            continue
        score_d = str(dd.splitlines()[-1])
        #final_data = i[9:-11] + score_d[27:-1]
        final_data = i.replace('ensemble_','').replace('.inp_output','').replace('_',' ') + score_d[27:-1]
        all_data.append(final_data)
        
    return all_data    
        




while(1):
    fid_awake = open("awake_me","w+")
    
    now = datetime.datetime.now()
    fid_awake.write(str(now))
    fid_awake.close()
    
    files=os.listdir("./")
    
    run_files = []
    for i in files:
        if i.endswith("runconfig"):
            run_files.append(i)
        
    for i in run_files:
        
        d = ss_importdata(i)
        input_dir = d[1]
        input_dir = input_dir.rstrip()
        if d[0].count("submit") == 1:            
            
            os.system("./tools/gen_script.sh "+ input_dir)
            os.system("./tools/condor_submit.sh "+ input_dir)
            
            print (input_dir)
            print ("submitted")
            
            ss_write_in_file(i, ["started",input_dir])
        elif d[0].count("started") == 1:
            cond_data = subprocess.Popen(["condor_q","sshanke3",], stdout=subprocess.PIPE)
            msg = cond_data.communicate()[0]
            #print (msg)
            count_ = 0
            for m in msg.splitlines():
                #print(str(m))
                if str(m).count(input_dir)>0:
                    count_=count_ +1
            
            if count_ == 0:
                ss_write_in_file(i, ["finished",input_dir])
                

        elif d[0].count("finished") == 1:            
            print("completed")
            data = read_energies(input_dir)
            
            ss_write_in_file(i, ["completed",input_dir]+data)
        

    
    time.sleep(30)
            
        

        
    #break
