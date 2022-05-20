#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 11:39:30 2020

@author: sudhanshu
"""
import os
import sys
import shutil
from string import *
from pyrosetta import *
from pyrosetta.teaching import *
from pyrosetta.rosetta.core.pose import *
import numpy as np
import pandas as pd
import subprocess
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import pylab as plb
from scipy.optimize import curve_fit
import scipy




class fit_with_vm():

    def __init__(self):
        self.x = []
        self.y = []
        self.gaus_num = 4
        self.max_eqn = 5
        self.third_val = 1
        self.dir = ""
        self.quality_check = "RSQR" #["RSQR","MSE"]
        self.user_name = ""
        self.result_dir = ""
        self.file_number = 0
    
    def ss_exist(self,nm,typ):
        out=False
        if os.path.exists(nm):
            if typ=='dir':
                if os.path.isdir(nm):
                    out = True
            elif typ == 'file':
                if os.path.isfile(nm):
                    out = True
                                    
        return out  
        
        
    def input_dir_name(self,dir_name):
        self.dir = dir_name
        
    def result_dir_name(self,rdir):
        self.result_dir = rdir
    
    def create_calc_done_tag(self):
        tag_file = self.result_dir + "/" +self.user_name + "_fitting_data_" +str(self.file_number)+".tag"
        fid = open(tag_file,"w+")
        fid.write("Hi!")
        fid.close()
        
    def check_calc_done_tag(self):
        tag_file = self.result_dir + "/" +self.user_name + "_fitting_data_" +str(self.file_number)+".tag"
        
        return(self.ss_exist(tag_file,'file'))
            
        
    def user_name_current(self,uname):
        self.user_name = uname
        
    def data_file_read(self,data_p,two_col=0):
        data= open(data_p,"r").readlines() 
        all_data=[]
        for i in data:
            d = i.split()
            all_data.append(d)     
        all_data=np.array(all_data).astype(float)
        if all_data.shape[1]==1:
            self.x = np.array(range(len(all_data[:,0])))*15
        else:
            self.x = all_data[:,0]
            
        data_np = all_data[:,-1]*627.509    
        data_np = data_np - np.min(data_np)   
        self.y = data_np    
        #return(x,data_np)
        return ("X and Y initiated")
    
    def data_file_read_server(self,data_p,two_col=0):
        data= open(data_p,"r").readlines() 
        all_data=[]
        counter = 0
        #print(data)
        for i in data:
            counter +=1
            if counter <3:
                continue
            d = i.split()
            all_data.append(d)     
        
        all_data=np.array(all_data).astype(float)
        #self.all_data = all_data
        all_data = all_data[all_data[:,0].argsort()]
        if all_data.shape[1]==1:
            self.x = np.array(range(len(all_data[:,0])))*15
        else:
            self.x = all_data[:,0]
            
        data_np = all_data[:,-1]*627.509    
        data_np = data_np - np.min(data_np)   
        self.y = data_np    
        
        self.file_number = int(data_p.split("/")[-1][2:-10])
        #return(x,data_np)
        return ("X and Y initiated")
    
    def data_file_read_server_auto(self, two_col=0):
        data_p = self.dir_name +"/" 
        data= open(data_p,"r").readlines() 
        all_data=[]
        counter = 0
        for i in data:
            counter +=1
            if counter <3:
                continue
            d = i.split()
            all_data.append(d)     
        
        all_data=np.array(all_data).astype(float)
        #self.all_data = all_data
        all_data = all_data[all_data[:,0].argsort()]
        if all_data.shape[1]==1:
            self.x = np.array(range(len(all_data[:,0])))*15
        else:
            self.x = all_data[:,0]
            
        data_np = all_data[:,-1]*627.509    
        data_np = data_np - np.min(data_np)   
        self.y = data_np    
        #return(x,data_np)
        return ("X and Y initiated")
    
    def gaus1(self, x, a,b,c):
        mu = b
        kappa = c
        #  {\color{DarkOrange}f(x) = \sum\limits_{i=1}^{n} {a_i. \frac {\exp [c_i. \cos(x - b_i) ]} {2 .\pi . i_0(c_i)}    }+ d}
        #return a*np.exp(( np.cos(np.deg2rad(x - mu)))/c) # THIS simplified form can also be used 
        return (a/(2. * np.pi * scipy.special.i0(kappa)))*np.exp((kappa * np.cos(np.deg2rad(x - mu))))  
       
    
    def cos_funct(self,x,a,b):
        return (a*np.cos(np.deg2rad(x - b))) 
        

    def eqn_x(self, x,*abc): #number of equations
        count = len(abc)//3
        val = 0
        for i in range(count): 
            val = val+ self.gaus1(x,abc[i*3],abc[i*3+1],abc[i*3+2])
        return val
    

    def mse(self, d1,d2):
        return sum((d1-d2)**2)/len(d1)
    
    def rsquare_data(self, d11,d22): #d11 : experminetal data X data
        
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
        #print(d1_max_pos)
        
        x_points = np.array([np.min(d1), np.max(d1)])
        y_points = np.array([y_pred[d1_min_pos], y_pred[d1_max_pos]])
        
        xy_out = [x_points, y_points]   
        
        
        return r_sq,slope,xy_out,~np.isnan(nan_vals)
        #plot method
        #plt.text(0.2,0.1,  "slope %5.2f" % model.coef_)
        #plt.plot(x,y_pred)

    def return_gaus_type(self, gaus_number):
        return self.eqn_x


    def run_gaus_for_start_point2(self, start_v,pltx=0, extra_d=0):
        x = self.x
        y = self.y
        third_val = self.third_val
        gaus_number = self.gaus_num
        curr_gaus = self.return_gaus_type(gaus_number) 
        p0 = np.ones(gaus_number*3)*start_v
        p0[1::3]=np.cumsum([360/gaus_number]*gaus_number)-(180/gaus_number)
        #p0[2::3]=10
        popt = pcov = np.array([])
        
        # FOR SPECIAL COMPLEX CASES
        # bound_l = (-500, -200,-np.inf)*gaus_number  #make third walue 100 for omega
        # bound_u = ( 500, 460, 60)*gaus_number
        
        bound_l = (-50, -200,0)*gaus_number  #make third walue 100 for omega
        bound_u = ( 500, 360, 6)*gaus_number
        bound = (bound_l,bound_u)
        
        try:
            popt,pcov = curve_fit(curr_gaus,x,y,maxfev=300000, p0=p0, bounds=bound)
            rv = self.rsquare_data(y,curr_gaus(x,*popt))[0]
            msev = self.mse(y,curr_gaus(x,*popt))
        except:
            print(" Not Converged!")
            rv=0
            msev=100000
        if msev != 100000:
            if (pltx > 0):
                fig = plt.figure()       
                self.plt_prop()
                x2=range(0,360+2,1)
                y2 = curr_gaus(x2,*popt)
                
                plt.plot(x2,y2,'r-',label='fit')
                plt.plot(x,y,'b+:',label='data') 
                plt.legend()
                plt.title('Fig. Fit for best MSE')
                plt.xlabel('Angle (°)')
                plt.ylabel('Energy Kcal/Mol')     
                plt.xticks(range(0,361,30))
                if (pltx ==2):
                    for i in range(gaus_number):
                        y3 = self.gaus1(x2,popt[i*3],popt[i*3 +1],popt[i*3 +2])
                        plt.plot(x2,y3)
    #            for i in range(len(popt)//3):
    #                print("a"+str(i+1)+ " = "+ "%11.3f" % popt[i*3],
    #                      ", b"+str(i+1)+ " = "+ "%11.3f" % popt[i*3+1],
    #                      ", c"+str(i+1)+ " = "+ "%11.3f" % popt[i*3+2])
                
                print("-"*20)
                a_ = popt[0::3]
                b_ = popt[1::3]
                c_ = popt[2::3]
                aa_=''
                bb_=''
                cc_=''
                for i in range(len(a_)):
                    aa_ = aa_ + ("%0.6e" % a_[i]) + "  " 
                    bb_ = bb_ + ("%0.6e" % b_[i]) + "  " 
                    cc_ = cc_ + ("%0.6e" % c_[i]) + "  " 
                aa_ = self.e_handle(aa_) + "# magnitudes of the distributions"
                bb_ = self.e_handle(bb_) + "# midpoints of the distributions"
                cc_ = self.e_handle(cc_) + "# twice the squares of the widths of the distributions"
                #dd_ = ("%10.7f" % (-(min(y2))+extra_d))+"   # intercept (coefficient of zeroth order term)"
                dd_ = ("%10.7f" % 0)+"   # intercept (coefficient of zeroth order term)"
                ee_ = " 1   # the type of equation 0:gaussian and 1:Von-Mises"
                result_file = self.result_dir + "/" + self.user_name+"_fitting_data_" +str(self.file_number)+".dat"
                fid = open(result_file,'w+')
                fid.write("a    "+ aa_ +"\n")
                fid.write("b    "+ bb_ +"\n")
                fid.write("c    "+ cc_ +"\n")
                fid.write("d    "+ dd_ +"\n")   
                fid.write("e    "+ ee_ +"\n")  
                fid.close()
                
                figure_file = self.result_dir + "/" +self.user_name + "_fitting_data_" +str(self.file_number)+".png"
                fig.savefig(figure_file, dpi=200, format="png")
                #print("a    "+ aa_)
                #print("b    "+ bb_)
                #print("c    "+ cc_)
                #print("c   ", popt[2::3])
                #print ("d    "+ dd_)
        return (popt, msev, rv)
    
    def find_best_set_of_equations(self, pltx=0):
        if self.check_calc_done_tag():
            #print("ALRWEADY DONE")
            return 1
        
        start_gaus = self.gaus_num
        max_rsq = 0
        best_gaus_num = start_gaus
        best_start_v = 0.5
         
        for i in range(start_gaus,self.max_eqn+1):
            self.gaus_num = i
            rsq1,start_v=self.find_best_fit2(1)
           #tmp1,mse1,rsq1=self.run_gaus_for_start_point2( start_v,0)
            if rsq1 > max_rsq:
                max_rsq = rsq1
                best_gaus_num = i
                best_start_v = start_v
            
                if rsq1>0.9999:
                    #print (best_start_v, rsq1)
                    break
        
        self.gaus_num = best_gaus_num
        out_final=self.run_gaus_for_start_point2( best_start_v ,pltx)
        self.gaus_num = start_gaus
        self.create_calc_done_tag()
        return(out_final)
                
                

    def e_handle(self,line):
        line = " "+line
        line=line.replace("e+00","    ")
        line= line.replace("e+0", "e")
        line = line.replace("e-0","e-")
        for i in range(5):
            line = line.replace("e"+str(i),"e"+str(i)+"  ")
        for i in range(5):
            line = line.replace("e-"+str(i)+"","e-"+str(i)+" ")
        line = line.replace(" -","-")
        return line
    
    def find_best_fit2(self ,returnv = 0):
        #data_p ="/home/sudhanshu/Desktop/projects/7_quantum/thp2x_b3lyp/trial_set/inps/chi_data"
        gaus_number = self.gaus_num
        curr_gaus = self.return_gaus_type(gaus_number)
        x = self.x
        y = self.y
        n = len(x)                          #the number of data
        rsq=[]
        for ii in range(1,10):
            i = ii*0.5
            if returnv == 0: 
                print("\rTrying " + str(i),end="")
            popt,msev,rv = self.run_gaus_for_start_point2(i )
            rsq.append([i,rv,msev])
            if msev<=0.00001:        
                break
        #print("\r----------------") 
        # Tried R-squared and MSE both, if energy difference is low, MSE works better.
        rsq = np.array(rsq)
        
        if self.quality_check == "MSE":
            min_mse = np.min(rsq[:,2])
            pos_min = np.where(rsq[:,2]==min_mse)
            val_min = rsq[pos_min[0][0],0]
            
            if returnv == 0: 
                print ("\nMin MSE:")
                popt,m,r = self.run_gaus_for_start_point2( val_min,1 )
                plt.title('Fig. Fit for best MSE')
                
                print ("Best MSE for " + str(gaus_number) + " gausian equation is " +
                       ("%6.4f" % min_mse) + ", start_point_v: " + str(val_min))
            if returnv == 1:  
                return(min_mse, val_min)
        elif self.quality_check == "RSQR":
        
            
            max_rsq = np.max(rsq[:,1])
            pos_max_rsq = np.where(rsq[:,1]==max_rsq)
            val_max_rsq = rsq[pos_max_rsq[0][0],0]
    
    
    #################
            if returnv == 0: 
                print ("\nBest r-squared:")
                popt,m,r = self.run_gaus_for_start_point2(val_max_rsq,1 )
                plt.title('Fig. Fit for best R-squared')    
                
                print ("Best R-square for " + str(gaus_number) + " gausian equation is " +
                       ("%6.4f" % max_rsq) + ", start_point_v: " + str(val_max_rsq))
            if returnv == 1:  
                return(max_rsq, val_max_rsq)
        
        
        
    
    def plt_prop(self):
        plt.xlabel("", fontsize = 15)
        plt.xticks(fontsize = 13)
        plt.grid()
        plt.ylabel(" ", fontsize = 15)
        plt.yticks(fontsize = 13)
        plt.title(" ", fontsize = 15)
    
    def all_data(pltx=1):
        main_dir = "/home/sudhanshu/Desktop/projects/7_quantum/pentose/"
        final_out = []
        for i in ['phi','psi']:
            if pltx == 1:
                plt.figure()
            d_angle = [ ]
            for j in [60,180,300]:
                file = main_dir+i+"/"+i+"_"+str(j)
                fid = open(file,'r')
                data = fid.readlines()
                fid.close()
                data2=[]
                for k in data:
                    data2.append(k.split())
                data2 = np.array(data2).astype(np.float)
                d_angle.append(data2)
                if pltx == 1:
                    dplt = data2[:,1]-min(data2[:,1])
                    plt.plot(data2[:,0],dplt*627.509,"o-")
            d_angle = np.array(d_angle)
            final_out.append(d_angle)
            av_plt=np.sum(d_angle,0)/3
            av_plt[:,1]= av_plt[:,1] - min(av_plt[:,1])
            if pltx == 1:
                plt.legend([60,180,300])
                plt.plot(data2[:,0],(av_plt*627.509)[:,1],"k--")
                plt.grid()
                plt.xticks(range(0,361,30))
                #print(data2)
        return np.array(final_out)

# data_file_psi_60 = "/home/sudhanshu/Desktop/projects/7_quantum/pentose/psi/psi_60"  
# server_file = "/home/sudhanshu/Desktop/projects/pipeline/my_web/mount_jazz/p_1.runconfig"
# data_file_psi_180 = "/home/sudhanshu/Desktop/projects/7_quantum/pentose/psi/psi_180" 
# gg = fit_with_vm()
# gg.data_file_read_server(server_file)
# #gg.data_file_read(data_file_psi_180)
# gg.run_gaus_for_start_point2(0.5,1);

   
