#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 15:06:56 2019

@author: sudhanshu
"""
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import pylab as plb
from scipy.optimize import curve_fit
import sys
import os
from src.clustering_parameter_functions import link_name_converter
#from scipy import asarray as ar,exp


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



class fit_with_gaussian_file_data:
    def __init__(self,in_file, gaus_num=6, ss_plot = 1):
        '''same code as fit_with_gaussian_ in 7_quantum
        only difference is it can be applied directly to a file'''
        self.x = []
        self.y = []
        self.print_val = 1
        self.write_to_file = 0
        self.gaus_num = gaus_num
        self.third_val = 200
        self.back_dir = "./.data_parameters"
        self.param_dir = "./params/"
        self.dir = ""
        self.link_name =""
        
        self.ss_plot = ss_plot
        if not ss_exist(self.back_dir,"dir"):
            os.mkdir(self.back_dir)
        
        if not ss_exist(self.param_dir,"dir"):
            os.mkdir(self.param_dir)
        self.in_file = in_file
        self.parameter_file_init = in_file.split("/")[-1].replace(".runconfig","")+"_"+str(self.gaus_num)+"_"  
        self.out_file_name = (self.param_dir + self.parameter_file_init[:-2] + "g" + self.parameter_file_init[-2:-1] + ".params").replace(".exphist","")      
        self.data_file_read(in_file)
        self.find_best_fit()
        
        
        
        
        
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

        #possible link name file:
        link_name_file = self.in_file.replace("exphist","links")        
        if os.path.exists(link_name_file):
            fid = open(link_name_file,"r")
            data = fid.readline()
            fid.close()
            
            self.link_name = data


        
        #return(x,data_np)
        return ("X and Y initiated")
    
    #def von_mises(x,a,b,c):
    #    return (a/(2*np.pi))*np.exp(c*(np.cos*(x-b)))

    def gaus1(self,x,a,b,c):
        return a*np.exp((-1*(x-b)**2)/c)    
    
    
    def gausx(self,x,*abc):
        count = len(abc)//3
        val = 0
        for i in range(count): 
            val = val+ self.gaus1(x,abc[i*3],abc[i*3+1],abc[i*3+2])
        return val
        
    def return_gaus_type(self,gaus_number):
        return self.gausx
      
    def mse(self,d1,d2):
        return sum((d1-d2)**2)/len(d1)
    
    def rsquare_data(self,d11,d22): #d11 : experminetal data X data
        
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
    
    
    
    def run_gaus_for_start_point(self, start_v,pltx=0, third_val=200,  extra_d=0):
        gaus_number = self.gaus_num
        x = self.x
        y = self.y
        third_val = self.third_val
        
        curr_gaus = self.return_gaus_type(gaus_number) 
        p0 = np.ones(gaus_number*3)*start_v
        p0[1::3]=np.array([360/gaus_number]*gaus_number)-(180/gaus_number) + min(x)
        #p0[2::3]=1
        popt = np.array([])
        pcov = np.array([])
        
        # |      _._
        # |     / | \
        # |    /  |  \
        # |   /   a   \
        # | _/    |    \_
        # <---b--->
        #    <----c---->
        #           a      b     c
        bound_l = (-200, -360, third_val)*gaus_number  #make third walue 100 for omega
        bound_u = ( 2000, 360, third_val+3000)*gaus_number   # first value 1000 for mirror phi
        bound=(bound_l,bound_u)
        
        # you may need to play with bound balues for faster and accurate calculation
        
        
        
        self_param_file = self.parameter_file_init  +str(start_v) + "_" + str(third_val) + "_" + str(extra_d)
        
        #print(self.parameter_file_init  +str(start_v) + "_" + str(third_val) + "_" + str(extra_d))
        
        param_file = self.back_dir+"/"+self_param_file

        try:
            
            if ss_exist( param_file, "file" ):
                fidx = open(param_file,"r")
                p_file_data = fidx.readlines()
                fidx.close()
                #print("from file" + param_file)
                popt = p_file_data[0].strip().split(",")
                #print(popt, len(popt))
                
                if len(popt) < 2:
                    popt = []
                else:
                    popt = np.array(popt).astype(float)
                #print(popt)
                msev = float(p_file_data[1])
                rv = float(p_file_data[2])
                #return ( popt,msev,rv )
            else:
                popt,pcov = curve_fit(curr_gaus,x,y,maxfev=300000, p0=p0, bounds=bound)
                rv = self.rsquare_data(y,curr_gaus(x,*popt))[0]
                msev = self.mse(y,curr_gaus(x,*popt))
        except:
            print(" Not Converged!",end="")
            rv=0
            msev=100000
        if msev != 100000:
            if (pltx == 1):
                plt.figure()  
                if self.ss_plot == 1:
                    self.plt_prop()
                x2=range(int(np.min([min(x),0])),int(np.min([min(x),0]))+362,1)
                y2 = curr_gaus(x2,*popt)
                
                if self.ss_plot == 1:
                    plt.plot(x,y,'b+:',label='data') 
                    plt.plot(x2,y2,'r-',label='fit')
                    
                    plt.legend()
                    plt.title('Fig. Fit for best MSE with'+ str(start_v))
                    plt.xlabel('Angle (°)')
                    plt.ylabel('Energy Kcal/Mol')     
                    plt.xticks(range(int(np.min([min(x),0])),int(np.min([min(x),0]))+363,30))
                
                
                if self.print_val == 1:
                
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
                    
                    if self.link_name == "":
                        use_link_name = "LINKID  " + link_name_converter(self.out_file_name.split("/")[-1])
                    else:
                        use_link_name = self.link_name
                        
                    
                    
                    if self.write_to_file == 0:
                        print("# PAGS Parameters.")
                        print("FUNCTION  GAUSSIAN")
                        print(use_link_name)
                        print("a    "+ aa_)
                        print("b    "+ bb_)
                        print("c    "+ cc_)
                        #print("c   ", popt[2::3])
                        print ("d    "+("%10.7f" % (-(min(y2))+extra_d))+"   # intercept (coefficient of zeroth order term)")
                        
                    else:
                        fid_out = open(self.out_file_name, "w+")
                        fid_out.write("# PAGS Parameters for linkage "+
                                      use_link_name+ "\n")
                        fid_out.write("FUNCTION  GAUSSIAN"+ "\n")
                        fid_out.write( use_link_name+ "\n")
                        fid_out.write("a    "+ aa_ + "\n")
                        fid_out.write("b    "+ bb_ + "\n")
                        fid_out.write("c    "+ cc_ + "\n")
                        fid_out.write("d    "+("%10.7f" % (-(min(y2))+extra_d))+"   # intercept (coefficient of zeroth order term)\n")
                        self.write_to_file = 0
                
                self.popt = popt
                self.msev = msev
                self.rv = rv
        #print ("fitting done")
        fidx = open(param_file,"w+")
        popt_str=""
        for i in popt:
            popt_str = popt_str + str(i) + ","
        fidx.write(popt_str[:-1]+"\n")
        
        fidx.write(str(msev)+"\n")
        fidx.write(str(rv)+"\n")   
        fidx.close()        
        
        return ( popt,msev,rv )
     
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
    
    def find_best_fit(self):
        #data_p ="/home/sudhanshu/Desktop/projects/7_quantum/thp2x_b3lyp/trial_set/inps/chi_data"
        gaus_number = self.gaus_num
        curr_gaus = self.return_gaus_type(gaus_number)
        x = self.x
        y = self.y
        n = len(x)                          #the number of data
        rsq2=[]
        t_val = [600,500,400,300,200,100]
        for k in t_val:
            self.third_val = k
            for ii in range(1,50):
                i = ii*10 +190
                print("\rTrying " + str(i),end="")
                popt,msev,rv = self.run_gaus_for_start_point( i )
                rsq2.append([i,rv,msev,k])
                if msev<=0.001:        
                    break
        print("\r----------------") 
        # Tried R-squared and MSE both, if energy difference is low, MSE works better.

        
        rsq2 = np.array(rsq2)
        self.rsq2 = rsq2
        counter = 0
        for i in t_val:
            counter += 1
            rsq = rsq2[np.where( rsq2[:,3]==i)[0],0:3]
            
            min_mse = np.min(rsq[:,2])
            pos_min = np.where(rsq[:,2]==min_mse)
            val_min = rsq[pos_min[0][0],0]
            
            
            #print(val_min)
            
            #print ("\nMin MSE:")
            self.third_val = i
            self.print_val = 0
            popt,m,r = self.run_gaus_for_start_point( val_min, 1 )
            if self.ss_plot == 1:
                # plt.title('Fig. Fit for best MSE with:'+ str(val_min) + " #:" + str(counter))
                plt.title('Fig. Fit for best MSE: #' + str(counter), fontsize=20)
                plt.show(block=False)
        
            
            #print ("Best MSE for " + str(gaus_number) + " gausian equation is " +
                 #  ("%6.4f" % min_mse) + ", start_point_v: " + str(val_min))
        if self.ss_plot == 1:    
            inp_sel = input("Input selection number: (0 [default] for no selection) ")
            if len(inp_sel.strip()) == 0:
                inp_sel = "0"
            inp_sel = int( inp_sel )-1
            if inp_sel == -1:
                return
            inp_sel = t_val[inp_sel]
            
            rsq = rsq2[np.where( rsq2[:,3]==inp_sel)[0],0:3]
                
            min_mse = np.min(rsq[:,2])
            pos_min = np.where(rsq[:,2]==min_mse)
            val_min = rsq[pos_min[0][0],0]
            self.third_val = inp_sel
            self.print_val = 1
            self.write_to_file = 1
            
                    
            popt,m,r = self.run_gaus_for_start_point( val_min, 1 )
        
            print ("# Best MSE for " + str(gaus_number) + " gausian equation is " +
                     ("%6.4f" % min_mse) + ", start_point_v: " + str(val_min))

    
    
    #################
        #print ("\nBest r-squared:")
        #self.third_val = k_max_rsq
        #popt,m,r = self.run_gaus_for_start_point(val_max_rsq,1 )
        #plt.title('Fig. Fit for best R-squared')    
        
        #print ("Best R-square for " + str(gaus_number) + " gausian equation is " +
         #      ("%6.4f" % max_rsq) + ", start_point_v: " + str(val_max_rsq))
        
    
        
    
    def plt_prop(self):
        plt.xlabel("", fontsize = 15)
        plt.xticks(fontsize = 13)
        plt.grid()
        plt.ylabel(" ", fontsize = 15)
        plt.yticks(fontsize = 13)
        plt.title(" ", fontsize = 20)
    
if __name__ == '__main__':
    
    len_argv = len(sys.argv)
    
    if len_argv == 1:
        print("Give histogram file for calculation.")
    elif len_argv ==2:    
        x=fit_with_gaussian_file_data(sys.argv[1])    
        if x.ss_plot==1:
            plt.show()
    elif len_argv ==3:    
        x=fit_with_gaussian_file_data(sys.argv[1], int(sys.argv[2]))
        # if x.ss_plot==1:
        #     plt.show()
    elif len_argv ==4:    
        x=fit_with_gaussian_file_data(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
        
        # if x.ss_plot==1:
        #     plt.show()
    else:
        print ("wrong input.")
