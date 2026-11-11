'''
    fixfield_ML_merge.py
    
    merge fixfield output for machine learning

'''

import os
import sys, getopt
import pandas as pd
import numpy as np
import re
import time
import scipy

m_pion = 139.57 #pion mass in MeV
m_muon = 105.66 #muon mass in MeV
m_part = 0
part_type = ''

def get_number(string):#only return the number in a string for sotrting
    #print(int(float(re.findall(r'\d+', string)[0])))
    return int(float(re.findall(r'\d+', string)[0]))

def writeFile(new_df, fileName, output_format):
    if output_format == "fixfield":
        new_df.to_csv(fileName, sep=' ', index=False)#, header=None)
    else:
        new_df.to_csv(fileName, sep=' ', index=False)#, header=None)

def main(output_lst):#names=['x','y','z','px','py','pz','xp','yp','E'] -> bdsim output format, columns=['x', 'xp', 'y', 'yp', 'z', 'E'] -> LhARA output format

    #fixfield format m0 [MeV*c^-2], q [unit charge], x [m], y [m], z [m], Ptot [MeV/c], ux [], uy [], uz []
    output_file_lst = output_lst[1:] #extract the output file names and exclude the argv length
    
    ML_input_file_lst = []
    ML_output_file_lst = []
    
    ML_input_lst = []
    ML_output_lst = []
    
    ML_input_particle_sample_file_name = "ML_test_particle_input"
    ML_output_particle_file_name = "ML_test_particle_out"
    
    for i in range(0,len(output_file_lst)):
        
        if ML_input_particle_sample_file_name in output_file_lst[i]:
            ML_input_file_lst.append(output_file_lst[i])
        
        elif ML_output_particle_file_name in output_file_lst[i]:
            ML_output_file_lst.append(output_file_lst[i])
            
    ML_input_file_lst = sorted(ML_input_file_lst, key=get_number)
    ML_output_file_lst = sorted(ML_output_file_lst, key=get_number)
    
    print("start printing input and output list")
    print("====================================================================")
    print(ML_input_file_lst)
    print(ML_output_file_lst)
    print("====================================================================")
        
    fixfield_headers =  ["m0 [MeV*c^-2]", "q [unit charge]", "x [m]", "y [m]", "z [m]", "Ptot [MeV/c]", "ux []", "uy []", "uz []"]
    ML_sample_headers =  ["rmin [m]", "rmax [m]", "zmax [m]", "Cell opening angle [deg]", "r0 [m]", "k_val",
                            "C1", "lambda0",  "mag0 B0 [T]",  "mag0 cent_theta[deg]", "mag0 length_theta[deg]", 
                            "mag1 B0 [T]",  "mag1 cent_theta[deg]", "mag1 length_theta[deg]", "mag2 B0 [T]",  "mag2 cent_theta[deg]", "mag2 length_theta[deg]", 
                            "m0 [MeV*c^-2]", "q [unit charge]", "x [m]", "y [m]", "z [m]", "Ptot [MeV/c]", "ux []", "uy []", "uz []"]
    
    print("input col num : "+str(len(ML_sample_headers))+", output col num : "+str(len(fixfield_headers)))
    
    starttime = time.time()
    
    for i in range(0,len(ML_input_file_lst)):
        
        df = pd.DataFrame()
        df = pd.read_csv(ML_input_file_lst[i], sep = "\s+", skiprows=1, names=ML_sample_headers) #input for ML contains all info of FFA magnets
        ML_input_lst.append(df)
        
    ML_input_merged_df = pd.concat(ML_input_lst, ignore_index=True)#ignore_index=True in order to assign new index to row without original index used and same index assigned
    ML_input_merged_df.reset_index(drop=True) #to give proper row index
    
    for i in range(0,len(ML_output_file_lst)):
        
        df = pd.DataFrame()
        df = pd.read_csv(ML_output_file_lst[i], sep = "\s+", skiprows=1, names=fixfield_headers) #output for ML same with fixfield particle format
        ML_output_lst.append(df)
        
    ML_output_merged_df = pd.concat(ML_output_lst, ignore_index=True)
    ML_output_merged_df_with_lost_part = pd.DataFrame(columns=fixfield_headers) #pd dataframe that contains output + lost particles
    
    
    input_leng = len(ML_input_merged_df["x [m]"].to_numpy()) #length of input - total number of particles at begining
    input_Ptot_lst = ML_input_merged_df["Ptot [MeV/c]"].to_numpy()
    output_Ptot_lst = ML_output_merged_df["Ptot [MeV/c]"].to_numpy()
    print(ML_input_merged_df["m0 [MeV*c^-2]"])
    a = 0
    for i in range(0,input_leng):
        if a<len(ML_output_merged_df["x [m]"].to_numpy()):
            Ptot_ratio = (input_Ptot_lst[i]-output_Ptot_lst[a])/input_Ptot_lst[i]
            if Ptot_ratio<1e-6: #when particle in output and input are the same particle
                ML_output_merged_df_with_lost_part.loc[i] = ML_output_merged_df.loc[a]#len(ML_output_merged_df_with_lost_part)
                a = a + 1
            else:
                lost_part = {"m0 [MeV*c^-2]":ML_input_merged_df["m0 [MeV*c^-2]"][int(i/3)],
                         "q [unit charge]":ML_input_merged_df["q [unit charge]"][i],
                         "x [m]":0,
                         "y [m]":0,
                         "z [m]":0,
                         "Ptot [MeV/c]":ML_input_merged_df["Ptot [MeV/c]"][i],
                         "ux []":0,
                         "uy []":0,
                         "uz []":0}
                ML_output_merged_df_with_lost_part.loc[i] = lost_part
        else: # when input particle is not found in output(if particle is lost) assign 0 to x, y, z, ux, uy, uz to indicate lost poarticle
            lost_part = {"m0 [MeV*c^-2]":ML_input_merged_df["m0 [MeV*c^-2]"][int(i/3)],
                         "q [unit charge]":ML_input_merged_df["q [unit charge]"][i],
                         "x [m]":0,
                         "y [m]":0,
                         "z [m]":0,
                         "Ptot [MeV/c]":ML_input_merged_df["Ptot [MeV/c]"][i],
                         "ux []":0,
                         "uy []":0,
                         "uz []":0}
            ML_output_merged_df_with_lost_part.loc[i] = lost_part
        
        if i%(input_leng/10) == 0:
            print(str(i)+" particles processed")
            
    

    timepass = time.time()-starttime
    print("total surviving particles : "+str(a))
    print("merging completed in "+str(timepass)+" sec. for "+str(len(ML_input_file_lst))+" files and total of "+str(len(ML_input_merged_df["m0 [MeV*c^-2]"]))+" particles")
    print("input part num : "+str(len(ML_input_merged_df["m0 [MeV*c^-2]"]))+", output part num : "+str(len(ML_output_merged_df_with_lost_part["m0 [MeV*c^-2]"])))
    
    out_name = "data/ML_input_sample/ML_input_merged.dat"
    
    print("writing to file : "+out_name)
    
    writeFile(ML_input_merged_df, out_name, "ML_input")
    
    out_name = "data/ML_output_sample/ML_output_merged.dat"
    
    print("writing to file : "+out_name)
    
    writeFile(ML_output_merged_df_with_lost_part, out_name, "fixfield")
    
if __name__ == "__main__" :
    main(sys.argv)