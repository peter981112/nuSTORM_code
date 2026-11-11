'''
    fixfield_read_FFA_def.py
    
    read .latt FFA def file

'''

import os
import sys
import pandas as pd
import numpy as np
import scipy

def read_file(filename): #return flux in format of list
    filein = open(filename, "r")            # Open file for reading
    inputd = filein.readlines()

    leng = len(inputd)
    print("the number of line of the txt file : " + str(filename) + " list is : " + str(leng))
    inputdata = []

    for i in range(0,leng):
        tokens = inputd[i].split()
        #print(str(i)+" th line start =============================")
        #print(tokens)
        token_data = []
        for j in range(0,len(tokens)):
            if '//' in tokens[j]:
                break
            else:
                token_data.append(tokens[j])
                
        #print(token_data)
        inputdata.append(token_data)
            
    return inputdata

input_FFA_def_file = sys.argv[1]
FFA_def_file_data = read_file(input_FFA_def_file)

FFA_type = FFA_def_file_data[8][1]

cell_open_ang = float(FFA_def_file_data[11][0])

FFA_mag_B0 = [float(FFA_def_file_data[14][2]), float(FFA_def_file_data[18][2]), float(FFA_def_file_data[22][2])]
C1_val = float(FFA_def_file_data[15][3])
lambda0_val = float(FFA_def_file_data[15][1])
cent_theta = [float(FFA_def_file_data[14][0]), float(FFA_def_file_data[18][0]), float(FFA_def_file_data[22][0])]
length_theta = [float(FFA_def_file_data[16][0]), float(FFA_def_file_data[20][0]), float(FFA_def_file_data[24][0])]
k_val = float(FFA_def_file_data[14][3])
r_0_val = float(FFA_def_file_data[14][1])
r_max_val = float(FFA_def_file_data[10][1])
r_min_val = float(FFA_def_file_data[10][0])
z_max_val = float(FFA_def_file_data[10][2])

FFA_out_string = ""
FFA_out_string = FFA_out_string+"FFA_type = '"+str(FFA_type)+"'\n\ncell_open_ang = "+str(cell_open_ang)+"\n\nFFA_mag_B0 = "+str(FFA_mag_B0)
FFA_out_string = FFA_out_string+"\nC1_val = "+str(C1_val)+"\nlambda0_val = "+str(lambda0_val)+"\ncent_theta = "+str(cent_theta)+"\nlength_theta = "+str(length_theta)
FFA_out_string = FFA_out_string+"\nk_val = "+str(k_val)+"\nr_0_val = "+str(r_0_val)
FFA_out_string = FFA_out_string+"\nr_max_val = "+str(r_max_val)+"\nr_min_val = "+str(r_min_val)+"\nz_max_val = "+str(z_max_val)

print("FFA data====================================")
print(FFA_out_string)
