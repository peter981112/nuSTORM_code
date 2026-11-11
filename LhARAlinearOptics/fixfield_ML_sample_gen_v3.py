'''
    fixfield_ML_sample_gen_v3.py
    
    generate input sample for ML machine learning 

'''

import os
import sys, getopt
import pandas as pd
import numpy as np
import re
import scipy

def get_number(string):#only return the number in a string for sotrting
    #print(int(float(re.findall(r'\d+', string)[0])))
    return int(float(re.findall(r'\d+', string)[0]))

def read_single_FFA_file(filename): #read data from a FFA def file with a single FFA and return a pandas dataframe
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

    FFA_def_file_data = inputdata

    FFA_type = FFA_def_file_data[8][1]
    print("reading FFA_type : "+str(FFA_type))

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
    
    FFA_data = {'B0' : FFA_mag_B0,
            'C1' : [C1_val, C1_val, C1_val],
            'lambda0' : [lambda0_val, lambda0_val, lambda0_val],
            'cent_theta' : cent_theta, #angle in degree
            'length_theta' : length_theta,
            'k' : [k_val, k_val, k_val],
            'r_0' : [r_0_val, r_0_val, r_0_val],
            'r_max' : [r_max_val, r_max_val, r_max_val],
            'r_min' : [r_min_val, r_min_val, r_min_val],
            'z_max' : [z_max_val, z_max_val, z_max_val]}
            
    FFA_df = pd.DataFrame(FFA_data)
    
    #return FFA_mag_B0, C1_val, lambda0_val, cent_theta, length_theta, k_val, r_0_val, r_max_val, r_min_val, z_max_val, cell_open_ang
    return FFA_df, cell_open_ang

def apply_random_bias_to_FFA_df(FFA_df, rand_bias):

    FFA_mag_B0 = [FFA_df['B0'][0]+rand_bias, FFA_df['B0'][1]-rand_bias, FFA_df['B0'][2]+rand_bias]
    C1_val = FFA_df['C1'][0]+rand_bias
    lambda0_val = FFA_df['lambda0'][0]+0.2*rand_bias
    cent_theta = [FFA_df['cent_theta'][0]+rand_bias, FFA_df['cent_theta'][1]+2*rand_bias, FFA_df['cent_theta'][2]+3*rand_bias]
    length_theta = [FFA_df['length_theta'][0]+1.2*rand_bias, FFA_df['length_theta'][1]+1.5*rand_bias, FFA_df['length_theta'][2]+1.2*rand_bias]
    k_val = FFA_df['k'][0]+2*rand_bias
    r_0_val = FFA_df['r_0'][0]+5*rand_bias
    r_max_val = FFA_df['r_max'][0]+5*rand_bias
    r_min_val = FFA_df['r_min'][0]+5*rand_bias
    z_max_val = FFA_df['z_max'][0]+0.1*rand_bias
    
    FFA_data = {'B0' : FFA_mag_B0,
                'C1' : [C1_val, C1_val, C1_val],
                'lambda0' : [lambda0_val, lambda0_val, lambda0_val],
                'cent_theta' : cent_theta, #angle in degree
                'length_theta' : length_theta,
                'k' : [k_val, k_val, k_val],
                'r_0' : [r_0_val, r_0_val, r_0_val],
                'r_max' : [r_max_val, r_max_val, r_max_val],
                'r_min' : [r_min_val, r_min_val, r_min_val],
                'z_max' : [z_max_val, z_max_val, z_max_val]}
                
    FFA_df_biased = pd.DataFrame(FFA_data)
    
    return FFA_df_biased

def writeFile(new_df, fileName, output_format):
    if output_format == "fixfield":
        new_df.to_csv(fileName, sep=' ', index=False)#, header=None)
    else:
        new_df.to_csv(fileName, sep=' ', index=False, header=None)

    return

def gen_FFA_def_string(FFA_df, cell_open_ang, FFA_type):

    FFA_def_string = "FFA def for ML sample\n\n1				//lattice periodicity [int]\n1				//number of elements [int]\n1				//number of cell types [int]\n\n\n\n"
    FFA_def_string = FFA_def_string + "1 "+FFA_type+" 3			// rclo=31.19248518\n0.005000			//\n"
    FFA_def_string = FFA_def_string + str(FFA_df['r_min'][0])+" "+ str(FFA_df['r_max'][0]) + " "+ str(FFA_df['z_max'][0]) + "	//collimators rmin [m], rmax [m], zmax [m]\n"
    FFA_def_string = FFA_def_string + str(cell_open_ang) + "			//Cell opening angle [deg]\n"
    FFA_def_string = FFA_def_string + "0.000000				//deltar [m]\n\n"
    
    #FFA_data_ML_sample format rmin [m], rmax [m], zmax [m] Cell opening angle [deg] r0 [m] k_val C1 lambda0  mag0(B0 [T]  cent_theta[deg]  length_theta[deg])... ...
    FFA_data_ML_sample = str(FFA_df['r_min'][0])+" "+ str(FFA_df['r_max'][0])+" "+str(FFA_df['z_max'][0])+" "+str(cell_open_ang)+" "+str(FFA_df['r_0'][0])
    FFA_data_ML_sample = FFA_data_ML_sample + " "+ str(FFA_df['k'][0])+" "+str(FFA_df['C1'][0])+" "+str(FFA_df['lambda0'][0])
    
    #Indev FFA format 
    #5.300000-mag center ang	34.052527-r_0 cnet rad	-3.100000-B0	7.252286-k val
    #-4.600000-opening ang(negative num)	0.700000-lambda0	-1.000000-extrapollation order	2.240000-C1
    #4.600000-opening ang(positive num)	0.700000-lambda0	-1.000000-extrapollation order	2.240000-C1
    for i in range(0,3):
        FFA_def_string = FFA_def_string + str(FFA_df['cent_theta'][i]) + "	"+str(FFA_df['r_0'][i])+"	"+str(FFA_df['B0'][i])+"	"+str(FFA_df['k'][i])+"\n"
        FFA_def_string = FFA_def_string + str(-FFA_df['length_theta'][i]/2) + "	"+str(FFA_df['lambda0'][i])+"	-1.000000	"+str(FFA_df['C1'][i])+"\n"
        FFA_def_string = FFA_def_string + str(FFA_df['length_theta'][i]/2) + "	"+str(FFA_df['lambda0'][i])+"	-1.000000	"+str(FFA_df['C1'][i])+"\n"
        FFA_def_string = FFA_def_string + "\n"
        FFA_data_ML_sample = FFA_data_ML_sample + " "+ str(FFA_df['B0'][i])+" "+str(FFA_df['cent_theta'][i])+" "+str(FFA_df['length_theta'][i])

    return FFA_def_string
    

def main(Beamline_type_num, part_per_set, output_dir):
    
    print("Beamline_type_num : "+str(Beamline_type_num)+", part_per_set : "+str(part_per_set)+", output_dir : "+str(output_dir))
    
    FFA_def_file_name = "FFA_ML_test"
    FFA_def_file_ext = ".latt"
    input_particle_file_name_prefix = "ML_test_particle_input"
    input_particle_file_ext = ".dat"
    ML_input_particle_sample_file_name= "ML_test_particle_input"
    FFA_input_data_file_name= "ML_test_FFA_input"
    ML_input_particle_sample_file_ext = ".txt"
    
    FFA_def_file_dir_header = "/home/wonjong/LhARAlinearOptics_workspace/FixField_test_dir/fixfield_simplified/data/ML_base_FFA_def_files/"
    FFA_def_file_header = "FFA_def_"
    base_beamline_type_num = len(os.listdir(FFA_def_file_dir_header))
    base_beamline_FFA_df_lst = []
    base_cell_open_ang_lst = []
    
    for Beamline_num in range(0,base_beamline_type_num): #read base FFA beamline data
        FFA_file_def_dir = FFA_def_file_dir_header+FFA_def_file_header+str(Beamline_num+1)+"/"
        FFA_file_lst = os.listdir(FFA_file_def_dir)
        FFA_file_lst = sorted(FFA_file_lst, key=get_number)
        base_beamline_FFA_df_lst.append([])
        base_cell_open_ang_lst.append([])
        
        for FFA_file_num in range(0,len(FFA_file_lst)):
            FFA_df, cell_open_ang = read_single_FFA_file(FFA_file_def_dir+FFA_file_lst[FFA_file_num])
            base_beamline_FFA_df_lst[Beamline_num].append(FFA_df)
            base_cell_open_ang_lst[Beamline_num].append(cell_open_ang)
        

    FFA_index = 0
    for Beamline_num in range(0,Beamline_type_num):
    
        beamline_type = int(Beamline_num%base_beamline_type_num) #represents the type of beamline in numbers
        FFA_num_in_beam = len(base_beamline_FFA_df_lst[beamline_type]) #number of FFA magnets in the beamline
        print("FFA_num_in_beam : "+str(FFA_num_in_beam))
        
        FFA_def_file_name_lst =[]
        FFA_input_data_file_name_lst = []
        FFA_df_lst = []
        cell_open_ang_lst = base_cell_open_ang_lst[beamline_type]
        FFA_def_string_lst = []
        
        input_particle_file_name = output_dir + input_particle_file_name_prefix + "-" + str(FFA_index) + input_particle_file_ext
        
        FFA_type = 'ffag-r-enge'
        
        for i in range(0,FFA_num_in_beam):
        
            FFA_def_file_name_lst.append(output_dir + FFA_def_file_name + "-"+str(i)+"-" + str(FFA_index) + FFA_def_file_ext)#str(FFA_num_1) + FFA_def_file_ext
            FFA_input_data_file_name_lst.append(output_dir + FFA_input_data_file_name + "-"+str(i)+"-" + str(FFA_index) + ML_input_particle_sample_file_ext)

            rand_bias = 0.00*(np.random.rand()-0.5)#put random bias to get different FFA samples for ML
            FFA_df_lst.append(apply_random_bias_to_FFA_df(base_beamline_FFA_df_lst[beamline_type][i], rand_bias))

            FFA_def_string_lst.append(gen_FFA_def_string(FFA_df_lst[i], cell_open_ang_lst[i], FFA_type))


        Ptot_res = 25
        Ptot_min = 2.5e3#^3 for GeV
        Ptot_max = 8.5e3
        Ptot_lst = np.linspace(Ptot_min, Ptot_max, Ptot_res)
        dp_val = Ptot_lst[1]-Ptot_lst[0]
        
        part_mass_res = 15
        mass_min = 100#100MeV
        mass_max = 200#200MeV
        mass_lst = np.linspace(mass_min, mass_max, part_mass_res)
        dm_val = mass_lst[1]-mass_lst[0]
        
        x_z_max = 0.2 #unit m
        
        u_x_z_max = 0.012
        
        total_part_num = part_per_set*Ptot_res*part_mass_res
        
        m0 = np.zeros(total_part_num)
        q = np.zeros(total_part_num)
        x = np.zeros(total_part_num)
        y = np.zeros(total_part_num)
        z = np.zeros(total_part_num)
        Ptot = np.zeros(total_part_num)
        ux = np.zeros(total_part_num)
        uy = np.zeros(total_part_num)
        uz = np.zeros(total_part_num)
        
        a = 0 #for index of np lists
        for m_index in range(0, len(mass_lst)):
            for Ptot_index in range(0, len(Ptot_lst)):
                for part in range(0, part_per_set):
                    Ptot[a] = Ptot_lst[Ptot_index] + dp_val*np.random.rand()/2 #add random P up to half of dP so P value is uniqe for each particle
                    m0[a] = mass_lst[m_index] + dm_val*np.random.rand()/2
                    q[a] = 1.0 #+1e charge
                    #unique P value will be used for particle identification during merging process
                    
                    x[a] = 2*(np.random.rand()-0.5)*x_z_max
                    #y[a] = np.zeros(total_part_num)
                    z[a] = 2*(np.random.rand()-0.5)*x_z_max
                    
                    ux[a] = 2*(np.random.rand()-0.5)*u_x_z_max
                    uy[a] = 1
                    uz[a] = 2*(np.random.rand()-0.5)*u_x_z_max
                    norm_factor = np.sqrt(ux[a]**2+uy[a]**2+uz[a]**2)
                    ux[a] = ux[a]/norm_factor
                    uy[a] = uy[a]/norm_factor
                    uz[a] = uz[a]/norm_factor
                    
                    a = a+1
            
        #FFA_data_ML_sample format rmin [m], rmax [m], zmax [m] Cell opening angle [deg] r0 [m] k_val C1 lambda0  mag0(B0 [T]  cent_theta[deg]  length_theta[deg])... ...
        FFA_sample_df_1 = pd.DataFrame()
        part_df = pd.DataFrame()
        fixfield_headers =  ["m0 [MeV*c^-2]", "q [unit charge]", "x [m]", "y [m]", "z [m]", "Ptot [MeV/c]", "ux []", "uy []", "uz []"]
        ML_sample_headers =  ["rmin [m]", "rmax [m]", "zmax [m]", "Cell opening angle [deg]", "r0 [m]", "k_val",
                            "C1", "lambda0",  "mag0 B0 [T]",  "mag0 cent_theta[deg]", "mag0 length_theta[deg]", 
                            "mag1 B0 [T]",  "mag1 cent_theta[deg]", "mag1 length_theta[deg]", "mag2 B0 [T]",  "mag2 cent_theta[deg]", "mag2 length_theta[deg]", 
                            "m0 [MeV*c^-2]", "q [unit charge]", "x [m]", "y [m]", "z [m]", "Ptot [MeV/c]", "ux []", "uy []", "uz []"]
        FFA_sample_headers =  ["rmin [m]", "rmax [m]", "zmax [m]", "Cell opening angle [deg]", "r0 [m]", "k_val",
                            "C1", "lambda0",  "mag0 B0 [T]",  "mag0 cent_theta[deg]", "mag0 length_theta[deg]", 
                            "mag1 B0 [T]",  "mag1 cent_theta[deg]", "mag1 length_theta[deg]", "mag2 B0 [T]",  "mag2 cent_theta[deg]", "mag2 length_theta[deg]"]

        #df = pd.read_csv(input_file, sep = "\s+", skiprows=1, names=['m0', 'q', 'x', 'y', 'z', 'p', 'ux', 'uy','uz'])#, header=None
        #calculate total momentum

        part_df_data = {"m0 [MeV*c^-2]" : m0,
                "q [unit charge]" : q,
                "x [m]" : x,
                "y [m]" : y,
                "z [m]" : z,
                "Ptot [MeV/c]" : Ptot,
                "ux []" : ux,
                "uy []" : uy,
                "uz []" : uz}
        
        part_df = pd.DataFrame(part_df_data)
        part_df.columns = fixfield_headers
        
        print("input_particle_file_name : "+input_particle_file_name+", ML_input_particle_sample_file_name"
                            +ML_input_particle_sample_file_name)
        
        out_form = "fixfield"
        writeFile(part_df, input_particle_file_name, out_form) #write particle def file so fixfield can run it
        
        FFA_sample_df_lst = []
        
        for i in range(0,FFA_num_in_beam):
            FFA_sample_df_lst.append(pd.DataFrame(part_df_data))
            
            FFA_sample_df_lst[i]["rmin [m]"] = FFA_df_lst[i]['r_min'][0]
            FFA_sample_df_lst[i]["rmax [m]"] = FFA_df_lst[i]['r_max'][0]
            FFA_sample_df_lst[i]["zmax [m]"] = FFA_df_lst[i]['z_max'][0]
            FFA_sample_df_lst[i]["Cell opening angle [deg]"] = cell_open_ang_lst[i]
            FFA_sample_df_lst[i]["r0 [m]"] = FFA_df_lst[i]['r_0'][0]
            FFA_sample_df_lst[i]["k_val"] = FFA_df_lst[i]['k'][0]
            FFA_sample_df_lst[i]["C1"] = FFA_df_lst[i]['C1'][0]
            FFA_sample_df_lst[i]["lambda0"] = FFA_df_lst[i]['lambda0'][0]
            
            for j in range(0,3):
                mag_name = "mag"+str(i)
                FFA_sample_df_lst[i][mag_name+" B0 [T]"] = FFA_df_lst[i]['B0'][j]
                FFA_sample_df_lst[i][mag_name+" cent_theta[deg]"] = FFA_df_lst[i]['cent_theta'][j]
                FFA_sample_df_lst[i][mag_name+" length_theta[deg]"] = FFA_df_lst[i]['length_theta'][j]
        
            FFA_sample_df_lst[i] = FFA_sample_df_lst[i].reindex(columns=FFA_sample_headers)

            print("FFA_def_file_name : "+str(FFA_def_file_name_lst[i]))
            
            text_file = open(FFA_def_file_name_lst[i], "w") #write 1st FFA def file
            text_file.write(FFA_def_string_lst[i])
            text_file.close()  

            writeFile(FFA_sample_df_lst[i], FFA_input_data_file_name_lst[i], out_form)#"ML_part") #write particle def file so fixfield can run it
   
            
        FFA_index = FFA_index + 1
    
if __name__ == "__main__" :
    print(sys.argv)
    if len(sys.argv) == 4:
        #main(3, 20, "data/ML_input_sample/")
        main(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3])#format of pion input file, format of input, format of output
    else:
        print("wrong input. Input format should be : Beamline_type_num, part_per_set, output_dir")