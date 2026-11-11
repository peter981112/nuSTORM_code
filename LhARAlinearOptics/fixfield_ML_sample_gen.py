'''
    fixfield_ML_sample_gen.py
    
    generate input sample for ML machine learning 

'''

import os
import sys, getopt
import pandas as pd
import numpy as np
import scipy

m_pion = 139.57 #pion mass in MeV
m_muon = 105.66 #muon mass in MeV
m_part = 0
part_type = ''

def writeFile(new_df, fileName, output_format):
    if output_format == "fixfield":
        new_df.to_csv(fileName, sep=' ', index=False)#, header=None)
    else:
        new_df.to_csv(fileName, sep=' ', index=False, header=None)

    return
def main(FFA_type_num, part_per_set, output_dir):#energy_res, part_mass_res#names=['x','y','z','px','py','pz','xp','yp','E'] -> bdsim output format, columns=['x', 'xp', 'y', 'yp', 'z', 'E'] -> LhARA output format
    
    print("FFA_type_num : "+str(FFA_type_num)+", part_per_set : "+str(part_per_set)+", output_dir : "+str(output_dir))
    for FFA_num in range(0,FFA_type_num):
    
        FFA_def_file_name = "FFA_ML_test"
        FFA_def_file_ext = ".latt"
        input_particle_file_name = "ML_test_particle"
        input_particle_file_ext = ".dat"
        ML_input_particle_sample_file_name= "ML_test_particle_input"
        ML_input_particle_sample_file_ext = ".txt"
        
        FFA_def_file_name = output_dir + FFA_def_file_name + "-" + str(FFA_num) + FFA_def_file_ext
        input_particle_file_name = output_dir + input_particle_file_name + "-" + str(FFA_num) + input_particle_file_ext
        ML_input_particle_sample_file_name = output_dir + ML_input_particle_sample_file_name + "-" + str(FFA_num) + ML_input_particle_sample_file_ext
        
        FFA_type = 'ffag-r-enge'

        cell_open_ang = 15.0
        
        rand_bias = 0.2*(np.random.rand()-0.5)#put random bias to get different FFA samples for ML

        FFA_mag_B0 = [-2.0+rand_bias, 2.576994-rand_bias, -2.0+rand_bias]
        C1_val = 2.24+rand_bias
        lambda0_val = 0.7+0.2*rand_bias
        cent_theta = [3.6+rand_bias, 7.5+2*rand_bias, 11.4+3*rand_bias]
        length_theta = [1.5+0.7*rand_bias, 1.65+0.8*rand_bias, 1.5+0.7*rand_bias]
        k_val = 24.968629+10*rand_bias
        r_0_val = 48.898169+10*rand_bias
        r_max_val = 50.098833+10*rand_bias
        r_min_val = 48.098833+10*rand_bias
        z_max_val = 0.3+0.1*rand_bias

        
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
        

        #theta_lst = np.linspace(df['cent_theta'][0]-df['length_theta'][0]/2, df['cent_theta'][2]+df['length_theta'][2]/2, theta_res)
        
        #fixfield format m0 [MeV*c^-2], q [unit charge], x [m], y [m], z [m], Ptot [MeV/c], ux [], uy [], uz []
        form = "fixfield"
        
        Ptot_res = 25
        Ptot_min = 3.5e3#^3 for GeV
        Ptot_max = 6.5e3
        Ptot_lst = np.linspace(Ptot_min, Ptot_max, Ptot_res)
        dp_val = Ptot_lst[1]-Ptot_lst[0]
        
        part_mass_res = 15
        mass_min = 100#100MeV
        mass_max = 200#200MeV
        mass_lst = np.linspace(mass_min, mass_max, part_mass_res)
        
        x_z_max = 0.2 #unit m
        
        u_x_z_max = 0.02
        
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
                    m0[a] = mass_lst[m_index]
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
        ML_sample_df = pd.DataFrame()
        part_df = pd.DataFrame()
        fixfield_headers =  ["m0 [MeV*c^-2]", "q [unit charge]", "x [m]", "y [m]", "z [m]", "Ptot [MeV/c]", "ux []", "uy []", "uz []"]
        ML_sample_headers =  ["rmin [m]", "rmax [m]", "zmax [m]", "Cell opening angle [deg]", "r0 [m]", "k_val",
                            "C1", "lambda0",  "mag0 B0 [T]",  "mag0 cent_theta[deg]", "mag0 length_theta[deg]", 
                            "mag1 B0 [T]",  "mag1 cent_theta[deg]", "mag1 length_theta[deg]", "mag2 B0 [T]",  "mag2 cent_theta[deg]", "mag2 length_theta[deg]", 
                            "m0 [MeV*c^-2]", "q [unit charge]", "x [m]", "y [m]", "z [m]", "Ptot [MeV/c]", "ux []", "uy []", "uz []"]

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
        ML_sample_df = pd.DataFrame(part_df_data)
        
        ML_sample_df["rmin [m]"] = FFA_df['r_min'][0]
        ML_sample_df["rmax [m]"] = FFA_df['r_max'][0]
        ML_sample_df["zmax [m]"] = FFA_df['z_max'][0]
        ML_sample_df["Cell opening angle [deg]"] = cell_open_ang
        ML_sample_df["r0 [m]"] = FFA_df['r_0'][0]
        ML_sample_df["k_val"] = FFA_df['k'][0]
        ML_sample_df["C1"] = FFA_df['C1'][0]
        ML_sample_df["lambda0"] = FFA_df['lambda0'][0]
        
        for i in range(0,3):
            mag_name = "mag"+str(i)
            ML_sample_df[mag_name+" B0 [T]"] = FFA_df['B0'][i]
            ML_sample_df[mag_name+" cent_theta[deg]"] = FFA_df['cent_theta'][i]
            ML_sample_df[mag_name+" length_theta[deg]"] = FFA_df['length_theta'][i]
            
        
        part_df.columns = fixfield_headers
        ML_sample_df = ML_sample_df.reindex(columns=ML_sample_headers)

        print("FFA_def_file_name : "+FFA_def_file_name+", input_particle_file_name : "+input_particle_file_name+", ML_input_particle_sample_file_name"+ML_input_particle_sample_file_name)
        
        text_file = open(FFA_def_file_name, "w") #write FFA def file
        text_file.write(FFA_def_string)
        text_file.close()  
        
        out_form = "fixfield"
        writeFile(part_df, input_particle_file_name, out_form) #write particle def file so fixfield can run it
        
        writeFile(ML_sample_df, ML_input_particle_sample_file_name, out_form)#"ML_part") #write particle def file so fixfield can run it
    
if __name__ == "__main__" :
    if len(sys.argv) == 4:
        #main(3, 20, "data/ML_input_sample/")
        main(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3])#format of pion input file, format of input, format of output
    else:
        print("wrong input. Input format should be : pion_input_file format_of_input format_of_output")