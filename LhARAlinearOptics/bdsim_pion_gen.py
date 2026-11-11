'''
    bdsim_pion_gen.py
    
    generate bdsim pion data

'''


import os
import sys, getopt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import traceSpace as tSp
import scipy as sp
import random 

def writeFile(new_df, fileName):

    new_df.to_csv(fileName, sep=' ', index=False, header=None)

    return

def gaussian(x,mu,sigma):
    return np.random.normal(0,sigma)

def getParabolic(x,p1):
    ran = random.random()
    a = np.array( [ 1., 0., (-3.*p1*p1), (2.*p1*p1*p1*(2.*ran - 1.)) ] )
    r = np.roots(a)
    isol = 0
    for ri in r:
        if not isinstance(ri, complex):
            if ri >= -p1:
                if ri <= p1:
                    isol += 1
                    p = ri
    if isol != 1:
        raise Exception("Simulation.getParabolic; p multiply defined")

    return p

def main(argv):#names=['x','y','z','px','py','pz','xp','yp','E'] -> bdsim output format, columns=['x', 'xp', 'y', 'yp', 'z', 'E'] -> LhARA output format
    
    nSamples = int(argv) #number of particles to generate
    print("number of particles to generate : "+str(nSamples))
    gaussian_FWHM = 2e-9 #full width half max for beam pulse in time units of ns
    gaussian_sigma = gaussian_FWHM/2.35 #convert FWHM to sigma of gaussian
    gaussian_sigma = gaussian_sigma*3e8 #convert time sigma to spacial spread
    parabolic_width = 1e-9
    parabolic_width = parabolic_width*3e8
    r = 5. #radius of beam in cm
    
    x_lst = np.zeros(nSamples)
    y_lst = np.zeros(nSamples)
    z_lst = np.zeros(nSamples)
    
    px_lst = np.zeros(nSamples)
    py_lst = np.zeros(nSamples)
    pz_lst = np.zeros(nSamples)
    
    xp_lst = np.zeros(nSamples)
    yp_lst = np.zeros(nSamples)
    E_lst = np.linspace(4.5,5.5,nSamples) #uniform distribution of E between 4.5-5.5GeV
    
    r_lst = r*np.random.random(nSamples)
    theta_lst = 2*np.pi*np.random.random(nSamples)
    x_lst = r_lst*np.cos(theta_lst)
    y_lst = r_lst*np.sin(theta_lst)
            
    z_lst = np.vectorize(getParabolic)(z_lst,parabolic_width)#give parabolic distribution to z
    
    pion_bdsim_df = pd.DataFrame({'x':x_lst,'y':y_lst,'z':z_lst,'px':px_lst,'py':py_lst,'pz':pz_lst,'xp':xp_lst,'yp':yp_lst,'E': E_lst})
    
    #pion_bdsim_df['z'] = np.vectorize(getParabolic)(pion_bdsim_df['z'],parabolic_width)
    
    columns=['x', 'xp', 'y', 'yp', 'z', 'E']
    
    pion_LhARA_df = pion_bdsim_df[columns]
    pion_LhARA_df['x'] = pion_LhARA_df['x']/1e2
    pion_LhARA_df['y'] = pion_LhARA_df['y']/1e2
    pion_LhARA_df['z'] = pion_LhARA_df['z']/1e2
    
    bdsim_out_name = '/home/wonjong/LhARAlinearOptics_workspace/BDSIM_test_dir/bdsim-scratch/KDE_convert/custom_test_bdsim_pions.txt'
    LhARA_out_name = '/home/wonjong/LhARAlinearOptics_workspace/BDSIM_test_dir/bdsim-scratch/KDE_convert/custom_test_LhARA_pions.txt'
    
    pion_bdsim_df.to_csv(bdsim_out_name, sep=' ', index=False, header=None)
    pion_LhARA_df.to_csv(LhARA_out_name, sep=' ', index=False, header=None)
    
if __name__ == "__main__" :
    main(sys.argv[1])
