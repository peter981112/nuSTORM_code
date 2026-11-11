'''
    bdsim_convert.py
    
    modified from Paul Kyberd's KDEdataFromFluka.py code to modify bdsim pion data

'''

import os
import sys, getopt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import traceSpace as tSp
import scipy as sp
import random 

class bdsim_convert:

    def __init__(self):
        print("instantiate")
 
    def beamGen_bdsim(self, data, columns, num_samples):
        bdsim_data = data[columns]
    # Compute Scott's factor for KDE bandwidth
        scott_factor = bdsim_data.shape[0]**(-1./(bdsim_data.shape[1]+4))
        scott_factor /= 10 # Scott's factor is too large for our purposes, so we reduce it

    # Estimate underlying distribution of the data using KDE
        bdsim = sp.stats.gaussian_kde(bdsim_data.T,bw_method=scott_factor)

    # Generate samples from the KDE distribution
        samples = bdsim.resample((num_samples,)).T

    # Create a new dataframe to store the samples
        bdsim_sample_df = pd.DataFrame(samples, columns=columns)
    
        return bdsim_sample_df


    def bdsim_plotting(self, df,new_df):
    # quick plotting to compare the original and the KDE generated data 
        plt.hist(df['x'][:10000], bins=60, range = [-15,15], alpha=0.5, label='Original')
        plt.hist(new_df['x'][:10000], bins=60, range = [-15,15], alpha=0.5, label='KDE')
        plt.legend(loc='upper right')
        plt.show()

        plt.hist(df['xp'][:10000], bins=60, range = [-0.05,0.05], alpha=0.5, label='Original')
        plt.hist(new_df['xp'][:10000], bins=60, range = [-0.05,0.05], alpha=0.5, label='KDE')
        plt.legend(loc='upper right')
        plt.show()

        plt.hist(df['p'][:10000], bins=60, range = [1.7,2.3], alpha=0.5, label='Original')
#        plt.hist(new_df['p'][:10000], bins=60, range = [1.7,2.3], alpha=0.5, label='KDE')
        plt.legend(loc='upper right')
        plt.show()

        plt.hist(df['p'][:10000], bins=60, range = [1.7,2.3], alpha=0.5, label='Original')
#        plt.hist(new_df['p'][:10000], bins=60, range = [1.7,2.3], alpha=0.5, label='KDE')
        plt.legend(loc='upper right')
        plt.show()

        plt.hist2d(df['x'][:10000],df['xp'][:10000],bins=[50,50], range = [[-15,15],[-0.05,0.05]])
        plt.title('Original data x-px phase space')
        plt.show()
        plt.hist2d(new_df['x'][:10000], new_df['xp'][:10000],bins=[50,50], range = [[-15,15],[-0.05,0.05]])
        plt.title('KDE data x-px phase space')
        plt.show()     
          
    def readFile(self, flukaDataFile):
    #load Fluka data into a pandas dataframe
        df = pd.read_csv(flukaDataFile, sep = "\s+", header=None, names=['x','y','z','px','py','pz','xp','yp','E'])
    #calculate total momentum
        df['p'] = (np.sqrt(df['px']**2+df['py']**2+df['pz']**2)).to_list()
        df['z'] = 0.000001

        return df

    def writeFile(self, new_df, fileName):

        new_df.to_csv(fileName, sep=' ', index=False, header=None)

        return
 
    def openGenFile(self, fileName):
        self._genFile = open(fileName, "r")

    def readGenFile(self):

        data = self._genFile.readline()
        data = data.strip()
        data = data.split()
        if not data: return None, None
        #   create a trace space vector, and alter the distance from (assumed) mm to m
        pion = tSp.traceSpace(0.0, float(data[0])/100.0, float(data[1])/100.0, -50.0, float(data[2]), float(data[3]))
        pPion = float(data[5])

        return pion, pPion
 
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

def main(argv):

    nSamples = 10000
    #   Instantiate the KDE class
    flukaDataFile = argv #'selected_pions_7.2GeV_2.0mm.txt'
    print("data source file for BDSIM : "+flukaDataFile)
    gaussian_FWHM = 2e-9 #full width half max for beam pulse in time units of ns
    gaussian_sigma = gaussian_FWHM/2.35 #convert FWHM to sigma of gaussian
    gaussian_sigma = gaussian_sigma*3e8 #convert time sigma to spacial spread
    parabolic_width = 1e-9
    parabolic_width = parabolic_width*3e8
    KDEdata = bdsim_convert()

    df = KDEdata.readFile(flukaDataFile)

    print ("df is \n", df)
    nSamples = len(df['x'].to_numpy())
    print("Converting "+str(nSamples)+" samples with bunch width of : "+str(parabolic_width/3e8)+" ns which is : "+str(parabolic_width)+" m")

    #num_samples is used set to the new sample size
    new_df = KDEdata.beamGen_bdsim(df, columns=['x','y','z','px','py','pz','xp','yp','E'], num_samples=nSamples)
    
    #new_df['z'] = np.vectorize(gaussian)(new_df['z'],0,gaussian_sigma)
    new_df['z'] = np.vectorize(getParabolic)(new_df['z'],parabolic_width)

    print (" new_df is \n", new_df)
    #quick plotting for visual inspection
    #KDEdata.bdsim_plotting(df,new_df)

    #save the new data to a text file
    genDataFile = "/home/wonjong/LhARAlinearOptics_workspace/BDSIM_test_dir/bdsim-scratch/KDE_convert/converted_bdsim_pions.txt"
    print("writing new data file")
    KDEdata.writeFile(new_df, genDataFile)

    #   Read the generated data file
    print("reading new data file")
    KDEdata.openGenFile(genDataFile)
    
    print("mean energy of pions : "+str(new_df.loc[:, 'E'].mean())+" GeV")
    print("ref part KE : "+str((new_df.loc[:, 'E'].mean()-0.13957061)*1000)+" MeV")

    #   Read the generated data file
    for i in range(nSamples):
        pion, pPion =  KDEdata.readGenFile()
        if pion == None: break
        if i%3000==0:
            print (str(i)+" th main: pion is ", pion)
            print (str(i)+" th main: pPion is ", pPion)


if __name__ == "__main__" :
    main(sys.argv[1])
