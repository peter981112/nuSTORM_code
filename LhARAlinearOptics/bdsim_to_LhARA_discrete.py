'''
    bdsim_to_LhARA_discrete.py
    
    generate bdsim pion data

'''

import os
import sys, getopt
import pandas as pd
import numpy as np
import scipy


def writeFile(new_df, fileName):

    new_df.to_csv(fileName, sep=' ', index=False, header=None)

    return

def main(argv):#names=['x','y','z','px','py','pz','xp','yp','E'] -> bdsim output format, columns=['x', 'xp', 'y', 'yp', 'z', 'E'] -> LhARA output format
    bdsim_file = argv
    df = pd.read_csv(bdsim_file, sep = "\s+", header=None, names=['x','y','z','px','py','pz','xp','yp','E'])
    #calculate total momentum
    df['p'] = (np.sqrt(df['px']**2+df['py']**2+df['pz']**2)).to_list()
    
    columns=['x', 'xp', 'y', 'yp', 'z', 'E']
    pion_LhARA_df = df[columns]
    pion_LhARA_df['x'] = pion_LhARA_df['x']/1e2
    pion_LhARA_df['y'] = pion_LhARA_df['y']/1e2
    pion_LhARA_df['z'] = pion_LhARA_df['z']/1e2
    
    LhARA_out_name = '/home/wonjong/LhARAlinearOptics_workspace/31-UserDirectory/KDE_convert/custom_test_LhARA_pions.txt'
    
    print("writing to output LhARA file : "+LhARA_out_name)
    
    writeFile(pion_LhARA_df, LhARA_out_name)
    
    
    
if __name__ == "__main__" :
    main(sys.argv[1])
