'''
    pion_converter_all_format.py
    
    convert pion between bdsimm, LhARA, and fixfield

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
if ('muon' in sys.argv[1].lower()) or ('mu' in sys.argv[1].lower()):
    m_part = m_muon
    part_type = 'muon'
    
elif ('pion' in sys.argv[1].lower()) or ('pi' in sys.argv[1].lower()):
    m_part = m_pion
    part_type = 'pion'
    
print("creating "+part_type+" file")
    
def writeFile(new_df, fileName, output_format):
    if output_format == "fixfield":
        new_df.to_csv(fileName, sep=' ', index=False)#, header=None)
    else:
        new_df.to_csv(fileName, sep=' ', index=False, header=None)

    return

def to_bdsim(input_df, input_format):
    bdsim_df = input_df
    columns = ['x','y','z','px','py','pz','xp','yp','E']
    
    if input_format == "lhara":#['x', 'xp', 'y', 'yp', 'z', 'E']
        bdsim_df['x'] = bdsim_df['x']*1e1
        bdsim_df['y'] = bdsim_df['y']*1e1
        bdsim_df['z'] = bdsim_df['z']*1e1
        
        p = bdsim_df['E'].to_numpy()
        p = np.sqrt(p**2-(m_part/1e3)**2)#enrgy to mom
        px = p*bdsim_df['xp'].to_numpy()
        py = p*bdsim_df['yp'].to_numpy()
        pz = np.sqrt(p**2-px**2-py**2)
        
        bdsim_df['px'] = px
        bdsim_df['py'] = py
        bdsim_df['pz'] = pz
        
    elif input_format == "bdsim":#['x','y','z','px','py','pz','xp','yp','E'])
        a = 0
        
    elif input_format == "fixfield":#=['m0', 'q', 'x', 'y', 'z', 'Ptot', 'ux', 'uy','uz'])
        bdsim_df['x'] = bdsim_df['x']*1e3
        bdsim_df['y'] = bdsim_df['y']*1e3
        bdsim_df['z'] = bdsim_df['z']*1e3
        
        p = bdsim_df['p'].to_numpy()
        xp = bdsim_df['ux'].to_numpy()
        yp = bdsim_df['uz'].to_numpy()
        uz = bdsim_df['uy'].to_numpy()#in fixfield y aixis is beam axis
        
        p = p/1e3#MeV to GeV
        E = np.sqrt(p**2+(m_part/1e3)**2)
        
        xp = xp/uz
        yp = yp/uz
        
        pz = p*uz
        px = p*xp
        py = p*yp
        
        bdsim_df['E'] = E
        bdsim_df['px'] = px
        bdsim_df['py'] = py
        bdsim_df['pz'] = pz
        bdsim_df['xp'] = xp
        bdsim_df['yp'] = yp
        
    bdsim_df = bdsim_df[columns]
    return bdsim_df

def convert_to(input_df, output_format):
    output_df = input_df
    LhARA_columns = ['x', 'xp', 'y', 'yp', 'z', 'E']
    bdsim_columns = ['x','y','z','px','py','pz','xp','yp','E']
    fixfield_columns = ['m0', 'q', 'x', 'y', 'z', 'p', 'ux', 'uy','uz']
    
    if output_format == "lhara":#['x', 'xp', 'y', 'yp', 'z', 'E']
        output_df['x'] = output_df['x']/1e1
        output_df['y'] = output_df['y']/1e1
        output_df['z'] = output_df['z']/1e1
        
        output_df = output_df[LhARA_columns]
        return output_df
        
    elif output_format == "bdsim":#['x','y','z','px','py','pz','xp','yp','E'])
        output_df = output_df[bdsim_columns]
        return output_df
        
    elif output_format == "fixfield":#=['m0', 'q', 'x', 'y', 'z', 'Ptot', 'ux', 'uy','uz'])
        output_df['x'] = output_df['x']/1e3
        output_df['z'] = output_df['y']/1e3 #for fixfield y is beam axis
        output_df['y'] = np.zeros(len(output_df['xp'].to_numpy()))#output_df['z']/1e3
               
        ux = output_df['xp'].to_numpy()
        uy = output_df['yp'].to_numpy()
        uz = np.zeros(len(ux))
        m = np.zeros(len(ux))
        q = np.zeros(len(ux))
        
        p = output_df['E'].to_numpy()
        p = np.sqrt(p**2-(m_part/1e3)**2)#enrgy to mom
        
        for i in range(0,len(ux)):
            uz[i] = 1
            m[i] = m_part
            q[i] = 1
            
        norm_factor = np.sqrt(ux**2+uy**2+uz**2)
        ux = ux/norm_factor
        uy = uy/norm_factor
        uz = uz/norm_factor
        
        output_df['m0'] = m
        output_df['q'] = q
        output_df['p'] = p
        output_df['ux'] = ux
        output_df['uy'] = uz
        output_df['uz'] = uy
        
        output_df = output_df[fixfield_columns]
        return output_df

def main(pion_input, input_format, output_format):#names=['x','y','z','px','py','pz','xp','yp','E'] -> bdsim output format, columns=['x', 'xp', 'y', 'yp', 'z', 'E'] -> LhARA output format
    #LhARA format ['x', 'xp', 'y', 'yp', 'z', 'E']
    #bdsim format ['x','y','z','px','py','pz','xp','yp','E']
    #LhARA and bdsim E in GeV
    #fixfield format m0 [MeV*c^-2], q [unit charge], x [m], y [m], z [m], Ptot [MeV/c], ux [], uy [], uz []
    input_file = pion_input
    in_form = input_format.lower()
    out_form = output_format.lower()
    df = pd.DataFrame()
    out_df = pd.DataFrame()
    fixfield_headers =  ["m0 [MeV*c^-2]", "q [unit charge]", "x [m]", "y [m]", "z [m]", "Ptot [MeV/c]", "ux []", "uy []", "uz []"]
    if in_form == "lhara":
        df = pd.read_csv(input_file, sep = "\s+", header=None, names=['x', 'xp', 'y', 'yp', 'z', 'E'])
    elif in_form == "bdsim":
        df = pd.read_csv(input_file, sep = "\s+", header=None, names=['x','y','z','px','py','pz','xp','yp','E'])
    elif in_form == "fixfield":
        df = pd.read_csv(input_file, sep = "\s+", skiprows=1, names=['m0', 'q', 'x', 'y', 'z', 'p', 'ux', 'uy','uz'])#, header=None
    #calculate total momentum
    bdsim_df = to_bdsim(df, in_form)
    out_df = convert_to(bdsim_df, out_form)
    
    if out_form == "fixfield":
        print("out_form = "+out_form)
        out_df.columns = fixfield_headers

    out_name = 'custom_test_'+part_type+'_'+in_form+'.txt'
    
    print("writing to output "+out_form+" file : "+out_name)
    
    writeFile(out_df, out_name, out_form)
    
    
if __name__ == "__main__" :
    if len(sys.argv) == 4:
        main(sys.argv[1], sys.argv[2], sys.argv[3])#format of pion input file, format of input, format of output
    else:
        print("wrong input. Input format should be : pion_input_file format_of_input format_of_output")
