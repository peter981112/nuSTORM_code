"""
convert bdsim and LhARA output ROOT and .dat file to input file format of each
"""
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import struct            as strct
import numpy             as np
import pandas as pd
import math              as mth
import scipy
import os
import io
import sys
import time

import Particle as Prtcl
import BeamLine as BL
import BeamLineElement  as BLE
import BeamIO   as bmIO

import ROOT
import pybdsim

def writeFile(new_df, fileName):

    new_df.to_csv(fileName, sep=' ', index=False, header=None)

    return

def main(input_format, part_type):
    in_form = input_format.lower()
    particle_type = part_type.lower()
    out_df = pd.DataFrame()
    LhARA_columns = ['x', 'xp', 'y', 'yp', 'z', 'E']
    bdsim_columns = ['x','y','z','px','py','pz','xp','yp','E']
    fixfield_columns = ['m0', 'q', 'x', 'y', 'z', 'p', 'ux', 'uy','uz']
    #=================================   LhARA output convert   =============================================================
    LhARA_inputfile = "/home/wonjong/LhARAlinearOptics_workspace/31-UserDirectory_2026_Jan/99-Scratch/custom_pion_flux_test.dat"
    beamlinesummaryfile   = None

    #nEvts = 1000
    nEvts = 74434

    #bdsim_inputfile = "custom_pion_test.root"
    #bdsim_inputfile = "z_spread_5GeV.root"
    bdsim_inputfile = "batch_test_reinstall.root"
    #bdsim_inputfile = "batch_test.root"
    if in_form == "lhara":
        HOMEPATH    = os.getenv('HOMEPATH')
        print("Start LhARA anal")
        print("         ----> HOMEPATH:", HOMEPATH)
            
        #.. File handling:
        print("         ----> Check input and beamlinesummary files:")

        #.. ----> Input file, if specified:
        if LhARA_inputfile != None and not os.path.isfile(LhARA_inputfile):
            LhARA_inputfile = os.path.join(HOMEPATH, LhARA_inputfile)
        if not os.path.isfile(LhARA_inputfile): 
            print("             ----> Input file does not exist.")
            print("                   Exit.")
            sys.exit(1)

        ibmIOr = bmIO.BeamIO(None, LhARA_inputfile)
        #ibmIOr = bmIO.BeamIO("99-Scratch", "Data4Tests.dat")
        print("             ----> Input file:", LhARA_inputfile)

        print("     <---- Initialisation complete.")

        print("     ----> Read data file:")

        EndOfFile = False
        iEvt = 0
        iCnt = 0
        Scl  = 10

        starttime = time.time()

        while not EndOfFile:
            EndOfFile = ibmIOr.readBeamDataRecord()
            if not EndOfFile:
                iEvt += 1
                if (iEvt % Scl) == 0:
                    print("         ----> Read event ", iEvt)
                    iCnt += 1
                    if iCnt == 10:
                        iCnt = 1
                        Scl  = Scl * 10
            if iEvt < 1:
                print(BL.BeamLine.getinstance())
                print(Prtcl.Particle.getParticleInstances()[iEvt-1])
            if nEvts != None and iEvt == nEvts:
                break
            
        print("     <----", iEvt, "events read")

        print(" <---- Data-file reading done.")
            
        print(" ----> Plot progression:")

        #Prtcl.Particle.plotTraceSpaceProgression()
        font = {'family': 'serif', \
                'color':  'darkred' \
                }

        plt.rcParams["figure.figsize"] = (7.5, 10.)

        nLoc_LhARA   = []
        xLoc_LhARA   = []
        xpLoc_LhARA  = []
        yLoc_LhARA   = []
        ypLoc_LhARA  = []
        ELoc_LhARA   = []
        ELab_LhARA   = []
        zLoc_LhARA   = []
        rLoc_LhARA   = []
        tLoc_LhARA   = []
        survive_part_num_LhARA   = []
        Scl    = []
        s_lst_LhARA   = []

        iSrc = BLE.Source.getinstances()[0]
        logE = False
        if iSrc.getMode() == 0:
            logE = True
                
        nPrtcl = 0
        iLoc_max = 0
        for iPrtcl in Prtcl.Particle.getinstances():#get iLoc_max the total number of 
            nPrtcl += 1
            iLoc = -1
            for iTrcSpc in iPrtcl.getTraceSpace():
                iLoc += 1      
            iLoc_max = iLoc
            if nPrtcl>2:
                break

        nPrtcl = 0  
        print("iLoc_max : "+str(iLoc_max))      
        for iPrtcl in Prtcl.Particle.getinstances():
            particleMASS = 139.57061 #pion mass in MeV
            nPrtcl += 1
            pion_life_time = 2.6033e-8
            if isinstance(iPrtcl, Prtcl.ReferenceParticle):
                iRefPrtcl = iPrtcl
                s_lst_LhARA = iRefPrtcl.gets()
                print("s_lst_LhARA leng : "+str(len(s_lst_LhARA)))
                print("s_lst_LhARA iLoc_maxth ele : "+str(s_lst_LhARA[iLoc_max]))
                continue
            iLoc = -1
            ds = 0
            
            for iTrcSpc in iPrtcl.getTraceSpace():
                iLoc += 1
                if iLoc == iLoc_max:#only consider the final element
                    
                    p0 = mth.sqrt(np.dot(iRefPrtcl.getPrOut()[iLoc][:3], \
                                            iRefPrtcl.getPrOut()[iLoc][:3]))
                    E0  = iRefPrtcl.getPrOut()[iLoc][3]
                    b0  = p0/E0
                    E   = E0 + iPrtcl.getTraceSpace()[iLoc][5] * p0
                    gamma = E/particleMASS
                    p   = mth.sqrt(E**2 - particleMASS**2)
                    b   = p/E
                    E  -= particleMASS
                    D   = mth.sqrt(1. + \
                                    2.*iPrtcl.getTraceSpace()[iLoc][5]/b0 +
                                    iPrtcl.getTraceSpace()[iLoc][5]**2)
                    
                    t = iPrtcl.getTraceSpace()[iLoc][4]*b0 / (b*scipy.constants.c)

                    eps = ( iPrtcl.getTraceSpace()[iLoc][1]**2 +    \
                            iPrtcl.getTraceSpace()[iLoc][3]**2  ) / \
                            (2.*D**2)

                    if iLoc_max>0:
                        ds = iPrtcl.gets()[iLoc] - iPrtcl.gets()[iLoc-1]
                    #no longer need decay processing as it is now included in LhARA code
                    #if np.random.random()>np.exp(-ds/(scipy.constants.c*pion_life_time*np.sqrt(gamma**2-1))):
                        #break        


                    
                    """
                    eps = (p - p0) / p0
                    LhARA_columns = ['x', 'xp', 'y', 'yp', 'z', 'E']
                    """

                    xLoc_LhARA.append(iPrtcl.getTraceSpace()[iLoc][0])
                    xpLoc_LhARA.append(iPrtcl.getTraceSpace()[iLoc][1])
                    yLoc_LhARA.append(iPrtcl.getTraceSpace()[iLoc][2])
                    ypLoc_LhARA.append(iPrtcl.getTraceSpace()[iLoc][3])
                    zLoc_LhARA.append(iPrtcl.getTraceSpace()[iLoc][4]*b0)
                    ELoc_LhARA.append(iPrtcl.getTraceSpace()[iLoc][5])
                    ELab_LhARA.append(E)
                    tLoc_LhARA.append(t)
                    Scl.append(eps)
                
                
        print(len(xLoc_LhARA))
        data = {'x' : xLoc_LhARA,
                'xp': xpLoc_LhARA,
                'y' : yLoc_LhARA,
                'yp': ypLoc_LhARA, 
                'z' : zLoc_LhARA, 
                'E' : ELab_LhARA}

        out_df = pd.DataFrame(data)

        out_name = 'LhARA_final_ele_output.txt'
            
        print("writing to output LhARA file : "+out_name)

        writeFile(out_df, out_name)
        
    elif in_form == "bdsim":
        pybdsim.Data.LoadROOTLibraries()
        ROOT.gSystem.Load("librebdsim")
        ROOT.gSystem.Load("libbdsimRootEvent")
        
        nLoc_bdsim   = []
        xLoc_bdsim   = []
        xpLoc_bdsim  = []
        yLoc_bdsim   = []
        ypLoc_bdsim  = []
        zLoc_bdsim   = []
        rLoc_bdsim   = []
        tLoc_bdsim   = []
        ELoc_bdsim   = []
        ELab_bdsim   = []
        Scl    = []
        p = 0
        px = 0
        py = 0
        pz = 0

        ele_name_lst = []

        d = pybdsim.DataPandas.BDSIMOutput(bdsim_inputfile)
        iLoc_max = len(d.get_sampler_names())

        d = pybdsim.Data.Load(bdsim_inputfile)
        print("using bdsim file : "+bdsim_inputfile)
        #d = pybdsim.Data.Load("batch_test.root")
        prim_data = pybdsim.Data.SamplerData(d,0)
        tot_part_count = pybdsim.Data.GetNEventsInBDSIMFile(bdsim_inputfile)#len(prim_data.data['partID']) #total number of pions at the begining
        print("bdsim with "+str(tot_part_count)+" particles")
        starttime = time.time()
        #iLoc_max = 50#133
        nPrtcl = 0

        particleMASS = 0.13957061 #pion mass in GeV
        nPrtcl += 1
        pion_life_time = 2.6033e-8 #pion lifetime in sec

        d_data = pybdsim.Data.SamplerData(d,iLoc_max)

        ele_name_lst = d_data.samplerName
        #if np.random.random()>np.exp(-ds/(scipy.constants.c*pion_life_time*np.sqrt(gamma**2-1))):
            #break        

        #if iLoc_max > (len(xLoc_bdsim)-1):
            #nLoc_bdsim = (iPrtcl.getLoc_bdsimation())
            
        #eps = (p - p0) / p0
        
        xLoc_bdsim = d_data.data['x']
        xpLoc_bdsim = d_data.data['xp']
        yLoc_bdsim = d_data.data['y']
        ypLoc_bdsim = d_data.data['yp']
        ELoc_bdsim = d_data.data['energy']
        ELab_bdsim = d_data.data['kineticEnergy']
        #Scl = (eps)
        
        x = np.array(d_data.data['x'])
        y = np.array(d_data.data['y'])
        
        zLoc_bdsim = (d_data.data['S'])
        rLoc_bdsim = ((x**2+y**2)**0.5)#append transverse distance from z axis
        #tLoc_bdsim = (t)
        
        curr_part_num = len(d_data.data['partID'])
        part_ID_lst = np.array(d_data.data['partID'])
        part_index_lst = []
        for i in range(0,curr_part_num): #remove all the non pions for current stage of analysis. Also remove all particles that escpaed the beam pipe with radius 0.5m, and too low energy particles
            #proton 2212 electron 11 positron -11 gamma or photon 22 neutron 2112 pion + 211 pion - -211 pion 0 111 muon - 13 muon + -13
            if (particle_type == 'pion') or (particle_type == 'pi'):
                if (part_ID_lst[i] == 211) and (rLoc_bdsim[i]<0.35) and (ELab_bdsim[i]>4.3):
                    part_index_lst.append(i)
            elif (particle_type == 'muon') or (particle_type == 'mu'):
                if (part_ID_lst[i] == -13) and (rLoc_bdsim[i]<0.35) and (ELab_bdsim[i]>4.3):
                    part_index_lst.append(i)
                
        #print("number of surviving part : "+str(len(part_index_lst)))
        """
        eps = (p - p0) / p0
        bdsim_columns = ['x','y','z','px','py','pz','xp','yp','E']
        """
        xLoc_bdsim = xLoc_bdsim[part_index_lst]
        xpLoc_bdsim = xpLoc_bdsim[part_index_lst]
        yLoc_bdsim = yLoc_bdsim[part_index_lst]
        ypLoc_bdsim = ypLoc_bdsim[part_index_lst]
        ELoc_bdsim = ELoc_bdsim[part_index_lst]
        ELab_bdsim = ELab_bdsim[part_index_lst]
        ELab_bdsim = ELab_bdsim * 1000 #convert GeV to MeV for comp with LhARA
        rLoc_bdsim = rLoc_bdsim[part_index_lst]
        
        p = ELab_bdsim
        p = np.sqrt(p**2-particleMASS**2)#enrgy to mom
        px = p*xpLoc_bdsim
        py = p*ypLoc_bdsim
        pz = np.sqrt(p**2-px**2-py**2)

        #zLoc_bdsim = [i for i in zLoc_bdsim if i != 0] #remove 0 s values
        zLoc_bdsim = zLoc_bdsim[part_index_lst]
        
        nLoc_bdsim = (len(xLoc_bdsim))

        timepass = time.time()-starttime
        print("number of serviving pions : "+str(nLoc_bdsim))
        print("bdsim data processing done in "+str(timepass)+" sec for "+str(tot_part_count)+" particles. Start writting")

        data = {'x' : xLoc_bdsim*1e3,#to convert to cm from m
                'y' : yLoc_bdsim*1e3,
                'z' : zLoc_bdsim*1e3,
                'px' : px,
                'py' : py,
                'pz' : pz,
                'xp' : xpLoc_bdsim,
                'yp' : ypLoc_bdsim,
                'E' : ELab_bdsim}

        out_df = pd.DataFrame(data)

        out_name = 'bdsim_final_ele_'+particle_type+'_output.txt'
            
        print("writing to output LhARA file : "+out_name)

        writeFile(out_df, out_name)




if __name__ == "__main__" :
    if len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])#format of pion input file
    else:
        print("wrong input. Input format should be : pion_input_file format_of_input format_of_output")