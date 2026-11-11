from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import struct            as strct
import numpy             as np
import math              as mth
import scipy
import os
import io
import sys
import time

import ROOT
import pybdsim

pybdsim.Data.LoadROOTLibraries()
ROOT.gSystem.Load("librebdsim")
ROOT.gSystem.Load("libbdsimRootEvent")

logE = False
font = {'family': 'serif', \
        'color':  'darkred' \
        }

plt.rcParams["figure.figsize"] = (7.5, 10.)

nLoc   = []
xLoc   = []
xpLoc  = []
yLoc   = []
ypLoc  = []
zLoc   = []
sLoc   = []
rLoc   = []
tLoc   = []
ELoc   = []
ELab   = []
Scl    = []

ele_name_lst = []

bdsim_file_name = "field_map_test_out.root"

#d = pybdsim.DataPandas.BDSIMOutput("custom_pion_test.root")
d = pybdsim.DataPandas.BDSIMOutput(bdsim_file_name)
iLoc_max = len(d.get_sampler_names())#50#132

d = pybdsim.Data.Load(bdsim_file_name)
#d = pybdsim.Data.Load("custom_pion_test.root")
#d = pybdsim.Data.Load("z_spread_5GeV.root")
#d = pybdsim.Data.Load("couple_part_test.root")
#d = pybdsim.Data.Load("batch_test_reinstall_skim.root")
#d = pybdsim.Data.Load("batch_test_reinstall.root")
#d = pybdsim.Data.Load("batch_test.root")
prim_data = pybdsim.Data.SamplerData(d,0)
tot_part_count = len(prim_data.data['partID']) #total number of pions at the begining
print("plotting with "+str(tot_part_count)+" particles")
starttime = time.time()
nPrtcl = 0

particleMASS = 0.13957061 #pion mass in GeV
nPrtcl += 1
pion_life_time = 2.6033e-8 #pion lifetime in sec

for iLoc in range(0, iLoc_max):
    
    d_data = pybdsim.Data.SamplerData(d,iLoc)

    ele_name_lst.append(d_data.samplerName)
    #if np.random.random()>np.exp(-ds/(scipy.constants.c*pion_life_time*np.sqrt(gamma**2-1))):
        #break        

    #if iLoc > (len(xLoc)-1):
        #nLoc.append(iPrtcl.getLocation()[iLoc])
        
    #eps = (p - p0) / p0
    
    xLoc.append(d_data.data['x'])
    xpLoc.append(d_data.data['xp'])
    yLoc.append(d_data.data['y'])
    ypLoc.append(d_data.data['yp'])
    ELoc.append(d_data.data['energy'])
    ELab.append(d_data.data['kineticEnergy'])
    #Scl.append(eps)
    
    x = np.array(d_data.data['x'])
    y = np.array(d_data.data['y'])
    
    zLoc.append(d_data.data['T'])
    sLoc.append(d_data.data['S'])
    rLoc.append((x**2+y**2)**0.5)#append transverse distance from z axis
    #tLoc[iLoc].append(t)
    
    curr_part_num = len(d_data.data['partID'])
    part_ID_lst = np.array(d_data.data['partID'])
    pion_index_lst = []
    for i in range(0,curr_part_num): #remove all the non pions for current stage of analysis. Also remove all particles that escpaed the beam pipe with radius 0.5m, and too low energy particles
        if (part_ID_lst[i] == 211) and (rLoc[iLoc][i]<0.5) and (ELab[iLoc][i]>4.3):
            pion_index_lst.append(i)
            
    xLoc[iLoc] = xLoc[iLoc][pion_index_lst]
    xpLoc[iLoc] = xpLoc[iLoc][pion_index_lst]
    yLoc[iLoc] = yLoc[iLoc][pion_index_lst]
    ypLoc[iLoc] = ypLoc[iLoc][pion_index_lst]
    ELoc[iLoc] = ELoc[iLoc][pion_index_lst]
    ELab[iLoc] = ELab[iLoc][pion_index_lst]
    rLoc[iLoc] = rLoc[iLoc][pion_index_lst]
    zLoc[iLoc] = zLoc[iLoc][pion_index_lst]

    zLoc[iLoc] = [i*0.3 for i in zLoc[iLoc] if i != 0] #remove 0 s values
    
    nLoc.append(len(xLoc[iLoc]))
    if iLoc%10 == 0:
        print(str(iLoc)+" th beam line element processed")
        print(str(len(xLoc[iLoc]))+" particles remaining.")
        
timepass = time.time()-starttime
time_per_ele = timepass/iLoc_max
print("number of serviving pions : "+str(nLoc))
print("data processing done in "+str(timepass)+" sec. "+str(time_per_ele)+" sec per beamline element for "+str(tot_part_count)+" particles. Start plotting")
starttime = time.time()
plotFILE = 'BDSIM_ParticleProgressionPlot.pdf'
with PdfPages(plotFILE) as pdf:
    for iLoc in range(0, iLoc_max):
        fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(6., 6.), \
                                layout="constrained")
        # add an artist, in this case a nice label in the middle...
        Ttl = str(iLoc)+" th beam line element : "+str(ele_name_lst[iLoc])#nLoc[iLoc]
        fig.suptitle(Ttl, fontdict=font)

        #axs[0, 0].set_title('x,y')
        axs[0, 0].hist2d(xLoc[iLoc], yLoc[iLoc], bins=100)
        axs[0, 0].set_xlabel('x (m)')
        axs[0, 0].set_ylabel('y (m)')
    
        #axs[0, 1].set_title('delta')
        
        axs[0, 1].hist2d(rLoc[iLoc], ELab[iLoc], bins=100)#2d hist of energy and r
        axs[0, 1].set_xlabel('r (m)')
        axs[0, 1].set_ylabel('Kinetic energy (GeV)')
    
        #axs[1, 0].set_title('x, xprime')
        axs[1, 0].hist2d(xLoc[iLoc], xpLoc[iLoc], bins=100)
        axs[1, 0].set_xlabel('x (m)')
        axs[1, 0].set_ylabel('xprime (m)')

        #axs[1, 1].set_title('y, yprime')
        axs[1, 1].hist2d(yLoc[iLoc], ypLoc[iLoc], bins=100)
        axs[1, 1].set_xlabel('y (m)')
        axs[1, 1].set_ylabel('yprime (m)')

        axs[2, 0].hist(ELab[iLoc], 100)
        axs[2, 0].set_yscale('linear')
        if logE:
            axs[2, 0].set_yscale('log')
        axs[2, 0].set_xlabel('Kinetic energy (GeV)')
        axs[2, 0].set_ylabel('Number')

        axs[2, 1].hist(zLoc[iLoc], 100)
        axs[2, 1].set_yscale('linear')
        if logE:
            axs[2, 1].set_yscale('log')
        axs[2, 1].set_xlabel('s (m)')
        axs[2, 1].set_ylabel('Number')


        pdf.savefig()
        plt.close()
        if iLoc%10 == 0:
            print(str(iLoc)+" th beam line element plotted")
            
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6., 6.), layout="constrained")
    fig.suptitle('number of surviving pions', fontdict=font)
    
    z_loc_bdsim_0th_comp = [item[0] for item in sLoc]
    axs.plot(z_loc_bdsim_0th_comp, nLoc, label='bdsim')
    axs.legend()
    axs.set_xlabel('s (m)')
    axs.set_ylabel('surviving pion number')
    axs.set_ylim([0, nLoc[0]*1.15])
    pdf.savefig()
    plt.close()
timepass = time.time()-starttime
print("Plotting done in "+str(timepass)+" sec.")