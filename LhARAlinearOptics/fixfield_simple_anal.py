"""
simple plotting for fixfield
just start and end
"""
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import struct            as strct
import numpy             as np
import pandas as pd
import scipy
import os
import io
import sys
import time

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

#fixfield_file_name = ["custom_test_pions.txt", "custum_test_out.dat"]
fixfield_file_name = ["custom_test_pion_bdsim.txt", "custum_test_out.dat"]
#fixfield_file_name = ["custom_test_muon_bdsim.txt", "custum_test_out.dat"]

print("input : "+fixfield_file_name[0])

iLoc_max = len(fixfield_file_name) #afterwards to be fixed to num of element

df = pd.read_csv(fixfield_file_name[0], sep = "\s+", skiprows=1, names=['m0', 'q', 'x', 'y', 'z', 'p', 'ux', 'uy','uz'])

tot_part_count = len(df['m0']) #total number of pions at the begining
print("plotting with "+str(tot_part_count)+" particles")

starttime = time.time()
nPrtcl = 0

particleMASS = 139.57061 #pion mass in MeV
nPrtcl += 1
pion_life_time = 2.6033e-8

ds = 0

for iLoc in range(iLoc_max):
    
    df = pd.read_csv(fixfield_file_name[iLoc], sep = "\s+", skiprows=1, names=['m0', 'q', 'x', 'y', 'z', 'p', 'ux', 'uy','uz'])
    
    p   = df['p'].to_numpy()/1e3 #MeV to GeV
    m   = df['m0'].to_numpy()/1e3 #MeV to GeV
    E   = np.sqrt(p**2+m**2)
    gamma = E/particleMASS
    b   = p/E
    E  -= m#kinetic energy
    
    xp = df['ux'].to_numpy()
    yp = df['uz'].to_numpy()
    uz = df['uy'].to_numpy()#in fixfield y aixis is beam axis
    
    xp = xp/uz
    yp = yp/uz

    
    t = df['z'].to_numpy() / (b*scipy.constants.c)


    #if iLoc>0:
        #ds = iPrtcl.gets()[iLoc] - iPrtcl.gets()[iLoc-1]
    #if np.random.random()>np.exp(-ds/(scipy.constants.c*pion_life_time*np.sqrt(gamma**2-1))):
        #break        

    if iLoc > (len(xLoc)-1):
        nLoc.append(len(df['m0'].to_numpy()))

    
    """
    eps = (p - p0) / p0
    """
    x = df['x'].to_numpy()
    y = df['z'].to_numpy()#for fixfield y is the beam axis
    
    
    xLoc.append(x)
    xpLoc.append(xp)
    yLoc.append(y)
    ypLoc.append(yp)
    zLoc.append(df['y'].to_numpy())
    rLoc.append((x**2+y**2)**0.5)#append transverse distance from z axis
    #ELoc.append(iPrtcl.getTraceSpace()[iLoc][5])
    ELab.append(E)
    tLoc.append(t)
    #Scl.append(eps)
        
timepass = time.time()-starttime
time_per_ele = timepass/iLoc_max
print("number of serviving pions : "+str(nLoc))
print("data processing done in "+str(timepass)+" sec. "+str(time_per_ele)+" sec per beamline element for "+str(tot_part_count)+" particles. Start plotting")
starttime = time.time()
plotFILE = 'fixfield_ParticleProgressionPlot.pdf'
with PdfPages(plotFILE) as pdf:
    for iLoc in range(0, iLoc_max):
        fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(6., 6.), \
                                layout="constrained")
        # add an artist, in this case a nice label in the middle...
        Ttl = str(iLoc)+" th beam line element : "#+str(ele_name_lst[iLoc])#nLoc[iLoc]
        if iLoc == 0:
            Ttl = Ttl + "source"
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
    
    sLoc = [zLoc[0],[420]]
    z_loc_bdsim_0th_comp = [item[0] for item in sLoc]
    axs.plot(z_loc_bdsim_0th_comp, nLoc, label='fixfield')
    axs.legend()
    axs.set_xlabel('s (m)')
    axs.set_ylabel('surviving pion number')
    axs.set_ylim([0, nLoc[0]*1.15])
    pdf.savefig()
    plt.close()
timepass = time.time()-starttime
print("Plotting done in "+str(timepass)+" sec.")