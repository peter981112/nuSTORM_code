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

import Particle as Prtcl
import BeamLine as BL
import BeamLineElement  as BLE
import BeamIO   as bmIO

import ROOT
import pybdsim

#=================================   LhARA plotting   =============================================================
LhARA_inputfile = "/home/wonjong/LhARAlinearOptics_workspace/31-UserDirectory/99-Scratch/custom_pion_flux_test.dat"
beamlinesummaryfile   = None

#nEvts = 100
nEvts = 74434

#bdsim_inputfile = "custom_pion_test.root"
#bdsim_inputfile = "z_spread_5GeV.root"
bdsim_inputfile = "batch_test_reinstall.root"
#bdsim_inputfile = "batch_test.root"

# Two distributions to be compared:
Hist_LhARA_begin_prod_st_xy  = ROOT.TH2D("Hist_LhARA_begin_prod_st_xy", "Hist_LhARA_begin_prod_st_xy",  100, -0.2, 0.2, 100, -0.2, 0.2)
root_hist_lst_LhARA = []
Hist_bdsim_begin_prod_st_xy  = ROOT.TH2D("Hist_bdsim_begin_prod_st_xy", "Hist_bdsim_begin_prod_st_xy",  100, -0.2, 0.2, 100, -0.2, 0.2)
root_hist_lst_bdsim = []
KS_prob_xy_binned = 0
KS_prob_xy_binned_lst = []
KS_prob_xy_unbinned = 0

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
for iPrtcl in Prtcl.Particle.getinstances():
    particleMASS = 139.57061 #pion mass in MeV
    nPrtcl += 1
    pion_life_time = 2.6033e-8
    if isinstance(iPrtcl, Prtcl.ReferenceParticle):
        iRefPrtcl = iPrtcl
        s_lst_LhARA = iRefPrtcl.gets()
        print(s_lst_LhARA)
        continue
    iLoc = -1
    ds = 0
    
    for iTrcSpc in iPrtcl.getTraceSpace():
        iLoc += 1
        
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

        if iLoc>0:
            ds = iPrtcl.gets()[iLoc] - iPrtcl.gets()[iLoc-1]
        if np.random.random()>np.exp(-ds/(scipy.constants.c*pion_life_time*np.sqrt(gamma**2-1))):
            break        

        if iLoc > (len(xLoc_LhARA)-1):
            nLoc_LhARA.append(iPrtcl.getLocation()[iLoc])
            xLoc_LhARA.append([])
            xpLoc_LhARA.append([])
            yLoc_LhARA.append([])
            ypLoc_LhARA.append([])
            ELoc_LhARA.append([])
            ELab_LhARA.append([])
            Scl.append([])
            zLoc_LhARA.append([])
            rLoc_LhARA.append([])
            tLoc_LhARA.append([])

        
        """
        eps = (p - p0) / p0
        """
        x = iPrtcl.getTraceSpace()[iLoc][0]
        y = iPrtcl.getTraceSpace()[iLoc][2]
        
        if iLoc == 2:
            Hist_LhARA_begin_prod_st_xy.Fill(x,y)
        
        if (iLoc>=len(root_hist_lst_LhARA)) and (len(root_hist_lst_LhARA)<200):#nPrtcl == 1:
            hist_name = "Hist_LhARA_xy_at_" + str(iLoc) + "_th_ele"
            root_hist_lst_LhARA.append(ROOT.TH2D(hist_name, hist_name,  100, -0.5, 0.5, 100, -0.5, 0.5))
            #print(len(root_hist_lst_LhARA))
        root_hist_lst_LhARA[iLoc].Fill(x,y)
        
        xLoc_LhARA[iLoc].append(iPrtcl.getTraceSpace()[iLoc][0])
        xpLoc_LhARA[iLoc].append(iPrtcl.getTraceSpace()[iLoc][1])
        yLoc_LhARA[iLoc].append(iPrtcl.getTraceSpace()[iLoc][2])
        ypLoc_LhARA[iLoc].append(iPrtcl.getTraceSpace()[iLoc][3])
        zLoc_LhARA[iLoc].append(iPrtcl.getTraceSpace()[iLoc][4]*b0)
        rLoc_LhARA[iLoc].append((x**2+y**2)**0.5)#append transverse distance from z axis
        ELoc_LhARA[iLoc].append(iPrtcl.getTraceSpace()[iLoc][5])
        ELab_LhARA[iLoc].append(E)
        tLoc_LhARA[iLoc].append(t)
        Scl[iLoc].append(eps)

emit_x_LhARA = []#RMS emittance of beam
emit_y_LhARA = []
beta_x_LhARA = []#twiss para
gamma_x_LhARA = []
beta_y_LhARA = []
gamma_y_LhARA = []
for i in range(0,len(xLoc_LhARA)):
    survive_part_num_LhARA.append(len(xLoc_LhARA[i]))
    
    x_np_arr = np.array(xLoc_LhARA[i])
    xp_np_arr = np.array(xpLoc_LhARA[i])
    y_np_arr = np.array(yLoc_LhARA[i])
    yp_np_arr = np.array(ypLoc_LhARA[i])
    
    x_bar = sum(x_np_arr)/len(x_np_arr)
    x_bar_sq = sum(x_np_arr**2)/len(x_np_arr)
    xp_bar = sum(xp_np_arr)/len(xp_np_arr)
    xp_bar_sq = sum(xp_np_arr**2)/len(xp_np_arr)
    x_xp_bar = sum(x_np_arr*xp_np_arr)/len(x_np_arr)
    
    y_bar = sum(y_np_arr)/len(y_np_arr)
    y_bar_sq = sum(y_np_arr**2)/len(y_np_arr)
    yp_bar = sum(yp_np_arr)/len(yp_np_arr)
    yp_bar_sq = sum(yp_np_arr**2)/len(yp_np_arr)
    y_xp_bar = sum(y_np_arr*yp_np_arr)/len(x_np_arr)
    
    emit_x = np.sqrt(x_bar_sq*xp_bar_sq - x_xp_bar**2)
    emit_y = np.sqrt(y_bar_sq*yp_bar_sq - y_xp_bar**2)
    
    emit_x_LhARA.append(emit_x)
    emit_y_LhARA.append(emit_y)
    beta_x_LhARA.append((x_bar_sq-x_bar**2)/emit_x)
    beta_y_LhARA.append((y_bar_sq-y_bar**2)/emit_y)
    gamma_x_LhARA.append((xp_bar_sq-xp_bar**2)/emit_x)
    gamma_y_LhARA.append((yp_bar_sq-yp_bar**2)/emit_y)

print("emit_x_LhARA : "+str(emit_x_LhARA))
timepass = time.time()-starttime
time_per_ele = timepass/iLoc
print("number of serviving pions : "+str(survive_part_num_LhARA))
print("LhARA data processing done in "+str(timepass)+" sec. "+str(time_per_ele)+" sec per beamline element for "+str(survive_part_num_LhARA[0])+" particles.")        
#============================================   bdsim plotting   ==============================================

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

emit_x_bdsim = []#RMS emittance of beam
emit_y_bdsim = []
beta_x_bdsim = []#twiss para
gamma_x_bdsim = []
beta_y_bdsim = []
gamma_y_bdsim = []

ele_name_lst = []

d = pybdsim.DataPandas.BDSIMOutput(bdsim_inputfile)
iLoc_max = len(d.get_sampler_names())+1

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

for iLoc in range(0, iLoc_max):
    
    d_data = pybdsim.Data.SamplerData(d,iLoc)

    ele_name_lst.append(d_data.samplerName)
    #if np.random.random()>np.exp(-ds/(scipy.constants.c*pion_life_time*np.sqrt(gamma**2-1))):
        #break        

    #if iLoc > (len(xLoc_bdsim)-1):
        #nLoc_bdsim.append(iPrtcl.getLoc_bdsimation()[iLoc])
        
    #eps = (p - p0) / p0
    
    xLoc_bdsim.append(d_data.data['x'])
    xpLoc_bdsim.append(d_data.data['xp'])
    yLoc_bdsim.append(d_data.data['y'])
    ypLoc_bdsim.append(d_data.data['yp'])
    ELoc_bdsim.append(d_data.data['energy'])
    ELab_bdsim.append(d_data.data['kineticEnergy'])
    #Scl.append(eps)
    
    x = np.array(d_data.data['x'])
    y = np.array(d_data.data['y'])
    
    if iLoc == 2:
        for i in range(0,len(x)):
            Hist_bdsim_begin_prod_st_xy.Fill(x[i],y[i])
            
    hist_name = "Hist_bdsim_xy_at_" + str(iLoc) + "_th_ele"
    print(hist_name)
    root_hist_lst_bdsim.append(ROOT.TH2D(hist_name, hist_name,  100, -0.5, 0.5, 100, -0.5, 0.5))
    for i in range(0,len(x)):
        root_hist_lst_bdsim[iLoc].Fill(x[i],y[i])
    
    zLoc_bdsim.append(d_data.data['S'])
    rLoc_bdsim.append((x**2+y**2)**0.5)#append transverse distance from z axis
    #tLoc_bdsim[iLoc].append(t)
    
    curr_part_num = len(d_data.data['partID'])
    part_ID_lst = np.array(d_data.data['partID'])
    pion_index_lst = []
    for i in range(0,curr_part_num): #remove all the non pions for current stage of analysis. Also remove all particles that escpaed the beam pipe with radius 0.5m, and too low energy particles
        if (part_ID_lst[i] == 211) and (rLoc_bdsim[iLoc][i]<0.35) and (ELab_bdsim[iLoc][i]>4.3):
            pion_index_lst.append(i)
            
    xLoc_bdsim[iLoc] = xLoc_bdsim[iLoc][pion_index_lst]
    xpLoc_bdsim[iLoc] = xpLoc_bdsim[iLoc][pion_index_lst]
    yLoc_bdsim[iLoc] = yLoc_bdsim[iLoc][pion_index_lst]
    ypLoc_bdsim[iLoc] = ypLoc_bdsim[iLoc][pion_index_lst]
    ELoc_bdsim[iLoc] = ELoc_bdsim[iLoc][pion_index_lst]
    ELab_bdsim[iLoc] = ELab_bdsim[iLoc][pion_index_lst]
    ELab_bdsim[iLoc] = ELab_bdsim[iLoc] * 1000 #convert GeV to MeV for comp with LhARA
    rLoc_bdsim[iLoc] = rLoc_bdsim[iLoc][pion_index_lst]

    zLoc_bdsim[iLoc] = [i for i in zLoc_bdsim[iLoc] if i != 0] #remove 0 s values
    
    x_np_arr = np.array(xLoc_bdsim[iLoc])
    xp_np_arr = np.array(xpLoc_bdsim[iLoc])
    y_np_arr = np.array(yLoc_bdsim[iLoc])
    yp_np_arr = np.array(ypLoc_bdsim[iLoc])
    
    x_bar = sum(x_np_arr)/len(x_np_arr)
    x_bar_sq = sum(x_np_arr**2)/len(x_np_arr)
    xp_bar = sum(xp_np_arr)/len(xp_np_arr)
    xp_bar_sq = sum(xp_np_arr**2)/len(xp_np_arr)
    x_xp_bar = sum(x_np_arr*xp_np_arr)/len(x_np_arr)
    
    y_bar = sum(y_np_arr)/len(y_np_arr)
    y_bar_sq = sum(y_np_arr**2)/len(y_np_arr)
    yp_bar = sum(yp_np_arr)/len(yp_np_arr)
    yp_bar_sq = sum(yp_np_arr**2)/len(yp_np_arr)
    y_xp_bar = sum(y_np_arr*yp_np_arr)/len(x_np_arr)
    
    emit_x = np.sqrt(x_bar_sq*xp_bar_sq - x_xp_bar**2)
    emit_y = np.sqrt(y_bar_sq*yp_bar_sq - y_xp_bar**2)
    
    emit_x_bdsim.append(emit_x)
    emit_y_bdsim.append(emit_y)
    beta_x_bdsim.append((x_bar_sq-x_bar**2)/emit_x)
    beta_y_bdsim.append((y_bar_sq-y_bar**2)/emit_y)
    gamma_x_bdsim.append((xp_bar_sq-xp_bar**2)/emit_x)
    gamma_y_bdsim.append((yp_bar_sq-yp_bar**2)/emit_y)
    
    nLoc_bdsim.append(len(xLoc_bdsim[iLoc]))
    if iLoc%10 == 0:
        print(str(iLoc)+" th beam line element processed")
        print(str(len(xLoc_bdsim[iLoc]))+" particles remaining.")
        #print(str(len(zLoc_bdsim[iLoc])))
        
timepass = time.time()-starttime
time_per_ele = timepass/iLoc_max
zLoc_bdsim[0] = np.zeros(1) #to avoide no element in 0th comp
print([item[0] for item in zLoc_bdsim])
print("emit_x_bdsim : "+str(emit_x_bdsim))
print("number of serviving pions : "+str(nLoc_bdsim))
print("bdsim data processing done in "+str(timepass)+" sec. "+str(time_per_ele)+" sec per beamline element for "+str(tot_part_count)+" particles. Start plotting")
starttime = time.time()

KS_prob_xy_binned = Hist_bdsim_begin_prod_st_xy.KolmogorovTest(Hist_LhARA_begin_prod_st_xy, "UU")
print("KS prob at iLoc = 2 is : "+str(KS_prob_xy_binned))

iLoc_LhARA = 0        
plotFILE = "bdsim-LhARA_comp_ParticleProgressionPlot_diff.pdf"
s_lst_servive = []
servive_ratio = []
with PdfPages(plotFILE) as pdf:
    
    for iLoc in range(0, iLoc_max):
        fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(6., 6.), \
                                layout="constrained")
        # add an artist, in this case a nice label in the middle...
        while True: #get beamline element index for LhARA based on if the s value is close enough
            if iLoc == 0:
                break
            if np.abs(s_lst_LhARA[iLoc_LhARA]-zLoc_bdsim[iLoc][0])<1e-2:
                break
            iLoc_LhARA = iLoc_LhARA + 1
            
        s_lst_servive.append(s_lst_LhARA[iLoc_LhARA])
        servive_ratio.append(nLoc_bdsim[iLoc]/survive_part_num_LhARA[iLoc_LhARA])
        
        KS_prob_xy_binned_lst.append(root_hist_lst_bdsim[iLoc].KolmogorovTest(root_hist_lst_LhARA[iLoc_LhARA], "UU"))#get KS prob
        
        Ttl = nLoc_LhARA[iLoc_LhARA] + " at s = " + str(s_lst_LhARA[iLoc_LhARA]) + " m"
        fig.suptitle(Ttl, fontdict=font)
        
        #axs[1, 0].set_title('x, xprime')
        hist_xy = axs[0, 0].hist2d(xLoc_LhARA[iLoc], yLoc_LhARA[iLoc], bins=100)
        axs[0, 0].title.set_text('LhARA x-y')
        axs[0, 0].set_xlabel('x (m)')
        axs[0, 0].set_ylabel('y (m)')

        #axs[1, 1].set_title('y, yprime')
        axs[0, 1].hist2d(rLoc_LhARA[iLoc], ELab_LhARA[iLoc], bins=100)#2d hist of energy and r
        axs[0, 1].title.set_text('LhARA r-KE')
        axs[0, 1].set_xlabel('r (m)')
        axs[0, 1].set_ylabel('Kinetic energy (MeV)')
        
        #axs[0, 0].set_title('x,y')
        axs[1, 0].hist2d(xLoc_bdsim[iLoc], yLoc_bdsim[iLoc], bins=100, range = [[hist_xy[1][0], hist_xy[1][-1]], [hist_xy[2][0], hist_xy[2][-1]]])#, range = [[xmin, xmax], [ymin, ymax]]
        #print(hist_temp)
        axs[1, 0].title.set_text('bdsim x-y')
        axs[1, 0].set_xlabel('x (m)')
        axs[1, 0].set_ylabel('y (m)')
    
        #axs[0, 1].set_title('delta')
        
        axs[1, 1].hist2d(rLoc_bdsim[iLoc], ELab_bdsim[iLoc], bins=100)#2d hist of energy and r
        axs[1, 1].title.set_text('bdsim r-KE')
        axs[1, 1].set_xlabel('r (m)')
        axs[1, 1].set_ylabel('Kinetic energy (MeV)')

        axs[2, 0].hist(ELab_bdsim[iLoc], 100, alpha=0.75, histtype='step', label='bdsim')
        axs[2, 0].hist(ELab_LhARA[iLoc_LhARA], 100, alpha=0.75, histtype='step', label='LhARA')
        axs[2, 0].set_yscale('linear')
        if logE:
            axs[2, 0].set_yscale('log')
        axs[2, 0].legend()
        axs[2, 0].set_xlabel('Kinetic energy (MeV)')
        axs[2, 0].set_ylabel('Number')

        #axs[2, 1].hist(zLoc_bdsim[iLoc], 100, alpha=0.75, label='bdsim')
        axs[2, 1].hist(zLoc_LhARA[iLoc_LhARA], 100, alpha=0.75, label='LhARA')
        axs[2, 1].set_yscale('linear')
        if logE:
            axs[2, 1].set_yscale('log')
        axs[2, 1].legend()
        axs[2, 1].set_xlabel('s (m)')
        axs[2, 1].set_ylabel('Number')
        
        if iLoc%10 == 0:
            print(str(iLoc)+" th beam line element plotted")

        pdf.savefig()
        plt.close()
        
    print("Pion survival ratio : LhARA - "+str(100*survive_part_num_LhARA[-1]/survive_part_num_LhARA[0])+" %, bdsim : "+str(100*nLoc_bdsim[-1]/nLoc_bdsim[0])+" %")
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6., 6.), layout="constrained")
    fig.suptitle('number of surviving pions', fontdict=font)
    
    axs.plot(s_lst_LhARA, survive_part_num_LhARA, label='LhARA')
    z_loc_bdsim_0th_comp = [item[0] for item in zLoc_bdsim]
    axs.plot(z_loc_bdsim_0th_comp, nLoc_bdsim, label='bdsim')
    axs.legend()
    axs.set_xlabel('s (m)')
    axs.set_ylabel('surviving pion number')
    axs.set_ylim([0, nLoc_bdsim[0]*1.15])
    pdf.savefig()
    plt.close()
    
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6., 6.), layout="constrained")
    fig.suptitle('ratio of number of surviving pions bdsim/LhARA', fontdict=font)
    
    axs.plot(s_lst_servive, servive_ratio, label='ratio')
    axs.legend()
    axs.set_xlabel('s (m)')
    axs.set_ylabel('ratio of surviving pion number')
    axs.set_ylim([0, 1.15])
    pdf.savefig()
    plt.close()
    
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6., 6.), layout="constrained")
    fig.suptitle('emit_x', fontdict=font)
    
    axs.plot(s_lst_LhARA, emit_x_LhARA, label='LhARA')
    axs.plot(z_loc_bdsim_0th_comp, emit_x_bdsim, label='bdsim')
    axs.legend()
    axs.set_xlabel('s (m)')
    axs.set_ylabel('emit_x')
    #axs.set_ylim([0, nLoc_bdsim[0]*1.15])
    pdf.savefig()
    plt.close()
    
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6., 6.), layout="constrained")
    fig.suptitle('emit_y', fontdict=font)
    
    axs.plot(s_lst_LhARA, emit_y_LhARA, label='LhARA')
    axs.plot(z_loc_bdsim_0th_comp, emit_y_bdsim, label='bdsim')
    axs.legend()
    axs.set_xlabel('s (m)')
    axs.set_ylabel('emit_y')
    #axs.set_ylim([0, nLoc_bdsim[0]*1.15])
    pdf.savefig()
    plt.close()
    
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6., 6.), layout="constrained")
    fig.suptitle('beta_x', fontdict=font)
    
    axs.plot(s_lst_LhARA, beta_x_LhARA, label='LhARA')
    axs.plot(z_loc_bdsim_0th_comp, beta_x_bdsim, label='bdsim')
    axs.legend()
    axs.set_xlabel('s (m)')
    axs.set_ylabel('beta_x')
    #axs.set_ylim([0, nLoc_bdsim[0]*1.15])
    pdf.savefig()
    plt.close()
    
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6., 6.), layout="constrained")
    fig.suptitle('gamma_x', fontdict=font)
    
    axs.plot(s_lst_LhARA, gamma_x_LhARA, label='LhARA')
    axs.plot(z_loc_bdsim_0th_comp, gamma_x_bdsim, label='bdsim')
    axs.legend()
    axs.set_xlabel('s (m)')
    axs.set_ylabel('gamma_x')
    #axs.set_ylim([0, nLoc_bdsim[0]*1.15])
    pdf.savefig()
    plt.close()
    
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6., 6.), layout="constrained")
    fig.suptitle('beta_y', fontdict=font)
    
    axs.plot(s_lst_LhARA, beta_y_LhARA, label='LhARA')
    axs.plot(z_loc_bdsim_0th_comp, beta_y_bdsim, label='bdsim')
    axs.legend()
    axs.set_xlabel('s (m)')
    axs.set_ylabel('beta_y')
    #axs.set_ylim([0, nLoc_bdsim[0]*1.15])
    pdf.savefig()
    plt.close()
    
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6., 6.), layout="constrained")
    fig.suptitle('gamma_y', fontdict=font)
    
    axs.plot(s_lst_LhARA, gamma_y_LhARA, label='LhARA')
    axs.plot(z_loc_bdsim_0th_comp, gamma_y_bdsim, label='bdsim')
    axs.legend()
    axs.set_xlabel('s (m)')
    axs.set_ylabel('gamma_y')
    #axs.set_ylim([0, nLoc_bdsim[0]*1.15])
    pdf.savefig()
    plt.close()
    
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6., 6.), layout="constrained")
    fig.suptitle('Kolmogorov probability', fontdict=font)
    
    axs.plot(z_loc_bdsim_0th_comp, KS_prob_xy_binned_lst, label='Kolmogorov probability')
    axs.legend()
    axs.set_xlabel('s (m)')
    axs.set_ylabel('Kolmogorov probability')
    #axs.set_ylim([0, nLoc_bdsim[0]*1.15])
    pdf.savefig()
    plt.close()

    


iLoc_LhARA = 0        
plotFILE = "bdsim-LhARA_comp_ParticleProgressionPlot_Phase_space_diff.pdf"
with PdfPages(plotFILE) as pdf:
    
    for iLoc in range(0, iLoc_max):
        fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(6., 6.), \
                                layout="constrained")
        # add an artist, in this case a nice label in the middle...
        while True:
            if iLoc == 0:
                break
            if np.abs(s_lst_LhARA[iLoc_LhARA]-zLoc_bdsim[iLoc][0])<1e-2:
                break
            iLoc_LhARA = iLoc_LhARA + 1
            
        Ttl = nLoc_LhARA[iLoc_LhARA] + " at s = " + str(s_lst_LhARA[iLoc_LhARA]) + " m"
        fig.suptitle(Ttl, fontdict=font)
        
        #axs[0, 0].set_title('x,y')
        axs[0, 0].hist(xLoc_bdsim[iLoc], 100, alpha=0.75, histtype='step', label='bdsim')
        axs[0, 0].hist(xLoc_LhARA[iLoc_LhARA], 100, alpha=0.75, histtype='step', label='LhARA')
        axs[0, 0].legend()
        axs[0, 0].title.set_text('bdsim and LhARA x')
        axs[0, 0].set_xlabel('x (m)')
        axs[0, 0].set_ylabel('N')
    
        #axs[0, 1].set_title('delta')
        
        axs[0, 1].hist(yLoc_bdsim[iLoc], 100, alpha=0.75, histtype='step', label='bdsim')
        axs[0, 1].hist(yLoc_LhARA[iLoc_LhARA], 100, alpha=0.75, histtype='step', label='LhARA')
        axs[0, 1].legend()
        axs[0, 1].title.set_text('bdsim and LhARA y')
        axs[0, 1].set_xlabel('y (m)')
        axs[0, 1].set_ylabel('N')
    
        #axs[1, 0].set_title('x, xprime')
        hist_xxp = axs[1, 0].hist2d(xLoc_bdsim[iLoc], xpLoc_bdsim[iLoc], bins=100)
        axs[1, 0].title.set_text('bdsim x-xp')
        axs[1, 0].set_xlabel('x (m)')
        axs[1, 0].set_ylabel('xprime (m)')

        #axs[1, 1].set_title('y, yprime')
        axs[1, 1].hist2d(xLoc_LhARA[iLoc_LhARA], xpLoc_LhARA[iLoc_LhARA], bins=100, range = [[hist_xxp[1][0], hist_xxp[1][-1]], [hist_xxp[2][0], hist_xxp[2][-1]]])#2d hist of energy and r
        axs[1, 1].title.set_text('LhARA x-xp')
        axs[1, 1].set_xlabel('x (m)')
        axs[1, 1].set_ylabel('xprime (m)')

        hist_yyp = axs[2, 0].hist2d(yLoc_bdsim[iLoc], ypLoc_bdsim[iLoc], bins=100)
        axs[2, 0].title.set_text('bdsim y-yp')
        axs[2, 0].set_xlabel('y (m)')
        axs[2, 0].set_ylabel('yprime (m)')

        #axs[1, 1].set_title('y, yprime')
        axs[2, 1].hist2d(yLoc_LhARA[iLoc_LhARA], ypLoc_LhARA[iLoc_LhARA], bins=100, range = [[hist_yyp[1][0], hist_yyp[1][-1]], [hist_yyp[2][0], hist_yyp[2][-1]]])#2d hist of energy and r
        axs[2, 1].title.set_text('LhARA y-yp')
        axs[2, 1].set_xlabel('y (m)')
        axs[2, 1].set_ylabel('yprime (m)')
        
        if iLoc%10 == 0:
            print(str(iLoc)+" th beam line element plotted")

        pdf.savefig()
        plt.close()
        
    print("Pion survival ratio : LhARA - "+str(100*survive_part_num_LhARA[-1]/survive_part_num_LhARA[0])+" %, bdsim : "+str(100*nLoc_bdsim[-1]/nLoc_bdsim[0])+" %")
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6., 6.), layout="constrained")
    fig.suptitle('number of surviving pions', fontdict=font)
    
    axs.plot(s_lst_LhARA, survive_part_num_LhARA, label='LhARA')
    z_loc_bdsim_0th_comp = [item[0] for item in zLoc_bdsim]
    axs.plot(z_loc_bdsim_0th_comp, nLoc_bdsim, label='bdsim')
    axs.legend()
    axs.set_xlabel('s (m)')
    axs.set_ylabel('surviving pion number')
    axs.set_ylim([0, nLoc_bdsim[0]*1.15])
    pdf.savefig()
    plt.close()
    
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6., 6.), layout="constrained")
    fig.suptitle('ratio of number of surviving pions bdsim/LhARA', fontdict=font)
    
    axs.plot(s_lst_servive, servive_ratio, label='ratio')
    axs.legend()
    axs.set_xlabel('s (m)')
    axs.set_ylabel('ratio of surviving pion number')
    axs.set_ylim([0, 1.15])
    pdf.savefig()
    plt.close()
    
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6., 6.), layout="constrained")
    fig.suptitle('emit_x', fontdict=font)
    
    axs.plot(s_lst_LhARA, emit_x_LhARA, label='LhARA')
    axs.plot(z_loc_bdsim_0th_comp, emit_x_bdsim, label='bdsim')
    axs.legend()
    axs.set_xlabel('s (m)')
    axs.set_ylabel('emit_x')
    #axs.set_ylim([0, nLoc_bdsim[0]*1.15])
    pdf.savefig()
    plt.close()
    
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6., 6.), layout="constrained")
    fig.suptitle('emit_y', fontdict=font)
    
    axs.plot(s_lst_LhARA, emit_y_LhARA, label='LhARA')
    axs.plot(z_loc_bdsim_0th_comp, emit_y_bdsim, label='bdsim')
    axs.legend()
    axs.set_xlabel('s (m)')
    axs.set_ylabel('emit_y')
    #axs.set_ylim([0, nLoc_bdsim[0]*1.15])
    pdf.savefig()
    plt.close()
    
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6., 6.), layout="constrained")
    fig.suptitle('beta_x', fontdict=font)
    
    axs.plot(s_lst_LhARA, beta_x_LhARA, label='LhARA')
    axs.plot(z_loc_bdsim_0th_comp, beta_x_bdsim, label='bdsim')
    axs.legend()
    axs.set_xlabel('s (m)')
    axs.set_ylabel('beta_x')
    #axs.set_ylim([0, nLoc_bdsim[0]*1.15])
    pdf.savefig()
    plt.close()
    
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6., 6.), layout="constrained")
    fig.suptitle('gamma_x', fontdict=font)
    
    axs.plot(s_lst_LhARA, gamma_x_LhARA, label='LhARA')
    axs.plot(z_loc_bdsim_0th_comp, gamma_x_bdsim, label='bdsim')
    axs.legend()
    axs.set_xlabel('s (m)')
    axs.set_ylabel('gamma_x')
    #axs.set_ylim([0, nLoc_bdsim[0]*1.15])
    pdf.savefig()
    plt.close()
    
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6., 6.), layout="constrained")
    fig.suptitle('beta_y', fontdict=font)
    
    axs.plot(s_lst_LhARA, beta_y_LhARA, label='LhARA')
    axs.plot(z_loc_bdsim_0th_comp, beta_y_bdsim, label='bdsim')
    axs.legend()
    axs.set_xlabel('s (m)')
    axs.set_ylabel('beta_y')
    #axs.set_ylim([0, nLoc_bdsim[0]*1.15])
    pdf.savefig()
    plt.close()
    
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6., 6.), layout="constrained")
    fig.suptitle('gamma_y', fontdict=font)
    
    axs.plot(s_lst_LhARA, gamma_y_LhARA, label='LhARA')
    axs.plot(z_loc_bdsim_0th_comp, gamma_y_bdsim, label='bdsim')
    axs.legend()
    axs.set_xlabel('s (m)')
    axs.set_ylabel('gamma_y')
    #axs.set_ylim([0, nLoc_bdsim[0]*1.15])
    pdf.savefig()
    plt.close()
    
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6., 6.), layout="constrained")
    fig.suptitle('Kolmogorov probability', fontdict=font)
    
    axs.plot(z_loc_bdsim_0th_comp, KS_prob_xy_binned_lst, label='Kolmogorov probability')
    axs.legend()
    axs.set_xlabel('s (m)')
    axs.set_ylabel('Kolmogorov probability')
    #axs.set_ylim([0, nLoc_bdsim[0]*1.15])
    pdf.savefig()
    plt.close()
        
        
print("KS_prob_xy_binned_lst - ")
print(KS_prob_xy_binned_lst)   

timepass = time.time()-starttime
print("Plotting done in "+str(timepass)+" sec.")