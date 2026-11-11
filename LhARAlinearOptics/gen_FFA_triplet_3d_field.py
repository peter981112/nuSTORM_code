import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybdsim
import time
from subprocess import check_call as _check_call

def init_mag_para():
    return 0

def get_gen_azim_ang(theta, eta, r, r0):
    return theta - np.tan(eta)*np.log(r/r0)

def Fringe_func_approx(r0, theta_end, theta_strat, lambda0, C1):
    return 1/(1+np.exp(C1*r0*(theta_strat-theta_end)/lambda0)) #approx fringe field drop with simple angle instead of generalized azimuthal angle
    

def main():
   
    data = {'B0' : [-1.51226, 2.33691, -1.51226],
    'C1' : [2.240, 2.240, 2.240],
    'lambda0' : [0.695, 0.695, 0.695],
    'cent_theta' : [2.352, 4.250, 6.148], #angle in degree
    'length_theta' : [1.7, 1.870, 1.7],
    'k' : [72.341, 72.341, 72.341],
    'r_0' : [98.0493, 98.0493, 98.0493],
    'r_max' : [99.049, 99.049, 99.049],
    'r_min' : [97.049, 97.049, 97.049],
    'z_max' : [0.500, 0.500, 0.500]}
    
    cell_open_ang = 8.5
    
    df = pd.DataFrame(data)
    
    res = 100 #resolution of field map
    theta_res = 100 #theta resolution of field map
    
    # generate x and y points along their own axes
    #(x,y) cart coor with origin as mag center
    x = np.linspace(df['r_min'][0]-df['r_0'][0], df['r_max'][0]-df['r_0'][0], res)
    y = np.linspace(-df['z_max'][0], df['z_max'][0], res)
    #theta_lst = np.linspace(df['cent_theta'][0]-df['length_theta'][0]/2, df['cent_theta'][2]+df['length_theta'][2]/2, theta_res)
    theta_lst = np.linspace(df['cent_theta'][1]-cell_open_ang/2, df['cent_theta'][1]+cell_open_ang/2, theta_res)
    
    # define functions for each field component
    def fx(x,y):
        return 0
    def fy(x,y):
        return 0
    def fz(x, y, theta, B0, C1, lambda0, cent_theta, length_theta, k, r0, r_max, r_min, z_max):#FFA field
        theta_rad = np.deg2rad(theta)
        theta_open = cent_theta-length_theta/2
        theta_close = cent_theta+length_theta/2
        theta_open_rad = np.deg2rad(theta_open)
        theta_close_rad = np.deg2rad(theta_close)
        r = np.sqrt(y**2+(x+r0)**2)
        z = y
        fringe_drop = Fringe_func_approx(r0, theta_rad, theta_open_rad, lambda0, C1)*Fringe_func_approx(r0, theta_close_rad, theta_rad, lambda0, C1)
        #print("fringe_drop - "+str(fringe_drop)+" r/r0 - "+str(r/r0))
        return B0*fringe_drop*(r/r0)**k
    
    data = []
    theta = df['cent_theta'][1]
    pars = [df.iloc[0].to_numpy(), df.iloc[1].to_numpy(), df.iloc[2].to_numpy()]#[B0, r_0, theta, cent_theta-length_theta/2, cent_theta+length_theta/2, lambda0, C1, k]
    # loop over and build up 3d lists of lists of lists
    count = 0
    starttime = time.time()
    for thetai in range(theta_res):
        z_at_i = df['r_0'][0]*thetai*cell_open_ang*np.pi/(theta_res*180)
        count = count + 1
        data_2d = []
        for xi in range(res):
            v = np.zeros(6*res).reshape((res,6))
            for yi in range(res):
                B_vec = np.zeros(3)
                for mag_index in range(3):
                    B_vec = B_vec + [fx(x[xi],y[yi]),fz(x[xi],y[yi],theta_lst[thetai],*pars[mag_index]),fy(x[xi],y[yi])]
                v[yi] = [x[xi],y[yi],z_at_i,*B_vec]
            data_2d.append(v)
        data.append(data_2d)
        if count%(theta_res/5) == 0:
            print("field map "+str(100*count/theta_res)+" % done")

    # convert to numpy array
    data = np.array(data)
    timepass = time.time()-starttime
    print("field map at theta of "+str(theta)+" degree")
    print("3D field map of "+str(res)+"x"+str(res)+"x"+str(theta_res)+" generated in "+str(timepass)+" sec. ")#+str(time_per_ele)+" sec per beamline element for "+str(tot_part_count)+" particles. Start plotting")
    
    # construct a BDSIM format field object and write it out
    f = pybdsim.Field.Field3D(data)
    f.Write('3dexample.dat')

    fd = pybdsim.Field.Field2D(data, doublePrecision=True)
    fd.Write('3dexample_dp.dat')

    # compress the result
    _check_call(['gzip', "-f", "3dexample.dat"])
    _check_call(['gzip', "-f", "3dexample_dp.dat"])
    f.Write('3dexample.dat') # write again to keep original
    
    Plot(data[int(theta_res/2)])
    
    center_B_at_O = []
    for theta_val in theta_lst:
        B_at_O = 0
        for mag_index in range(3):
            B_at_O = B_at_O + fz(0,0,theta_val,*pars[mag_index])
        center_B_at_O.append(B_at_O)
        
    fig, ax = plt.subplots()
    ax.plot(df['r_0'][0]*np.deg2rad(theta_lst), center_B_at_O)
    ax.set(xlabel ='s (m)', ylabel='B_z (T)')
    fig.suptitle('B_z at center along FFA magnet')
    fig.savefig('3dexample_B_along_theta.png',dpi=400)
    
    return data

def Plot(array):

    a = array #shortcut
    #print(a)
    mag = np.sqrt(a[:,:,3]**2 + a[:,:,4]**2 + a[:,:,5]**2)
    print(mag)

    fig, ax = plt.subplots()

    p = ax.pcolor(a[:,:,0], a[:,:,1], mag, cmap=plt.cm.RdBu, vmin=abs(mag).min(), vmax=abs(mag).max())
    cb = fig.colorbar(p, ax=ax)
    cb.set_label('B(T)')

    #plt.figure()
    #plt.quiver(a[:,:,0], a[:,:,1], a[:,:,2], a[:,:,3], mag, cmap=plt.cm.magma, pivot='mid', scale=None)
    ax.set(xlabel ='X (m)', ylabel='Y (m)')
    fig.suptitle('2D Example Data')
    #plt.tight_layout()
    #plt.axes().set_aspect('equal', 'datalim')
    fig.savefig('2dexample.png',dpi=600)
    
if __name__ == "__main__":
    main()
