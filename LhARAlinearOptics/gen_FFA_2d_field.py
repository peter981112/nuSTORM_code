import matplotlib.pyplot as plt
import numpy as np
import pybdsim
import time
from subprocess import check_call as _check_call

def get_gen_azim_ang(theta, eta, r, r0):
    return theta - np.tan(eta)*np.log(r/r0)

def Fringe_func_approx(r0, theta_end, theta_strat, lambda0, C1):
    return 1/(1+np.exp(C1*r0*(theta_strat-theta_end)/lambda0)) #approx fringe field drop with simple angle instead of generalized azimuthal angle
    

def main():
   
    B0 = -1.50830
    C1 = 2.240
    lambda0 = 0.693
    cent_theta = 7.040 #angle in degree
    length_theta = 1.948
    k = 56.643
    r_0 = 92.1871
    r_max = 93.200
    r_min = 91.200
    z_max = 0.500

    res = 500 #resolution of field map
    
    # generate x and y points along their own axes
    #(x,y) cart coor with origin as mag center
    x = np.linspace(r_min-r_0, r_max-r_0, res)
    y = np.linspace(-z_max, z_max, res)
    theta_lst = np.linspace(cent_theta-length_theta/2, cent_theta+length_theta/2, res)
    
    # define functions for each field component
    def fx(x,y):
        return 0
    def fy(x,y):
        theta_rad = np.deg2rad(theta)
        theta_open_rad = np.deg2rad(theta_open)
        theta_close_rad = np.deg2rad(theta_close)
        r = np.sqrt(y**2+(x+r0)**2)
        z = y
        fringe_drop = Fringe_func_approx(r0, theta_rad, theta_open_rad, lambda0, C1)*Fringe_func_approx(r0, theta_close_rad, theta_rad, lambda0, C1)
        #print("fringe_drop - "+str(fringe_drop)+" r/r0 - "+str(r/r0))
        return B0*fringe_drop*(r/r0)**k
    def fz(x, y, B0, r0, theta, theta_open, theta_close, lambda0, C1, k):#FFA field
        return 0
    
    data = []
    theta = cent_theta
    pars = [B0, r_0, theta, cent_theta-length_theta/2, cent_theta+length_theta/2, lambda0, C1, k]
    # loop over and build up 3d lists of lists of lists
    count = 0
    starttime = time.time()
    for xi in x:
        count = count + 1
        v = []
        for yi in y:
            v.append([xi,yi,fx(xi,yi),fy(xi,yi),fz(xi,yi,*pars)])
        data.append(v)
        if count%(res/5) == 0:
            print("field map "+str(100*count/res)+" % done")

    # convert to numpy array
    data = np.array(data)
    timepass = time.time()-starttime
    print("field map at theta of "+str(theta)+" degree")
    print("2D field map of "+str(res)+"x"+str(res)+" generated in "+str(timepass)+" sec. ")#+str(time_per_ele)+" sec per beamline element for "+str(tot_part_count)+" particles. Start plotting")
    
    # construct a BDSIM format field object and write it out
    f = pybdsim.Field.Field2D(data)
    f.Write('2dexample.dat')

    fd = pybdsim.Field.Field2D(data, doublePrecision=True)
    fd.Write('2dexample_dp.dat')

    # compress the result
    _check_call(['gzip', "-f", "2dexample.dat"])
    _check_call(['gzip', "-f", "2dexample_dp.dat"])
    f.Write('2dexample.dat') # write again to keep original
    
    Plot(data)
    
    center_B_at_O = []
    for theta_val in theta_lst:
        pars = [B0, r_0, theta_val, cent_theta-length_theta/2, cent_theta+length_theta/2, lambda0, C1, k]
        B_at_O = fz(0,0,*pars)
        center_B_at_O.append(B_at_O)
        
    fig, ax = plt.subplots()
    ax.plot(r_0*np.deg2rad(theta_lst), center_B_at_O)
    ax.set(xlabel ='s (m)', ylabel='B_z (T)')
    fig.suptitle('B_z at center along FFA magnet')
    fig.savefig('2dexample_B_along_theta.png',dpi=400)
    
    return data

def Plot(array):

    a = array #shortcut
    #print(a)
    mag = np.sqrt(a[:,:,2]**2 + a[:,:,3]**2 + a[:,:,4]**2)
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
