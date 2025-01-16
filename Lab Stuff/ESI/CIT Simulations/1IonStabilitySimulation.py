import numpy as np
import matplotlib.pyplot as plt
import CIT_Tools as CIT_T
import scipy.special as special


omega_RF = 2* np.pi*300*1e3
z0 = 9*1e-3
r0 = 1e-2
xi = 2.36*1e-2

CylinderTrap = CIT_T.CIT(omega_RF,xi,z0,r0)

DC_V = 0
RF_V = 400



r_init = r0 * np.random.rand()
z_init = z0 * np.random.rand()


def l_n(n,trap):
    return np.pi * (2*n + 1) / (2*trap.z0)

def sum_j_r(r,z,j,trap):
    n = 0
    sum = 0
    while n < j:
        sum+= (-1)**n * special.iv(1,r*trap.pn(n)) / np.i0(trap.r0 * trap.pn(n)) * np.cos(trap.pn(n) * z)
        


def sum_j_z(r,z,j,trap):
    n = 0
    sum = 0
    while n < j:
        sum += (-1)**(n+1) * np.i0(r*trap.pn(n)) / np.i0(trap.r0 * trap.pn(n)) * np.sin(trap.pn(n) * z)

def dx2dt2(r,z,t,DC,RF,trap,m,q):
    
    
    dr2 = -2 / (trap.z0 * m/q) * (DC + RF*np.cos(trap.omega_RF * t)) * sum_j_r(r,z,10,trap)

    dz2 = -2 / (trap.z0 * m/q) * (DC + RF*np.cos(trap.omeag_RF * t)) * sum_j_z(r,z,10,trap)
    
    return np.array([dr2,dz2]) 


def VelocityVerlet(r_init,vr_init,z_init,vz_init,dt,start,stop,DC,RF,trap,m,q):
    currentZ = z_init
    currentVZ = vz_init

    currentR = r_init
    currentVR = vr_init
    now = start
    
    Z = []
    R = []
    times = []
    while now < stop:

        Z.append(currentZ)
        R.append(currentR)
        times.append(now)


        nextZ = currentZ + currentVZ * dt + 0.5 * dx2dt2(currentR,currentZ,now,DC,RF,trap,m,q)[1] * dt**2
        nextVZ = currentVZ + (dx2dt2(currentR,currentZ,now,DC,RF,trap,m,q)[1] + dx2dt2(currentR,currentZ,now+dt,DC,RF,trap,m,q)[1])/2 *dt

        nextR = currentR + currentVR * dt + 0.5 * dx2dt2(currentR,currentZ,now,DC,RF,trap,m,q)[0] * dt**2
        nextVR = currentVZ + (dx2dt2(currentR,currentZ,now,DC,RF,trap,m,q)[0] + dx2dt2(currentR,currentZ,now+dt,DC,RF,trap,m,q)[0])/2 *dt

        now += dt
    return Z, R, times

dt = 1e-7
ZS,RS,times = VelocityVerlet(0.1*1e-3,0,-0.1*1e-3,0,dt,0,3*1e-5,0,400,CylinderTrap,443 * 1.66*1e-27, 1.6*1e-19)
