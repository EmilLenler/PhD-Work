import numpy as np

def csch(x):
    return 1/np.sinh(x)

max_n = 10
T = 0.5*1e-3 #1mK
hbar = 1*1e-34#hbar
kb = 1.380*1e-23 #Boltzmann constant
omega_i = 272*2*np.pi*1e3
omega_o = 772*2*np.pi*1e3
beta = 1/(kb*T)
E0 = 1/2*hbar*omega_o
rho_LIST = []
Z_real = csch(beta*hbar*omega_o/2)/2*csch(beta*hbar*omega_i/2)/2



for i in range(max_n):
    i_phase_factor = np.exp(-beta*(hbar*omega_i*(i+1/2)))
    for j in range(max_n):
        o_phase_factor = np.exp(-beta*(hbar*omega_o*(j+1/2)))
        rho_LIST.append(i_phase_factor*o_phase_factor)
rho_array = np.array(rho_LIST,dtype = np.complex64)
Z_N = np.sum(rho_array)
rho_NORMAL = rho_array/Z_N

N = len(rho_NORMAL)
rho_MAT = np.zeros((N,N))


# for n in range(rho_NORMAL):
#     if n + N > 