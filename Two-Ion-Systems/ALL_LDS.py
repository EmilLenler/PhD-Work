import numpy as np
import TwoIons as TI
import matplotlib.pyplot as plt

plt.rc('font',size = 14)

m_ba = 135*1.66*1e-27
charge_ba = 1*1.6*1e-19
q_set = 0.4
a_set = -0.01
m_por = 9000*1.66*1e-27
charge_por = 24*1.6*1e-19
trap = TI.Trap(2*np.pi*5.2*1e6,0.248,3/4*2.7*1e-3,3/4*3.5*1e-3)
q_set = 0.4
a_set = -0.01

por_masses = [(i+1)/12*m_por for i in range(12)]
por_charges = [(i+1)/12*charge_por for i in range(12)]
Por_Systems = [TI.two_ion_system(m_ba,charge_ba,por_m,por_c,trap) for (por_m,por_c) in zip(por_masses,por_charges)]


# freq_inz = []
# freq_outz = []
# freq_inr = []
# freq_outr = []
# vec_inz = []
# vec_outz = []
# vec_inr = []
# vec_outr = []
k = 2*np.pi/(1762*1e-9)
LD_inz = []
LD_outz = []
LD_inr = []
LD_outr = []
for sys in Por_Systems:
    freqs, vecs = sys.vibrational_freqs_and_modes(a_set,q_set)
    # freq_inz.append(freqs[0])
    # freq_outz.append(freqs[1])
    # freq_inr.append(freqs[2])
    # freq_outr.append(freqs[3])
    # vec_inz.append(vecs[0])
    # vec_outz.append(vecs[1])
    # vec_inr.append(vecs[2])
    # vec_outr.append(vecs[3])
    ### Calculate LDs on barium
    print('a_outz/sqrt(w_outz) = ',vecs[1][0]/np.sqrt(freqs[1]))
    print('a_inr/sqrt(w_inr) = ',vecs[2][0]/np.sqrt(freqs[2]))
    LD_inz.append( TI.dicke_param(freqs[0],k,sys.m1,vecs[0][0]))
    LD_outz.append(TI.dicke_param(freqs[1],k,sys.m1,vecs[1][0]))
    LD_inr.append(TI.dicke_param(freqs[2],k,sys.m1,vecs[2][0]))
    LD_outr.append(TI.dicke_param(freqs[3],k,sys.m1,vecs[3][0]))
fig,ax = plt.subplots(figsize = (10,10))
N_s = np.arange(0,12)+1
#ax.scatter(N_s,LD_inz,label = 'In phase axial')
ax.scatter(N_s,LD_outz,label = 'Out of phase axial',color = 'blue')
ax.scatter(N_s,LD_inr,label =' In phase radial',color = 'red')
#ax.scatter(N_s,LD_outr,label = 'Out of phase radial')
ax.legend()
ax.set_xticks(np.arange(0,13))
#ax.set_ylim(0,0.03)
ax.set_ylabel(r'Lamb-Dicke parameter $\eta$')
ax.set_xlabel('Number of porypherin rings')
plt.show()