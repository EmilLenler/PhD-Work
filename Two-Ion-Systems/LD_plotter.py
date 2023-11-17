import numpy as np
import TwoIons as TI
import matplotlib.pyplot as plt
plt.rc('font',size = 14)
r0s = np.linspace(1,3.5,100)*1e-3
m_ba = 135*1.66*1e-27
charge_ba = 1*1.6*1e-19
q_set = 0.4
a_set = -0.01
m_por = 9000*1.66*1e-27
charge_por = 24*1.6*1e-19

width = 10
height = 15

Trap1 = TI.Trap(2*np.pi*5.2*1e6,0.248,3/4*2.7*1e-3,3/4*3.5*1e-3)
por_masses = [(i+1)/12*m_por for i in range(12)]
por_charges = [(i+1)/12*charge_por for i in range(12)]
Por_Systems = [TI.two_ion_system(m_ba,charge_ba,por_m,por_c,Trap1) for (por_m,por_c) in zip(por_masses,por_charges)]
fig,ax = plt.subplots(figsize=(width,height))

Number_as = 10
LD_s = np.zeros((Number_as,np.size(Por_Systems)))
for j,a in enumerate(np.linspace(a_set*1e-2*1/2,a_set*1e-1/2,Number_as)):
        ba_freqs = np.array([TISys.secular_frequencies(a,q_set)[0] for TISys in Por_Systems])/(2*np.pi*1e3)
        freqs = np.array([TISys.vibrational_freqs_and_modes(a,q_set)[0][0] for TISys in Por_Systems])
        amplitudes = np.array([TISys.vibrational_freqs_and_modes(a,q_set)[1][0][1] for TISys in Por_Systems])
        LD_s[j] = np.array([TI.dicke_param(freq,2*np.pi/(550*1e-9),sys.m2,amp)for freq,sys,amp in zip(freqs,Por_Systems,amplitudes)])
        ax.scatter(np.arange(0,12)+1,LD_s[j]**2,label = r'barium frequency = {0} kHz'.format(np.round(ba_freqs[j],1)))
ax.legend()
print(LD_s)
ax.set_ylim(0,0.06)
ax.set_xticks(np.arange(0,13))
ax.set_xlabel('Number of porypherin rings')
ax.set_ylabel(r'$\eta^2$')
plt.show()

