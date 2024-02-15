import numpy as np
import TwoIons as TI
import matplotlib.pyplot as plt
plt.rc('font',size = 14)
r0s = np.linspace(1,3.5,100)*1e-3
m_ba = 135*1.66*1e-27
charge_ba = 1*1.6*1e-19
q_set = 0.4
a_set = -0.01
V_RF_our = (2*np.pi*5.2*1e6)**2*q_set*m_ba*(3/4*3.5*1e-3)**2/(2*charge_ba)
V_DC_our = -a_set*m_ba*(3/4*2.7*1e-3)**2*(2*np.pi*5.2*1e6)**2/(4*charge_ba*0.248)
print('Our RF Voltage is', V_RF_our,'Volts')
print('Our DC Voltage is ', V_DC_our,'Volts')

omega_rfs = np.sqrt((2*charge_ba*V_RF_our)/(m_ba*q_set*r0s**2))

z0s = np.sqrt((-4*charge_ba*0.248*V_DC_our)/(m_ba*a_set*omega_rfs**2))
a_check = (-4*charge_ba*0.248*V_DC_our)/(m_ba*z0s**2*omega_rfs**2)
print('Checking a',a_check)
q_check = (2*charge_ba*V_RF_our)/(m_ba*r0s**2*omega_rfs**2)
print('q check',q_check)

#We've checked that our system parameters are functioning. Now to make a bunch of systems.
Traps = [TI.Trap(w_rf,0.248,z0,r0) for w_rf,z0,r0 in zip(omega_rfs,z0s,r0s)]
m_por = 9000*1.66*1e-27
charge_por = 24*1.6*1e-19

Systems = [TI.two_ion_system(m_ba,charge_ba,m_por,charge_por,trap) for trap in Traps]
w_s= [sys.vibrational_freqs_and_modes(a_set,q_set)[0] for sys in Systems]
w_ins = np.array([w[0] for w in w_s])
v_s = [sys.vibrational_freqs_and_modes(a_set,q_set)[1] for sys in Systems]
v_ins = np.array([v[0] for v in v_s])
def r0_to_rf_MHz_divide_2pi(r0):
    rf = np.sqrt(2*charge_ba*V_RF_our/(m_ba*q_set*(r0*1e-3)**2))*1e-6/(2*np.pi)
    return np.round(rf,2)
alphas = np.array([1,2,3,4,5])
fig,ax = plt.subplots(figsize = (10,15))
for alpha in alphas:
    ax.plot(r0s*1e3,w_ins/(2*np.pi*1e3*alpha),label = r'$\alpha = ${0}'.format(alpha))
ax.legend()
ax.set_xlabel(r'$r_0$ / mm')
ax.set_ylabel(r'$\Omega_c$ / kHz')
ax2 = ax.twiny()
ax1ticks = ax.get_xticks()
ax2ticks = ax1ticks
ax2.set_xticks(ax2ticks)
ax2.set_xticklabels(r0_to_rf_MHz_divide_2pi(ax2ticks))
ax2.set_xbound(ax.get_xbound())
ax2.set_xlabel('RF Frequency / MHz')
fig.suptitle(r'Carrier Rabi frequency for a given $ \alpha = \frac{\delta}{\Omega_c}$, assuming 12 porypherin rings')



fig,ax = plt.subplots(figsize = (10,15))
for alpha in alphas:
    dickes = np.array([TI.dicke_param(w,2*np.pi/(1762*1e-9),m_ba,v[0]) for w,v in zip(w_ins,v_ins)])
    print(np.shape(dickes))
    ax.plot(r0s*1e3,w_ins/(2*np.pi*1e3*alpha)*dickes,label = r'$\alpha = ${0}'.format(alpha))
ax.set_ylim(0,1.8)
ax.legend()
ax.set_xlabel(r'$r_0$ / mm')
ax.set_ylabel(r'$\eta\Omega_c$ / kHz')
ax2 = ax.twiny()
ax1ticks = ax.get_xticks()
ax2ticks = ax1ticks
ax2.set_xticks(ax2ticks)
ax2.set_xticklabels(r0_to_rf_MHz_divide_2pi(ax2ticks))
ax2.set_xbound(ax.get_xbound())
ax2.set_xlabel('RF Frequency / MHz')


plt.show()[