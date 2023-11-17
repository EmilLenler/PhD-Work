import numpy as np
import TwoIons as TI
import matplotlib.pyplot as plt

plt.rc('font', size = 16)
Trap1 =TI.Trap(2*np.pi*5.2*1e6,0.248,3/4*2.7*1e-3,3/4*3.5*1e-3)

m_ba = 135*1.66e-27
m_por = 9000*1.66e-27
charge_ba = 1*1.6e-19
charge_por = 24*1.6e-19
mu = m_por/m_ba
rho = charge_por/charge_ba
System1 = TI.two_ion_system(m_ba,charge_ba,m_por,charge_por,Trap1)


fig,ax = plt.subplots()
#ax.plot(System1.qrange,System1.lower_curve(System1.qrange), color = 'k')
#ax.plot(System1.qrange,System1.higher_curve(System1.qrange),color = 'k')
#ax.fill_between(System1.qrange,System1.lower_curve(System1.qrange),System1.higher_curve(System1.qrange),color = 'tab:blue')
#ax.plot(rho/mu*System1.qrange,System1.lower_curve(System1.qrange),color = 'red')
#ax.plot(rho/mu*System1.qrange,System1.higher_curve(System1.qrange),color = 'red')
#ax.set_ylim(-0.20,0)

BothAxis,Both,Ion1,Ion2 = System1.Master_Stability()
ax.scatter(Both[1],Both[0],color = 'green',s = 2)
ax.scatter(Ion1[1],Ion1[0],color = 'tab:blue', s = 2)
ax.scatter(Ion2[1],Ion2[0],color = 'tab:red', s = 2)   
ax.scatter(BothAxis[1],BothAxis[0],color = 'grey', s = 2)
#print(Ion2[0])
#print(Ion2[1])
plt.show()
