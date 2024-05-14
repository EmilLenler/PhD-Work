import numpy as np
import TwoIons as TI
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches 
from matplotlib import ticker as mtick


plt.rc('font', size = 22)
Trap1 =TI.Trap(2*np.pi*5.2*1e6,0.248,3/4*2.7*1e-3,3/4*3.5*1e-3)

m_ba = 135*1.66e-27
m_por = 9000*1.66e-27
charge_ba = 1*1.6e-19
charge_por = 24*1.6e-19
mu = m_por/m_ba
rho = charge_por/charge_ba
System1 = TI.two_ion_system(m_ba,charge_ba,m_por,charge_por,Trap1)

def Ba_to_PolyPor(x):
    return rho/mu*x
def Polypor_to_Ba(x):
    return mu/rho*x
fig,ax = plt.subplots(figsize = (15,10))
#ax.plot(System1.qrange,System1.lower_curve(System1.qrange), color = 'k')
#ax.plot(System1.qrange,System1.higher_curve(System1.qrange),color = 'k')
#ax.fill_between(System1.qrange,System1.lower_curve(System1.qrange),System1.higher_curve(System1.qrange),color = 'tab:blue')
#ax.plot(rho/mu*System1.qrange,System1.lower_curve(System1.qrange),color = 'red')
#ax.plot(rho/mu*System1.qrange,System1.higher_curve(System1.qrange),color = 'red')
#ax.set_ylim(-0.20,0)

BothAxis,Both,Ion1,Ion2 = System1.Master_Stability()
ax.scatter(Both[1],Both[0],s = 1,color = 'tab:blue',label = 'Ba + Polyporypher')
ax.scatter(Ion1[1],Ion1[0], s = 1,color = 'tab:red',label = 'Single ion stability for ion 1 only')
#ax.scatter(Ion2[1],Ion2[0],color = 'tab:red', s = 2)   
ax.scatter(BothAxis[1],BothAxis[0], s = 1,color = 'tab:green',label = 'Single ion stability for both ions, as well as alignment along axis')
#ax.legend()

blue_patch = mpatches.Patch(color='tab:blue', label=r'Both ions stable, but not aligned along $z$-axis') 
red_patch = mpatches.Patch(color='tab:red', label=r'Only Ba is stable')
green_patch = mpatches.Patch(color='tab:green', label=r'Both ions stable, aligned along $z$-axis') 
plt.legend(handles = [green_patch,blue_patch,red_patch])
ax.set_ylabel(r'$a$ parameter for Barium')
ax.set_xlabel(r'$q$ parameter for Barium')
secax_x = ax.secondary_xaxis('top',functions=(Ba_to_PolyPor,Polypor_to_Ba))
secax_y = ax.secondary_yaxis('right',functions=(Ba_to_PolyPor,Polypor_to_Ba))
secax_x.set_xlabel(r'$q$ parameter for polyporphyrin')
secax_y.set_ylabel(r'$a$ parameter for polyporphyrin')
x_ticks = ax.get_xticks()
y_ticks = ax.get_yticks()
print(x_ticks)
secax_x.set_xticks(Ba_to_PolyPor(x_ticks))
secax_y.set_yticks(Ba_to_PolyPor(y_ticks))
#print(Ion2[0])
#print(Ion2[1])
#Perhaps make nice curves instead.
# # BothACurve = []
# # for q in Both[1]:
# #     a = np.min(np.array(Both[0])[np.array(Both[1]) == q])
# #     BothACurve.append(a)
# # fig, ax = plt.subplots()
# # ax.plot(np.array(Both[1]),np.array(BothACurve)[np.array])
plt.show()
