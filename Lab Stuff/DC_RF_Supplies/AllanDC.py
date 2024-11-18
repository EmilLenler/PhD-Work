import matplotlib.pyplot as plt
import numpy as np
import allan_variance as av


data = np.loadtxt(r'C:\Users\au581149\PhD-Work\Lab Stuff\DC_RF_Supplies\scope_log_20240906_cit_dc.txt',skiprows = 1,delimiter = ',')
t0 = data[0,0]
data[:,0]-=t0
print(data[:,0])
ts = data[:,0]
ts =ts[data[:,0]>10]
Vs = data[:,1][data[:,0]>10]
fig,ax = plt.subplots()
ax.set_xlabel('Time / s')
ax.set_ylabel('DC Voltage (Monitor) / V')
ax.plot(ts,Vs)
plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)


tau,av = av.compute_avar(Vs,(np.diff(ts))[0])
figav,axav = plt.subplots()
axav.plot(tau,av)
axav.set_xscale('log')
axav.set_yscale('log')
axav.set_xlabel(r'$\tau$ / s')
axav.set_ylabel(r'$\sigma_A$')

maxdiff = np.max(Vs)-np.min(Vs)
std = np.std(Vs)
print('Max difference between two datapoints is ', maxdiff, 'V')
print('With respect to the mean this is a ', maxdiff/np.mean(Vs)*100, '% deviation')
print('STD is ',std, 'V')
print('With respect to the mean this is a ', std/np.mean(Vs)*100, '% deviation')



plt.show()
