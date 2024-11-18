import numpy as np
import allan_variance as av
import matplotlib.pyplot as plt
import pandas as pd

data = np.loadtxt(r'C:\Users\au581149\PhD-Work\Lab Stuff\DC_RF_Supplies\scope_log_20240904_RF_cit.txt',skiprows=1)
figData,axData = plt.subplots()
axData.set_ylabel('Monitor peak-to-peak voltage / V')
times = data[:,0]-data[0,0]
# axData.plot(times[times<=1000],2.88*data[:,1][times<=1000],label = r'Monitor (heating up)')
axData.set_xlabel('Time / s')
plt.plot(times,data[:,2]*1.44,label = 'c3')
mask = times>1000
# axData.plot(times[mask],2.88*data[:,1][mask],label = r'Monitor (after heating)')
vPPs = 2.88*data[:,1]
maximumDeviation = np.max(vPPs[mask])-np.min(vPPs[mask])
print('Max difference in Vpp after heating is: ',maximumDeviation, ' V')
print('With respect to the mean of the Vpp data after heating, this is a deviation of ', maximumDeviation/np.mean(vPPs[mask])*100, ' %')
print('STD in Vpp after heating is: ', np.std(vPPs[mask]))
print('With respect to the mean of the Vpp data after heating, this is a deviation of ', np.std(vPPs[mask])/np.mean(vPPs[mask])*100, ' %')
axData.legend()






tau1,av1 = av.compute_avar(data[:,1][mask],np.diff(data[:,0])[0])
tau2,av2 = av.compute_avar(data[:,2][mask],np.diff(data[:,0])[0])
fig1,ax = plt.subplots()
ax.plot(tau1,np.sqrt(av1)/(np.sqrt(av1[0])),label = 'RF Monitor')
ax.plot(tau2,np.sqrt(av2)/(np.sqrt(av2[0])),label = 'Function Generator')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$\tau$')
ax.set_ylabel(r'$\sigma_{Allan}$')

ax.legend()
fig1.suptitle('Stability of RF supply')

plt.show()