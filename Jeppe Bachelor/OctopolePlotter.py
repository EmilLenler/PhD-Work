import numpy as np
import matplotlib.pyplot as plt


Delay_250 = np.array([0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2])
Count_250 = np.array([0,5,51,92,93,69,72,49,41,54,41,35,13])

Delay_450 = np.array([0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.25,1.5,1.75,2])
Count_450 = np.array([0,7,80,160,119,124,87,78,74,75,55,50,27,28,23,15])
Delay_1k = np.array([0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.25,1.5,1.75,2])
Count_1k = np.array([0,10,117,215,177,185,166,135,108,130/(1.15),127/(1.15),101/(1.15),65/(1.15),58/1.15,42/1.15,41/1.15]) #Note these have been normalized to DC current as current went up to 1150 at the 0.8ms mark-


N_pulses = 50

#plt.plot(Delay_250,Count_250)
fig_250,ax_250 = plt.subplots()
ax_250.errorbar(Delay_250,Count_250/N_pulses,np.sqrt(Count_250)/N_pulses,ls = '--', marker = "o")
ax_250.set_xlabel('Delay between detection and octopole opening (ms)')
ax_250.set_ylabel('Counts per pulse')
ax_250.set_title('DC Current: 250 c/s')

fig_450,ax_450 =plt.subplots()
ax_450.set_title('DC Current: 450 c/s')


ax_450.errorbar(Delay_450,Count_450/50,np.sqrt(Count_450)/50,ls = '--', marker = 'o')
ax_450.set_xlabel('Delay between detection and octopole opening (ms)')
ax_450.set_ylabel('Counts per pulse')

fig_1k,ax_1k = plt.subplots()
ax_1k.set_title('DC Current: 1000 c/s')
ax_1k.errorbar(Delay_1k,Count_1k/50,np.sqrt(Count_450)/50,ls = '--', marker = 'o')
ax_1k.set_xlabel('Delay between detection and octopole opening (ms)')
ax_1k.set_ylabel('Counts per pulse')

fig_comp,ax_comp = plt.subplots()
Norm_250 = Count_250/np.max(Count_250)
Norm_450 = Count_450/np.max(Count_450)
Norm_1k = Count_1k/np.max(Count_1k)
ax_comp.errorbar(Delay_250,Norm_250,np.sqrt(Count_250)/np.max(Count_250),ls = '--',marker = 'o', label = '250 c/s DC')
ax_comp.errorbar(Delay_450,Norm_450,np.sqrt(Count_450)/np.max(Count_450),ls = '--',marker = 'o', label = '450 c/s DC')
ax_comp.errorbar(Delay_1k,Norm_1k,np.sqrt(Count_1k)/np.max(Count_1k),ls = '--',marker = 'o', label = '1000 c/s DC')
ax_comp.legend()
ax_comp.set_xlabel('Delay between detection and octopole opening (ms)')
ax_comp.set_ylabel('Detector signal normalized to peak signal')
ax_comp.set_title('Comparison of ion pulse shape')
plt.show()
