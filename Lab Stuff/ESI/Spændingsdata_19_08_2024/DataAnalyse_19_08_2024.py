from PrettyPlot import makePrettyPlot
import numpy as np
import matplotlib.pyplot as plt


HcapData = np.loadtxt('HCAP.txt',skiprows=1)
TLData = np.loadtxt('TubeLens.txt',skiprows=1)
SkimmerData = np.loadtxt('Skimmer.txt',skiprows=1)


HCFig,HCax = plt.subplots(2,1)
HCFig.suptitle('Detector signal vs. Heated Capillary voltage')
HCax[0].scatter(HcapData[:,0],HcapData[:,2], label = 'Chann. 1')
HCax[0].scatter(HcapData[:,0],HcapData[:,1], label = 'Chann. 2')


HCax[1].scatter(HcapData[:,0],HcapData[:,2]/np.max(HcapData[:,2]), label = 'Chann. 1')
HCax[1].scatter(HcapData[:,0],HcapData[:,1]/np.max(HcapData[:,1]), label = 'Chann. 2')
HCax[1].set_xlabel('Voltage / V')
HCax[0].set_ylabel('Counts acc. over 10 seconds')
HCax[1].set_ylabel('Normalized signal')
HCax[0].legend()
HCax[1].legend()
HCax[0].set_ylim(0,250000)
HCax[1].set_ylim(0,1.1)





SkimFig, SkimAx = plt.subplots(2,1)
SkimFig.suptitle('Detector signal vs. skimmre voltage')
SkimAx[0].scatter(SkimmerData[:,0],SkimmerData[:,1], label = 'Chann. 1')
SkimAx[0].scatter(SkimmerData[:,0],SkimmerData[:,2], label = 'Chann. 2')


SkimAx[1].scatter(SkimmerData[:,0],SkimmerData[:,1]/np.max(SkimmerData[:,1]), label = 'Chann. 1')
SkimAx[1].scatter(SkimmerData[:,0],SkimmerData[:,2]/np.max(SkimmerData[:,2]), label = 'Chann. 2')
SkimAx[1].set_xlabel('Voltage / V')
SkimAx[0].set_ylabel('Counts acc. over 10 seconds')
SkimAx[1].set_ylabel('Normalized signal')
SkimAx[0].legend()
SkimAx[1].legend()
SkimAx[0].set_ylim(0,250000)
SkimAx[1].set_ylim(0,1.1)



TLFig, TLAx = plt.subplots(2,1)
TLFig.suptitle('Detector signal vs. tube lens voltage')
TLAx[0].scatter(TLData[:,0],TLData[:,1], label = 'Chann. 1')
TLAx[0].scatter(TLData[:,0],TLData[:,2], label = 'Chann. 2')


TLAx[1].scatter(TLData[:,0],TLData[:,1]/np.max(TLData[:,1]), label = 'Chann. 1')
TLAx[1].scatter(TLData[:,0],TLData[:,2]/np.max(TLData[:,2]), label = 'Chann. 2')
TLAx[1].set_xlabel('Voltage / V')
TLAx[0].set_ylabel('Counts acc. over 10 seconds')
TLAx[1].set_ylabel('Normalized signal')
TLAx[0].legend()
TLAx[1].legend()
TLAx[0].set_ylim(0,250000)
TLAx[1].set_ylim(0,1.1)

plt.show()