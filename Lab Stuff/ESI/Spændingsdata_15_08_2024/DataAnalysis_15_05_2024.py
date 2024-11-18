import matplotlib.pyplot as plt
from PrettyPlot import makePrettyPlot
import numpy as np

HCData = np.loadtxt(r'C:\Users\au581149\PhD-Work\Lab Stuff\ESI\Spændingsdata_15_08_2024\Heated Capillary.txt', skiprows = 2)
Oct1Data = np.loadtxt(r'C:\Users\au581149\PhD-Work\Lab Stuff\ESI\Spændingsdata_15_08_2024\Octopole 1.txt', skiprows = 2)
Oct2Data = np.loadtxt(r'C:\Users\au581149\PhD-Work\Lab Stuff\ESI\Spændingsdata_15_08_2024\Octopole 2.txt', skiprows = 2)
SkimData = np.loadtxt(r'Lab Stuff\ESI\Spændingsdata_15_08_2024\Skimmer.txt', skiprows = 2)
TubeData = np.loadtxt(r'Lab Stuff\ESI\Spændingsdata_15_08_2024\Tube Lens.txt', skiprows  =2)

figHC,axHC = makePrettyPlot(HCData,sqrtN = True)
axHC.set_title('2nd Channeltron Signal dependance on capillary voltage')
axHC.set_ylim(0,75000)

figOct1,axOct1 = makePrettyPlot(Oct1Data,sqrtN = True)
axOct1.set_title('2nd Channeltron Signal dependance on 1st octopole voltage')
axOct1.set_yscale('log')
axOct1.axvline(31, ls = '--', label = 'PT Entrance Voltage (31V)')
axOct1.legend()

figOct2,axOct2 = makePrettyPlot(Oct2Data,sqrtN = True)
axOct2.set_title('2nd Channeltron Signal dependance on 2nd octopole voltage')
axOct2.axvline(31.5,ls = '--', label = '1st Octopole DC (31.5V). Assumed energy of ions.')
axOct2.legend()

figSkim,axSkim = makePrettyPlot(SkimData,sqrtN = True)
axSkim.set_title('2nd Channeltron Signal dependance on skimmer voltage')
axSkim.axvline(120,ls = '--', label = 'Tube lens voltage (120V)')
axSkim.axvline(31.5, ls = '--', c = 'k', label = '1st Octopole DC (31.5V)')
axSkim.legend()

figTube,axTube = makePrettyPlot(TubeData,sqrtN = True)
axTube.set_title('2nd Channeltron Signal dependance on tube lens voltage')


# fig,(ax1,ax2) = plt.subplots(2,1,sharex = True)
# ax1.scatter(Oct1Data[:,0],Oct1Data[:,1])
# ax2.scatter(Oct1Data[:,0],Oct1Data[:,1])
# ax1.set_ylim(1000000,4000000)
# ax2.set_ylim(0,7*1e4)
# ax1.spines.bottom.set_visible(False)
# ax2.spines.top.set_visible(False)
# ax1.xaxis.tick_top()
# ax1.tick_params(labeltop=False)  # don't put tick labels at the top
# ax2.xaxis.tick_bottom()
# d = .5  # proportion of vertical to horizontal extent of the slanted line
# kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
#               linestyle="none", color='k', mec='k', mew=1, clip_on=False)
# ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
# ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)
# ax1.ticklabel_format(style = 'plain')
# fig.suptitle('2nd Channeltron Signal dependance on 1st octopole voltage')
# ax2.set_ylabel('Accumulated ion counts over 10 seconds', position = (0,0.5))



for axe in [axHC,axOct1,axOct2,axSkim,axTube]:
    axe.set_xlabel('Voltage / V')
    axe.set_ylabel('Accumulated ion counts over 10 seconds')



plt.show()