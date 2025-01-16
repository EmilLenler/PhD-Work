import numpy as np
import matplotlib.pyplot as plt




class Experiment:
    def __init__(self,GateWidth1,GateWidth2):
        
        self.ArmTime = 4
        
        self.OctDelay1 = np.random.rand()
        self.OctWidth1 = np.random.rand()
        self.OctDelay2 = np.random.rand()
        self.OctWidth2 = np.random.rand()
        self.OctHi = np.random.rand()
        self.OctLo = np.random.rand()


        self.EntrDelay1 = np.random.rand()
        self.EntrWidth1 = np.random.rand()
        self.EntrDelay2 = np.random.rand()
        self.EntrWidth2 = np.random.rand()
        self.EntrHi = np.random.rand()
        self.EntrLo = np.random.rand()

        self.CylinderDCDelay1 = np.random.rand()
        self.CylinderDCWidth1 = np.random.rand()
        self.CylinderDCDelay2 = np.random.rand()
        self.CylinderDCWidth2 = np.random.rand()
        self.CylinderDCHi = np.random.rand()
        self.CylinderDCLo = np.random.rand()


        self.CylinderRFDelay1 = np.random.rand()
        self.CylinderRFWidth1 = np.random.rand()
        self.CylinderRFDelay2 = np.random.rand()
        self.CylinderRFWidth2 = np.random.rand()
        self.CylinderRFHi = np.random.rand()
        self.CylinderRFLo = 0

        self.ExitDelay1 = np.random.rand()
        self.ExitWidth1 = np.random.rand()
        self.ExitDelay2 = np.random.rand()
        self.ExitWidth2 = np.random.rand()
        self.ExitHi = np.random.rand()
        self.ExitLo = np.random.rand()

        self.GateDelay1 = np.random.rand()
        self.GateWidth1 = GateWidth1
        self.GateDelay2 = np.random.rand()
        self.GateWidth2 = GateWidth2
        self.GateHi = 1
        self.GateLo = 0
        


    def PlotTimings(self):
        Prefixes = ['Oct','Entr','CylinderDC','CylinderRF','Exit','Gate']
        fig,ax = plt.subplots(6,1,figsize = (10,10))
        Time = np.linspace(0,self.ArmTime,10000)

        for j,prefix in enumerate(Prefixes):
            componentDict = {key: val for key,val in self.__dict__.items() if key.startswith(prefix)}
            D1 = [val for key,val in componentDict.items() if key.endswith('Delay1')][0]
            W1 = [val for key,val in componentDict.items() if key.endswith('Width1')][0]
            D2 = [val for key,val in componentDict.items() if key.endswith('Delay2')][0]
            W2 = [val for key,val in componentDict.items() if key.endswith('Width2')][0]
            
            HiVal = [val for key,val in componentDict.items() if key.endswith('Hi')][0]
            LoVal = [val for key,val in componentDict.items() if key.endswith('Lo')][0]
            LineX = [0,D1,D1,(D1+W1),(D1+W1),(D1+W1+D2),(D1+W1+D2),(D1+W1+D2+W2),(D1+W1+D2+W2),Time[-1]]#,W1,D2,W2,Time[-1]]
            LineY = [LoVal,LoVal,HiVal,HiVal,LoVal,LoVal,HiVal,HiVal,LoVal,LoVal]
            ax[j].plot(LineX,LineY)
            ax[j].set_xlim(0,Time[-1])
            ax[j].set_ylabel(prefix)
            # ax[j].set_ylim(0,1.5*HiVal)

        ax[-1].set_xlabel('Time / UNIT')
        plt.show() 
        return fig,ax

exp1 = Experiment(np.random.rand(), np.random.rand())
exp1.PlotTimings()
