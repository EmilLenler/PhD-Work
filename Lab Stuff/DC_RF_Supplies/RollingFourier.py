import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp


DataFrame = pd.read_csv(r'C:\Users\au581149\PhD-Work\Lab Stuff\DC_RF_Supplies\ScopeCIT_RF.csv',skiprows = 3)
DataFrame= DataFrame.to_numpy()
print(DataFrame)
V = DataFrame[:,1]
t = DataFrame[:,0]
print('Sample rate is ', 1/np.diff(t)[0])
fig,ax = plt.subplots()
ax.plot(t,V)


fig2,ax2=plt.subplots()
#mostly just following the stfft example from scipy documentation
T_x = (np.diff(t)[0]) # Sample rate
N = len(t)
print(T_x)
window = np.blackman(np.floor(5/(300000)/T_x)) #Blackman window over 5 RF cycles
print(len(window))
SFT = sp.signal.ShortTimeFFT(window,hop = 10,fs = 1/T_x)

Sx = SFT.stft(V)
def f(q):
    return q/(len(window)*T_x)
print(Sx.shape)

freqs = []
for j in range(Sx.shape[0]):
    freqs.append(f(j))
freqs = np.array(freqs)
ax2.plot(freqs,np.abs(Sx[:,5]))
plt.show()


# plt.show()