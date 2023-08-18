import numpy as np
import TwoIons as TI
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from Utils import gssmax
import numpy.linalg as LA

Trap1 = TI.Trap(2*np.pi*5.2*1e6,0.248,3/4*2.7*1e-3,3/4*3.5*1e-3)
m_ba = 135*1.66*1e-27
charge_ba = 1*1.6*1e-19
a = -0.00074
q = 0.2
print('Secular Frequencies: ',Trap1.Trap_secular_frequencies(-0.00074,0.2)/(2*np.pi)*1e-3)
System_Ba = TI.two_ion_system(m_ba,charge_ba,m_ba,charge_ba,Trap1)
m_por = 9000*1.66*1e-27
charge_por = 24*1.6*1e-19

por_masses = [(i+1)/12*m_por for i in range(12)]
por_charges = [(i+1)/12*charge_por for i in range(12)]
Por_Systems = [TI.two_ion_system(m_ba,charge_ba,por_m,por_c,Trap1) for (por_m,por_c) in zip(por_masses,por_charges)]


frequencies = []
vectors = []
for System in Por_Systems:
    fs,vecs = System.vibrational_freqs_and_modes(a,q)
    frequencies.append(fs[0])
    vectors.append(vecs[0])
frequencies= np.array(frequencies)
print('2 Ion Frequencies', frequencies*1e-3/(2*np.pi))
# print(vectors)

k = 2*np.pi/(1762*1e-9) #wave vector for 1762nm
dickes = np.array([TI.dicke_param(freq,k,m_ba,v[0]) for freq,v in zip(frequencies,vectors)])
print('Lamb Dickes are', dickes)
fig,ax = plt.subplots(figsize = (10,15))
print('Plotting Dicke')
ax.scatter(np.arange(0,12)+1,dickes)
ax.set_xlabel('Number of porypherin rings')
ax.set_ylabel('Lamb-dicke parameter for out-of-phase axial motion on Ba ion')
fig.suptitle('Lamb-Dicke parameters for two-ion system of barium and porypherine')
ax.set_xticks(np.arange(0,12)+1)

fig,ax = plt.subplots(figsize = (10,15))
print('Plotting Pulse durations')
ax.scatter(np.arange(0,12)+1,np.pi/(dickes*1*1e3*2*np.pi))
ax.set_ylabel(r'$\pi$ pulse duration / s')
ax.set_xlabel('Number of porypherin rings')
ax.set_xticks(np.arange(0,12)+1)
fig.suptitle(r'$\pi$ pulse durations for 1kHz carrier Rabi frequency')
fig,ax = plt.subplots(figsize = (10,15))
print('1ms pulse')
ax.scatter(np.arange(0,12)+1,np.pi/(1e-3*dickes*2*np.pi)*1e-3)
ax.set_xlabel('Number of porypherin rings')
ax.set_ylabel('Carrier rabi frequency / kHz')
fig.suptitle(r'Required carrier rabi frequency to achieve 1ms $\pi$ pulse on sideband')
fig,ax = plt.subplots(figsize = (10,15))
print('10ms pulse')
ax.scatter(np.arange(0,12)+1,np.pi/(10*1e-3*dickes*2*np.pi)*1e-3)
fig.suptitle(r'Required carrier rabi frequency to achieve 10ms $\pi$ pulse on sideband')
ax.set_xticks(np.arange(0,12)+1)
ax.set_xlabel('Number of porypherines')
ax.set_ylabel('Rabi frequency / kHz')
state = np.array([0,1,0,0],dtype = np.cdouble)
ExcState = np.array([0,0,0,1],dtype = np.cdouble)
dt = 1e-11
ts = np.arange(0,1e-4,dt)
hbar = 1.05*1e-34
state1s = np.zeros(ts.size,dtype = np.clongdouble)
state2s = np.zeros(ts.size,dtype = np.clongdouble)
state3s = np.zeros(ts.size,dtype = np.clongdouble)
state4s = np.zeros(ts.size,dtype = np.clongdouble)





t_end = 1e-3
def TimeEv(t,phi):
   return 1/(1j*hbar)*np.matmul(TI.Dicke_Hamilton(35*1e3,-frequencies[2],frequencies[2],t,m_ba,vectors[2][0],k),phi)
UnExcitedEvolution = solve_ivp(TimeEv,[0,t_end],state,t_eval = np.linspace(0,t_end,10000),max_step = 1e-8)
ExcitedEvolution = solve_ivp(TimeEv,[0,t_end],np.array([0,0,0,1],dtype = np.clongdouble),t_eval=np.linspace(0,t_end,10000),max_step = 1e-8)
print(UnExcitedEvolution.y)
fig,ax = plt.subplots(2,1,figsize = (10,15))
print('UneXc Evo')
ax[0].plot(UnExcitedEvolution.t/(1e-3),UnExcitedEvolution.y[0]*np.conjugate(UnExcitedEvolution.y[0]),label = r'$\vert 0,e\rangle$')
ax[0].plot(UnExcitedEvolution.t/(1e-3),UnExcitedEvolution.y[1]*np.conjugate(UnExcitedEvolution.y[1]),label = r'$\vert 0,g\rangle$')
ax[0].plot(UnExcitedEvolution.t/(1e-3),UnExcitedEvolution.y[2]*np.conjugate(UnExcitedEvolution.y[2]),label = r'$\vert 1,e\rangle$')
ax[0].plot(UnExcitedEvolution.t/(1e-3),UnExcitedEvolution.y[3]*np.conjugate(UnExcitedEvolution.y[3]),label = r'$\vert 1,g\rangle$')
ax[0].plot(UnExcitedEvolution.t/(1e-3),np.abs(UnExcitedEvolution.y[0])**2+np.abs(UnExcitedEvolution.y[1])**2+np.abs(UnExcitedEvolution.y[2])**2+np.abs(UnExcitedEvolution.y[3])**2,color = 'k')
ax[0].legend()
ax[0].set_title(r'Initialization in ground state $\vert 0,g\rangle$ ')
print(np.max(np.abs(np.imag(UnExcitedEvolution.y[0]*np.conjugate(UnExcitedEvolution.y[0])))))
print(np.max(np.abs(np.imag(UnExcitedEvolution.y[1]*np.conjugate(UnExcitedEvolution.y[1])))))
print(np.max(np.abs(np.imag(UnExcitedEvolution.y[2]*np.conjugate(UnExcitedEvolution.y[2])))))
print(np.max(np.abs(np.imag(UnExcitedEvolution.y[3]*np.conjugate(UnExcitedEvolution.y[3])))))
print('ExcEvo')
ax[1].plot(ExcitedEvolution.t/(1e-3),ExcitedEvolution.y[0]*np.conjugate(ExcitedEvolution.y[0]),label = r'$\vert 0,e\rangle$')
ax[1].plot(ExcitedEvolution.t/(1e-3),ExcitedEvolution.y[1]*np.conjugate(ExcitedEvolution.y[1]),label = r'$\vert 0,g\rangle$')
ax[1].plot(ExcitedEvolution.t/(1e-3),ExcitedEvolution.y[2]*np.conjugate(ExcitedEvolution.y[2]),label = r'$\vert 1,e\rangle$')
ax[1].plot(ExcitedEvolution.t/(1e-3),ExcitedEvolution.y[3]*np.conjugate(ExcitedEvolution.y[3]),label = r'$\vert 1,g\rangle$')
ax[1].plot(ExcitedEvolution.t/(1e-3),np.abs(ExcitedEvolution.y[0])**2+np.abs(ExcitedEvolution.y[1])**2+np.abs(ExcitedEvolution.y[2])**2+np.abs(ExcitedEvolution.y[3])**2,color = 'k')

ax[1].legend()
ax[1].set_xlabel('Time / ms')
ax[0].set_ylabel('State Population')
ax[1].set_ylabel('State Population')
ax[1].set_title(r'Initialization in 1st blue sideband $\vert 1,g\rangle$')


fig,ax = plt.subplots(2,1,figsize = (10,15))
ax[0].plot(UnExcitedEvolution.t/1e-3,np.abs(UnExcitedEvolution.y[0])**2+np.abs(UnExcitedEvolution.y[2])**2)
ax[1].plot(UnExcitedEvolution.t/1e-3,np.abs(ExcitedEvolution.y[0])**2+np.abs(ExcitedEvolution.y[2])**2)
ax[0].set_title(r'Initialization in ground state $\vert 0,g\rangle$ ')
ax[1].set_title(r'Initialization in 1st blue sideband $\vert 1,g\rangle$')

for axe in ax:
    axe.set_ylim(0,1.05)
    axe.set_ylabel('Excited state population')
for axe in ax:
    axe.axhline(1,ls = '--',color = 'k')
ax[1].set_xlabel('Time / ms')
# eigvals,eigvecs = LA.eig(TI.Dicke_Hamilton(2*np.pi*6.61*1e3,-frequencies[2]-2*np.pi*0*1e3,frequencies[2],0,m_ba,vectors[2][0],k))
# print('Eigenvalues are:',eigvals)
# print('Eigenvectors are',eigvecs)
# print(eigvecs[:,0])
# print(eigvecs[:,1])
# print(eigvecs[:,2])
# print(eigvecs[:,3])
# def flopmax(delta):
#     def FTimeEv(t,phi):
#         return 1/(1j*hbar)*np.matmul(TI.Dicke_Hamilton(2*np.pi*3.5*1e3,-frequencies[2],frequencies[2]-2*np.pi*1e3*delta,t,m_ba,vectors[2][0],k),phi)
#     flopsolve = solve_ivp(FTimeEv,[0,1e-3],np.array([0,0,0,1],dtype = np.clongdouble),t_eval=np.linspace(0,1e-3,10000),max_step = 1e-8)
# #     return np.max(np.abs(flopsolve.y[0])**2)
# # goodflop = gssmax(flopmax,1,100,n_runs = 10)
# print(' Bonus detuning is -2pi',goodflop,'kHz')
# # t_end_L = 1e-1
# # def BigTimeEv(t,phi):
# #     return 1/(1j*hbar)*np.matmul(TI.Dicke_Hamilton(2*np.pi*6.61*1e3,-frequencies[-1],frequencies[-1],t,m_ba,vectors[-1][0],k),phi)
# BigUnExcited = solve_ivp(BigTimeEv,[0,t_end_L],state,t_eval=np.linspace(0,t_end_L,10000),max_step = 1e-8)
# BigExcited = solve_ivp(BigTimeEv,[0,t_end_L],ExcState,t_eval=np.linspace(0,t_end_L,10000),max_step = 1e-8)
# fig,ax = plt.subplots(2,1,figsize = (10,15))
# fig.suptitle('12 Porypherin System')
# labels = [r'$\vert 0,e\rangle$',r'$\vert 0,g\rangle$',r'$\vert 1,e\rangle$',r'$\vert 1,g\rangle$']
# for j,y_s in enumerate(BigUnExcited.y):
#   ax[0].plot(BigUnExcited.t,np.abs(y_s)**2,label = labels[j])
# for j,y_s in enumerate(BigExcited.y):
#     ax[1].plot(BigExcited.t,np.abs(y_s)**2,label = labels[j])
# ax[0].plot(BigUnExcited.t,np.abs(BigUnExcited.y[0])**2+np.abs(BigUnExcited.y[1])**2+np.abs(BigUnExcited.y[2])**2+np.abs(BigUnExcited.y[3])**2,ls = '--',color = 'k')
# ax[1].plot(BigExcited.t,np.abs(BigExcited.y[0])**2+np.abs(BigExcited.y[1])**2+np.abs(BigExcited.y[2])**2+np.abs(BigExcited.y[3])**2,ls = '--',color = 'k')


# for axe in ax:
#     axe.legend()
# def BigLongExcitation(t,phi):
#      return 1/(1j*hbar)*np.matmul(TI.Dicke_Hamilton(6.5*1e3,-frequencies[-1],frequencies[-1],t,m_ba,vectors[-1][0],k),phi)
# LongUnexcited = solve_ivp(BigLongExcitation,[0,0.01],state,t_eval = np.arange(0,0.01,1e-7))
# LongExcited = solve_ivp(BigLongExcitation,[0,0.01],ExcState,t_eval = np.arange(0,0.01,1e-7))



# fig,ax = plt.subplots(2,1,figsize = (10,15))
# for j,y_s in enumerate(LongUnexcited.y):
#     ax[0].plot(LongUnexcited.t,np.abs(y_s)**2,label = labels[j])
# for j,y_s in enumerate(LongExcited.y):
#     ax[1].plot(LongExcited.t,np.abs(y_s)**2,label = labels[j])
# for axe in ax:
#     axe.legend()

# print(1/(1j*hbar)*np.matmul(TI.Dicke_Hamilton(2*np.pi*6.5*1e3,-frequencies[-1],frequencies[-1],0,m_ba,vectors[-1][0],k),ExcState)*1e-8)
# print(1/(1j*hbar)*TI.Dicke_Hamilton(2*np.pi*6.5*1e3,-frequencies[-1],frequencies[-1],0,m_ba,vectors[-1][0],k))
# fig,ax = plt.subplots(figsize = (10,15))
#ax.plot(LongExcited.t,np.exp(-1j*frequencies[-1]*LongExcited.t)*np.exp(1j*frequencies[-1]*LongExcited.t)-1)

# def checker(t,phi):
#     return np.array(,dtype = np.com)

ts = np.linspace(0,1e-3,10000)
GS_pop = np.ones(20)*1/(20)
Shelv_Probs = [Por_Systems[2].shelving_prob(GS_pop,32*1e3,t,frequencies[2],k,m_ba,vecs[2][0],) for t in ts]
fig,ax = plt.subplots(figsize = (10,15))
ax.plot(ts,Shelv_Probs)
plt.show()
