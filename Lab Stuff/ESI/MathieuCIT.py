import numpy as np
import matplotlib.pyplot as plt

def  b1Curve(q):
    return 1-q-q**2/8 + q**3/64 -q**4/1536 - 11*q**5/36864 + 49*q**8/589824 - 55*q**7/9437184 - 83*q**8/35389440

def a0Curve(q):
    return -q**2/2 +7*q**4/128 -29*q**6/2304 + 68687*q**8/18874368
qs = np.linspace(0,1.5,100000)
fig,ax = plt.subplots()
ax.set_ylim(-0.7,0.2)
ax.set_xlim(0,1.6)
#Determine intersections!
upperq = qs[np.argmin(np.abs(b1Curve(qs)+2*a0Curve(1/2*qs)))]
lowerq = qs[np.argmin(np.abs(a0Curve(qs)+2*b1Curve(1/2*qs)))]
rightq = qs[np.argmin(np.abs(-2*b1Curve(1/2*qs)-b1Curve(qs)))]

print(lowerq)
print(upperq)
print(rightq)
e = 1.6*1e-19
m = 478*1.66*1e-27
omega = 2*np.pi*300*1e3
chi2 = 2.36*1e-4
Z = 1
def RFToq(V):
    return 8*V*e*Z/(m*chi2*omega**2)
def qToRF(q):
    return q*m*chi2*omega**2/(8*e*Z)

def DCtoA(V):
    return -16*e*Z*V/(m*chi2*omega**2)
def AtoDC(a):
    return a*m*chi2*omega**2/(-16*e*Z)





def upperLine(qs):
    res = []
    for q in qs:
        if q<=upperq:
            res.append(-2*a0Curve(1/2*q))
        if q>upperq:
            if q>=rightq:
                break
            else:
                res.append(b1Curve(q))

    return np.array(res)
def lowerLine(qs):
    res = []
    for q in qs:
        if q<=lowerq:
            res.append(a0Curve(q))
        if q>lowerq:
            if q>=rightq:
                break
            else:
                res.append(-2*b1Curve(1/2*q))
    return np.array(res)

ax.plot(qs[qs<rightq],upperLine(qs), c = 'r',label = 'm = 478')
ax.plot(qs[qs<rightq],lowerLine(qs), c = 'r')
ax.plot(480/479*qs[qs<rightq],480/479*upperLine(qs), c = 'g',label = 'm = 479')
ax.plot(480/479*qs[qs<rightq],480/479*lowerLine(qs), c = 'g')
ax.plot(480/478*qs[qs<rightq],480/478*upperLine(qs), c = 'b',label = 'm = 480')
ax.plot(480/478*qs[qs<rightq],480/478*lowerLine(qs), c = 'b')
ax.set_xlabel(r'$q_{z,478}$')
ax.set_ylabel(r'$a_{z,478}$')
secaxX = ax.secondary_xaxis('top',functions = (qToRF,RFToq))
secaxY = ax.secondary_yaxis('right',functions = (AtoDC,DCtoA))
secaxX.set_xlabel('RF Voltage / V')
secaxY.set_ylabel('DC Voltage / V')

y1 = np.max(upperLine(qs))
x1 = qs[np.argmax(upperLine(qs))]

y2 = np.max(480/478*upperLine(qs))
x2 = 480/478*qs[np.argmax(480/478*upperLine(qs))]

slope100 = (y1-y2)/(x1-x2)
offset = y1-slope100*x1

def scanLine(q,adj):
    return slope100*q*(1+adj)+offset
fig2,ax2 = plt.subplots()
ax2.plot(qs[qs<rightq],upperLine(qs), c = 'r',label = 'm = 478')
ax2.plot(qs[qs<rightq],lowerLine(qs), c = 'r')


ax2.plot(480/479*qs[qs<rightq],480/479*upperLine(qs), c = 'g',label = 'm = 479')
ax2.plot(480/479*qs[qs<rightq],480/479*lowerLine(qs), c = 'g')

ax2.plot(480/478*qs[qs<rightq],480/478*upperLine(qs), c = 'b',label = 'm = 480')
ax2.plot(480/478*qs[qs<rightq],480/478*lowerLine(qs), c = 'b')

# ax2.plot(qs,scanLine(qs,-0.001),c = 'k',ls = '--',label = r'99.9% slope')

ax2.set_xlim(0.77,0.788)
ax2.set_ylim(0.145,0.15155)
secax2X = ax2.secondary_xaxis('top',functions = (qToRF,RFToq))
secax2Y = ax2.secondary_yaxis('right',functions = (AtoDC,DCtoA))
secax2X.set_xlabel('RF Voltage / V')
secax2Y.set_ylabel('DC Voltage / V')
ax2.set_xlabel(r'$q_{z,478}$')
ax2.set_ylabel(r'$a_{z,478}$')

DC1 = (AtoDC(y1))
RF1 = (qToRF(x1))

DC2 = AtoDC(y2)
RF2 = qToRF(x2)

print(DC1,RF1)
print(DC2,RF2)

print('Slope is ')
print(0.999*(DC1-DC2)/(RF1-RF2), 'V_DC / V_RF')
print('Offset is ', AtoDC(offset))


fig3,ax3 = plt.subplots()
ax3.plot(qs[qs<rightq],upperLine(qs), c = 'r',label = 'm = 478')
ax3.plot(qs[qs<rightq],lowerLine(qs), c = 'r')


ax3.plot(480/479*qs[qs<rightq],480/479*upperLine(qs), c = 'g',label = 'm = 479')
ax3.plot(480/479*qs[qs<rightq],480/479*lowerLine(qs), c = 'g')

ax3.plot(480/478*qs[qs<rightq],480/478*upperLine(qs), c = 'b',label = 'm = 480')
ax3.plot(480/478*qs[qs<rightq],480/478*lowerLine(qs), c = 'b')
ax3.set_xlim(0,0.2)
ax3.set_ylim(-0.011,0.011)
ax3.set_xlabel(r'$q_{z,478}$')
ax3.set_ylabel(r'$a_{z,478}$')
secax3x= ax3.secondary_xaxis('top',functions = (qToRF,RFToq))
secax3y = ax3.secondary_yaxis('right',functions = (AtoDC,DCtoA))
secax3x.set_xlabel('RF Voltage / V')
secax3y.set_ylabel('DC Voltage / V')
ax.legend()
ax2.legend()
ax3.legend()

fig4,ax4 = plt.subplots()
ax4.plot(qs[qs<upperq],(-2*a0Curve(1/2*qs)/(-2*480/479*a0Curve(1/2*479/480*qs)))[qs<upperq], c = 'k')



plt.show()


