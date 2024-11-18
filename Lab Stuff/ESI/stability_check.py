import numpy as np
import matplotlib.pyplot as plt


def lower_curve(q):
    #Returns the lower-stability line for a given q- parameter
    #Variables:
    # q - array or single q parameter value
    return -1/2*q**2+7/128*q**4-29/2304*q**6+69697/18874368*q**8
def higher_curve(q):
    #Returns the upper stability curve for a given q - parameter
    #Variables:
    # q - array or single q parameter value
    return 1-q-1/8*q**2+1/64*q**3-1/1536*q**4-11/36864*q**5


qs = np.linspace(0,1,1000)
plt.plot(qs,lower_curve(qs))
plt.plot(qs,higher_curve(qs))
plt.show()