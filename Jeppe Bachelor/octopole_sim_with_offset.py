import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

L = 0.133 #meter, Octopole length
d = 0 #meter, Output depth

N = 100 #Sections
dx = L/N #meter, Length intervals
output_x = int(np.floor(d/dx)) #How many intervals the output depth fills
dv = 1
dt = dx / dv

k = 1.38 * 10**(-23) #Boltzmanns Constant
T = 300 #Kelvin, Temperature
m = 8*10**(-25) #kg, Mass

#Function for offsetting the velocity
def offset(x):
    a = 0
    b = 0
    return a * x + b

#Boltzmann probability density function
def boltzmann(v, x):
    sigma2 = (k*T)/m
    return np.sqrt(2/(sigma2*np.pi)) * np.exp(-1/(2*sigma2)*(v - offset(x))**2)

#Probability for particle going from section j to i
#"Short" is when particle goes directly to i,
#"Long" is when particle reflects at x = 0, and then reaches i
def probability(i, j):
    #Boltzmann at the x-value for section j
    func = lambda v: boltzmann(v, (j + 1/2) * dx)
    
    #If in output section
    if j == N - output_x:
        if i == N - output_x:
            return 1
        return 0
    #Probabilty of going to the output section
    elif i == N - output_x:
        short = 1/2 * quad(func, (i - j - 1/2) * dv, np.inf)[0]
        long = 1/2 * -quad(func, -(i + j + 1/2) * dv, -np.inf)[0]
    #Probabilty fo going to section i
    else:
        short = 1/2 * quad(func, (i - j - 1/2) * dv, (i - j + 1/2) * dv)[0]
        long = 1/2 * -quad(func, -(i + j + 1/2) * dv, -(i + j + 3/2) * dv)[0]
    return short + long

#Create markov chain
matrix = []
for j in range(N + 1 - output_x):
    row = [probability(i, j) for i in range(N + 1 - output_x)]
    matrix.append(row)
    
#Initial density vector
density_vector = [1/N for _ in range(N - output_x)] + [output_x/N]
print(density_vector)
output = [density_vector[-1]]

while density_vector[-1] < 0.999:
    density_vector = np.dot(density_vector, matrix)
    output.append(density_vector[-1])
    #print("{:.3f}".format(output[-1]))

output_change = [output[i + 1] - output[i] for i in range(len(output) - 1)]

#Plots
fig, ax = plt.subplots()
ax.grid()
ax.plot([1000 * dt * i for i in range(len(output))], output)
ax.set_xlabel("Time (ms)")
ax.set_ylabel("Count Percentage")

fig, ax = plt.subplots()
ax.grid()
ax.plot([1000 * dt * i for i in range(len(output_change))], output_change)
ax.set_xlabel("Time (ms)")
ax.set_ylabel("Count Rate Percentage")
plt.show()