import numpy as np
from numba import njit, float64, int64, boolean
import numba as nb
k_B = 1.381e-23 # J/K

@njit(float64[:](float64[:]))
def numba_erf(s):
    z2 = np.abs(s)
    t = 1 / (1 + 0.32759109962 * z2)
    res = (    - 1.061405429 ) * t
    res = (res + 1.453152027 ) * t
    res = (res - 1.421413741 ) * t
    res = (res + 0.2844966736) * t
    res =((res - 0.254829592 ) * t) * np.exp(-z2*z2)
    res = res + 1
    return res*np.sign(s)

@njit(nb.types.Tuple((float64[:], float64[:], float64[:], int64))(float64[:], float64[:], float64[:], int64, float64, boolean[:], int64, float64, float64, float64, float64, float64))
def handle_collisions(vxs, vys, vzs, N_sim, t_step, not_out_of_bounds, T, mass, mass_gas, reduced_mass, 
                      buffer_gas_n_density, cross_section_ion_buffer_gas):
    
    mean_gas_speed = np.sqrt(8*k_B*T / (np.pi*mass_gas))
    median_gas_speed = np.sqrt(2*k_B*T / mass_gas)

    vs = np.sqrt(vxs**2 + vys**2 + vzs**2)
    s = vs / median_gas_speed
    relative_speed = mean_gas_speed*((s+1/(2*s))*np.sqrt(np.pi)/2*numba_erf(s) + 1/2 * np.exp(-s**2))

    collision_prob = 1 - np.exp(- relative_speed * t_step * buffer_gas_n_density * cross_section_ion_buffer_gas)
    random_ns = np.random.random(N_sim)
    should_collide = (random_ns < collision_prob) * not_out_of_bounds

    # Calculate changes in velocity based on statistical description - Delahaye 2019
    std_x = np.sqrt(1 + 2*k_B*T / (mass_gas*vxs**2))
    std_y = np.sqrt(1 + 2*k_B*T / (mass_gas*vys**2))
    std_z = np.sqrt(1 + 2*k_B*T / (mass_gas*vzs**2))

    epsilon_x = np.empty(N_sim)
    epsilon_y = np.empty(N_sim)
    epsilon_z = np.empty(N_sim)
    for i in range(N_sim):
        epsilon_x[i] = np.random.normal(0, std_x[i])
        epsilon_y[i] = np.random.normal(0, std_y[i])
        epsilon_z[i] = np.random.normal(0, std_z[i])

    delta_vxs = (-reduced_mass / mass * vxs * (1+epsilon_x)) * should_collide
    delta_vys = (-reduced_mass / mass * vys * (1+epsilon_y)) * should_collide
    delta_vzs = (-reduced_mass / mass * vzs * (1+epsilon_z)) * should_collide
    return delta_vxs, delta_vys, delta_vzs, np.sum(should_collide)