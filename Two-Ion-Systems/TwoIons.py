
import numpy as np
from scipy.optimize import fsolve
import math

class Trap:
  #Class that holds different constants relevant for the trap.
  def __init__(self,omega_RF,kappa,z0,r0):
    #Initialization of the class
    #Variables:
    # omega_rf - RF frequency of trap
    # kappa - geometric constant of the trap. This is usually obtained from simulations.
    # z0 half the distance between endcap electrodes
    # r0 half the distance between diametrically opposing rods in the trap.
    self.omega_RF = omega_RF
    self.kappa = kappa
    self.z0 = z0
    self.r0 = r0
  def Trap_secular_frequencies(self,a,q):
    #Returns the single-ion frequencies [omega_z,omega_r]
    #Variables:
    # a - a parameter, descrbies the strength of the endcap voltage
    # q - q parameter, describes the strength of the RF voltage
    return np.array([self.omega_RF*np.sqrt(np.abs(a/2)),self.omega_RF/2*np.sqrt(q**2/2+a)])




class two_ion_system:
  # Class holding many of the methods for a system of two particular ions in a given trap
  def __init__(self,m1,q1,m2,q2,trap,*,qrange=np.linspace(0,0.9,1000)):
    #Initializes the two-ion system class.

    #Variables:
    # m1 - mass of ion 1
    # m2 - mass of ion 2
    # q1 - charge of ion 1
    # q2 - charge of ion 2
    # trap - trap object that the ions sit in, made by using the Trap class in this code
    # qrange - range of q parameters used for some of the methdos within.
    self.m1 = m1
    self.m2 = m2
    self.q1 = q1
    self.q2 = q2
    self.trap = trap
    self.qrange = qrange
    self.omega_RF = self.trap.omega_RF
    self.stable_q = None
    self.stable_a = None
  def xi_parameter(self,dicke_param,n_start,n_end):
    # xi function as defined in https://arxiv.org/abs/2004.02959
    if n_start == 0:
      return 0
    else:
      return dicke_param*np.sqrt(math.factorial(np.max([n_start,n_end]))/(math.factorial(np.min([n_start,n_end]))))

  def lower_curve(self,q):
    #Returns the lower-stability line for a given q- parameter
    #Variables:
    # q - array or single q parameter value
    return -1/2*q**2+7/128*q**4-29/2304*q**6+69697/18874368*q**8
  def higher_curve(self,q):
    #Returns the upper stability curve for a given q - parameter
    #Variables:
    # q - array or single q parameter value
    return 1-q-1/8*q**2+1/64*q**3-1/1536*q**4-11/36864*q**5

  def Master_Stability(self): #Create the stable area for a given q range
      lower_as = self.lower_curve(self.qrange)
      higher_as = self.higher_curve(self.qrange)
      lowest_a = min(min(lower_as),min(higher_as))
      a_list = np.linspace(0,lowest_a,1000)
      mu = self.m2/self.m1
      rho = self.q2/self.q1
      Both_Stable_a = []
      Both_Stable_q = []
      Ion1_Stab_List_a = []
      Ion1_Stab_List_q = []
      Ion2_Stab_List_a = []
      Ion2_Stab_List_q = []
      Both_Stable_plus_axis_a = []
      Both_Stable_plus_axis_q = []
      for q in self.qrange:
        for a in a_list:
          ion1_stable = False
          ion2_stable = False
          along_axis = False #Reset bools
          omega1z,omega2z,omega1r,omega2r = self.secular_frequencies(a,q)
          zeta_z = self.m1*self.m2*omega1z**2*omega2z**2/(self.m1*omega1z**2+self.m2*omega2z**2)
          zeta_r = self.m1*self.m2*omega1r**2*omega2r**2/(self.m1*omega1r**2+self.m2*omega2r**2)
          #if a>self.lower_curve(q) and rho/mu*a>self.lower_curve(rho/mu*q) and a<self.higher_curve(q) and rho/mu*a<self.higher_curve(rho/mu*q):# and zeta_z<zeta_r:
          if a>self.lower_curve(q) and a<self.higher_curve(q):
            ion1_stable = True
          if rho/mu*a>self.lower_curve(rho/mu*q) and rho/mu*a<self.higher_curve(rho/mu*q):
            ion2_stable = True
          if zeta_z<zeta_r:
            along_axis = True
          if ion2_stable and ion1_stable and along_axis:
            Both_Stable_plus_axis_a.append(a)
            Both_Stable_plus_axis_q.append(q)
          elif ion2_stable and ion1_stable:
            Both_Stable_a.append(a)
            Both_Stable_q.append(q)
          elif ion2_stable:
            Ion2_Stab_List_a.append(a)
            Ion2_Stab_List_q.append(q)
          elif ion1_stable:
            Ion1_Stab_List_a.append(a)
            Ion1_Stab_List_q.append(q)

      return [Both_Stable_plus_axis_a,Both_Stable_plus_axis_q],[Both_Stable_a,Both_Stable_q],[Ion1_Stab_List_a,Ion1_Stab_List_q],[Ion2_Stab_List_a,Ion2_Stab_List_q]
  def create_stable_area(self): 
    #Create the stable two-ion stability area using the q's from self.qrange
    # The stability range is calculated as the area that fulfills stability of both ions if they were alone, while the two-ion system is still trapped along the trap axis
    lower_as = self.lower_curve(self.qrange)
    higher_as = self.higher_curve(self.qrange)
    lowest_a = min(min(lower_as),min(higher_as))
    a_list = np.linspace(0,lowest_a,1000)
    mu = self.m2/self.m1
    rho = self.q2/self.q1
    self.stable_q = []
    self.stable_a = []
    # print(self.qrange)
    # print(a_list)
    for q in self.qrange:
      for a in a_list:
        omega1z,omega2z,omega1r,omega2r = self.secular_frequencies(a,q)
        zeta_z = self.m1*self.m2*omega1z**2*omega2z**2/(self.m1*omega1z**2+self.m2*omega2z**2)
        zeta_r = self.m1*self.m2*omega1r**2*omega2r**2/(self.m1*omega1r**2+self.m2*omega2r**2)
        if a>self.lower_curve(q) and rho/mu*a>self.lower_curve(rho/mu*q) and a<self.higher_curve(q) and rho/mu*a<self.higher_curve(rho/mu*q) and zeta_z<zeta_r:
          self.stable_q.append(q)
          self.stable_a.append(a)
    return None
  def secular_frequencies(self,a,q): 
    #Calculate the secular frequencies of the ions. i.e. if they were alone in the trap.
    #Variables: 
    # a - a parameter
    # q - q parameter

    #Returns: 
    #omega1z - axial frequency of ion 1
    #omega2z - axial frequency of ion 2
    #omega1r - radial frequency of ion 1
    #omega2r - radial frequency of ion 2
    rho = self.q2/self.q1
    mu = self.m2/self.m1
    a2 = rho/mu *a
    q2 = rho/mu*q

    omega1z = self.omega_RF * np.sqrt(np.abs(a)/2)
    omega2z = self.omega_RF * np.sqrt(np.abs(a2)/2)
    omega1r = self.omega_RF / 2 * np.sqrt(q**2 / 2 + a)
    omega2r = self.omega_RF / 2 * np.sqrt(q2**2 / 2 + a2)
    return omega1z,omega2z,omega1r,omega2r

  def vibrational_freqs_and_modes(self,a,q): 
      #Calculate eigenmodes and eigenfrequencies of axial and radial direction.
      # Variables:
      # a- a parameter
      # q - q parameter

      #Returns:
      # in_phase_z - frequency of axial in-phase mode in Hz
      # out_of_phase_z - frequency of the axial out of ohagse mode in Hz
      # in_phase_r - frequency of the in phase radial mode in Hz
      # out_of_phase_r - frequency of the out of phase radial mode in Hz
      #vec_in_phase_z - normalized in-phase axial eigenmode
      #vec_out_of_phase_z - normalizeed out-of-phase axial eigenmode
      #vec_in_phase_r - normalized in-phase radial eigenmode
      #vec_out_of_phase_r - normalized out-of-phase radial eigenmode
      omega1z,omega2z,omega1r,omega2r = self.secular_frequencies(a,q)
      mu = self.m2/self.m1
      rho = self.q2/self.q1
      K11 = omega1z**2*(1+2/(1+1/rho))
      K22 = omega2z**2*(1+2/(1+rho))
      K12 = -2*omega1z**2/(np.sqrt(mu)*(1+1/rho))
      K33 = omega1r**2-omega1z**2/(1+1/rho)
      K44 = omega2r**2-omega1z**2/(mu*(1+1/rho))
      K34 = -1/2*K12

      in_phase_z = np.sqrt((K11+K22-np.sqrt((K11-K22)**2+4*K12**2))/2)
      out_of_phase_z = np.sqrt((K11+K22+np.sqrt((K11-K22)**2+4*K12**2))/2)
      in_phase_r = np.sqrt((K33+K44+np.sqrt((K33-K44)**2+4*K34**2))/2)
      out_of_phase_r = np.sqrt((K33+K44-np.sqrt((K33-K44)**2+4*K34**2))/2)

      a_plus_normalization = 1/np.sqrt(1+(K12/(out_of_phase_z**2-K22))**2)
      a_plus_1 = a_plus_normalization
      a_plus_2 = a_plus_normalization*(K12/(out_of_phase_z**2-K22))
      vec_out_of_phase_z = [a_plus_1,a_plus_2]
      a_minus_normalization = 1/np.sqrt(1+(K12/(in_phase_z**2-K22))**2)
      a_minus_1 = a_minus_normalization
      a_minus_2 = a_minus_normalization*(K12/(in_phase_z**2-K22))
      vec_in_phase_z = [a_minus_1,a_minus_2]

      b_plus_normalization = 1/np.sqrt(1+(K34/(in_phase_r**2-K44))**2)
      b_plus_1 = b_plus_normalization
      b_plus_2 = b_plus_normalization*(K34/(in_phase_r**2-K44))
      vec_in_phase_r = [b_plus_1,b_plus_2]
      b_minus_normalization = 1/np.sqrt(1+(K34/(out_of_phase_r**2-K44))**2)
      b_minus_1 = b_minus_normalization
      b_minus_2 = b_minus_normalization*(K34/(out_of_phase_r**2-K44))
      vec_out_of_phase_r = [b_minus_1,b_minus_2]

      return [in_phase_z,out_of_phase_z,in_phase_r,out_of_phase_r],[vec_in_phase_z,vec_out_of_phase_z,vec_in_phase_r,vec_out_of_phase_r]

  def create_resonance_line_1oz_2or(self):
      #Calculates resonance line for the parametric resonance of 1 out-of-phase axial mode to 2 out of phagse radial. i.e. where omega_oz = 2 omega_or
      #Returns:
      # resonant_as - array of a- values where resonance occurs, each corresponds to the same index of q-value from self.qrange
      resonant_as = []
      for q in self.qrange:
        def helper_oo(a):
          return self.vibrational_freqs_and_modes(a,q)[0][1]-2*self.vibrational_freqs_and_modes(a,q)[0][3]
        resonant_a = fsolve(helper_oo,self.lower_curve(q)/10)[0]
        resonant_as.append(resonant_a)
      return resonant_as
  def create_resonance_line_1iz_2or(self):
    #Calculates resonance line for the parametric resonance of 1 in-phase axial mode to 2 out of phagse radial. i.e. where omega_oz = 2 omega_or
    #Returns:
    # resonant_as - array of a- values where resonance occurs, each corresponds to the same index of q-value from self.qrange
    resonant_as = []
    for q in self.qrange:
      def helper_io(a):
        return self.vibrational_freqs_and_modes(a,q)[0][0]-2*self.vibrational_freqs_and_modes(a,q)[0][3]
      resonant_a = fsolve(helper_io,self.lower_curve(q)/10)[0]
      resonant_as.append(resonant_a)
    return resonant_as
  def transfer_rate_1oz_2or(self,a,q):
      #Returns the transfer rate of parametric transfer from 2 out of phase radial phonons to 1 out of phase axial phonon, or vice versa.
      #Variables:
      # a - a parameter
      # q - q parameter
      eps0 = 8.85*1e-12 #Vacuum Permittivity
      hbar = 1.05*1e-34 # hbar
      freqs = self.vibrational_freqs_and_modes(a,q)[0]
      vecs = self.vibrational_freqs_and_modes(a,q)[1]
      #print(np.array(freqs)*1e-3/(2*np.pi))
      #print(vecs)
      mu = self.m2/self.m1
      rho = self.q2/self.q1
      g_z = np.sqrt(hbar/(2*self.m1*freqs[1]))
      g_r = np.sqrt(hbar/(2*self.m1*freqs[3]))
      #print(np.array(freqs)*1e-3/(2*np.pi))
      #print(g_z,g_r)
      z1_eq = (self.q1*self.q2/(4*np.pi*eps0*self.m1*(self.omega_RF*np.sqrt(abs(a)/2))**2*(1+1/rho)**2))**(1/3)
      gamma = 3*self.q1*self.q2/(8*np.pi*eps0*(z1_eq)**4*(1+1/rho)**4)
      a_plus = vecs[1]
      b_minus = vecs[3]
      coeff = 3*np.sqrt(2)*gamma*g_z*g_r**2/6*(a_plus[0]*b_minus[0]**2+1/mu*a_plus[0]*b_minus[1]**2-1/np.sqrt(mu)*a_plus[1]*b_minus[0]**2-1/mu**(3/2)*a_plus[1]*b_minus[1]**2-2/np.sqrt(mu)*a_plus[0]*b_minus[0]*b_minus[1]+2/mu*a_plus[1]*b_minus[0]*b_minus[1])
      coeff_hbar = coeff/hbar
      #print(coeff)
      return coeff_hbar
  def z1_eq(self,a):
    # Calculates and returns the equilibrium position of ion 1 in the trap. From thos the equilibrium position of ion 2 is simple to calculate.
    eps0 = 8.85*1e-12 #Vacuum Permittivity
    mu = self.m2/self.m1
    rho = self.q2/self.q1
    return (self.q1*self.q2/(4*np.pi*eps0*self.m1*(self.omega_RF*np.sqrt(abs(a)/2))**2*(1+1/rho)**2))**(1/3)
  def transfer_rate_1iz_2or(self,a,q):
      #Returns the transfer rate of parametric transfer from 2 out of phase radial phonons to 1 in phase axial phonon, or vice versa.
      #Variables:
      # a - a parameter
      # q - q parameter
      eps0 = 8.85*1e-12 #Vacuum Permittivity
      hbar = 1.05*1e-34 # hbar
      freqs = self.vibrational_freqs_and_modes(a,q)[0]
      vecs = self.vibrational_freqs_and_modes(a,q)[1]
      #print(np.array(freqs)*1e-3/(2*np.pi))
      #print(vecs)
      mu = self.m2/self.m1
      rho = self.q2/self.q1
      g_z = np.sqrt(hbar/(2*self.m1*freqs[0]))
      g_r = np.sqrt(hbar/(2*self.m1*freqs[3]))
      #print(np.array(freqs)*1e-3/(2*np.pi))
      #print(g_z,g_r)
      z1_eq = (self.q1*self.q2/(4*np.pi*eps0*self.m1*(self.omega_RF*np.sqrt(abs(a)/2))**2*(1+1/rho)**2))**(1/3)
      gamma = 3*self.q1*self.q2/(8*np.pi*eps0*(z1_eq)**4*(1+1/rho)**4)
      a_plus = vecs[0]
      b_minus = vecs[3]
      coeff = 3*np.sqrt(2)*gamma*g_z*g_r**2/6*(a_plus[0]*b_minus[0]**2+1/mu*a_plus[0]*b_minus[1]**2-1/np.sqrt(mu)*a_plus[1]*b_minus[0]**2-1/mu**(3/2)*a_plus[1]*b_minus[1]**2-2/np.sqrt(mu)*a_plus[0]*b_minus[0]*b_minus[1]+2/mu*a_plus[1]*b_minus[0]*b_minus[1])
      coeff_hbar = coeff/hbar
      #print(coeff)
      return coeff_hbar
