import numpy as np

class CIT:
    def __init__(self,omega_RF,xi,z0,r0):
        #Initalzie the CIT class,
        #Parameters:
        # omega_RF -- RF frequency of trap
        # xi -- geometric parameter used for calculating the a and q parameters (see https://pubs.acs.org/doi/10.1021/jasms.3c00349)
        # z0 -- distance to endcap electrodes from center of trap
        # r0 -- distance to cylinder walls from center of trap

        self.omega_RF = omega_RF
        self.xi = xi
        self.z0 = z0
        self.r0 = r0
    
    def pn(self,n):
        #Determine pn used for calculating the potential. Refer to https://doi.org/10.1016/0020-7381(77)80034-X for more information
        #Parameters:
        # n -- integer

        return np.pi*(2*n+1)/(2*self.z0)

    def CITPotential(self,r,z,V0,tol = 1e-3):
        #Return the CIT potential at (r,z)
        # r -- float with the radial coordinate at which potential is to be calcualted
        # z -- float with the axial coordinate at which potential is to be calculated
        # V0 -- float equal to potential difference between the cylinder wall and endcaps
        # tol -- relative tolerance at which to stop summing, when calculating the potential


        #Since potential is an infinite sum, we loop until potential converges
        p0 = self.pn(0)
        InitVal = 4 * (1/np.pi) * (np.i0(p0*r)/np.i0(p0*self.r0))*np.cos(p0*z)

        #Initalize the CurrValues and NextValues array, 
        CurrValues = [InitVal]
        NextValues = [InitVal]
        
        #Determine the potential using only n = 0
        CurrentPotential = np.sum(CurrValues)


        p1 = self.pn(1)
        NextVal = 4 * (-1)/(3*np.pi) * (np.i0(p1*r)/np.i0(p1*self.r0)) * np.cos(p1*z)
        NextValues.append(NextVal)

        #Determine potential using n = 0 --> 1
        NextPotential = np.sum(NextValues)

        #Start looping, increasing the number of terms in sum until a relative tolerance of tol is reached
        n=2
        while n<11:#np.abs((CurrentPotential-NextPotential)/CurrentPotential)>tol:
            CurrValues.append(NextVal)
            NextVal = 4*(-1)**n/((2*n+1)*np.pi) * (np.i0(self.pn(n)*r)/np.i0(self.pn(n)*self.r0)) * np.cos(self.pn(n)*z)
            
            NextValues.append(NextVal)
            CurrentPotential = np.sum(CurrValues)
            NextPotential = np.sum(NextValues)

            n+=1
        return NextPotential*V0

    def Term(self,r,z,n):
        #Return the nth term of the CIT potential at (r,z)
        # r -- float with the radial coordinate at which potential is to be calcualted
        # z -- float with the axial coordinate at which potential is to be calculated
        # n -- integer with the term number

        return 4*(-1)**n/((2*n+1)*np.pi) * (np.i0(self.pn(n)*r)/np.i0(self.pn(n)*self.r0)) * np.cos(self.pn(n)*z)
    def CITPotentialV(self,r,z,V0,nmax,tol = 1e-3):
        #Return the CIT potential at (r,z)
        # r -- float with the radial coordinate at which potential is to be calcualted
        # z -- float with the axial coordinate at which potential is to be calculated
        # V0 -- float equal to potential difference between the cylinder wall and endcaps
        # nmax -- integer with the maximum number of terms to be summed

        ns = np.array(range(nmax))
        terms = self.Term(r,z,ns)
        return np.sum(terms)*V0
