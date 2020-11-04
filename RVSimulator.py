import numpy as np
import matplotlib.pyplot as plt

### Constants and conversions

jupiter_mass_to_kg = 1.898e27
solar_mass_to_kg = 1.989e30
days_to_seconds = 86400.
au_to_meter = 1.496e11
G = 6.67408e-11 #m3 kg-1 s-2

### Kepler Solver based on http://alpheratz.net/dynamics/twobody/KeplerIterations_summary.pdf

def keplerstart(e,M):
    t34=e**2.
    t35=e*t34
    t33=np.cos(M)
    result=M+(-0.5*t35+e+(t34+1.5*t33*t35)*t33)*np.sin(M)
    return result 

def eps3(e,M,x):
    t1=np.cos(x)
    t2=-1.+e*t1
    t3=np.sin(x)
    t4=e*t3
    t5=-1.*x+t4+M
    t6=t5/(0.5*t5*t4/t2+t2)
    result=t5/((0.5*t3-1/6*t1*t6)*e*t6+t2)
    return result

def keplersolve(e,M,precision):
    Mnorm=np.mod(M,2.*np.pi)
    E0=keplerstart(e,Mnorm)
    dE=precision+1.
    count=0
    while dE > precision:
        #print count
        E=E0-eps3(e,Mnorm,E0)
        dE=np.absolute(E-E0)
        E0=E
        count=count+1
        if count == 100:
            print('Failed to converge!')
            return 0
    return E


class RadialVelocity:

    def __init__(self, pl_semi, pl_period, ecc, incl, longitude, pl_msini, st_mass, st_dist):

        self.pl_semi = pl_semi
        self.pl_period = pl_period * days_to_seconds 
        self.ecc = ecc
        self.incl = incl
        self.longitude = longitude
        self.pl_msini = pl_msini
        self.st_mass = st_mass
        self.st_dist = st_dist

    
    def Kepler(self):

        arr_size = self.pl_period
        self.day_arr=np.arange(arr_size+1,dtype=np.float)        #these are the hours on orbit counting from periastron
        M_array=np.arange(arr_size+1,dtype=np.float)        # mean anomaly array
        M_array=2.*np.pi/self.pl_period*M_array             #this is an array for the fraction of the orbit based on mean anomaly
        self.E_array=np.arange(arr_size+1,dtype=np.float)   #this is the array for the eccentric anomaly 


        for t in range(len(self.E_array)):
            result=keplersolve(ecc, M_array[t], 1e-14)
            self.E_array[t]=result


        #Array of the true anomaly (based on equation from HM Schmid's Diploma Thesis)
        self.phi_array= np.arctan(np.tan(self.E_array/2.)/(np.sqrt((1-self.ecc)/(1+self.ecc))))*2.

    def RadialVelocityCurve(self):

        K_term1=((self.pl_msini*jupiter_mass_to_kg)/(self.pl_msini*jupiter_mass_to_kg+self.st_mass*solar_mass_to_kg))
        K_term2=(self.pl_semi*au_to_meter*np.sin(self.incl/180.*np.pi))/(np.sqrt(1.-self.ecc*self.ecc))
        K=K_term1*K_term2*np.sqrt(G)

        #Radial velocity for each hour / each value of phi 
        self.v_r=K*(np.cos((self.longitude/180.*np.pi)+self.phi_array)+self.ecc*np.cos((self.longitude/180.*np.pi)))

        #Compute the apparent separation and the best hours to observe the planet
        #These are the indices for the fastest and slowest velocities
        index1=np.where(self.v_r == np.max(self.v_r))
        best1=int(index1[0])

        index2=np.where(self.v_r==np.min(self.v_r))
        best2=int(index2[0])

        #these are the best hours on orbit:
        print('The best hours to image the planet (after perihelion) are:', self.day_arr[best1], self.day_arr[best2])
        #print phi_array[best1]*180/np.pi, phi_array[best2]*180/np.pi

        #These are the corresponding separations:
        r_best1_au=a*(1-ecc*np.cos(self.E_array[best1]))
        r_best2_au=a*(1-ecc*np.cos(self.E_array[best2]))

        r_best1_arcsec=r_best1_au/distance
        r_best2_arcsec=r_best2_au/distance

        print('The corresponding separations [in AU] are:', r_best1_au,r_best2_au)
        print('The projected separations [in arcsec] are:',r_best1_arcsec,r_best2_arcsec)


        plt.figure()
        plt.plot(self.day_arr/self.pl_period, self.v_r)
        plt.title('Radial Velocity Curve')
        plt.xlabel('Phase')
        plt.ylabel('Radial velocity [m/s]')
        plt.show()




