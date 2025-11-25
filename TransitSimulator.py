import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
from scipy import interpolate
import numpy.random as rand
# Constants

jup_rad = 69911000 # in meter
sun_rad = 696340000 # in meter
au = 149597870691 # in meter
days_to_minutes = 24*60
jupiter_mass_to_kg = 1.898e27
solar_mass_to_kg = 1.989e30
days_to_seconds = 86400.
au_to_meter = 1.496e11
G = 6.67408e-11 #m3 kg-1 s-2



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
        if count == 1000:
            print('Failed to converge!')
            return 0
    return E







class System:
    def __init__(self, name, pl_mass, pl_radius, pl_period, st_mass, st_radius, incl = 90.0, ecc = 0.0, long=0.0,debug=False):

        self.name = name
        self.pl_mass = pl_mass # in Jupiter masses
        self.pl_radius = pl_radius # in Jupiter radii
        self.pl_period = pl_period # in days
        self.st_mass = st_mass # in solar masses
        self.st_radius = st_radius # in solar radii
        self.ecc = ecc # number
        self.incl = incl # in radians
        self.long = long # in degrees
        self.longitude = self.long
        self.pl_semi = ( (pl_period * days_to_seconds)**2 * G * st_mass * solar_mass_to_kg / (4*np.pi**2)) ** (1/3) / au_to_meter
        self.pl_msini = self.pl_mass * np.sin(incl)
        if debug:
            print('a :', self.pl_semi)
            print('msini :',self.pl_msini)
            print('P :',self.pl_period)

        self.TransitParams()
    
    
    def TransitParams(self,debug=False):

        # Transit depth
        self.delta = (self.pl_radius * jup_rad / (self.st_radius * sun_rad))**2
        
        # Total transit time
        self.tT = self.pl_period * days_to_minutes * (self.st_radius * sun_rad) / (np.pi * self.pl_semi * au) * np.sqrt((1 + (self.pl_radius*jup_rad) / (self.st_radius*sun_rad))**2 - (((self.pl_semi * au) / (self.st_radius * sun_rad)) *np.cos(self.incl))**2)
        
        # Full transit time
        self.tF = np.sqrt(self.tT**2 - (4*(self.pl_period*days_to_minutes)**2*np.sqrt(self.delta)*(self.st_radius*sun_rad)**2) / (np.pi**2* (self.pl_semi * au)**2))
        
        if debug:
            print('Delta: ',self.delta)
            print('Total transit time: ',self.tT)
            print('Full transit time: ',self.tF)
    #def Image(self):
    
    
    def TransitModel(self,baseline=1.0,Nsteps=10000,plot=True,fsize=6):
        
        # theoretical transit light curve
    
        t = np.linspace(1/Nsteps, int(self.tT /2. + self.tT*baseline/2),Nsteps)
        f = np.ones(len(t))
        a = -1/(self.tT/2 - self.tF/2)
        b = -1*a*self.tT/2

        W = np.clip(a*t+b,0,1)
        f -= W*self.delta

        self.T = np.concatenate([np.flip(-t),t]) # In seconds.
        self.F = np.concatenate([np.flip(f),f])
        
        if plot:               
            plt.figure(figsize=(fsize*4/3,fsize))
            plt.title(f'Modelled transit Light Curve of {self.name}', size=16)
            plt.plot(self.T/60,self.F, 'C0', label='Theoretical Transit Curve')
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)
            plt.xlabel('Time [h]', size=20)
            plt.ylabel('Normalised flux', size=20)
            plt.ylim((1-1.2*self.delta, 1 + 0.2*self.delta))
            plt.show()


    def TransitData(self,c,SNR,plot=True,fsize=6):
        self.TransitModel(plot=False)
        TD = np.arange(np.min(self.T),np.max(self.T),c) # Note that T is already in minutes.
        FD = np.interp(TD,self.T,self.F)
        FDn = FD + rand.normal(loc=0.0,scale=SNR,size=len(FD))
        if plot:               
            plt.figure(figsize=(fsize*4/3,fsize))
            plt.title(f'Modelled transit Light Curve of {self.name}', size=16)
            plt.plot(self.T/60,self.F, 'C0', label='Theoretical Transit Curve')
            plt.errorbar(TD/60,FDn,fmt='.',yerr=SNR,color='black', label = 'Simulated observations')
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)
            plt.xlabel('Time [h]', size=20)
            plt.ylabel('Normalised flux', size=20)
            plt.ylim((1-1.2*self.delta-3*SNR, 1 + 0.2*self.delta +3*SNR))
            plt.legend()
            plt.show()

    
    def Kepler(self):
        self.M_array = np.linspace(0.0,2.0*np.pi,200)
        self.E_array = self.M_array*0.0
        for t in range(len(self.M_array)):
            result=keplersolve(self.ecc, self.M_array[t], 1e-14)
            self.E_array[t]=result


        beta = self.ecc / (1+ np.sqrt(1-self.ecc**2))

        self.phi_array = self.E_array + 2 * np.arctan((beta * np.sin(self.E_array))/(1-beta*np.cos(self.E_array)))
        # self.phi_array= np.arctan(np.tan(self.E_array/2.)/(np.sqrt((1-self.ecc)/(1+self.ecc))))*2.


    def RadialVelocityModel(self,fsize=6,debug=False,plot=True):
        self.Kepler()

        if debug:
            plt.figure(figsize=(fsize*4/3,fsize))
            plt.plot(self.M_array,self.E_array,label='E anomaly')
            plt.plot(self.M_array,self.phi_array,label='T anomaly')
            plt.legend()
            plt.show()
        
        K = np.sqrt(G*self.st_mass*solar_mass_to_kg / (self.pl_semi*au_to_meter)) * (self.pl_msini *jupiter_mass_to_kg) / (self.st_mass * solar_mass_to_kg) 
        self.K = K
        if debug:
            print('G',G)
            print('Mstar',self.st_mass)
            print('Msol',solar_mass_to_kg)
            print('a',self.pl_semi)
            print('AU',au_to_meter)
            print('msini',self.pl_msini)
            print('Mjup',jupiter_mass_to_kg)
            print('K',self.K)
        #Radial velocity for each hour / each value of phi 
        self.v_r=K*(np.cos((self.longitude/180.*np.pi)+self.phi_array)+self.ecc*np.cos((self.longitude/180.*np.pi)))


        if plot:
            plt.figure(figsize=(fsize*4/3,fsize))
            plt.plot(self.M_array/(2*np.pi), self.v_r)
            plt.title(f'Radial Velocity Curve of {self.name}', size=16)
            plt.xlabel(r'Orbital phase', size=20)
            plt.ylabel('Radial velocity [m/s]', size=20)
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)
            plt.show()


    def RadialVelocityData(self,N_points,SNR,plot=True,fsize=6):
        self.RadialVelocityModel(plot=False)
        TD = np.linspace(np.min(self.M_array),np.max(self.M_array),int(N_points)) # In seconds.
        FD = np.interp(TD,self.M_array,self.v_r)
        FDn = FD + rand.normal(loc=0.0,scale=SNR,size=len(FD))
        if plot:               
            plt.figure(figsize=(fsize*4/3,fsize))
            plt.title(f'Modelled radial velocity curve of {self.name}', size=16)
            plt.plot(self.M_array / (2*np.pi),self.v_r, 'C0', label='Theoretical RV Curve')
            plt.errorbar(TD / (2*np.pi),FDn,fmt='.',yerr=SNR,color='black', label = 'Simulated observations')
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)
            plt.xlabel('Orbital phase', size=20)
            plt.ylabel('Radial velocity [m/s]', size=20)
            plt.ylim((-1.1*self.K-3*SNR, 1.1*self.K +3*SNR))
            plt.legend()
            plt.show()
