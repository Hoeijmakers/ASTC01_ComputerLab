import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
from scipy import interpolate

# Constants

jup_rad = 69911000 # in meter
sun_rad = 696340000 # in meter
au = 149597870691 # in meter
days_to_hours = 24*60


class System:
    def __init__(self, name, pl_mass, pl_radius, pl_period, pl_semi, st_mass, st_radius, ecc, incl, long):
        self.name = name
        self.pl_mass = pl_mass # in Jupiter masses
        self.pl_radius = pl_radius # in Jupiter radii
        self.pl_period = pl_period # in days
        self.pl_semi = pl_semi # in AU
        self.st_mass = st_mass # in solar masses
        self.st_radius = st_radius # in solar radii
        self.ecc = ecc # number
        self.incl = incl # in degrees
        self.long = long # in degrees

    
    def TransitParams(self):
        
        # Transit depth
        self.delta = (self.pl_radius * jup_rad / (self.st_radius * sun_rad))**2
        
        # Total transit time
        self.tT = self.pl_period * days_to_hours * (self.st_radius * sun_rad) / (np.pi * self.pl_semi * au) * np.sqrt((1 + (self.pl_radius*jup_rad) / (self.st_radius*sun_rad))**2 - (((self.pl_semi * au) / (self.st_radius * sun_rad)) *np.cos(self.incl))**2)
        
        # Full transit time
        self.tF = np.sqrt(self.tT**2 - (4*(self.pl_period*days_to_hours)**2*np.sqrt(self.delta)*(self.st_radius*sun_rad)**2) / (np.pi**2* (self.pl_semi * au)**2))
        
        # Transit parameter
        #self.b = 
    
    
    #def Image(self):
    
    
    def Transit(self):
        
        # theoretical transit light curve
    
        t = np.linspace(int(-self.tT /2. - self.tT/5.), int(self.tT /2. + self.tT/5.),int(self.tT + 2*self.tT/5.))
        normalised_brightness = np.ones(len(t))
        dip = normalised_brightness - self.delta
        x_1 = [-self.tT / 2., -self.tF/2.]
        y_1 =[1, 1-self.delta]
        x_2 = [self.tF/2., self.tT / 2.]
        y_2 =[1-self.delta, 1]

        f_int1 = sc.interpolate.interp1d(x_1, y_1)
        f_int2 = sc.interpolate.interp1d(x_2, y_2)
        x1_new = np.linspace(x_1[0], x_1[1], 50)
        x2_new = np.linspace(x_2[0], x_2[1], 50)

        y1_new = f_int1(x1_new)
        y2_new = f_int2(x2_new)

        

        ### Creating some 'Data'

        sigma = 0.002

        dip_data = dip + np.random.rand(len(dip))/100. 
        dip_data1 = dip - np.random.rand(len(dip))/100. 
        bright1 = normalised_brightness + np.random.rand(len(normalised_brightness) )/100.
        bright2 = normalised_brightness - np.random.rand(len(normalised_brightness) )/100.
        y1_data1 = y1_new + np.random.rand(len(y1_new) )/90.
        y1_data2 = y1_new - np.random.rand(len(y1_new) )/90.
        y2_data1 = y2_new + np.random.rand(len(y2_new) )/90.
        y2_data2 = y2_new - np.random.rand(len(y2_new) )/90.
        
                
        plt.figure(figsize=(20,20))
        plt.plot(t[np.where(t < -self.tT / 2.)], normalised_brightness[np.where(t < -self.tT / 2.)], 'b--', label='Theoretical Transit Curve')
        plt.plot(t[np.where(t > self.tT / 2.)], normalised_brightness[np.where(t > self.tT / 2.)], 'b--')
        plt.plot(t[np.where((t > -self.tF / 2.) * (t < self.tF / 2.))], dip[np.where((t > -self.tF / 2.) * (t < self.tF / 2.))], 'b--')

        plt.plot(t[np.where((t > -self.tF / 2.) * (t < self.tF / 2.))], dip_data[np.where((t > -self.tF / 2.) * (t < self.tF / 2.))],'g.')
        plt.plot(t[np.where((t > -self.tF / 2.) * (t < self.tF / 2.))], dip_data1[np.where((t > -self.tF / 2.) * (t < self.tF / 2.))],'g.')
        plt.plot(t[np.where(t < -self.tT / 2.)], bright1[np.where(t < -self.tT / 2.)], 'g.', label='Data')
        plt.plot(t[np.where(t < -self.tT / 2.)], bright2[np.where(t < -self.tT / 2.)], 'g.')
        plt.plot(t[np.where(t > self.tT / 2.)], bright1[np.where(t > self.tT / 2.)], 'g.')
        plt.plot(t[np.where(t > self.tT / 2.)], bright2[np.where(t > self.tT / 2.)], 'g.')
        plt.plot(x1_new, y1_data1, 'g.')
        plt.plot(x1_new, y1_data2, 'g.')
        plt.plot(x2_new, y2_data1, 'g.')
        plt.plot(x2_new, y2_data2, 'g.')
   




        plt.plot(x_1, y_1, 'bo--')
        plt.plot(x_2, y_2, 'bo--')
        plt.legend()
        
        plt.ylim((dip[0]-0.02, normalised_brightness[0] + 0.02))
        plt.show()
        
        
    
    #def TransitMeasurement(self):


#sys1 = System('WASP-121b', 1.183, 1.865, 3.34, 0.0558, 1.45, 1.5, 0,  np.pi / 2., 0)

#sys1.TransitParams()
#sys1.Transit()

#print(sys1.delta)
