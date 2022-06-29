import numpy as np
import matplotlib.pyplot as plt
import re
from scipy.optimize import curve_fit

def media_movel(y,n):
    y = np.array(y)
    y_new = np.empty(len(y))
    for i in np.arange(len(y)-n):
        y_new[i] = np.mean(y[i:i+n+1])
    for j in np.arange(len(y)-n,len(y)):
        y_new[j] = np.mean([y[j]]+list(y_new[j-n:j])+list(y[j:])) 
    return y_new

class DATA:
    def __init__(self,material,group,aperture,power):
        PATH = 'Data/'
        file = f'zscan_{material}_g{group}_{aperture}_{power}mW/'
        self.power = power*1e-3
        self.material = material
        self.aperture = aperture
        
        self.positions = np.load(open(PATH+file+file[:-1]+'_dl_pos.npy','rb'))
        self.positions = np.abs(self.positions-np.max(self.positions))*1e-3 #meters

        self.power_points = np.load(open(PATH+file+file[:-1]+'_pwr.npy','rb')) #Watts
        self.log_file = open(PATH+file+file[:-1]+'_log.csv','r').read()
        values = ['start_pos','end_pos','step','exposure_time','avg_Power','aperture','no_lens']
        self.log = dict(zip(values,[float(n) for n in re.findall(r'\d+| \d+\.\d+|\.\d+',self.log_file)]))

        print(f'Material: {material} Aperture: {aperture} Power: {power}')
    
    def plot_power(self,title=None):
        title = f'Power Plot\n{self.material}-{self.aperture},{self.power}W' if title is None else title
        fig,ax = plt.subplots(constrained_layout=True)
        ax.set_title(title,fontsize=16)
        ax.plot(self.positions,self.power_points)
        return fig,ax
        
    def plot_normalized_trasmittance(self,title=None):
        title = f'$\Delta T(z)$\n{self.material}-{self.aperture},{self.power}W\n T(Z>>Z0) = {self.T[-1]}' if title is None else title
        fig,ax = plt.subplots( constrained_layout=True)
        ax.set_title(title,fontsize=14)
        ax.set_xlabel('$z$ $(m)$',fontsize=13)
        ax.set_ylabel('$\Delta T(z)$',fontsize=13)
        ax.plot(self.positions,self.DT)
        return fig,ax

class CLOSE(DATA):
    def __init__(self, material, group, aperture, power):
        super().__init__(material, group, aperture, power)

        self.T = self.power_points/self.power #Trasmitance
        self.DT = self.T/self.T[-1] #T/T(Z>>Z0)

        self.Tp = np.max(self.DT)
        self.Tv = np.min(self.DT)
        self.iTp = np.where(self.DT == self.Tp)[0][0]
        self.iTv = np.where(self.DT == self.Tv)[0][0]
        print(f'Peak: {self.Tp}\nValley: {self.Tv}')
        self.DTpv = np.abs(self.Tp-self.Tv)

        p_min = np.min([self.iTv,self.iTp])
        p_max = np.max([self.iTv,self.iTp])
        self.Zfocal_plane = self.positions[p_min]+(self.positions[p_max]-self.positions[p_min])/2
        self.iZfocal_plane = np.argmin(np.abs(self.DT[p_min:p_max]-1))+p_min
        self.positions -= self.Zfocal_plane

        self.z0 = np.abs(self.positions[p_max]-self.positions[p_min])/1.7 #meters

        print(f'Position of the focal plane: {self.Zfocal_plane} m, {self.positions[self.iZfocal_plane]} m')
        print(f'Rayleigh range: {self.z0} m')

    def set_S(self,Topen):
        self.S = self.T[-1]/Topen
        print(f'S: {self.S}')
    
    def set_expValues(self,wave_lenght,L_material):
        self.wave_lenght = wave_lenght #meters
        self.L = L_material #meters

        self.w0 = np.sqrt(self.wave_lenght*self.z0/np.pi) #meters^2
        self.I0 = 2*self.power/(np.pi*self.w0**2) #W/m^2
        print(f'W0: {self.w0} m^2\nI0: {self.I0} W/m^2')

    def calc_n2(self):
        self.DPhi = self.DTpv/(0.406*(1-self.S)**0.27)
        self.n2 = (self.wave_lenght/(2*np.pi))*self.DPhi/(self.I0*self.L)
        print(f'DPhi: {self.DPhi}')
        print(f'n2: {self.n2} m^2/W')


class OPEN(DATA):
    def __init__(self, material, group, aperture, power):
        super().__init__(material, group, aperture, power)
        self.noisy_power_points = np.copy(self.power_points)
        self.power_points = media_movel(self.power_points,7)
        self.T = self.power_points/self.power #Trasmitance
        self.DT = self.T/self.T[-1] #T/T(Z>>Z0)

    def set_focalPlane(self,Zfolcal_plane):
        self.Zfocal_plane = Zfolcal_plane
        self.positions -= self.Zfocal_plane
    
    def plot_noisy_filtered_points(self):
        fig,ax = plt.subplots(constrained_layout=True)
        ax.plot(self.positions,self.noisy_power_points,label='noisy')
        ax.plot(self.positions,self.power_points,label='filtered')
        ax.legend()
        title = f'$P(z)$\n{self.material}-{self.aperture},{self.power}W'
        ax.set_title(title,fontsize=14)
        ax.set_xlabel('$z$ $(m)$',fontsize=13)
        ax.set_ylabel('$P(z)$',fontsize=13)
        return fig,ax

    def set_expValues(self,z0,wave_lenght,L,I0):
        self.z0 = z0
        self.wave_lenght = wave_lenght
        self.L = L
        self.I0 = I0
    
    def curve(self,z,qo):
        return qo/(2*np.sqrt(2))*1/(1+z**2/self.z0**2)+1
    
    def plot_fitcurve(self):
        q0,dq0 = curve_fit(self.curve,self.positions,self.DT)
        self.q0 = q0[0]
        self.beta = self.q0/(self.I0*self.L)
        print(f'qo: {self.q0}')
        print(f'b: {self.beta}')
        fig,ax = plt.subplots(constrained_layout=True)
        ax.scatter(self.positions,self.DT,color='gray',marker='+',label='Experimental points')
        p = self.positions
        ax.plot(p,self.curve(p,self.q0),label='Fitting',color='red')
        ax.legend()
        title = f'$\Delta T(z)$\n{self.material}-{self.aperture},{self.power}W'
        ax.set_title(title,fontsize=14)
        ax.set_xlabel('$z$ $(m)$',fontsize=13)
        ax.set_ylabel('$\Delta T(z)$',fontsize=13)
        return fig,ax

class EXPERIMENT:
    def __init__(self,material,group,power):
        self.material = material
        self.group = group
        self.power = power
        self.close = CLOSE(material,group,'close',power)
        self.open = OPEN(material,group,'open',power)
        self.open.set_focalPlane(self.close.Zfocal_plane)
        self.close.set_S(self.open.T[-1])
    
    def set_experiment_values(self,wave_lenght,L_material):
        self.wave_lenght = wave_lenght
        self.L = L_material
        self.close.set_expValues(wave_lenght,L_material)
        self.open.set_expValues(self.close.z0,wave_lenght,L_material,self.close.I0)
    
    def get_values(self):
        return [self.material,self.group,self.power,self.close.S,self.close.z0,self.close.w0, self.close.I0, self.close.DPhi, self.close.n2]


    