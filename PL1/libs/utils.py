import numpy as np
import matplotlib.pyplot as plt
import re

class DATA:
    def __init__(self,material,group,aperture,power):
        PATH = 'Data/'
        file = f'zscan_{material}_300322_g{group}_{aperture}_{power}mW/'
        self.power = power*1e-3
        self.material = material
        self.aperture = aperture
        
        self.positions = np.load(open(PATH+file+file[:-1]+'_dl_pos.npy','rb'))
        self.positions = np.abs(self.positions-np.max(self.positions))*1e-3 #meters

        self.power_points = np.load(open(PATH+file+file[:-1]+'_pwr.npy','rb')) #Watts
        self.log_file = open(PATH+file+file[:-1]+'_log.csv','r').read()
        values = ['start_pos','end_pos','step','exposure_time','avg_Power','aperture','no_lens']
        self.log = dict(zip(values,[float(n) for n in re.findall(r'\d+| \d+\.\d+|\.\d+',self.log_file)]))

        self.T = self.power_points/self.power #Trasmitance
        self.DT = self.T/self.T[-1] #T/T(Z>>Z0)

        print(f'Material: {material} Aperture: {aperture} Power: {power}')
    
    def plot_power(self,title=None):
        title = f'Power Plot\n{self.material}-{self.aperture},{self.power}W' if title is None else title
        fig,ax = plt.subplots()
        ax.set_title(title,fontsize=16)
        ax.plot(self.positions,self.power_points)

        
    def plot_normalized_trasmittance(self,title=None):
        title = f'$\Delta T(Z)$\n{self.material}-{self.aperture},{self.power}W\n T(Z>>Z0) = {self.T[-1]}' if title is None else title
        fig,ax = plt.subplots()
        ax.set_title(title,fontsize=16)
        ax.plot(self.positions,self.DT)

class CLOSE(DATA):
    def __init__(self, material, group, aperture, power):
        super().__init__(material, group, aperture, power)

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

    def set_focalPlane(self,Zfolcal_plane):
        self.Zfocal_plane = Zfolcal_plane
        self.positions -= self.Zfocal_plane
    

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


    