import numpy as np
import matplotlib.pyplot as plt
import re

class DATA:
    def __init__(self,material,group,aperture,power,ra):
        PATH = 'Data/'
        file = f'zscan_{material}_300322_g{group}_{aperture}_{power}mW/'
        self.power = power*1e-3
        self.material = material
        self.aperture = aperture
        self.ra = ra
        self.positions = np.load(open(PATH+file+file[:-1]+'_dl_pos.npy','rb'))
        self.positions = np.abs(self.positions-np.max(self.positions))*1e-3

        self.power_points = np.load(open(PATH+file+file[:-1]+'_pwr.npy','rb'))
        self.log_file = open(PATH+file+file[:-1]+'_log.csv','r').read()
        values = ['start_pos','end_pos','step','exposure_time','avg_Power','aperture','no_lens']
        self.log = dict(zip(values,[float(n) for n in re.findall(r'\d+| \d+\.\d+|\.\d+',self.log_file)]))
        A = (ra)**2*np.pi
        self.I = self.power_points/A
        self.I0 = self.power/A

        self.T = self.I/self.I0 #Trasmitance
        self.DT = self.T/self.T[-1] #T/T(Z>>Z0)

        print(f'Material: {material} Aperture: {aperture} Power: {power}')
        print(f'\tI0 = {self.I0}')
    
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
    def __init__(self, material, group, aperture, power, ra):
        super().__init__(material, group, aperture, power, ra)

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
        self.I0_focalPlane = self.I[self.iZfocal_plane]
        self.positions -= self.Zfocal_plane

        self.z0 = np.abs(self.positions[p_max]-self.positions[p_min])/1.7

        print(f'Position of the focal plane: {self.Zfocal_plane}, {self.positions[self.iZfocal_plane]}')
        print(f'I0 focal plane: {self.I0_focalPlane}')
        print(f'Rayleigh range: {self.z0}')

    def set_S(self,Topen):
        self.S = self.T[-1]/Topen
        print(f'S: {self.S}')
    
    def set_expValues(self,wave_lenght,L_material):
        self.wave_lenght = wave_lenght
        self.L = L_material

    def calc_n2(self):
        self.DPhi = self.DTpv/(0.406*(1-self.S)**0.27)
        self.n2 = (self.wave_lenght/(2*np.pi))*self.DPhi/(self.I0_focalPlane*self.L)
        print(f'DPhi: {self.DPhi}')
        print(f'n2: {self.n2}')


class OPEN(DATA):
    def __init__(self, material, group, aperture, power, ra):
        super().__init__(material, group, aperture, power, ra)

    def set_focalPlane(self,Zfolcal_plane):
        self.Zfocal_plane = Zfolcal_plane
        self.positions -= self.Zfocal_plane
    

class EXPERIMENT:
    def __init__(self,material,group,power,r=[2e-3,12e-3]):
        self.material = material
        self.group = group
        self.power = power
        self.close = CLOSE(material,group,'close',power,r[0])
        self.open = OPEN(material,group,'open',power,r[1])
        self.open.set_focalPlane(self.close.Zfocal_plane)
        self.close.set_S(self.open.T[-1])
    
    def set_experiment_values(self,wave_lenght,L_material):
        self.wave_lenght = wave_lenght
        self.L = L_material
        self.close.set_expValues(wave_lenght,L_material)


    