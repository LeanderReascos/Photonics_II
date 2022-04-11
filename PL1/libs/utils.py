import numpy as np
import matplotlib.pyplot as plt
import re

class EXPERIMENT:
    def __init__(self,material,group,aperture,power):
        PATH = 'Data/'
        file = f'zscan_{material}_300322_g{group}_{aperture}_{power}mW/'
        self.power = power*1e-3
        self.material = material
        self.aperture = aperture
        self.positions = np.load(open(PATH+file+file[:-1]+'_dl_pos.npy','rb'))

        self.power_points = np.load(open(PATH+file+file[:-1]+'_pwr.npy','rb'))
        self.log_file = open(PATH+file+file[:-1]+'_log.csv','r').read()
        values = ['start_pos','end_pos','step','exposure_time','avg_Power','aperture','no_lens']
        self.log = dict(zip(values,[float(n) for n in re.findall(r'\d+| \d+\.\d+|\.\d+',self.log_file)]))
        A = (self.log['aperture']*1e-3)**2*np.pi
        self.I = self.power_points/A
        self.I0 = self.power/A

        self.T = self.I/self.I0 #Trasmitance
        self.DT = self.T/self.T[-1] #T/T(Z>>Z0)

        print(f'Material: {material} Aperture: {aperture} Power: {power}')
        print(f'\tI0 = {self.I0}')

    def correct_positions(self):
        #Z=0 := Max distance
        self.positions = np.abs(self.positions-np.max(self.positions))
        p_min = np.min([self.iTv,self.iTp])
        p_max = np.max([self.iTv,self.iTp])
        self.Zfocal_plane = self.positions[p_min]+(self.positions[p_max]-self.positions[p_min])/2
        self.positions -= self.Zfocal_plane

        self.z0 = np.abs(self.positions[p_max]-self.positions[p_min])/1.7

        print(f'Position of the focal plane: {self.Zfocal_plane}')
        print(f'Rayleigh range: {self.z0}')


    def get_DTpv(self):
        self.Tp = np.max(self.DT)
        self.Tv = np.min(self.DT)
        self.iTp = np.where(self.DT == self.Tp)[0][0]
        self.iTv = np.where(self.DT == self.Tv)[0][0]
        print(f'Peak: {self.Tp}\nValley: {self.Tv}')
        return np.abs(self.Tp-self.Tv)
    

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
    