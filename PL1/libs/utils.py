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
        self.I = self.power_points/A*self.log['exposure_time']
        self.I0 = self.power/A*self.log['exposure_time']

        print(f'Material: {material} Aperture: {aperture} Power: {power}')
        print(f'\tI0 = {self.I0}')

    def plot_power(self,title=None):
        title = f'{self.material}-{self.aperture},{self.power}mW' if title is None else title
        fig,ax = plt.subplots()
        ax.set_title(title,fontsize=16)
        ax.plot(self.positions,self.power_points)
        
