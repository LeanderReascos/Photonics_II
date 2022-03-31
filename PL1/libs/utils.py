import numpy as np
import matplotlib.pyplot as plt
import re

class EXPERIMENT:
    def __init__(self,material,group,type,power):
        PATH = '../Data/'
        file = f'zscan_{material}_300322_g{group}_{type}_{power}mW/'
        self.power = power
        self.positions = np.load(open(PATH+file+file[:-1]+'_dl_pos.npy','rb'))
        self.power_points = np.load(open(PATH+file+file[:-1]+'_pwr.npy','rb'))
        self.log_file = open(PATH+file+file[:-1]+'_log.csv','r').read()
        values = ['start_pos','end_pos','step','exposure_time','avg_Power','aperture','no_lens']
        self.log = dict(zip(values,[float(n) for n in re.findall(r'\d+| \d+\.\d+|\.\d+',self.log_file)]))

    def plot_power(self):
        fig,ax = plt.subplots()
        ax.plot(self.positions,self.power_points)
        