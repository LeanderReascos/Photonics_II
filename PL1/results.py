from sympy import arg
from libs.utils import EXPERIMENT
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

wave_lenght = 780e-9
L = {'CG': 1e-3,  'WINDOW': 1e-2}
path = 'graphics/'

if __name__ == '__main__':
    args = sys.argv[1:]
    option = args[0]
    args = args[1:]
    material = args[0]
    power = int(args[1])
    groups = [int(c) for c in args[2]]
    if option.upper() == 'N2':
        results = []
        for g in groups:
            print(f'\n*** MATERIAL: {material} GROUP: {g} ***')
            EX = EXPERIMENT(material.lower(),g,power)
            EX.set_experiment_values(wave_lenght,L[material.upper()])
            EX.close.calc_n2()
            results.append(EX.get_values())
            figc,axc = EX.close.plot_normalized_trasmittance()
            figo,axo = EX.open.plot_normalized_trasmittance()
            figc.savefig(path+f'g{g}_{power}mW_{material}_close.pdf',transparent=True)
            figo.savefig(path+f'g{g}_{power}mW_{material}_open.pdf',transparent=True)

        df = pd.DataFrame(results,columns=['material','group','power','S','z0', 'w0','I0', 'DPhi', 'n2'])
        df.to_csv(f'results/{material}_{power}.csv')
    
    if option.upper() == 'NOISY':
        for g in groups:
            print(f'\n*** MATERIAL: {material} GROUP: {g} ***')
            EX = EXPERIMENT(material.lower(),g,power)
            fig,ax = EX.open.plot_noisy_filtered_points()
            plt.show()
        
    if option.upper() == 'BETA':
        results = []
        for g in groups:
            print(f'\n*** MATERIAL: {material} GROUP: {g} ***')
            EX = EXPERIMENT(material.lower(),g,power)
            EX.set_experiment_values(wave_lenght,L[material.upper()])
            fig,ax = EX.open.plot_fitcurve()
            fig.savefig(path+f'g{g}_{power}mW_{material}_open_beta.pdf',transparent=True)
            results.append([EX.open.q0,EX.open.beta])

        df = pd.DataFrame(results,columns=['q0','beta'])
        df.to_csv(f'results/{material}_{power}_beta.csv')