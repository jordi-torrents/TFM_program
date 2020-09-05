import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

directory='results/'

subdirectories = os.listdir(directory)
fig_pol, axs_pol = plt.subplots(2)
fig_GNF, axs_GNF = plt.subplots(1)
fig_BaC, axs_BaC = plt.subplots(1)
inset_ax = axs_GNF.inset_axes([0.62, 0.1, 0.34, 0.25])

if 'VicsekM' in subdirectories:
    measurements = os.listdir(directory+'VicsekM')
    if 'polarization' in measurements:
        current_dir=directory+'VicsekM/polarization/'
        data = np.empty((0,3))
        for file in os.listdir(current_dir):
            file_data = pd.read_csv(current_dir+file, sep='\s+',comment='#', header=None).values
            data = np.append(data,[[float(file[-6:-4])/100, np.mean(file_data), np.var(file_data)]], axis=0)
        axs_pol[0].plot(data[:,0], data[:,1],'.-', label='VicsekM', lw=0.8)
        axs_pol[1].plot(data[:,0], data[:,2],'.-', lw=0.5)
    if 'GNF' in measurements:
        current_dir=directory+'VicsekM/GNF/'
        for file in os.listdir(current_dir):
            data = pd.read_csv(current_dir+file, sep='\s+',comment='#', header=None).values.T
            axs_GNF.plot(np.mean(data, axis=1)[5:-70],np.std(data, axis=1)[5:-70],'.-', label='VicsekM '+'$\eta=0$.'+file[3:5], ms=3, lw=1)
            inset_ax.plot(np.mean(data, axis=1),np.std(data, axis=1),'-', label='VicsekM '+'$\eta=0$.'+file[3:5], ms=3, lw=1)


if '_LevyM_' in subdirectories:
    measurements = os.listdir(directory+'_LevyM_')
    if 'polarization' in measurements:
        current_dir=directory+'_LevyM_/polarization/'
        data = np.empty((0,3))
        for file in os.listdir(current_dir):
            file_data = pd.read_csv(current_dir+file, sep='\s+',comment='#', header=None).values
            data = np.append(data,[[float(file[-6:-4])/100, np.mean(file_data), np.var(file_data)]], axis=0)
        axs_pol[0].plot(data[:,0], data[:,1],'.-', label='_LevyM_', lw=0.8)
        axs_pol[1].plot(data[:,0], data[:,2],'.-', lw=0.5)

if 'bur-coa' in subdirectories:
    measurements = os.listdir(directory+'bur-coa')
    if 'speed' in measurements and 'polarization' in measurements:
        current_dir = directory+'_LevyM_/'
        noise_value = np.random.choice(os.listdir(current_dir+'speed/'))
        print('choosen noise value for bur-coa figure:',noise_value)
        speed, std_speed = pd.read_csv(current_dir+'speed/'+noise_value, sep='\s+', comment='#', header=None).values.T
        polarization = pd.read_csv(current_dir+'polarization/'+noise_value, comment='#', header=None).values
        tau_max=200

        norm_polarization=polarization-np.mean(polarization)
        norm_speed=speed-np.mean(speed)
        norm=np.sqrt(np.sum(norm_polarization[tau_max:-tau_max]**2))*np.sqrt(np.sum(norm_speed[tau_max:-tau_max]**2))

        taus=np.array(range(-100,tau_max))
        adjust=np.empty(len(taus))
        for j in range(len(taus)):
            tau=taus[j]
            adjust[j]=np.sum(norm_polarization[tau_max+tau:-tau_max+tau]*norm_speed[tau_max:-tau_max])/norm
        tau = taus[np.argmax(adjust)]

        print('from NCC={:.3f} to {:.3f} with tau={}'.format(np.sum(norm_polarization[tau_max:-tau_max]*norm_speed[tau_max:-tau_max])/norm,np.max(adjust),tau))
        print('mea_speed = {:.3f} (s.d. {:.4f}) and mean_polarization={:.3f} (s.d. {:.4f})'.format(np.mean(speed),np.std(speed),np.mean(polarization),np.std(polarization)))

        norm_speed=norm_speed[:-tau]
        norm_polarization=norm_polarization[tau:]

        axs_BaC.plot(norm_speed)
        axs_BaC.plot(norm_polarization)






inset_ax.set(xscale='log',yscale='log', xticks=(10,1000))
axs_GNF.set(xscale='log', yscale='log', title=directory)
axs_GNF.legend()
fig_GNF.tight_layout()
fig_GNF.savefig(directory+'GNF.png', dpi=300)


axs_pol[0].legend()
axs_pol[0].set(ylabel=r'$\varphi$', title=directory)
axs_pol[1].set(xlabel=r'$\eta$', ylabel=r'$\sigma^2$')
fig_pol.tight_layout()
fig_pol.savefig(directory+'polarization.png', dpi=300)


fig_BaC.tight_layout()
fig_BaC.savefig(directory+'BaC.png', dpi=300)
