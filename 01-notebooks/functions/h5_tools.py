import h5py
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

def load_data(file, key):
    """Load dataset of a h5 file and return the entry as np.array
    """

    with h5py.File(file, 'r') as f:
        data = f[key]
        data = np.asarray(data)
            
    return data


def visit_func(name, node):
    """Return all groups and datasets name and shapes of h5 file called name
    Use it as in the following:
    with h5py.File(filename, 'a') as f: f.visititems(visit_func)
    """
    if isinstance(node, h5py.Group):
        print("group", node.name)
    elif isinstance(node, h5py.Dataset):
        if (node.dtype == 'object') :
            print (node.name, 'is an object Dataset')
        else:
            print('\t', node.name, node.shape)
    else:
        print(node.name, 'is an unknown type')


def elog_selected_entries(filename, entries=["measurement folder","Temperature, K","transmission (%)","Absorbers 80um","Absorbers 100um","Short comment","comment","Exposure time (s)","number of frames",]):
    """Return the elog entries given, taken from the google drive sheet of the beamtime
    """
    with h5py.File(filename, 'r') as f: 
        for entry in entries:
            print("{0:20} {1}".format(entry,  f[f'elog/{entry}/'].asstr()[()]) )
            
def elog_selected_entries_dict(filename, entries=["measurement folder","Temperature, K","transmission (%)","Absorbers 80um","Absorbers 100um","Short comment","comment","Exposure time (s)","number of frames",]):
    """Return a dictionary with the elog entries given, taken from the google drive sheet of the beamtime
    """
    d = {}
    with h5py.File(filename, 'r') as f: 
        for entry in entries:
            d[entry] = f[f'elog/{entry}/'].asstr()[()]
    return d

def plot_ttc(filename, rep=0, vs=(0,0.15)):
    """Plots the twotime correlation function for different momentum transfers"""
    ttc_qv = load_data(filename, '/xpcs/ttcs/q')
    time = load_data(filename, '/xpcs/ttcs/times')
    
    ttcs = np.empty(shape=(len(ttc_qv), len(time), len(time)))
    with h5py.File(filename, 'r') as f:
        ttcs = np.asarray(f['/xpcs/ttcs/ttc_rep_qs_avg'])[rep,:,:,:]
        
    f, axs = plt.subplots(2,5,figsize=(10,4), tight_layout=True)
    for i in range(len(ttc_qv)):
        ax = axs.ravel()[i]
        norm = np.mean(ttcs[i,:300,-300:])
        ax.imshow(ttcs[i,:,:]-norm, origin='lower', cmap='turbo', vmin=vs[0], vmax=vs[1], extent=(time[0],time[-1])*2 )
        ax.set_xlabel('$t_1 (s)$')
        ax.set_title(f'{ttc_qv[i]:.2f} 1/nm', fontsize=8)
        ax.set_ylabel('$t_2 (s)$')
        divider = make_axes_locatable(ax)
        cb1 = mpl.colorbar.ColorbarBase(divider.append_axes('right', size='5%', pad=0.05), cmap=mpl.cm.turbo,  orientation='vertical', norm=mpl.colors.Normalize(vmin=0,vmax=.15));
