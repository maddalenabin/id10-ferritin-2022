import h5py
import numpy as np


def load_data(file, key):
    """Load dataset of a h5 file and return the entry as np.array
    """

    with h5py.File(file, 'r') as f:
        data = f[key]
        data = np.asarray(data)
    
    return data


def visit_func(name, node):
    '''Return all groups and datasets name and shapes of h5 file called name
    Use it as in the following:
    with h5py.File(filename, 'a') as f: f.visititems(visit_func)
    '''
    if isinstance(node, h5py.Group):
        print("group", node.name)
    elif isinstance(node, h5py.Dataset):
        if (node.dtype == 'object') :
            print (node.name, 'is an object Dataset')
        else:
            print('\t', node.name, node.shape)
    else:
        print(node.name, 'is an unknown type')

# def elog_entries():
#     entries = ['elog/Absorbers 100um/', 'elog/Absorbers 80um/', 'elog/Exposure time (s)/', 'elog/LN flow meter/', 'elog/Mesh (spots x lines)/', 'elog/No of spots/', 'elog/Sample/', 'elog/Sample no. (label)/', 'elog/Short comment/', 'elog/Temperature, C/', 'elog/Temperature, K/', 'elog/Total exposure time (s)/', 'elog/comment/', 'elog/measurement folder/', 'elog/number of frames/', 'elog/position y, mm/', 'elog/position z, mm/', 'elog/reps per spot/', 'elog/scan number/', 'elog/transmission (%)/']    return entries
#     return entries

# def saxs_entries():
#     entries = ['saxs/I/', 'saxs/I_reps/', 'saxs/q/', 'saxs/scans/']
#     return entries

# def xpcs_entries():
#     g2s = ['xpcs/g2s/delay/', 'xpcs/g2s/g2s/', 'xpcs/g2s/q/']
#     ttcs = ['xpcs/ttcs/q/', 'xpcs/ttcs/times/', 'xpcs/ttcs/ttc_avg_int/', 'xpcs/ttcs/ttc_rep_qs_avg/']
#     return g2s + ttcs