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