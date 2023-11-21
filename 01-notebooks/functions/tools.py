import numpy as np
import pandas as pd
import pickle


def D_coeff(T, eta, Rh=6.6e-9):
    """Calculate the diffusion coefficient in nm^2/us from Stokes-Einsten equation
    Args
        T: temperature in K
        eta: viscosity of the solvent in N*s/m^2 
        Rh: hydrodynamic radius in m. Default value is for ferritin
    """
    if T == 243:
        T = 243
    kb = 1.380649e-23 # m2 kg s-2 K-1
    D = kb*T / (6*np.pi*eta*Rh) # m^2 / s
    
    return D#*1e12

def viscosity():
    """Viscosity values in Ns/m2 of water/glycerol mixture with 50w%, 50v%, 60v% glycerol.
       http://www.met.reading.ac.uk/~sws04cdw/viscosity_calc.html
    """
    eta_50w_gly = {300: 0.0046742, 290: 0.0067715, 280: 0.010448, 270: 0.017490, 260: 0.032565, 250: 0.069890, 240: 0.18269, 230: 0.63521, 220: 3.4169, 210: 37.846}
    eta_50v_gly = {300: 0.0063573, 290: 0.0094968, 280: 0.015197, 270: 0.026578, 260: 0.052163, 250: 0.11934, 240: 0.33719, 230: 1.2893, 220: 7.7852, 210: 98.826}
    eta_60v_gly = {300: 0.011691, 290: 0.018606, 280: 0.032089, 270: 0.061359, 260: 0.13400, 250: 0.34857, 240: 1.1499, 230: 5.3002, 220: 40.015, 210: 657.37}
    # eta c1: 51v%, 57w% if c0=6.5mg/ml (after adding gly), c1=31, c2=44, c4=64 mg/ml
    eta_c1 = {300: 0.0067262, 290: 0.010106, 280: 0.016282, 270: 0.028708, 260: 0.056900, 250: 0.13172, 240: 0.37757, 230: 1.4692, 220: 9.0605, 210: 117.89}
    # eta c2: 54v%, 60w%
    eta_c2 = {300: 0.0080120, 290: 0.012256, 280: 0.020171, 270: 0.036484, 260: 0.074554, 250: 0.17907, 240: 0.53670, 230: 2.2044, 220: 14.508, 210: 203.60}
    # eta c3: 59v%, 65w%
    eta_c3 = {300: 0.010947, 290: 0.017301, 280: 0.029595, 270: 0.056042, 260: 0.12098, 250: 0.31037, 240: 1.0070, 230: 4.5503, 220: 33.555, 210: 536.70}

    etas = pd.DataFrame()
    etas['temp'] = list(eta_50w_gly.keys())

    dbs = [eta_50w_gly, eta_50v_gly, eta_60v_gly, eta_c1, eta_c2, eta_c3]
    labels = ['50w%', '50v%', '60v%', 'c1', 'c2', 'c3']

    for db,label in zip(dbs, labels):
        if all( etas['temp'].values == np.array(list(db.keys())) ):
            etas[label] = list(db.values())
        else: 
            print(label, " didn't work")

    return etas


def update_D_coeff(T, Tr, D, dD, c, run):
    """Updates the diffusion coefficient calculated from the data
    """
    DiffC = pd.read_pickle("/cfs/home/mabi3848/id10-ferritin-2022/03-source/diffusion_coefficient")
    exists = DiffC['run'].str.contains(run)
    
    if exists:
        pass
    else:
        new_row = {'temperature': T, 'transmission': Tr, 'D': D, 'dD': dD, 'c': c, 'run': run}
        DiffC2 = DiffC.append(new_row, ignore_index=True)
        DiffC2.to_pickle("/cfs/home/mabi3848/id10-ferritin-2022/03-source/diffusion_coefficient")

    return


def load_npz_data(filename):
    """Loads data from npz file and return a dictionary
    """
    file = np.load(filename)
    data = {key: file[key] for key in file.files}

    return data