import numpy as np
import pandas as pd
import pickle


def D_coeff(T, eta, Rh=7.3e-9):
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


def viscosity_glywater(T,weight_percent=0,volume_percent=0,waterVol=0,glycerolVol=0):

    #Variables ----------------
    #T				#temperature (degrees Celcius)
    #waterVol 		#volume of water required (ml)
    #glycerolVol 	#volume of Glycerol used (ml)

    #Densities ----------------
    glycerolDen = (1273.3-0.6121*T)/1000 			#Density of Glycerol (g/cm3)
    waterDen = (1-(((abs(T-4))/622)**1.7)) 	#Density of water (g/cm3)

    if (waterVol==0)|(glycerolVol==0):

        if weight_percent!=0:
            glycerolVol=weight_percent/glycerolDen
            waterVol=1-glycerolVol

        if volume_percent!=0:
            glycerolVol=volume_percent
            waterVol=1-glycerolVol

    #Fraction cacluator ----------------

    glycerolMass=glycerolDen*glycerolVol
    waterMass=waterDen*waterVol
    totalMass=glycerolMass+waterMass
    mass_fraction=glycerolMass/totalMass
    vol_fraction= glycerolVol/(glycerolVol+waterVol)

    #print ("Mass fraction of mixture =", round(mass_fraction,5))
    #print ("Volume fraction of mixture =", round(vol_fraction,5))


    #Density calculator ----------------

    ##Andreas Volk polynomial method
    contraction_av = 1-(3.520E-8*((mass_fraction*100)))**3+(1.027E-6*((mass_fraction*100)))**2+2.5E-4*(mass_fraction*100)-1.691E-4
    contraction = 1+contraction_av/100

    ## Distorted sine approximation method
    #contraction_pc = 1.1*math.pow(math.sin(numpy.radians(math.pow(mass_fraction,1.3)*180)),0.85)
    #contraction = 1 + contraction_pc/100

    density_mix=(glycerolDen*vol_fraction+waterDen*(1-vol_fraction))*contraction

    #print ("Density of mixture =",round(density_mix,5),"g/cm3")


    #Viscosity calcualtor ----------------

    glycerolVisc=0.001*12100*np.exp((-1233+T)*T/(9900+70*T))
    waterVisc=0.001*1.790*np.exp((-1230-T)*T/(36100+360*T))

    a=0.705-0.0017*T
    b=(4.9+0.036*T)*np.power(a,2.5)
    alpha=1-mass_fraction+(a*b*mass_fraction*(1-mass_fraction))/(a*mass_fraction+b*(1-mass_fraction))
    A=np.log(waterVisc/glycerolVisc)

    viscosity_mix=glycerolVisc*np.exp(A*alpha)

    #print ("Viscosity of mixture =",round(viscosity_mix,5), "Ns/m2")
    return viscosity_mix
        ############################################################################

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


def conc_v_w_mol():
    """ Values of protein concentration, and gly fraction in v%, w% and mol%"""
    print("""
    c (mg/ml) 	 v% 	 w% 	 mol% 	 
    102         48.9 	 54.8 	 19.2
    135 	50.6 	 56.4 	 20.2
    314 	60.5 	 65.9 	 27.5
    SU13        55.0     60.7    23.2
    DLS         44.2     50.0    16.4
    """)
    return 