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
        glycerolDen_RT = (1273.3-0.6121*23)/1000 			#Density of Glycerol (g/cm3)
        waterDen_RT = (1-(((abs(23-4))/622)**1.7)) 	#Density of water (g/cm3)

        if weight_percent!=0:
            glycerolVol=weight_percent/glycerolDen_RT
            waterVol=(1-weight_percent)/waterDen_RT

        if volume_percent!=0:
            glycerolVol=volume_percent
            waterVol=1-glycerolVol

    #Fraction cacluator ----------------

    glycerolMass=glycerolDen*glycerolVol
    waterMass=waterDen*waterVol
    totalMass=glycerolMass+waterMass
    mass_fraction=glycerolMass/totalMass
    vol_fraction= glycerolVol/(glycerolVol+waterVol)

    print ("Mass fraction of mixture =", mass_fraction)
    print ("Volume fraction of mixture =", vol_fraction)


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

def new_viscosity():
    """ Viscosity from Fig 8, w=0.5, in Gozalez et al, J. Chem. Eng. Data 2011, 56, 1397–1406 """
    values = np.array([[372.9775207588715, -0.04301949799556977],
            [363.0182366027012, 0.013839677824416358],
            [303.0807927478778, 0.6183220416236305],
            [293.09546354534103, 0.7770083877441973],
            [283.0565712105636, 0.9432164445010943],
            [273.01343054445266, 1.157962842639872],
            [268.05850991201703, 1.2737738365139235],
            [263.0880689937548, 1.384495895376037],
            [260.0399446329281, 1.5005816934266483],
            [253.12847528580392, 1.679480055324712],
            [250.01458718290178, 1.8149804312662385],
            [240.02407031949036, 2.2234595590917317]])


    values = np.array([
            [372.9775207588715, -0.04301949799556977],
            [363.0182366027012, 0.013839677824416358],
            [303.0807927478778, 0.6183220416236305],
            [293.09546354534103, 0.7770083877441973],
            [283.0565712105636, 0.9432164445010943],
            [273.01343054445266, 1.157962842639872],
            [268.05850991201703, 1.2737738365139235],
            [263.0880689937548, 1.384495895376037],
            [253.12847528580392, 1.679480055324712],
            [260.0399446329281, 1.5005816934266483],
            [250.01458718290178, 1.8149804312662385],
            [240.02407031949036, 2.2234595590917317],
            [312.7858174759948, 0.47708955402380016],
            [323.0665871676021, 0.36775458188723664],
            [332.7152521232539, 0.26273639560393397],
            [343.1277841313972, 0.1827599506490475],
            [352.92131932036307, 0.10738436825754659],
            [277.82320341439174, 1.0394885242387013],
            [245.0194781899867, 2.0102566163127116],
            [234.63165695329354, 2.464881350381542],
            [230.05713276240007, 2.7127374370077986],
            [226.0864907244321, 2.932569398328127],
            [222.9073725894211, 3.1336333302342934]])

    return values


def new_viscosity_60():
    """ Viscosity from Fig 8, w=0.6, in Gozalez et al, J. Chem. Eng. Data 2011, 56, 1397–1406 """

    values = np.array([
                [373.1349992875986, 0.15225264370243913],
                [363.2653573629451, 0.2078408964197217],
                [333.32987306986377, 0.448342526191962],
                [353.40180955227834, 0.2835972293703012],
                [312.300113849873, 0.6988087972722777],
                [302.2318337190154, 0.8465338676348124],
                [292.5631853943155, 1.0183524383199052],
                [282.8506258004293, 1.2155786332515424],
                [272.3970144683216, 1.4599341912338364],
                [262.4755164939777, 1.7489294259991335],
                [253.83148314389217, 2.056217000455638],
                [242.7624410735338, 2.495141024334294],
                [237.61105474146547, 2.753010758575824],
                [234.0482693317978, 2.9411685402933667],
                [229.5971266065617, 3.2017785508946175],
                [225.32289077615073, 3.4694743648639035],
                [248.4967780117918, 2.253124061534731],
    ])
    return values


def new_viscosity_70():
    """ Viscosity from Fig 8, w=0.7, in Gozalez et al, J. Chem. Eng. Data 2011, 56, 1397–1406 """

    values = np.array([
        [373.0846175867325, 0.3152230559644161],
        [363.06795597519294, 0.3816084033905761],
        [353.20958556910796, 0.4741107656264221],
        [343.33719095213047, 0.574898213058916],
        [333.1268526822417, 0.686219825316194],
        [323.25948917231085, 0.8106724693488814],
        [313.299574255668, 0.9612511937375755],
        [303.1408467588993, 1.1340598739157637],
        [293.0201250832786, 1.339521415730656],
        [282.9775532704489, 1.5880331766716502],
        [273.0577614325029, 1.8571674072113433],
        [265.60113671574834, 2.113031326218516],
        [261.1286921807556, 2.2761147355459093],
        [256.8943447867446, 2.4341280320582444],
        [249.5183016345202, 2.7498888397475993],
        [245.77848061877376, 2.9443242048146603],
        [241.62193719894424, 3.1545748275086654],
        [239.0591331952836, 3.3031187488186244],
        [236.52932486867041, 3.453524815756936],
        [234.15496726832507, 3.587256046448975],
        [231.53459901487255, 3.7595900216561824],
        [229.3253141330217, 3.917828308545388]
    ])
    return values