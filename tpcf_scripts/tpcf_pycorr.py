import os
import tempfile
import numpy as np
from matplotlib import pyplot as plt
#from pancakes.utilities.cosmology import Cosmology
from pycorr import (TwoPointCorrelationFunction,
    LandySzalayTwoPointEstimator, setup_logging,
    )

# To activate logging
setup_logging()

#cosmo = Cosmology(omega_m=0.301)

load_previous_result = False

zmin = 0.4
zmax = 1.25

#if load_previous_result:
    #prev = LandySzalayTwoPointEstimator.load(
  #      '/cosma/home/analyse/epaillas/data/ds_boss/patchy/'
   #     f'v2/pycorr/Patchy_NGC_0001_z{zmin}-{zmax}_xi_s.npy'
#    )
 #   R1R2 = prev.R1R2

# read randoms
#fname = f"/cosma/home/analyse/epaillas/data/ds_boss/patchy/NGC/npy/"\
 #     f"Patchy_Randoms_x50_NGC_z{zmin}-{zmax}_sky.npy"

#data = np.load(fname)
#data_x = data[:, 0]
#data_y = data[:, 1]
#data_z = data[:, 2]
#data_dis = cosmo.ComovingDistance(data_z)
#randoms_positions1 = [data_x, data_y, data_z]
#randoms_weights1 = (data[:, 5] * data[:, 6])/(1 + 10000 * data[:, 3])


for nmock in range(1, 2):
    fname = f"/project/rrg-wperciva/tsfraser/Revolver/"\
            f"refAbacus_c000_ph006_z0.6-1.05_ELGs_mock.dat"

    fname_void_zbov = "/project/rrg-wperciva/tsfraser/Revolver/zobov-Voids_cat.txt"
    data = np.genfromtxt(fname_void_zbov)
    print(data[:,0].min(),data[:,1].max())
    data_x = data[:, 1]
    data_y = data[:, 2]
    data_z = data[:, 3]

    data_mocks = np.genfromtxt(fname)
    #data_dis = cosmo.ComovingDistance(data_z)
    data_positions1 = [data_x, data_y, data_z]
    #data_weights1 = #data[:, 3]

    edges = np.linspace(1e-9, 150, 151)


    data_mock_x = data_mocks[:,0]+1000
    data_mock_y = data_mocks[:,1]+1000
    data_mock_z = data_mocks[:,2]+1000

    data_positions2 = [data_mock_x,data_mock_y,data_mock_z]

    if not load_previous_result:
        if nmock == 1:
            R1R2 = None
        else:
            R1R2 = result.R1R2

    result = TwoPointCorrelationFunction(\
        's', edges, data_positions1=data_positions1,data_positions2 = data_positions2,\
        engine='corrfunc', nthreads=56,boxsize = 2000)
    print(type(result))
    output_fname = f"/project/rrg-wperciva/abacus_share/"\
        f"test_abacus_periodic_ELG.npy"
    print(result.name)
    result.save(output_fname)

