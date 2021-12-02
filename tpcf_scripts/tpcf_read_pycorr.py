import numpy as np
import matplotlib.pyplot as plt
import os
import astropy
from astropy.io import fits
import pycorr
import argparse


def run_pycorr():
    return None

def write_data():
    return None



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description= 'Runs TPCF analysis')
    parser.add_argument('--galaxy',help='directory of galaxy mocks')
    parser.add_argument('--void',help ='directory of void catalog')
    parser.add_argument('--galaxy-randoms',help='directory for galaxy randoms', default ='')
    parser.add_argument('--void_randoms',help = 'directory for void randoms', default ='')

    args =parser.parse_args()

    galaxy_path = args.galaxy
    void_path = args.void

    galaxy_rand_path =args.galaxy_randoms
    void_rand_path = args.void_randoms

    if len(galaxy_rand_path) > 0 and len(void_rand_path)< 1: 
        mode =2  #  IS THERE A NAME FOR THIS CASE?
    elif len(void_rand_path)>0 and len(galaxy_rand_path)< 1:
        mode =3 # PARTIAL AUTOCORR?

    elif len(void_rand_path)>0 and len(galaxy_rand_path)>0:
        mode = 4 # AUTOCORR
    else:
        mode =1 # PERIODIC BOX

    gal_fname = 'MOD_frac__0.005_cutsky_ELG_z1.100_AbacusSummit_base_c000_ph000.fits'
    void_fname = 'CutSkyvoxel-Voids_cat.txt'
    gal_rand_fname = 'MOD_frac__0.005_cutsky_ELG_random_S1000_1X.fits'
    void_rand_fname =''

    # OPENING THE FILES
    with fits.open(galaxy_path+gal_fname) as hdul:
        gal_data= hdul[1].data
    with fits.open(void_path+void_fname) as hdul:
        void_path =hdul[1].data
    if mode == 2:
        with fits.open(galaxy_rand_path+gal_rand_fname) as hdul:
            gal_rand_data =hdul[1].data
    if mode == 3:
        with fits.open(void_rand_path+void_rand_fname) as hdul:
            void_rand_data =hdul[1].data
    
    if mode == 4:
        with fits.open()

    # CURRENTLY ASSUMING ABACUS_FIRST_GEN_MOCKS





        
