import numpy as np
import matplotlib.pyplot as plt
import os
import astropy
from astropy.cosmology import Planck13 as cosmo
from astropy.io import fits
from pycorr (TwoPointCorrelationFunction,
    LandySzalayTwoPointEstimator, setup_logging,
    )
import argparse

def generate_positions(data,list_of_cols):
    '''
    GIVEN FITS DATA, RETURNS APPROPRIATE POSITIONS

    INPUTS: data -- fits data object
            list_of_cols -- length 3 list, could be RA,DEC,Z or X,Y,Z for sims

    OUTPUTS: positions -- 3 column object of positions.
    
    '''
    data_1 = data[list_of_cols[0]]
    data_2 = data[list_of_cols[1]]
    data_3 = data[list_of_cols[2]]

    positions = [data_1,data_2,data_3] 

    return positions

def compute_edges():
    return 1.

def run_pycorr(mode,positions,boxsize,edges):

    data_positions1 = positions[0]
    data_positions2 = positions[1]
    rand_positions1 = None
    rand_positions2 =None
    if mode ==2 :
        rand_positions1 = positions[2]
        rand_positions2 = None
    elif mode ==3:
        rand_positions1 = positions[2]
        rand_positions2= None
    elif mode == 4:
        rand_positions1 = positions[2]
        rand_positions2 = positions[3]

    
    results =TwoPointCorrelationFunction('s', edges, data_positions1=data_positions1,\
         data_positions2 = data_positions2,\
        random_positions1= rand_positions1,random_positions2 =rand_positions2,\
         nthreads=56,boxsize = boxsize)


    return(results)

def write_data(result,galaxy_fname):
    result.save(galaxy_fname+'TPCF.npy')
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
        void_data =hdul[1].data
    if mode == 2:
        with fits.open(galaxy_rand_path+gal_rand_fname) as hdul:
            gal_rand_data =hdul[1].data
    if mode == 3:
        with fits.open(void_rand_path+void_rand_fname) as hdul:
            void_rand_data =hdul[1].data
    
    if mode == 4:
        with fits.open(void_rand_path+void_rand_fname) as hdul1:
            void_rand_data =hdul1[1].data

        with fits.open(galaxy_rand_path+gal_rand_fname) as hdul2:
            galaxy_rand_data= hdul2[1].data

    # CURRENTLY ASSUMING ABACUS_FIRST_GEN_MOCKS


    boxsize =2000
    gal_pos = generate_positions(gal_data,['RA','DEC','Z'])
    void_pos = generate_positions(void_data,['RA','DEC','Z'])
    pos_list = []
    pos_list.append(gal_pos)
    pos_list.append(void_pos)
    if mode == 2:
        gal_rand_pos = generate_positions(gal_rand_data,['RA','DEC','Z'])
        pos_list.append(gal_rand_pos)
    if mode == 3:
        void_rand_pos = generate_positions(void_rand_data,['RA','DEC','Z'])
        pos_list.append(void_rand_pos)
    if mode == 4:

        gal_rand_pos = generate_positions(gal_rand_data,['RA','DEC','Z'])
        void_rand_pos = generate_positions(void_rand_data,['RA','DEC','Z'])
        pos_list.append(gal_rand_pos,void_rand_pos)
        
    edges = compute_edges()
    results = run_pycorr(mode,pos_list,boxsize,edges)
    write_data(results,gal_fname)