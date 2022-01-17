
from astropy.io import fits
import numpy as np
from pancakes.pancakes.utilities.cosmology import Cosmology
from pycorr import (TwoPointCorrelationFunction,
    LandySzalayTwoPointEstimator, setup_logging,
    TwoPointCounter,AnalyticTwoPointCounter
    )
import time
import matplotlib.pyplot as plt
cosmo = Cosmology(omega_m=0.315)
load_previous_result = False

# To activate logging
#setup_logging()

# read shifted random positions and weights
#fname = 'boss_randoms_rec_IFFTP_RecSym.fits'
#with fits.open(fname) as hdul:
#    data = hdul[1].data
#data_ra = data['RA_REC'].astype(np.float64)
#data_dec = data['DEC_REC'].astype(np.float64)
#data_z = data['Z_REC'].astype(np.float64)
#data_dis = cosmo.ComovingDistance(data_z)
#srandoms_positions = [data_ra, data_dec, data_dis]
#srandoms_weights = np.ones(data_ra.shape[0])#data['WEIGHT_FKP'].astype(np.float64)

def txt_extract(fname,cols,skip_lines,type,shift):
    data = np.genfromtxt(fname,skip_header =skip_lines)

    if type == 'box': # USES X,Y,Z COORDS

        data_x = data[:,cols[0]].astype(np.float64)
        data_y = data[:,cols[1]].astype(np.float64)
        data_z =data[:,cols[2]].astype(np.float64)

        data_positions = [data_x,data_y,data_z]
        data_w = np.ones(data_x.shape[0])

        return data_positions,data_w

    else: # Real sim, weights are what you set them as, using RA,DEC, Z
        data_ra = data[:,cols[0]].astype(np.float64)
        data_dec = data[:,cols[1]].astype(np.float64)
        data_z =data[:,cols[2]].astype(np.float64) 

        if shift == True:

            data_dis = cosmo.ComovingDistance(data_z)

            data_positions = [data_ra,data_dec,data_dis]

        else:
            data_positions = [data_ra,data_dec,data_z]
        data_w = np.ones(data_ra.shape[0])

        return data_positions, data_w

def fits_extract(fname,cols,type,shift,zmin,zmax):
    with fits.open(fname,ignore_missing_end=True) as hdul:
        data = hdul[1].data

    if type == 'box':

        data_x = data[cols[0]]#.astype(np.float64)
        data_y = data[cols[1]]#.astype(np.float64)
        data_z = data[cols[2]]#.astype(np.float64)
#
        data_positions = [data_x,data_y,data_z]
        data_w = np.ones(data_x.shape[0])
        fits.close()
        return data_positions,data_w

    else:

        
        data_ra = data[cols[0]].astype(np.float64)
        data_dec = data[cols[1]].astype(np.float64)
        data_z = data[cols[2]].astype(np.float64)

        data_ra = data_ra[(data_z>zmin)&(data_z<zmax)]
        data_dec = data_dec[(data_z>zmin)&(data_z<zmax)]
        data_z =  data_z[(data_z>zmin)&(data_z<zmax)]

        print(data_ra.shape)
        print(fname,'file')
        print(data_ra[0],'input?')
        if shift == True:


            data_dis = cosmo.ComovingDistance(data_z) # Double check this


            data_positions = [data_ra,data_dec,data_dis]
            data_w = np.ones(data_ra.shape[0])
        else:
            data_positions = [data_ra,data_dec,data_z]
            data_w = np.ones(data_ra.shape[0])

        hdul.close()

        return data_positions,data_w



# read shifted galaxy positions and weights
fname = f"/project/rrg-wperciva/tsfraser/ELGMocks/CubicBox/ELG/z1.100/AbacusSummit_base_c000_ph000/"\
            f"ELG_frac_001snap16_ph000.gcat.txt"


fname_vox = '/project/rrg-wperciva/tsfraser/ELGvoids/CubicBox/voxel-Voids_cat.txt'

gal_positions,gal_weights =  txt_extract(fname,cols = np.array([0,1,2]),skip_lines =0,type ='notsurvey',shift =False)

vox_positions,_ = txt_extract(fname_vox,cols = np.array([1,2,3]),skip_lines =2 , type ='notsurvey', shift = False)

X = gal_positions[0]
Y= gal_positions[1]
Z =gal_positions[2]

vox_X = vox_positions[0]
vox_Y = vox_positions[1]
vox_Z = vox_positions[2]
print(len(X))

f,ax = plt.subplots(1,3,figsize =(30,10))

ax[0].set_title('CubicBox: ELG')
ax[0].set_xlabel('X [Mpc]')
ax[0].set_ylabel('Y [Mpc]')
ax[0].scatter(X[::1000],Y[::1000], color = 'blue', label = 'galaxies')
ax[0].scatter(vox_X[::100],vox_Y[::100],color = 'red',label = 'voids')
ax[0].legend()
ax[0].grid(True)

#plt.show()
ax[1].set_title('CubicBox: ELG')
ax[1].set_xlabel('Y [Mpc]')
ax[1].set_ylabel('Z [Mpc]')
ax[1].scatter(Y[::1000],Z[::1000],color = 'blue', label = 'galaxies')
ax[1].scatter(vox_Y[::100],vox_Z[::100],color = 'red',label = 'voids')
ax[1].legend()
ax[1].grid(True)
#plt.show()
ax[2].set_title('CubicBox: ELG')
ax[2].set_xlabel('X [Mpc]')
ax[2].set_ylabel('Z [Mpc]')
ax[2].scatter(X[::1000],Z[::1000],color = 'blue',label ='galaxies')
ax[2].scatter(vox_X[::100],vox_Z[::100],color = 'red',label = 'voids')
ax[2].legend()
ax[2].grid(True)
plt.savefig('ELG_box.pdf')