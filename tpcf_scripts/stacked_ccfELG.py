
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
    with fits.open(fname) as hdul:
        data = hdul[1].data

    if type == 'box':

        data_x = data[cols[0]].astype(np.float64)
        data_y = data[cols[1]].astype(np.float64)
        data_z = data[cols[2]].astype(np.float64)

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

        if shift == True:


            data_dis = cosmo.ComovingDistance(data_z) # Double check this


            data_positions = [data_ra,data_dec,data_dis]

        else:
            data_positions = [data_ra,data_dec,data_z]
        data_w = np.ones(data_ra.shape[0])

        fits.close()

        return data_positions,data_w



# read shifted galaxy positions and weights
fname_gals = f"/project/rrg-wperciva/tsfraser/ELGMocks/CubicBox/ELG/z1.100/AbacusSummit_base_c000_ph000/"\
            f"ELG_frac_001snap16_ph000.gcat.txt"

fname_randoms = "/project/rrg-wperciva/tsfraser/Revolver/rands_nz10_sv3y5_LRG_prefiltered.fits"


zobovfname = '/project/rrg-wperciva/tsfraser/CutSky/stacked_run2_lrgy5sv3zobov-Voids_cat.txt'
voxelfname = '/project/rrg-wperciva/tsfraser/LRGvoids/CutSky/stacked_run2_lrgy5sv3voxel-Voids_cat.txt'


gal_positions,gal_weights =  fits_extract(fname_gals,cols = ['RA','DEC','Z'],type ='survey',shift =True,zmin =0.65,zmax=0.95)
rand_positions,rand_weights = fits_extract(fname_randoms,cols = ['RA','DEC','Z'],type = 'survey',shift =True,zmin=0.65,zmax=0.95)

cols = np.array([0,1,2])
skip_lines = 2
type = 'survey'
shift = False
vox_positions,vox_weights = txt_extract(voxelfname,cols,skip_lines,type,shift)
zob_positions,zob_weights = txt_extract(zobovfname,cols,skip_lines,type,shift)




#zmin = 1.05
#zmax = 1.15
#idx = (data_z >zmin) & (data_z<zmax)
#data_dis = cosmo.ComovingDistance(data_z)
#print(data_dis[idx].shape,'GAL REDSHIFTS')

#data['WEIGHT_FKP'].astype(np.float64) * \
    #data['WEIGHT_SYSTOT'].astype(np.float64)  * \

#print(data_positions.shape,data_weights.shape,'GAL SIZES AND WEIGHT SIZES')
    #(data['WEIGHT_CP'].astype(np.float64) + \
    #data['WEIGHT_NOZ'].astype(np.float64) - 1)


#print(void_dis.shape,'VOID REDSHIFTS')
# read original random positions and weights
#fname = '/project/rrg-wperciva/tsfraser/ELGMocks/CutSky/ELG/z1.100/MOD_frac__caty5sv3_cutsky_ELG_random_S1000_1X.fits'

#with fits.open(fname) as hdul:
#    data = hdul[1].data
#data_z = data['Z'].astype(np.float64)
#idx = (data_z > zmin) & (data_z < zmax)
#data_ra = data['RA'].astype(np.float64)[idx]
#data_dec = data['DEC'].astype(np.float64)[idx]
#data_dis = cosmo.ComovingDistance(data_z)[idx]
#print(data_dis.shape,'RANDOM REDSHIFTS')
#randoms_positions = [data_ra, data_dec, data_dis]
#randoms_weights = np.ones(data_ra.shape[0])#data['WEIGHT_FKP'].astype(np.float64)[idx]


edges = np.linspace(1e-9, 150, 101)
boxsize =2000
print('running zobov 2PCF')
start = time.time()
res_zob =TwoPointCorrelationFunction('s',edges,data_positions1 = gal_positions , \
    data_weights1 = gal_weights, data_positions2 = zob_positions, data_weights2 = zob_weights,\
         randoms_position1 = rand_positions,randoms_weights1 =rand_weights,engine='corrfunc',boxsize =boxsize,nthreads =56) 
end = time.time()
print('finished zobov 2PCF, {:.3f}'.format(end-start))




print('running voxel 2PCF')
start = time.time()
res_vox = TwoPointCorrelationFunction('s',edges,data_positions1 = gal_positions, \
    data_weights1 = gal_weights , data_positions2 = vox_positions,data_weights2 = vox_weights,\
        randoms_position1 = rand_positions, randoms_weights1 = rand_weights,engine='corrfunc',boxsize =boxsize,nthreads =56)
end = time.time()
print('finished zobov 2PCF, {:.3f}'.format(end-start))


ax = plt.gca()
ax.plot(res_vox.sep, res_vox.corr)
ax.plot(res_zob.sep, res_zob.corr)
ax.legend()
ax.set_xlabel('$r$')
ax.set_ylabel(r'$r^{2}\xi(r)$')
ax.grid(True)
plt.savefig('LRG_stack.pdf')
plt.show()
