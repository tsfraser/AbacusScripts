
from astropy.io import fits
import numpy as np
from pancakes.pancakes.utilities.cosmology import Cosmology
from pycorr import (TwoPointCorrelationFunction,
    LandySzalayTwoPointEstimator, setup_logging,
    TwoPointCounter,AnalyticTwoPointCounter,NaturalTwoPointEstimator
    )
import time
import matplotlib.pyplot as plt
from contrast import contrast
from contrast.contrast import box
from contrast.contrast.box import correlation_functions
from contrast.contrast.box.correlation_functions import tpcf_r


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

def fits_extract(fname,cols,type,shift):
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

        if shift == True:


            data_dis = cosmo.ComovingDistance(data_z) # Double check this


            data_positions = [data_ra,data_dec,data_dis]

        else:
            data_positions = [data_ra,data_dec,data_z]
        data_w = np.ones(data_ra.shape[0])

        fits.close()

        return data_positions,data_w


def plot_tpcf(res,title,isr2,ismultiple,void_labels,fname):
    if isr2 is True:
        if ismultiple is True:
            res_list = []

            plt.figure(1,figsize =(8,6))
            plt.title(title,fontsize=16)
            plt.xlabel('$r$')
            plt.ylabel(r'$r^{2}\xi(r)$')
            for i in range(len(res)):
                res_list.append(res[i])

                plt.plot(res[i].sep,res[i].corr*res[i].sep**2,lw=3,label = void_labels[i])

            plt.legend(fontsize=12)
            plt.savefig(fname)
            plt.clf()
        else:
            plt.figure(1,figsize =(8,6))
            plt.title(title,fontsize=16)
            plt.xlabel('$r$')
            plt.ylabel(r'$r^{2}\xi(r)$')
            plt.plot(res.sep,res.corr*res.sep**2,lw=3,label = 'CCF')
            plt.legend(fontsize=12)
            plt.savefig(fname)
            plt.clf()
    else:
        if ismultiple is True:
            res_list = []

            plt.figure(1,figsize =(8,6))
            plt.title(title,fontsize=16)
            plt.xlabel('$r$')
            plt.ylabel(r'$\xi(r)$')
            for i in range(len(res)):
                res_list.append(res[i])

                plt.plot(res[i].sep,res[i].corr,lw=3,label = void_labels[i])

            plt.legend(fontsize=12)
            plt.savefig(fname)
            plt.clf()
        else:
            plt.figure(1,figsize =(8,6))
            plt.title(title,fontsize=16)
            plt.xlabel('$r$')
            plt.ylabel(r'$\xi(r)$')
            plt.plot(res.sep,res.corr,lw=3,label = 'CCF')
            plt.legend(fontsize=12)
            plt.savefig(fname)
            plt.clf()


    return None


def plot_tpcf_manual(term,sep,title,isr2,ismultiple,void_labels,fname):
    '''
    What you use if you are not using a Pycorr TPCF object
    
    '''
    if isr2 is True:
        if ismultiple is True:
            res_list = []

            plt.figure(1,figsize =(8,6))
            plt.title(title,fontsize=16)
            plt.xlabel('$r$')
            plt.ylabel(r'$r^{2}\xi(r)$')
            for i in range(len(term)):
                res_list.append(term)

                plt.plot(sep,res_list[i]*sep**2,lw=3,label = void_labels[i])

            plt.legend(fontsize=12)
            plt.savefig(fname)
            plt.clf()
        else:
            plt.figure(1,figsize =(8,6))
            plt.title(title,fontsize=16)
            plt.xlabel('$r$')
            plt.ylabel(r'$r^{2}\xi(r)$')
            plt.plot(sep,term*sep**2,lw=3,label = 'CCF')
            plt.legend(fontsize=12)
            plt.savefig(fname)
            plt.clf()
    else:
        if ismultiple is True:
            res_list = []

            plt.figure(1,figsize =(8,6))
            plt.title(title,fontsize=16)
            plt.xlabel('$r$')
            plt.ylabel(r'$\xi(r)$')
            for i in range(len(term)):
                res_list.append(term[i])

                plt.plot(sep,res_list[i],lw=3,label = void_labels[i])

            plt.legend(fontsize=12)
            plt.savefig(fname)
            plt.clf()
        else:
            plt.figure(1,figsize =(8,6))
            plt.title(title,fontsize=16)
            plt.xlabel('$r$')
            plt.ylabel(r'$\xi(r)$')
            plt.plot(sep,term,lw=3,label = 'CCF')
            plt.legend(fontsize=12)
            plt.savefig(fname)
            plt.clf()


    return None



if __name__ =='__main__':


    N_bins = 100
    edges = np.linspace(1e-9, 150, N_bins+1)
    boxsize =2000
# read shifted galaxy positions and weights
    fname = f"/project/rrg-wperciva/tsfraser/ELGMocks/CubicBox/ELG/z1.100/AbacusSummit_base_c000_ph000/"\
            f"ELG_frac_001snap16_ph000.gcat.txt"
    #data  = np.genfromtxt(fname)
    data_pos,data_weights = txt_extract(fname,[0,1,2],skip_lines=0,type='box',shift =False)
    data_x = data_pos[0][::100]
    data_y = data_pos[1][::100]
    data_z = data_pos[2][::100]
    data_pos = [data_x,data_y,data_z]
    data_weights = np.ones(data_x.shape[0])
#data_x = data[::100,0].astype(np.float64)
#data_y = data[::100,1].astype(np.float64)
#data_z = data[::100,2].astype(np.float64)
    print(data_x.shape,'Number of galaxies')



#zmin = 1.05
#zmax = 1.15
#idx = (data_z >zmin) & (data_z<zmax)
#data_dis = cosmo.ComovingDistance(data_z)
#print(data_dis[idx].shape,'GAL REDSHIFTS')


#data_positions = [data_x, data_y, data_z]
#data_weights = np.ones(data_x.shape[0])


#data['WEIGHT_FKP'].astype(np.float64) * \
    #data['WEIGHT_SYSTOT'].astype(np.float64)  * \

#print(data_positions.shape,data_weights.shape,'GAL SIZES AND WEIGHT SIZES')
    #(data['WEIGHT_CP'].astype(np.float64) + \
    #data['WEIGHT_NOZ'].astype(np.float64) - 1)

    zobovfname = '/project/rrg-wperciva/tsfraser/ELGvoids/CubicBox/zobov-Voids_cat.txt'
    zobovdata = np.genfromtxt(zobovfname,skip_header=2)
#    data = hdul[1].data
#
#zmin = 1.05
#zmax = 1.15

#void_z = data[:,3].astype(np.float64)
#zidx = (void_z >zmin) & (void_z<zmax)


    void_x = zobovdata[:,1].astype(np.float64)#[zidx]#data['RA'].astype(np.float64)
    void_y = zobovdata[:,2].astype(np.float64)#[zidx]#data['DEC'].astype(np.float64)
    void_z = zobovdata[:,3].astype(np.float64)



#[zidx]#data['Z'].astype(np.float64)
#void_dis = cosmo.ComovingDistance(void_z)


    zobov_void_positions = [void_x, void_y, void_z]
    void_weights = np.ones(void_x.shape[0])
    print(void_x.shape,'number of zobov voids')


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



    print('running zobov 2PCF')
    start = time.time()
    res_zob =TwoPointCorrelationFunction('s',edges,data_positions1 = data_pos , data_weights1 = data_weights, engine='corrfunc',boxsize =boxsize,nthreads =56) 
    end = time.time()
    print('finished zobov 2PCF, {:.3f}'.format(end-start))


    voxelfname = '/project/rrg-wperciva/tsfraser/ELGvoids/CubicBox/voxel-Voids_cat.txt'
    voxeldata = np.genfromtxt(voxelfname,skip_header=2)

    void_x = voxeldata[:,1].astype(np.float64)
    void_y = voxeldata[:,2].astype(np.float64)
    void_z = voxeldata[:,3].astype(np.float64)

    voxel_void_positions = [void_x, void_y, void_z]
    void_weights = np.ones(void_x.shape[0])
    print(void_x.shape,'number of voxel voids')


    print('running voxel 2PCF')
    start = time.time()
    res_vox = TwoPointCorrelationFunction('s',edges,data_positions1 = data_pos,data_weights1 = data_weights , engine='corrfunc',boxsize =boxsize,nthreads =56)
    end = time.time()
    print('finished zobov 2PCF, {:.3f}'.format(end-start))



    plot_tpcf([res_zob,res_vox],'ELG voids, box sims',True,True,['zobov','voxel'],'r2autoVoidCCF_ELG.pdf')

    plot_tpcf([res_zob,res_vox],'ELG voids, box sims',False,True,['zobov','voxel'],'autoVoidCCF_ELG.pdf')
    
    ax = plt.gca()
    ax.plot(res_vox.sep, res_vox.corr*res_vox.sep**2)
    ax.plot(res_zob.sep, res_zob.corr*res_zob.sep**2)
    ax.legend()
    ax.set_xlabel('$r$')
    ax.set_ylabel(r'$r^{2}\xi(r)$')
    ax.grid(True)
    plt.savefig('autoVoidELGCCF.pdf')
#plt.show()
    plt.clf()
    D1D2 =  TwoPointCounter('s',edges,positions1 = data_pos, positions2 = voxel_void_positions, weights1 = data_weights,\
        weights2= void_weights,position_type= 'xyz',engine ='corrfunc',\
        boxsize =boxsize, nthreads =56 )

    R1R2 = AnalyticTwoPointCounter('s',edges,boxsize= boxsize)

    result2 = NaturalTwoPointEstimator(D1D2=D1D2, R1R2=R1R2)
    plt.figure(1)
    plt.plot(result2.sep,result2.corr,'k-')
    plt.xlabel('$r$')
    plt.ylabel(r'$\xi(r)$')
    plt.grid(True)
    plt.savefig('MANUAL_autoVoidELGCCF.pdf')
    plt.clf()

    plot_tpcf(result2,'Manual ELG Void box sim CCF',True,False,['voids'],'MAN_auto_void_CCF_elg.pdf')


    plt.figure(2)
    plt.plot(result2.sep,result2.corr*result2.sep**2,'k-')
    plt.xlabel('$r$')
    plt.ylabel(r'$r^{2}\xi(r)$')
    plt.grid(True)
    plt.savefig('2_MANUAL_autoVoidELGCCF.pdf')

    sep =  result2.sep
    D1D2 = D1D2.normalized_wcounts()
    R1R2 = R1R2.normalized_wcounts()
    term = D1D2/R1R2 -1 


    plt.figure(3)
    plt.plot(sep,term,'k-')
    plt.xlabel('$r$')
    plt.ylabel(r'$\xi(r)$')
    plt.grid(True)
    plt.savefig('SUPERMANUAL_autoVoidELGCCF.pdf')
    plt.clf()

    plot_tpcf_manual(term,sep,'Manual TPCF, ELG box voids',False,False,['void'],'supermanual_auto_void_ELGCCF.pdf')


    delta_r = tpcf_r(data_pos, edges, box_size=boxsize,postions2= voxel_void_positions,nthreads= 56)

    plt.plot(sep, delta_r)
    plt.savefig('EPaillas.pdf')
