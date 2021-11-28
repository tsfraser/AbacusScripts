import numpy as np
import matplotlib.pyplot as plt
import astropy
from pycorr import NaturalTwoPointEstimator

dir =''

nomen = 'test_abacus_periodic_ELG.npy'
fname = dir+nomen

def read_file(fname):
    data = NaturalTwoPointEstimator.load(fname)
    sep = data.sep
    xi = data.corr

    return sep,xi

def plot_corr(sep,corr,title,label):
    plt.figure(1,figsize = (10,8))
    plt.title(title,fontsize = 18)
    if label is None:
        plt.plot(sep,corr, color = 'red')
    else:
        plt.plot(sep,corr, color = 'red', label = label)
    plt.xlabel('Separation, s', fontsize=  14)
    plt.ylabel('Xi', fontsize = 14)
    plt.legend(loc='best')
    plt.savefig(title+label+'.png')
    plt.show()

if __name__ == '__main__':

    #data = NaturalTwoPointEstimator.load(fname)
    sep,corr = read_file(fname)
    plot_corr(sep,corr,'ELG_old_HOD_Abacus_base','ELG')




  
