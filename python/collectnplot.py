import pylab as pl
import sys
from matplotlib import rc
rc('font',**{'family':'serif'})
"""
    mcs = (int)1e5;
    n_spins = 2;
    initial_temp = 1;
    final_temp = 1;
    temp_step = 1e-3;
"""

def colplot_data(L, mc, T, spin):
    # collects data from file and plots
    probname = "Prob_L"+str(L)+"_mc"+str(mc)+"_T"+str(T)+"_spin"+str(spin)+".dat"
    Enername = "Energyprob_L"+str(L)+"_mc"+str(mc)+"_T"+str(T)+"_spin"+str(spin)+".dat"

    mcs, E, Esqu, absM, Msqu, accepted = pl.loadtxt('../data/'+probname, usecols=(0,1,2,3,4,5), unpack=True)
    
    pl.figure()
    pl.plot(mcs, E, label=r'Energy')
    pl.plot(mcs, Esqu, label=r'Energy$^2$')
    pl.plot(mcs, absM, label=r'|Magnetic Moment|')
    pl.plot(mcs, Msqu, label=r'(Magnetic Moment)$^2$')
    pl.legend()
    pl.show()

if __name__ == '__main__':
    L = int(sys.argv[1])
    mc = int(sys.argv[2])
    T = int(sys.argv[3])
    spin = str(sys.argv[4])

    colplot_data(L, mc, T, spin)
