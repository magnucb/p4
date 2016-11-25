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

def colplot_4c(spin):
    # collects data from file and plots
    L = 20
    mc = int(1e5)
    temp = 100
    temps = [100, 120, 140, 160, 179, 199, 219, 240]
    for temp in temps:
        probname = "Prob_L"+str(L)+"_mc"+str(mc)+"_T"+str(temp)+"_spin"+str(spin)
        # Enername = "Energyprob_L"+str(L)+"_mc"+str(mc)+"_T"+str(T)+"_spin"+str(spin)
        mcs, E, Esqu, absM, Msqu = pl.loadtxt('../data/4c/'+probname+".dat", usecols=(0,1,2,3,4), unpack=True)

        pl.figure()
        pl.plot(mcs, E/E[-1], label=r'<E>')
        pl.plot(mcs, Esqu/Esqu[-1], label=r'<E$^2$>')
        pl.plot(mcs, absM/absM[-1], label=r'<|M|>')
        pl.plot(mcs, Msqu/Msqu[-1], label=r'<M$^2$>')
        pl.xlim([0, 14e3])
        pl.title("Values normalized with their final values, T=%g" % (temp/100.))
        pl.xlabel("Monte Carlo Cycles")
        pl.ylabel("Normalized values")
        pl.legend(loc="best")
        pl.savefig("../figs/4c/"+probname+".png")

    # for the sake of accepted configurations

    pl.figure()
    acceptlist = []
    for temp in temps:
        probname = "Prob_L"+str(L)+"_mc"+str(mc)+"_T"+str(temp)+"_spin"+str(spin)
        mcs, accepted = pl.loadtxt('../data/4c/'+probname+".dat", usecols=(0,5), unpack=True)
        acceptlist.append(accepted[-1])
        pl.plot(mcs, accepted/accepted[-1], label=r'T=%g' % (temp/100.))
        # pl.xlim([0, 14e3])
        pl.title("Accepted spin configurations for all temperatures")
        pl.xlabel("Monte Carlo Cycles")
        pl.ylabel("Accepted spin configurations")
        pl.legend(loc="best")
    
    pl.savefig("../figs/4c/acceptedspins_"+spin+".png")

    # for the sake of accepted configurations vs temp.

    return temps, acceptlist

def FourC():
    accepts = []
    spinconfigs = ["up", "random"]

    for spin in spinconfigs:
        temps, acceptlist = colplot_4c(spin)
        accepts.append(acceptlist)
    temps = pl.array(temps)/100.
    pl.figure()
    pl.plot(temps, accepts[0], label="In. spin "+spinconfigs[0])
    pl.plot(temps, accepts[1], label="In. spin "+spinconfigs[1])
    pl.title("Accepted spin configurations versus all temperatures")
    pl.xlabel("Temperature variation")
    pl.ylabel("Accepted spin configurations")
    pl.xlim([1., 2.40])
    # pl.yscale("log")
    pl.legend(loc="best")
    pl.savefig("../figs/4c/acceptedspins_vsTemps.png")


def FourD():
    # collects data from file and plots
    L = 20
    mc = int(1e5)
    temps = [100, 240]
    spinconfigs = ["up", "random"]
    most_often = {}

    for spin in spinconfigs:
        
        pl.figure()
        for temp in temps:

            Enername = "Energyprob_L"+str(L)+"_mc"+str(mc)+"_T"+str(temp)+"_spin"+str(spin)
            energies, variance = pl.loadtxt('../data/4c/'+Enername+".dat", usecols=(0,1), unpack=True)
            pl.hist(energies, normed=0, bins=100, histtype="step", label="Temp=%s" % temp)
            hist, bins = pl.histogram(energies, bins=len(pl.unique(energies)))
            E = (bins[:-1])[pl.argmax(hist)] + 0.5*(bins[1]-bins[0])
            most_often[spin+" "+str(temp)] = E, max(hist), variance[-1]

        pl.title("Energy occurrence histogram for spin %s" % spin)
        pl.xlabel("Occurring energies")
        pl.ylabel("Count of energy")
        pl.xlim([-820, -350])
        pl.legend(loc="best")
        pl.savefig("../figs/4d/probabilityhistogram_%s.png" % spin)
    for i, j in most_often.iteritems():
        print i, " energy:", j[0], "\n---          count:", j[1]
        print "        Prob of state: %g " %(j[1]/87000.)
        print "             Variance: %g " %(j[2])

if __name__ == '__main__':

    # FourC()
    FourD()