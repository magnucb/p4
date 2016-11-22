import numpy as np
import pylab as pl

"""
def makeLattice(N=2):
    # makes lattice matrix
    lattice = -1*pl.ones((N,N))
    for i in range(N):
        for j in range(N):
            lattice[i,j] = lattice[i,j]**np.randint(2)

    return lattice

def _Energy(lattice):
    # finds energy of the lattice
    tot_en = 0.
    rows, cols = pl.array(shape(lattice)).astype(int)
    for j in range(rows):
        for i in range(cols):
            try:
                tot_en -= lattice[j,i]*lattice[j,i-1]
            except IndexError:
                continue
            try:
                tot_en -= lattice[j,i]*lattice[j,i+1]
            except IndexError:
                continue
            try:
                tot_en -= lattice[j,i]*lattice[j-1,i]
            except IndexError:
                continue
            try:
                tot_en -= lattice[j,i]*lattice[j+1,i]
            except IndexError:
                continue
    return tot_en

def _Magnetism(lattice):
    # finds magnetic moment of the lattice
    tot_mag = sum(lattice)
    return tot_mag

def create_dataset(Nset=1000, Ncrystal=2):
    # makes a dataset from which to make a partition function
    en_mag = pl.zeros((Nset, 2))

    for j in pl.arange(Nset):
        lattice     = makeLattice(Ncrystal)
        en_mag[j]   = _Energy(lattice), _Magnetism(lattice)
    return en_mag
"""

def partfunc(energies):
    # makes a partition function for the current microstate
    Z = 0
    for i in range(len(energies)):
        Z += pl.exp(-energies[i]*beta)
    return Z

def en_expect_pow(energies, Z, power):
    # calculates expected energies
    E = 0
    for i in range(len(energies)):
        E += (energies[i]**power)*pl.exp(-energies[i]*beta)
    return E/Z

def magn_expect_pow(magnetic, energies, Z, power):
    # calculates expected magnetic moment
    M = 0
    for i in range(len(energies)):
        M += (abs(magnetic[i]**power))*pl.exp(-energies[i]*beta)
    return M/Z

def spec_heat(EnPow2_exp, En_expPow2, k_B, T):
    # calculates specific heat of the system
    return (EnPow2_exp - En_expPow2)/(k_B*(T**2))

def magn_mom(MagPow2_exp, Mag_expPow2, k_B, T):
    # calculates magnetic moment of the system
    return (MagPow2_exp - Mag_expPow2)/(k_B*T)

T = 1.
k_B = 1.
beta = 1./(k_B*T)
L = 2
tot_no_microstates = L**4 # 16
J = 1.

###
### analytical
###

energies = pl.zeros(tot_no_microstates)
magnetic = pl.zeros(tot_no_microstates)
energies[0], energies[1:9], energies[9:11], energies[11:15], energies[15] = -8*J, 0, 8*J, 0, -8*J
magnetic[0], magnetic[1:5], magnetic[5:11], magnetic[11:15], magnetic[15] = 4, 2, 0, -2, -4

Z = partfunc(energies)
en_exp      = en_expect_pow(energies, Z, power=1)
EnPow2_exp  = en_expect_pow(energies, Z, power=2)
magn_exp    = magn_expect_pow(magnetic, energies, Z, power=1)
MagPow2_exp = magn_expect_pow(magnetic, energies, Z, power=2)

# print "Z      :", Z
print "<E>    :", en_exp
# print "<E**2> :", EnPow2_exp
# print "<E>**2 :", en_exp**2
print "<M>    :", magn_exp
# print "<M**2> :", MagPow2_exp
# print "<M>**2 :", magn_exp**2
print "Cv     :", spec_heat(EnPow2_exp, en_exp**2, k_B, T)
print "X      :", magn_mom(MagPow2_exp, magn_exp**2, k_B, T)
