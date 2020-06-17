import math as m
from matplotlib import pyplot as plt
import numpy as np
import csv 
from scipy.optimize import fsolve

# phi = volume fraction
# v = volume of molecule
# y = dipole moment

INV_BETA = (1 / 298 * 1.38 * 10 ** (-23)) ** -1 # by default rtp

def entropy_v1 (v, phi):
    return (phi / v * m.log (phi))   

def dipole_diff_v1 (v, phi, y, net_dip):
    return (0.5 * phi / v) * ((y / net_dip * m.log(1 + net_dip)) - m.log(1 + y))

# maps phiA to free energy of mixing
def energy_mix_map (va, vb, ya, yb, chemA, chemB): 
    acc_y = []
    acc_x = []
    for i in range (1, 1000): # for volume of mixing, vol fractions must > 0 
        phiA = 0.001 * i
        phiB = 1 - phiA
        net_dip = (phiA * ya) + (phiB * yb)
        subA = dipole_diff_v1 (va, phiA, ya, net_dip) + entropy_v1 (va, phiA)
        subB = dipole_diff_v1 (vb, phiB, yb, net_dip) + entropy_v1 (vb, phiB)
        acc_x.append (phiA) 
        acc_y.append (subA + subB)
    plt.plot(acc_x, acc_y, color = "red")
    plt.title('free energy of mixing between ' + chemA + " and " + chemB)
    plt.xlabel('volume fraction of ' + chemA)
    plt.ylabel('free energy of mixing')
    plt.grid()
    plt.show()

# example shown below 
# energy_mix_map (134.05, 30.0, 2.84, 29.8, "trichloromethane", "water")
# energy_mix_map (67.23, 30.0, 14.27, 29.8, "methanol", "water")


# determines whether 2 liquids are miscible
def is_miscible_v1 (va, vb, ya, yb): 
    acc_y = []
    acc_x = []
    for i in range (1, 100): # for volume of mixing, vol fractions must > 0 
        phiA = 0.01 * i
        phiB = 1 - phiA
        net_dip = (phiA * ya) + (phiB * yb)
        subA = dipole_diff_v1 (va, phiA, ya, net_dip) + entropy_v1 (va, phiA)
        subB = dipole_diff_v1 (vb, phiB, yb, net_dip) + entropy_v1 (vb, phiB)
        acc_x.append (phiA) 
        acc_y.append (subA + subB)
    first_deri = np.gradient (acc_y)
    second_deri = np.gradient (first_deri) # skipping first_deri, since it's not needed
    for y in second_deri: 
        if y < 0:
            return False
    return True
            
# creates list of coordinates of vb and yb that are miscible and not miscible
def miscible_map (va, ya, substance):
    x_mis = []
    y_mis = []
    x_immis = []
    y_immis = []
    for vb in range (1,400):
        for yb in range (1,400):
            is_miscible = is_miscible_v1 (va, vb, ya, yb)
            if is_miscible: # is_miscible returns a boolean
                x_mis.append(vb)
                y_mis.append(yb)
            else: 
                x_immis.append(vb)
                y_immis.append(yb)
    # small circle size to minimize overlap
    plt.scatter(x_mis, y_mis, s = 0.2, label = "miscible", color = "#ff5e6c")
    plt.scatter(x_immis, y_immis, s = 0.2, label = "immiscible", color = "#feb300")
    plt.xlabel('volume')
    plt.ylabel('dipole moment')
    plt.title('miscibility map of ' + substance)
    plt.legend()
    plt.show()

# miscible_map, except now it also displays a few additional chemicals 
# also now computes in terms of molar volume and dielectric constant

def mole_to_molar (x):
    z =  x * 6.02 / 10
    return z

def dip_to_dielec (y):
    a = (3 * y) * (2 * (y ** 2) + 3 * y + 9 ) / ((y + 3) ** 2) 
    z = 0.25 * (1 + a + m.sqrt(9 + 2 * a + a ** 2))
    return z


def dielec_to_y(eps):
    """
    Convert dielectirc constant to the value of y
    """
    
    def solve_y(y):
        return dip_to_dielec (y)-eps
    
    y_guess = (eps-1)*(2*eps+1)/eps/3
    y = fsolve(solve_y, y_guess)[0]
    
    return y

def molecular_v(M,rho):
    """
    M = molar mass in g/mol
    rho = density in g/mL
    
    return molecular volume in Angstrom^3
    """
    
    v = M/rho*10/6.02
    
    return v


def miscible_map_v1 (va, ya, substance, index): #index is related to .csv file
    x_mis = []
    y_mis = []
    x_immis = []
    y_immis = []
    for vb in range (1, 500):
        for yb in range (1, 250):
            yb = 0.125 * yb
            is_miscible = is_miscible_v1 (va, vb, ya, yb)
            if is_miscible: # is_miscible returns a boolean
                x_mis.append(vb)
                y_mis.append(yb)
            else: 
                x_immis.append(vb)
                y_immis.append(yb)
    # convert x-axis from molecular volume to molar volume
    for i in range(len(x_mis)):
        x_mis[i] = mole_to_molar(x_mis[i])
        y_mis[i] = dip_to_dielec(y_mis[i])
    for i in range(len(x_immis)):
        x_immis[i] = mole_to_molar(x_immis[i])
        y_immis[i] = dip_to_dielec(y_immis[i])
    # small circle size to minimize overlap
    plt.scatter(x_mis, y_mis, s = 1, label = "miscible", color = "#E6E8C3")
    plt.scatter(x_immis, y_immis, s = 1, label = "immiscible", color = "#B8F09A")
    plt.xlabel('molar volume (cm^3)')
    plt.ylabel('dielectric constant')
    plt.title('miscibility map of ' + substance)
    real_xmis = []
    real_ymis = []
    real_ximmis = []
    real_yimmis = []
    with open ('data.csv', 'r', encoding = "unicode-escape") as csv_file:
        csv_reader = csv.reader(csv_file)
        next(csv_reader)
        for line in csv_reader:
            if line[index] == '1':
                real_xmis.append(float(line[6]))
                real_ymis.append(float(line[1]))
            elif line[index] == '0':
                real_ximmis.append(float(line[6]))
                real_yimmis.append(float(line[1]))
            else:
                pass
    plt.scatter(real_xmis, real_ymis, s = 5, label = "real miscible", color = "#000000", marker = "o")
    plt.scatter(real_ximmis, real_yimmis, s = 5, label = "real immiscible", color = "#FF0000", marker = "s")
    plt.legend(loc = 'upper right')
    plt.show()

# execute one at a time
# miscible_map_v1 (30, 29.8, "water", 7)
# miscible_map_v1 (67.23, 14.27, "methanol", 8)
miscible_map_v1 (molecular_v(84.16,0.779), dielec_to_y(2.0243), "cyclohexane", 9)



