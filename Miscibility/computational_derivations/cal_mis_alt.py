import math as m
from matplotlib import pyplot as plt
import numpy as np
import csv 

def cal_entropy (v, phi):
    """
    v is volume per molecule in angstrom
    """
    return (phi / v * m.log (phi))   

def molecular_v(M,rho):
    """
    M = molar mass in g/mol
    rho = density in g/mL
    
    return molecular volume in Angstrom^3
    """
    
    v = M/rho*10/6.02
    
    return v


# temperature is in Kelvins
def free_energy_mixing (solA, solB, phiA, phiB, temp):
    """
    compute entropy of mixing. However, expression here is in 1/ m^3.
    There is a need to multiply 10 ** -30 in the end to convert this 
    into 1/ angstrom^3 for unit consistency
    """
    Beta = 1 / (temp * 1.38 * (10 ** -23))

    # converting from megapascals to pascals 
    def solSI_diff_squared (solA, solB):
        return ((solA - solB) * 1000) ** 2 

    energy_mix_angstrom = Beta * phiA * phiB * solSI_diff_squared (solA, solB) * (10 ** -30)
    return energy_mix_angstrom

def energy_mix_map (solA, solB, volA, volB, chemA, chemB, temp):
    acc_x = []
    acc_y = []
    for i in range (1, 1000):
        phiA = 0.001 * i 
        phiB = 1 - phiA
        energy_mix = free_energy_mixing (solA, solB, phiA, phiB, temp) + cal_entropy (volA, phiA) + cal_entropy (volB, phiB)
        acc_y.append (energy_mix)
        acc_x.append (phiA)
    plt.plot(acc_x, acc_y, color = "red")
    plt.title('free energy of mixing between ' + chemA + " and " + chemB)
    plt.xlabel('volume fraction of ' + chemA)
    plt.ylabel('free energy of mixing')
    plt.grid()
    plt.show()

# energy_mix_map (30.5, 48.0 , 67.23, 30.0, "methanol", "water", 298)
# energy_mix_map (18.7, 48.0 , molecular_v(84.16,0.779), 30.0, "cyclohexane", "water", 298)
# energy_mix_map (30.5, 18.7, 67.23, molecular_v(84.16,0.779), "methanol", "cyclohexane" 298)


def is_miscible_v1 (solA, solB, volA, volB, temp): 
    acc_y = []
    acc_x = []
    for i in range (1, 1000): # for volume of mixing, vol fractions must > 0 
        phiA = 0.001 * i
        phiB = 1 - phiA
        energy_mix = free_energy_mixing (solA, solB, phiA, phiB, temp) + cal_entropy (volA, phiA) + cal_entropy (volB, phiB)
        acc_x.append (phiA) 
        acc_y.append (energy_mix)
    first_deri = np.gradient (acc_y)
    second_deri = np.gradient (first_deri) # skipping first_deri, since it's not needed
    for y in second_deri: 
        if y < 0:
            return False
    return True

def mole_to_molar (x):
    z =  x * 6.02 / 10
    return z

def miscible_map (solA, volA, substance, temp):
    x_mis = []
    y_mis = []
    x_immis = []
    y_immis = []
    for i in range (100,300):
        solB = 0.1 * i 
        for j in range (1,400):
            volB = j 
            is_miscible = is_miscible_v1(solA, solB, volA, volB, temp)
            if is_miscible: # is_miscible returns a boolean
                x_mis.append(volB)
                y_mis.append(solB)
            else: 
                x_immis.append(volB)
                y_immis.append(solB)
    # convert to molar volume
    for i in range(len(x_mis)):
        x_mis[i] = mole_to_molar(x_mis[i])
    for i in range(len(x_immis)):
        x_immis[i] = mole_to_molar(x_immis[i])
    # small circle size to minimize overlap
    plt.scatter(x_mis, y_mis, s = 0.2, label = "miscible", color = "#ff5e6c")
    plt.scatter(x_immis, y_immis, s = 0.2, label = "immiscible", color = "#feb300")
    plt.xlabel('molar volume (ml/mol)')
    plt.ylabel('solubility coefficient')
    plt.title('miscibility map of ' + substance)
    plt.legend()
    plt.show()

def miscible_map_v1 (solA, volA, substance, temp, index):
    x_mis = []
    y_mis = []
    x_immis = []
    y_immis = []
    for i in range (100,500):
        solB = 0.1 * i 
        for j in range (1,400):
            volB = j 
            is_miscible = is_miscible_v1(solA, solB, volA, volB, temp)
            if is_miscible: # is_miscible returns a boolean
                x_mis.append(volB)
                y_mis.append(solB)
            else: 
                x_immis.append(volB)
                y_immis.append(solB)
    # convert to molar volume
    for i in range(len(x_mis)):
        x_mis[i] = mole_to_molar(x_mis[i])
    for i in range(len(x_immis)):
        x_immis[i] = mole_to_molar(x_immis[i])
    # small circle size to minimize overlap
    plt.scatter(x_mis, y_mis, s = 0.2, label = "miscible", color = "#ff5e6c")
    plt.scatter(x_immis, y_immis, s = 0.2, label = "immiscible", color = "#feb300")
    plt.xlabel('molar volume (ml/mol)')
    plt.ylabel('solubility coefficient')
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
                real_ymis.append(float(line[4]))
            elif line[index] == '0':
                real_ximmis.append(float(line[6]))
                real_yimmis.append(float(line[4]))
            else:
                pass
    plt.scatter(real_xmis, real_ymis, s = 5, label = "real miscible", color = "#000000", marker = "o")
    plt.scatter(real_ximmis, real_yimmis, s = 5, label = "real immiscible", color = "#FF0000", marker = "s")
    plt.legend(loc = 'upper right')
    plt.show()



# miscible_map_v1 (30.5, 67.23, "methanol", 298, 8)
miscible_map_v1 (48.0, 30.0, "water", 298, 7)
# miscible_map_v1 (16.6, molecular_v(84.16,0.779), "cyclohexane", 298, 9)