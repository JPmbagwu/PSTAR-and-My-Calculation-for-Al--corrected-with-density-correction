#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 13:59:06 2024

@author: johnpaulmbagwu
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Constants
M0_c = 938.272  # MeV/c^2 (rest mass of proton)
z = 1  # Charge of proton

# Function to calculate the beta value
def beta(T):
    b = (1 - (1 / ((T/M0_c) + 1))**2)**(1/2)
    return b

# Function to calculate the mass stopping power with density correction
def dedx(T, I, Z_A, rho):
    be = beta(T)
    SP = 0.3071 * Z_A * z**2 * (1/be**2) * (13.8373 + np.log(be**2 / (1 - be**2)) - be**2 - np.log(I))
    SP_corr = SP * (rho / 2.7)**0.15  # Sternheimer density correction
    return SP_corr

# Read the TXT data file containing information about Aluminum Oxide
df = pd.read_csv('Aluminum.txt', delimiter=' ')

# Extract the stopping power and energy columns from the dataframe
stopping_power = df['TotalStp.Pow']
energy = df['KineticEnergy']

# Calculate the stopping power values using the dedx function with density correction
density_aluminum = 2.7  # g/cm^3 for aluminum
Z_A_Al = 13/26.98  # Atomic number-to-mass ratio for aluminum
SP_values_corr = dedx(energy, 166, Z_A_Al, density_aluminum)

# Plotting
plt.figure(figsize=(10, 6))

# Plotting both calculated and experimental stopping power
plt.plot(energy, SP_values_corr, label='My Calculation (with density correction)')
plt.plot(energy, stopping_power, label='PSTAR Data')

plt.title('Aluminum')
plt.xlabel('Energy (MeV)')
plt.ylabel('Stopping Power (MeV/(g/cm^2))')

# Set the x and y-axis scales to logarithmic for better visualization
plt.xscale('log')
plt.yscale('log')

# Set the x-axis and y-axis limits for better visibility of data points
plt.xlim(1e-2, 1e4)

plt.grid(True)
plt.legend()

plt.show()
