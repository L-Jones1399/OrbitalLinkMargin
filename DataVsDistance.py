# Data vs Distance Plots
# Author: Jackie Scholes
# Date: 1/11/2024
# Description: Various functions relating to link budgeting, primarily as a function of distance for lunar applications

import matplotlib.pyplot as plt
import numpy as np
import math
from BigLinkBudget import *

# Defining constants

c = 3E8 # Speed of light in metres per second
h = 6.625E-34 # Planck's constant in J/(Hz*photon)

# Defining input parameters

wavelength = 1.553E-6 # Wavelength of laser light in metres
bitrate = 1E9 # Target transmission bitrate (bits per second)

powers = np.array([0.5, 1, 3]) # Input values for power
distances = range(360000000, 405000000, 100) # Input values for distance [min, max, step]
margins = [] # Empty list for link margin values

time_slot = 1/bitrate # Length of the OOK time slot in seconds
photon_energy = h * c / wavelength # Energy of photon in joules

def main():

    main_plot()
    aperture_plot()
    photons_per_pulse()
    photons_per_bit()
    irradiance_plot()
    loss_Tx_vs_Rx()

def main_plot():

    # This function plots the link margin in dBm as a function of distance for two different scenarios
    # I'm pretty sure this is modelling the lunar gateway NRHO orbit (i forgot lol), but this can easily be adapted for LEO, etc.

    powers = np.array([0.04, 0.39]) # Input values for laser power in watts. First case un-amplified, second case using EDFA
    distances = [3000000, 70000000, 380000000] # Input values for distance
    margins = [] # Empty list for link margin values

    for power in powers:

        for distance in distances:

            if distance == 3000000 or 70000000:
                rx_aperture = 0.3 # Receiver aperture in metres
            else:
                rx_aperture = 0.7 # Receiver aperture in metres

            margins.append(link_budget(laser_power=power, wavelength=wavelength, tx_aperture=0.1, rx_aperture=rx_aperture, orbit_altitude=distance, elevation_angle=45, tx_transmission=1, rx_transmission=1, rx_obscuration=0, rx_smf_il=0))

            # Now print margin as a function of distance for each power
        
        plt.plot(distances, margins, label=power)
        plt.legend(title="Power (W)")
        plt.xlabel("Distance (m)")
        plt.ylabel("Link Margin (dBm)")
        plt.grid()
        plt.savefig("test.png", dpi=500)

        margins = [] 

    plt.show()

def aperture_plot():

    # This code plots link budgets as a function of transmitter and receiver apertures for some assumed inputs.

    power = 1 # Input value for transmitter laser power in watts
    distance = 400000000 # Input value for distance in metres (Earth to moon)
    apertures = [x / 100 for x in range(10, 50, 5)] # Input values for apertures
    margins = [] # Empty list for link margin values

    for D_Tx in apertures:

        for D_Rx in apertures:

            margins.append(link_budget(laser_power=power, wavelength=wavelength, tx_aperture=D_Tx, rx_aperture=D_Rx, orbit_altitude=distance, elevation_angle=45, tx_transmission=1, rx_transmission=1, rx_obscuration=0, rx_smf_il=0))

            # Now print margin as a function of distance for each power
        
        plt.plot(apertures, margins, label=power)
        plt.legend(apertures, title="Transmitter Size (m)")
        plt.xlabel("Receiver Size (m)")
        plt.ylabel("Link Margin (dBm)")
        plt.grid()
        plt.savefig("Apertures.png", dpi=500)

        margins = [] 

    plt.show()

def photons_per_pulse():

    # This code determines the number of photons that will arrive at the receiver per pulse as a function of distance and laser power
    # This is a useful metric, since some detectors require a certain number of photons per pulse for successful readout

    print("Photon energy is", photon_energy)
    margins = []
    D_Tx = D_Rx = 0.2 # Transmitter/receiver diameter in metres

    for power in powers:

        for distance in distances:

            p_av = convert_from_dBm(link_budget(laser_power=power, wavelength=wavelength, tx_aperture=D_Tx, rx_aperture=D_Rx, orbit_altitude=distance, elevation_angle=45, tx_transmission=1, rx_transmission=1, rx_obscuration=0, rx_smf_il=0)) # Average power at receiver in watts
            print("P_av is", p_av, "W")
            margins.append(2 * p_av * time_slot / photon_energy) # Number of photons per pulse
        
        plt.plot(distances, margins, label=power)
        plt.legend(title="Power (W)")
        plt.xlabel("Distance (m)")
        plt.ylabel("Photons per Pulse")
        plt.grid()
        plt.savefig("test.png", dpi=500)

        margins = [] 

    plt.show()

def photons_per_bit():

    # This code determines the number of photons that will arrive at the detector per bit as a function of transceiver diameter and laser power
    # The red line indicates the minimum number of photons for (I BELIEVE) heterodyne receiver-based OOK @ 4Gbps according to Hemmati (again i forgot lol)

    print("Photon energy is", photon_energy)
    margins = []
    power = 0.5
    distance = 400000000
    apertures = [x / 1000 for x in range(100, 250, 5)] # Input values for apertures

    for D_Tx in apertures:


        p_av = convert_from_dBm(link_budget(laser_power=power, wavelength=wavelength, tx_aperture=D_Tx, rx_aperture=D_Tx, orbit_altitude=distance, elevation_angle=45, tx_transmission=1, rx_transmission=1, rx_obscuration=0, rx_smf_il=0)) # Average power at receiver in watts
        print("P_av is", p_av, "W")
        margins.append(p_av * time_slot / photon_energy) # Number of photons per bit

        
    plt.plot(apertures, margins, label=power)
    plt.axhline(y=36 + convert_from_dB(6.9), color="red")
    plt.title("Required Photons per Bit vs Transceiver Diameter")
    plt.legend(title="Power (W)")
    plt.xlabel("Transceiver Aperture (m)")
    plt.ylabel("Photons per Bit")
    plt.grid(axis='x')
    plt.savefig("test.png", dpi=500)

    margins = [] 

    plt.show()

def irradiance_plot():

    # This code determines the irradiance at the receiver as a function of distance and laser power.
    # Once again, a useful metric to determine whether you are going to get signal readout.

    irradiances = [] # Initialising empty list

    D_Tx = D_Rx = 0.2 # Transceiver diameter in metres

    for power in powers:

        for distance in distances:

            p_av = convert_from_dBm(link_budget(laser_power=power, wavelength=wavelength, tx_aperture=D_Tx, rx_aperture=D_Rx, orbit_altitude=distance, elevation_angle=45, tx_transmission=1, rx_transmission=1, rx_obscuration=0, rx_smf_il=0)) # Average power at receiver in watts
            area_Rx = math.pi * (0.25/2)**2 # Area of Rx antenna
            irradiance = p_av / area_Rx
            irradiances.append(irradiance) # Number of photons per pulse
        
        plt.plot(distances, irradiances, label=power)
        plt.legend(title="Power (W)")
        plt.xlabel("Distance (m)")
        plt.ylabel("Irradiance (W/m^2)")
        plt.grid()
        plt.savefig("test.png", dpi=500)

        irradiances = [] 

    plt.show()

def loss_Tx_vs_Rx():

    # This thing is fucking awesome
    # It's a contour plot of the expected losses induced (ONLY BY FREE SPACE PATH LOSS) for given transmitter and receiver diagrams
    # I forgot I made this:D i like it

    # Setting up constants
    wavelength = 1.550E-6 # Laser wavelength in metres
    distance = 400000000 # Distance to the moon in metres

    # Setting up receiver and transmitter diameters
    d_tx = [x / 100 for x in range(10, 100, 1)]
    d_rx = [x / 100 for x in range(10, 100, 1)]

    D_TX, D_RX = np.meshgrid(d_tx, d_rx)

    # Initialising empty list for losses
    total_losses_dB = np.zeros_like(D_TX)

    for ix, i in enumerate(d_tx):

        for jx, j in enumerate(d_rx):

            fspl = (wavelength/(4 * math.pi * distance))**2 # Free space path loss
        
            theta_pe = 2E-6 # Pointing error angle
            beam_waist = i/(2*1.21) # Laser beam waist
            theta_d = (wavelength / (math.pi * beam_waist)) # Beam divergence angle

            pointing_loss = math.e**((-2 * theta_pe**2)/(theta_d**2)) # Pointing error induced loss

            # Convert losses to dB
            fspl_dB = convert_to_dB(fspl)
            pointing_loss_dB = convert_to_dB(pointing_loss)
            inherent_loss_dB = -6 # Tx and Rx losses, inhernetly 6dB
            coupling_loss_dB = -4 # Estimate of coupling loss in dB

            # Gains
            gain_tx = ((math.pi * i)/wavelength)**2 * 0.81
            gain_rx = ((math.pi * j)/wavelength)**2

            # Convert gains to dB
            gain_tx_dB = convert_to_dB(gain_tx)
            gain_rx_dB = convert_to_dB(gain_rx)

            print("fspl:", fspl_dB)
            print("pointing loss:", pointing_loss_dB)
            print("inherent loss:", inherent_loss_dB)
            print("Tx gain:", gain_tx_dB)
            print("Rx gain:", gain_rx_dB)

            total_losses_dB[ix, jx] = fspl_dB + pointing_loss_dB + inherent_loss_dB + coupling_loss_dB + gain_rx_dB + gain_tx_dB

    contour_plot = plt.contour(D_TX, D_RX, total_losses_dB, levels=[-100, -95, -90, -85])
    plt.clabel(contour_plot, inline=True, fontsize=8, fmt="%.1f dB")
    plt.xlabel("Transmitter Diameter (m)")
    plt.ylabel("Receiver Diameter (m)")
    plt.title("Total Losses in dB")
    plt.show()

main()