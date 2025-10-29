# Title: Big Link Budget (BLB)
# Author: Jackie Scholes
# Date: 26/11/2024
# Description: Program that does a real proper link budget based on an independent input config file

# Notes
# [S] denotes equation sourced from the spreadsheet
# This is assuming an Earth-to-space configuration. Certain parameters like elevation angle are assuming Earth-like atmospheric fading.

import math
import sys
import numpy as np
#from corelib import *

# Temporary input parameter names

#Name: convert_to_dBm
#Def: Converting laser power in watts to dBm
#Input: Laser power
#Output: Laser power in decibel milliwatts
def convert_to_dBm(laser_power):
    #The 1000 is to convert watts to milliwatts
    dbm = 10*math.log10((laser_power*1000)/1)
    return dbm

#For range loss, converts watts to decibel watts
def convert_to_dbW(laser_power):
    dbW = 10*math.log10((laser_power))
    #print("Laser power: ", dbW)
    return dbW

#laser_power = convert_to_dBm(0.5) # Laser power in watts
# laser_power = 0.5 # Laser power in watts
# wavelength = 1.54E-06 # Wavelength in metres
# tx_aperture = 0.1 # Transmitter aperture diameter in metres


# rx_aperture = 0.8 # Receiver aperture diameter in metres
# orbit_altitude = 1000374.2319998629 # Orbital altitude above the Earth in metres
'''
If elevation angle = 0, that is horizon. If angle = 90, that is straight overhead
elevation_angle = h, h = 90 - zenith
'''
# elevation_angle = 0 # Elevation angle from the Earth in degrees
# tx_transmission = 1 # Optical transmission ratio of the transmitter (unitless)
# rx_transmission = 0.8 # Optical transmission ratio of the receiver (unitless)
# rx_obscuration = 0.42 # Obscuration ratio of the receiver (unitless)
# rx_smf_il = -4 # Receiver SMF IL (dB) (look up later)
# rx_switch_il = -1.2 # Receiver switch IL (dB) (look up later)
# smf_coupling = -5 # Quick and dirty dB loss figure for SMF coupling

# Constants

#c = 3E8 # Speed of light in metres per second
c = 299792458 #New speed of light in m/s
planck = 6.62607015E-34 # Planck's constant in J/Hz
radius_earth = 6371*1000 # Radius of Earth in kilometres
radius_moon = 1737 # Radius of the Moon in kilometres
zenith = 0.857 # Zenith transmittance ratio (unitless) @1550nm (might need to research to make this variable)

# Input parameters (to eventually be passed from input file)

tx_planet = "earth"

def link_budget(laser_power, wavelength, tx_aperture, rx_aperture, orbit_altitude, elevation_angle, tx_transmission, rx_transmission, rx_obscuration, rx_smf_il, smf_coupling):
    h = 14/1000 #height above mean sea level (km)
    laser_power = convert_to_dBm(laser_power)
    beam_divergence, tx_gain = fnc_tx_gain(wavelength, tx_aperture)
    fspl = fnc_fspl(wavelength, elevation_angle, orbit_altitude)
    los_transmission = fnc_los_transmission(wavelength, h,elevation_angle)
    rx_gain = fnc_rx_gain(rx_aperture, wavelength)
    rx_obscuration_loss = fnc_rx_obscuration_loss(rx_obscuration)
    rx_absorption = fnc_rx_absorption(rx_transmission)

    
    link_budget_output = laser_power + tx_gain + fspl + rx_gain + rx_obscuration_loss + rx_absorption + rx_smf_il + smf_coupling

    #link_budget_output = laser_power + tx_gain + fspl + los_transmission + rx_gain + rx_obscuration_loss + rx_absorption + rx_smf_il + smf_coupling
    #print("LINK BUDGET", link_budget_output, "dB")
    return(link_budget_output)

#This is for intersatellite
def link_budget2(laser_power, wavelength, tx_aperture, rx_aperture, link_distance, tx_transmission, rx_transmission, rx_obscuration, rx_smf_il, smf_coupling):

    laser_power = convert_to_dBm(laser_power)
    beam_divergence, tx_gain = fnc_tx_gain(wavelength, tx_aperture)
    fspl = fnc_fspl2(wavelength, link_distance)
    #los_transmission = fnc_los_transmission(elevation_angle)
    rx_gain = fnc_rx_gain(rx_aperture, wavelength)
    rx_obscuration_loss = fnc_rx_obscuration_loss(rx_obscuration)
    rx_absorption = fnc_rx_absorption(rx_transmission)


    link_budget_output = laser_power + tx_gain + fspl + rx_gain + rx_obscuration_loss + rx_absorption + rx_smf_il + smf_coupling
    #print("LINK BUDGET", link_budget_output, "dB")
    return(link_budget_output)

#for range loss in db/m^2
def link_budget3(laser_power, wavelength, tx_aperture, rx_aperture, orbit_altitude, elevation_angle, tx_transmission, rx_transmission, rx_obscuration, rx_smf_il, smf_coupling):
    h = 2286/1000 #height above mean sea level (km)
    laser_power2 = convert_to_dbW(laser_power)
    beam_divergence, tx_gain = fnc_tx_gain(wavelength, tx_aperture)
    fspl = fnc_fspl3(laser_power, elevation_angle, orbit_altitude)
    los_transmission = fnc_los_transmission(wavelength, h,elevation_angle)
    #rx_gain = fnc_rx_gain(rx_aperture, wavelength)
    #rx_obscuration_loss = fnc_rx_obscuration_loss(rx_obscuration)
    #rx_absorption = fnc_rx_absorption(rx_transmission)

    
    link_budget_output = laser_power2 + tx_gain + fspl + rx_smf_il + smf_coupling

    #link_budget_output = laser_power + tx_gain + fspl + los_transmission + rx_gain + rx_obscuration_loss + rx_absorption + rx_smf_il + smf_coupling
    #print("LINK BUDGET", link_budget_output, "dB")
    return(link_budget_output)

def fnc_tx_gain(wavelength, tx_aperture):

    # [S] Antenna gain equation based on Gaussian beam divergence
    '''
    https://www.globalspec.com/reference/72958/203279/laser-beam-divergence
    https://en.wikipedia.org/wiki/Gaussian_beam#Beam_waist - 2*theta
    https://www.edmundoptics.com.au/knowledge-center/application-notes/lasers/gaussian-beam-propagation/ - diagram is wrong, should
    be capital theta
    '''

    '''
    Should TX aperture be /2 for radius instead of diameter??
    '''
    beam_divergence = 2 * wavelength / (math.pi * tx_aperture)
    tx_gain = 32/(beam_divergence**2)
    tx_gain = convert_to_dB(tx_gain)

    dviergence_angle_gain = 16/((380*10**(-6))**2)
    dviergence_angle_gain = convert_to_dB(dviergence_angle_gain)
    #print("tx_gain:", dviergence_angle_gain)
    #return beam_divergence, dviergence_angle_gain
    return beam_divergence, tx_gain

def fnc_link_distance(elevation_angle, orbit_altitude):

    # [S] Calculating link distance based on several input parameters

    # Determine radius of input celestial body
    if tx_planet == "earth":
        radius = radius_earth
    elif tx_planet == "moon":
        radius == radius_moon
    else:
        print("ERROR: Unsupported planet of origin")
        sys.exit()

    link_distance = -1 * radius * math.sin(math.radians(elevation_angle)) + math.sqrt(radius**2 * (math.sin(math.radians(elevation_angle)))**2 + orbit_altitude**2 + 2 * orbit_altitude * radius)
    #print("link_distance:", link_distance)
    #if elevation_angle ==180 or elevation_angle == 0:
        #print("link distance: ", link_distance)
    return link_distance


def fnc_fspl(wavelength, elevation_angle, orbit_altitude):

    link_distance = fnc_link_distance(elevation_angle, orbit_altitude) - 2286

    # [S] Isotropic free space path loss equation
    fspl = ((4 * math.pi * link_distance) / wavelength)**2
    fspl = -1 * convert_to_dB(fspl)
    #print("fspl:", fspl)
    return fspl

#range loss in units of dB/m^2
def fnc_fspl3(laser_power, elevation_angle, orbit_altitude):

    link_distance = fnc_link_distance(elevation_angle, orbit_altitude) - 2286
    #link_distance = 930*10**3

    # [S] Isotropic free space path loss equation
    fspl = ((4 * math.pi * link_distance**2) / laser_power)
    fspl = -1 * convert_to_dB(fspl)
    print("range loss:", fspl)
    return fspl

#This is for intersatellite links
#Change: Link_distance instead of orbit altitude - has to be in km
def fnc_fspl2(wavelength, link_distance):

    #link_distance = fnc_link_distance(elevation_angle, orbit_altitude)

    # [S] Isotropic free space path loss equation
    fspl = ((4 * math.pi * link_distance) / wavelength)**2
    fspl = -1 * convert_to_dB(fspl)
    #print("fspl:", fspl)
    return fspl

def fnc_relative_airmass(elevation_angle):

    # [S] Relative airmass equation
    theta = (elevation_angle + (244 / (165 + 47 * (elevation_angle**1.1)))) * math.pi / 180
    relative_airmass = 1 / math.sin(theta)
    #print("relative_airmass", relative_airmass)
    return relative_airmass

def fnc_los_transmission(wavelength, h, elevation_angle):
    wavelength = wavelength*10**6 #convert wavelength from m to um
    #h = height above mean sea level (km)
    a = 0.000487*wavelength**3 - 0.002237*wavelength**2+0.003864*wavelength-0.004442
    b = -0.00573*wavelength**3 + 0.02639*wavelength**2 - 0.04552*wavelength + 0.05164
    c = 0.02565*wavelength**3 - .1191*wavelength**2 + 0.20385*wavelength - .216
    d = -0.0638*wavelength**3 + 0.3034*wavelength**2 - 0.5083*wavelength + 0.425

    tau = a*h**3 + b*h**2 + c*h + d
    data2 = []
    dataliang = []
    As = -(4.3429*tau)/(math.sin(math.radians(elevation_angle)))

    # [S] Line-of-sight atmospheric transmission equation
    los_transmission = math.e**(fnc_relative_airmass(elevation_angle) * math.log(zenith))
    los_transmission = convert_to_dB(los_transmission)
    #print("los_transmission:", los_transmission)
    #print("ITU: ", As)
    return As

def fnc_rx_gain(rx_aperture, wavelength):

    # [S] Receiver antenna gain equation
    rx_gain = ((math.pi * rx_aperture) / wavelength)**2
    rx_gain = convert_to_dB(rx_gain)
    #print("rx_gain:", rx_gain)
    return rx_gain

def fnc_photon_energy(wavelength):

    # [S] Photon energy equation
    photon_energy = planck * c / wavelength
    return photon_energy

def fnc_rx_obscuration_loss(rx_obscuration):

    # [S] Receiver obscuration loss equation
    rx_obscuration_loss =  convert_to_dB(1 - rx_obscuration**2)
    #print("rx_obscuration_loss", rx_obscuration_loss)
    return rx_obscuration_loss

def fnc_rx_absorption(rx_transmission):

    # [S] Receiver absorption loss equation
    rx_absorption = convert_to_dB(rx_transmission)
    #print("rx_absorption:", rx_absorption)
    return rx_absorption



#Name: convert to dB
#Def: Converting the input gain into dB
#Input: Gain
#Output: Gain in dB
def convert_to_dB(input):
    dB = 10*math.log10(input)
    return dB

#Name: Convert from dbm
#Def converts the input from dbm into watts
#Input: dbm
#Output: watts
def convert_from_dBm(input):
    watt = (10**(input/10))/1000
    return watt

#I think it converts dbm back into gain? idk really
def convert_from_dB(input):
    output = 10**(input/10)
    return output


#link_budget(laser_power, wavelength, tx_aperture, rx_aperture, orbit_altitude, elevation_angle, tx_transmission, rx_transmission, rx_obscuration, rx_smf_il, smf_coupling)