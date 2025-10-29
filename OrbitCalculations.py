import scipy.optimize
from BigLinkBudget import *
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib.animation as animation
import scipy
from brokenaxes import brokenaxes

laser_power = 0.5 # Laser power in watts
wavelength = 1.54E-06 # Wavelength in metres
tx_aperture = 0.1 # Transmitter aperture diameter in metres
rx_aperture = 0.8 # Receiver aperture diameter in metres
orbit_altitude = 450000 # Orbital altitude above the Earth in metres
elevation_angle = 45 # Elevation angle from the Earth in degrees
tx_transmission = 1 # Optical transmission ratio of the transmitter (unitless)
rx_transmission = 0.8 # Optical transmission ratio of the receiver (unitless)
rx_obscuration = 0.42 # Obscuration ratio of the receiver (unitless)
rx_smf_il = -4 # Receiver SMF IL (dB) (look up later)
rx_switch_il = -1.2 # Receiver switch IL (dB) (look up later)
smf_coupling = -5 # Quick and dirty dB loss figure for SMF coupling

#Orbital parameters
radius = 6371000 #radius of earth in m
mu = 398600.435507*10E9 #Planetary mass of Earth in m^3/s^2. is equal to GM
eccentricity = 0.0002287
satellite_lowest_altitude = 413.866*1000 #lowest altitude of a satellite in m
satellite_lowest_altitude_km = satellite_lowest_altitude/1000
radius_perigee = radius+satellite_lowest_altitude
radius_perigee_km = radius_perigee/1000
print(radius_perigee_km)
true_anomaly = np.arange(0, 361, 0.5).tolist() #true anomaly of satellite in degrees

#Ground station location coordinates (km)
Ground_Station_Coords2 = np.array([4504.977303, 4504.977303])
Ground_Station_Coords =  np.array([-6371.0, 0]) #np.array([6371.0, 0])
Ground_Station_Coords3 = np.array([-6371.0, 0])
Ground_Station_x = Ground_Station_Coords[0]
Ground_Station_Y = Ground_Station_Coords[1]



#orbit calculations
def orbit_calculations(true_anomaly):
    angular_momentum = math.sqrt(radius_perigee*mu*(1+eccentricity))

    #Calcualting distance between earth centre and satellite (in m)
    distance = ((angular_momentum**2)/mu)*1/(1+eccentricity*math.cos(math.radians(true_anomaly)))
    apogee_distance = ((angular_momentum**2)/mu)*1/(1-eccentricity)
    altitude_apogee = apogee_distance-radius
    calc_sat_distance = distance-radius
    #print(calc_sat_distance)
    return distance, calc_sat_distance, apogee_distance

#Name: GroundStationAngle
#Purpose: To calculate the angle between the ground station (fixed point on Earth surface) to the satellite, as well as distance
#Input:
#Output:
def GroundStationAngle(true_anomaly, Ground_Station_Coords):
    #print("Ground station angle")
    Earth_Coords = np.array([0,0])
    distance, calc_sat_distance, apogee_distance = orbit_calculations(true_anomaly)
    altitude = distance-radius
    
    Sat_Coords_x = math.cos(math.radians(true_anomaly))*distance/1000
    Sat_Coords_y = math.sin(math.radians(true_anomaly))*distance/1000
    Sat_Coords = np.array([Sat_Coords_x, Sat_Coords_y])
    #Sat_Coords2 = np.array([Sat_Coords_y, Sat_Coords_x])

    #Create vectors
    Earth_Sat_Vector = np.subtract(Sat_Coords, Earth_Coords)
    Station_Sat_Vector = np.subtract(Sat_Coords,Ground_Station_Coords)

    #Straight line distance between Sat and Station
    Dist_Station_Sat = np.linalg.norm(Station_Sat_Vector)

    # if true_anomaly == 180 or true_anomaly==0:
    #     print("Dist_Station_Sat: ", Dist_Station_Sat)

    # Dot_Prod = np.dot(Earth_Sat_Vector, Station_Sat_Vector)
    # #Angle between two vectors, earth-sat and ground-sat
    # Angle_bw = math.acos(Dot_Prod/(np.linalg.norm(Earth_Sat_Vector)*Dist_Station_Sat))

    #Angle between satellite and station, as a function of the triangle
    # if true_anomaly > 180:
    #     degree_Value = 360-true_anomaly
    # else:
    #     degree_Value = true_anomaly
    # Angle_Sat_Ground = math.pi-Angle_bw-math.radians(degree_Value)
    # #Actual angle outside of triangle
    # # if true_anomaly > 180:
    # #     Real_Angle = math.degrees(math.pi + Angle_Sat_Ground)
    # # else:
    # #     Real_Angle = math.degrees(math.pi - Angle_Sat_Ground)
    # Real_Angle = (math.pi) - Angle_Sat_Ground
    # Real_Angle_rad = math.radians(Real_Angle)
    # #Angle from ground to satellite is measured from zenith, being 90 degrees. So we have to change it from 0 degrees to 90
    # Real_Angle2 = math.degrees(math.pi/2 - Real_Angle)

    '''
    Try 2
    '''
    Dot_Prod = np.dot(Station_Sat_Vector, Ground_Station_Coords)
    Angle_bw = math.acos(Dot_Prod/(np.linalg.norm(Station_Sat_Vector)*np.linalg.norm(Ground_Station_Coords)))
    Real_Angle = math.pi/2 - abs(Angle_bw)
    Real_Angle2 = math.degrees(Real_Angle)
    #print("REAL ANGLE 2: ", Real_Angle2)
    #print("Altitude at: ", true_anomaly,"is: ", altitude, "Real angle is: ", Real_Angle2, "dist ground-sat is: ", Dist_Station_Sat)
    
    #print(math.degrees(Angle_Sat_Ground))
    #print(str(true_anomaly) + "   "+str(Real_Angle))
    return Dist_Station_Sat, Real_Angle2, Sat_Coords, altitude

'''
Iterate through the orbit and calc distance for each degree of true anomaly
'''
def single_orbit_calc(Ground_Station_Coords):
    dB_values_orbit = []
    distance_satgnd = []
    distance_earth_sat = []
    elevation_array = []
    angle_array = []
    atmos_array = []
    for i, angle in enumerate((true_anomaly)):
        '''
        Now calculating from the ground station, using GroundStationAngle for the distance and angle
        '''
        distance, altitude, apogee_distance = orbit_calculations(angle)
        Dist_Station_Ground, Angle_Station_Ground, Sat_Coords, altitude = GroundStationAngle(angle, Ground_Station_Coords)
        
        distance_earth_sat.append(distance)
        elevation_array.append(Angle_Station_Ground)
        #print("Calc Angle: ", Angle_Station_Ground)

        '''
        Adjust the ground station angle (elevtation), to only range between 0 and 90. Currently it operates between -90 and 90
        '''
        if Angle_Station_Ground >= 5:
            angle_array.append(Angle_Station_Ground)
            distance_satgnd.append(Dist_Station_Ground)
            
            
            
            #don't need to convert to radians as the link budget code does that
            laser_dB = link_budget(laser_power, wavelength, tx_aperture, rx_aperture, altitude, Angle_Station_Ground, tx_transmission, rx_transmission, rx_obscuration, rx_smf_il, smf_coupling)

            sds = fnc_los_transmission(wavelength, 14/1000, Angle_Station_Ground)
            atmos_array.append(sds)
            print(altitude, ",", Angle_Station_Ground, ",", sds, ",", laser_dB)
            #laser_dB = link_budget(laser_power, wavelength, tx_aperture, rx_aperture, Dist_Station_Ground*1000, Angle_Station_Ground, tx_transmission, rx_transmission, rx_obscuration, rx_smf_il)
            #laser_dB = link_budget(laser_power, wavelength, tx_aperture, rx_aperture, altitude, elevation_angle, tx_transmission, rx_transmission, rx_obscuration, rx_smf_il)
        else:  
            laser_dB = None
            angle_array.append(None)
            distance_satgnd.append(None)
            atmos_array.append(None)
            #laser_dB = link_budget(laser_power, wavelength, tx_aperture, rx_aperture, altitude, Angle_Station_Ground, tx_transmission, rx_transmission, rx_obscuration, rx_smf_il, smf_coupling)


        dB_values_orbit.append(laser_dB)
        #print(str(angle) + "   "+str(laser_dB))
    
    #print(dB_values_orbit[0])
    #print("Altitude at apogee is (m): " + str(apogee_distance-radius))
    fig, ax1 = plt.subplots()
    ax1.plot(true_anomaly, dB_values_orbit)
    ax1.set_xlabel("True anomaly (degrees)")
    ax1.set_ylabel("Link margin (dB)")
    #plt.plot(true_anomaly, dB_values_orbit)
    #plt.plot(true_anomaly, distance_satgnd)
    #plt.ylabel("Link margin (dB)")
    #plt.xlabel("True anomaly (degrees)")
    #plt.title("eccentricity = " + str(eccentricity)+" , lowest alt = " + str(satellite_lowest_altitude/1000) + " km")
    plt.savefig("Eccentric_noatmos.png")
    #ax2 = ax1.twinx()

    # Plot data on the second y-axis (right)
    #ax2.plot(true_anomaly, distance_satgnd, 'r-')
    #ax2.set_ylabel('distance', color='r')
    #ax2.tick_params(axis='y', labelcolor='r')
    #fig.tight_layout() 

    plt.show()

    
    #plt.ylabel("Distance (km)")
    #plt.xlabel("True anomaly (degrees)")
    #plt.title("eccentricity = " + str(eccentricity)+" , lowest alt = " + str(satellite_lowest_altitude/1000) + " km")
    #plt.savefig("eccen="+str(eccentricity)+"_lowalt="+str(satellite_lowest_altitude/1000)+"_fixed.png")
    #plt.show()

    #plotting the distance vs angle for the ground station
    #Shows a parabola, but not a true one. Makes sense given 
    # plt.plot(true_anomaly, distance_satgnd)
    # plt.show()
    # v1 = np.gradient(distance_satgnd, true_anomaly)
    # a1 = np.gradient(v1, true_anomaly)
    # plt.figure(figsize=(10, 6))
    # #plt.plot(true_anomaly, distance_earth_sat, label='Distance (y)')
    # plt.plot(true_anomaly, v1, label='Velocity (v)', linestyle='--')
    # plt.plot(true_anomaly, a1, label='Acceleration (a)', linestyle=':')
    # plt.xlabel('Time (x)')
    # plt.ylabel('Value')
    # plt.legend()
    # plt.grid(True)
    # plt.title('Velocity and Acceleration vs True anomaly Sat-Gnd')
    # plt.show()

    # plt.plot(true_anomaly, elevation_array)
    # plt.show()
    # v2 = np.gradient(elevation_array, true_anomaly)
    # a2 = np.gradient(v2, true_anomaly)
    # plt.figure(figsize=(10, 6))
    # #plt.plot(true_anomaly, distance_earth_sat, label='Distance (y)')
    # plt.plot(true_anomaly, v2, label='Velocity (v)', linestyle='--')
    # plt.plot(true_anomaly, a2, label='Acceleration (a)', linestyle=':')
    # plt.xlabel('Time (x)')
    # plt.ylabel('Value')
    # plt.legend()
    # plt.grid(True)
    # plt.title('Velocity and Acceleration of elevation angle vs True anomaly ')
    # plt.show()

    

    # #plotting the distance vs angle for earth centre
    # v = np.gradient(distance_earth_sat, true_anomaly)
    # a = np.gradient(v, true_anomaly)
    # plt.figure(figsize=(10, 6))
    # #plt.plot(true_anomaly, distance_earth_sat, label='Distance (y)')
    # plt.plot(true_anomaly, v, label='Velocity (v)', linestyle='--')
    # plt.plot(true_anomaly, a, label='Acceleration (a)', linestyle=':')
    # plt.xlabel('Time (x)')
    # plt.ylabel('Value')
    # plt.legend()
    # plt.grid(True)
    # plt.title('Velocity and Acceleration vs True anomaly')
    # plt.show()

    # plt.plot(true_anomaly, distance_earth_sat)
    # plt.show()
    return dB_values_orbit


def flat_range_calc():
    orbit_range = np.arange(400*1000, 36000*1000, 1000).tolist()
    dB_values = []
    for i, value in enumerate((orbit_range)):
        laser_dB = link_budget(laser_power, wavelength, tx_aperture, rx_aperture, value, elevation_angle, tx_transmission, rx_transmission, rx_obscuration, rx_smf_il, smf_coupling)
        dB_values.append(laser_dB)
    #print(orbit_range)

    #print(dB_values[0])
    orbit_range_km = [x/1000 for x in orbit_range]
    plt.plot(orbit_range_km, dB_values)
    plt.title("Link margin vs Distance for " + str(convert_from_dBm(laser_power)) + "W laser")
    plt.ylabel("Link margin (dB)")
    plt.xlabel("Distance (km)")
    plt.show()

#laser_dB = link_budget(laser_power, wavelength, tx_aperture, rx_aperture, orbit_altitude, elevation_angle, tx_transmission, rx_transmission, rx_obscuration, rx_smf_il)
db_values = single_orbit_calc(Ground_Station_Coords)
db_values45 = single_orbit_calc(Ground_Station_Coords2)
db_valuesopp = single_orbit_calc(Ground_Station_Coords3)
new_list = [item for item in db_values if item is not None]
bax = brokenaxes(xlims=((0, 15), (345, 360)), hspace=0.05)
bax.plot(true_anomaly, db_values, label = "Station position A")
bax.set_xlabel('True anomaly (degrees)')
bax.set_ylabel('Link margin (dB)')
#plt.plot(new_list, label = "Station position A")
#plt.legend()
bax.legend()
plt.savefig("ISS_Broken.png")
#plt.show()
#plt.plot(true_anomaly, db_values45, label = "Station position B")
#plt.show()
#plt.plot(true_anomaly, db_valuesopp, label = "Station position C")

#plt.ylabel("Link margin (dB)")
#plt.xlabel("True anomaly (degrees)")
#plt.savefig("2D_ISS_OriginGS.png")

# plt.legend()
# #plt.title("eccentricity = " + str(eccentricity)+" , lowest alt = " + str(satellite_lowest_altitude/1000) + " km")
# plt.savefig("Eccentric_All_Gnd.png")
plt.show()
# flat_range_calc()

# laser_dB_mars = link_budget(laser_power, wavelength, tx_aperture, rx_aperture, 225000000*1000, elevation_angle, tx_transmission, rx_transmission, rx_obscuration, rx_smf_il)
# print("At Mars: " + str(laser_dB_mars))

# laser_dB_moon = link_budget(laser_power, wavelength, tx_aperture, rx_aperture, 385000*1000, elevation_angle, tx_transmission, rx_transmission, rx_obscuration, rx_smf_il)
# print("At Moon: " + str(laser_dB_moon))
satellite_trail_x = []
satellite_trail_y = []
for angle in true_anomaly:

    distance, calc_sat_distance, apogee_distance = orbit_calculations(angle)
    x = math.cos(math.radians(angle))*distance/1000
    satellite_trail_x.append(x)
    y = math.sin(math.radians(angle))*distance/1000
    satellite_trail_y.append(y)


def update(frame):
    #Updating the position of the satellite and adding the value to the trail seen in the animation
    # x,y = satellite_circle.center
    # distance, calc_sat_distance, apogee_distance = orbit_calculations(frame)
    # x = math.cos(math.radians(frame))*distance/1000
    # satellite_trail_x.append(x)
    # y = math.sin(math.radians(frame))*distance/1000
    # satellite_trail_y.append(y)

    #Updating
    satellite_circle.center = (satellite_trail_x[frame],satellite_trail_y[frame])
    
    #Updating the trailing line
    line2.set_xdata(satellite_trail_x[:frame])
    line2.set_ydata(satellite_trail_y[:frame])

    #Updating the satellite to ground station line
    Dist_Station_Ground, Angle_Station_Ground, Sat_Coords, altitude = GroundStationAngle(frame, Ground_Station_Coords)
    Station_Sat_Line.set_data([Ground_Station_x, Sat_Coords[0]], [Ground_Station_Y, Sat_Coords[1]])
    #Station_Sat_Line = ax.plot(Ground_Station_Coords, Sat_Coords)   
    

    #Updating the link budget annotation
    if db_values[frame] == None:
        db_value = 0
    else: db_value = db_values[frame]

    annotation.set_text("Link budget (dB): "+str(round(db_value,3)))

    return (satellite_circle,line2, annotation, Station_Sat_Line)

#Setting the initial values
distance, calc_sat_distance, apogee_distance = orbit_calculations(0)
apogee_dist_km = apogee_distance/1000
print("Apogee altitude (km): ", apogee_dist_km-(radius/1000))

# Dist_Station_Sat, Real_Angle, Sat_Coords = GroundStationAngle(1)
# print("DISTANCE: ", Dist_Station_Sat)
# print("REAL ANGLE: ", math.degrees(Real_Angle))

#Setting Earth at centre with radius in km
EarthCircle = plt.Circle((0, 0), 6371, color='royalblue')
#Initialising trail


#Setting the initial position of the satellite
init_sat_x = 6371+satellite_lowest_altitude/1000
satellite_trail_x.append(init_sat_x)
init_sat_y = 0
satellite_trail_y.append(init_sat_y)
satellite_circle = plt.Circle((init_sat_x, init_sat_y), 300, color='r')

station_circle = plt.Circle((Ground_Station_Coords[0], Ground_Station_Coords[1]), 200, color = "magenta")
station_circle2 = plt.Circle((Ground_Station_Coords2[0], Ground_Station_Coords2[1]), 500, color = "magenta")
station_circle3 = plt.Circle((Ground_Station_Coords3[0], Ground_Station_Coords3[1]), 500, color = "magenta")






#Initial position of the satellite. Doing it twice just because idk if it would work without one or other
def init():
    satellite_circle.center = (6371+satellite_lowest_altitude/1000, 0)
    #Actually drawing the satellite
    ax.add_patch(satellite_circle)
    return satellite_circle,



fig = plt.figure(figsize=(10, 6))

#Setting axis limits
ax = plt.axes(xlim=(1.5*min(satellite_trail_x), 1.5*max(satellite_trail_x)), ylim=(1.5*min(satellite_trail_y), 1.5*max(satellite_trail_y)))
#ax = plt.axes()
#Adding both the Earth and line in. No idea why the [0] is needed but it is
ax.add_patch(EarthCircle)
ax.add_patch(station_circle)
# plt.text(Ground_Station_Coords[0]+300, Ground_Station_Coords[1]+300, "A")
# ax.add_patch(station_circle2)
# plt.text(Ground_Station_Coords2[0]+300, Ground_Station_Coords2[1]+300, "B")
# ax.add_patch(station_circle3)
# plt.text(Ground_Station_Coords3[0]-1500, Ground_Station_Coords3[1]+300, "C")
line2 = ax.plot(satellite_trail_x[0], satellite_trail_y[0], lw = 2, linestyle = '--')[0]

#Create the line that connects the station to the satellite
Dist_Station_Ground, Angle_Station_Ground, Sat_Coords, altitude = GroundStationAngle(0, Ground_Station_Coords)
Station_Sat_Line = ax.plot([Ground_Station_x, Ground_Station_x], [Ground_Station_Y, Ground_Station_Y], lw = 2)[0]
#Station_Sat_Line = ax.plot([Ground_Station_x, Sat_Coords[0]], [Ground_Station_Y, Sat_Coords[1]])[0]

#Writing the link budget numbers to the top left of the figure
annotation = ax.annotate(
    '', 
    (0, 1),
    xytext=(4, -4),
    xycoords='axes fraction',
    textcoords='offset points',
    fontweight='bold',
    color='white',
    backgroundcolor='k',
    ha='left', va='top')

#Make the plot square for circle to appear circle
ax.set_aspect("equal", adjustable='box')


#Animation()
ani = animation.FuncAnimation(fig, update, init_func=init, frames=360, interval=20, blit=True)
#update(360 - 5)
plt.ylabel("Y (km)")
plt.xlabel("X (km)")
plt.show()
print("saving animation!")
FFwriter = animation.FFMpegWriter(fps=10)
#ani.save(filename="eccen="+str(eccentricity)+"_lowalt="+str(satellite_lowest_altitude/1000)+"_fixed.mp4", writer=FFwriter)
#ani.save(filename="eccen="+str(eccentricity)+"_lowalt="+str(satellite_lowest_altitude/1000)+"_fixed.gif", writer="pillow")
print("Done!")
plt.close()
