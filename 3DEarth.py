#code taken from https://gist.github.com/ravi4ram/a3b22fb69d966c6dfd0b2d6638bae7e2 

# Plot 3D orbit from TLE data

import numpy as np
import math
from sgp4 import api
from sgp4.conveniences import sat_epoch_datetime
import matplotlib.animation as animation
from BigLinkBudget import *

import matplotlib.pyplot as plt

laser_power = convert_to_dBm(0.5) # Laser power in watts
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

# TLE
nm = '2023-132A     '          
sat1_l1 = '1 25544U 98067A   25194.82082520  .00010000  00000+0  17963-3 0  9997'
sat1_l2 = '2 25544  51.6335 175.0526 0002287   7.9305 352.1726 15.50527943519299'

#TLE for second satellite
nm = '2023-132A     '          
sat2_l1 = '1 25544U 98067A   25212.99608426  .00012187  00000+0  21918-3 0  9998'
sat2_l2 = '2 25544  51.6364  84.9763 0002287 136.8880 197.5582 15.50259437522105'

Station_Coords = [6378.135,0,0]
Station_Coords = [-6354.282219769004, 551.0914877540154, -1.3122033207603634]

# get satellite info
def get_satellite(l1, l2):
    global mu, r, a
    satellite = api.Satrec.twoline2rv(l1, l2, api.WGS72)

    # constants
    mu = satellite.mu           # Earth’s gravitational parameter (km³/s²)
    r = satellite.radiusearthkm # Radius of the earth (km).
    print(r)
    # orbital parameters
    a = satellite.a * r
    # Altitude of the satellite at perigee 
    peri = satellite.altp * r
    # Altitude of the satellite at apogee
    apo  = satellite.alta * r
    print('a X e        : ', round(apo, 3), round(peri, 3) )   
    # Inclination
    i_deg = np.degrees(satellite.inclo)
    print('Inclination  : ', i_deg)
    # Eccentricity
    e = satellite.ecco
    print('Eccentricity : ', e)
    epoch = sat_epoch_datetime(satellite)
    print('Epoch        : ', epoch)
    #
    return satellite, epoch

# propagation starts at epoch and ends at full orbit time period
# function returns an array of datetime from start to end
def get_propagation_times(epoch):
    # DeprecationWarning: parsing timezone aware datetimes is deprecated;
    epoch = epoch.replace(tzinfo=None)
    
    # Epoch (start time)
    start = np.datetime64(epoch)

    # Period (end time)
    t = 2 * np.pi * (a**3/mu)**0.5
    # end time for full orbit    
    end = t
    
    # time steps
    dt = 100

    # create time array
    time_end = np.arange(0.0, end, dt)
    time_arr = start + time_end.astype('timedelta64[s]')
    # 
    return time_arr

# return state vectors for each datetime in the time array
def get_state_vectors(satellite, time_arr):
    position = []; velocity = [];
    for j in time_arr.tolist():
        jd, fr = api.jday(j.year, j.month, j.day, j.hour, j.minute, j.second)
        e, p, v = satellite.sgp4(jd, fr)

        #If e=/= that means there has been error
        if e==0:
            position.append(p); velocity.append(v)
        else:
            print(p); print(v)

    # slice into columns
    pos, vel   = list(zip(*position)), list(zip(*velocity))
    X, Y, Z    = np.array(pos[0]), np.array(pos[1]), np.array(pos[2])
    VX, VY, VZ = np.array(vel[0]), np.array(vel[1]), np.array(vel[2])
    state_vectors = [X, Y, Z, VX, VY, VZ]
    #
    print(X[0],Y[0],Z[0])
    return state_vectors

# plot
# Set 3D plot axes to equal scale. 
# Required since `ax.axis('equal')` and `ax.set_aspect('equal')` don't work on 3D.
def set_axes_equal_3d(ax: plt.Axes):
    """	
    https://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
    """
    # set axis limits
    def _set_axes_radius(ax, origin, radius):
        x, y, z = origin
        ax.set_xlim3d([x - radius, x + radius])
        ax.set_ylim3d([y - radius, y + radius])
        ax.set_zlim3d([z - radius, z + radius])
        return

    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])
    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    _set_axes_radius(ax, origin, radius)
    return



# plot orbit
def plot_xyz(state_vectors, r, Animation, state_vectors2):
    global fig, ax, orbit, satellite, Earth
    global X, Y, Z
    global u, v
    global Station_Sat_Line

    global X2, Y2, Z2
    global orbit2
    
    X, Y, Z = state_vectors[0], state_vectors[1], state_vectors[2]
    X2, Y2, Z2 = state_vectors2[0], state_vectors2[1], state_vectors2[2]
    
    fig = plt.figure()
    #ax = plt3.Axes3D(fig)
    ax = plt.axes(projection='3d')
    #ax.set_facecolor('black')

    # set labels
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    # set correct aspect ratio
    ax.set_box_aspect([1,1,1])
    set_axes_equal_3d(ax)
    # set view
    azm=150; ele=36; #39
    ax.view_init(elev=ele, azim=azm)
    # set limit
    size = 1.0
    #limit = max(max(X[0]), max(Y[0]), max(Z[0]))
    limit = max(max(X), max(Y), max(Z))
    ax.set_xlim(size*limit, -size*limit)
    ax.set_ylim(size*limit, -size*limit)
    ax.set_zlim(size*limit, -size*limit) 

    # earth
    #ax.scatter(0, 0, 0, marker='o', c='deepskyblue', s=3*r, alpha=0.5)

    u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:50j]
    x = r*np.cos(u)*np.sin(v)
    y = r*np.sin(u)*np.sin(v)
    z = r*np.cos(v)
    Earth = ax.plot_surface(x, y, z, alpha = 0.3)

    # Station
    ax.scatter(Station_Coords[0], Station_Coords[1], Station_Coords[2], c='r', s=100)

    #ax.scatter(-6354.282219769004, 551.0914877540154, -1.3122033207603634, c='r', s=100)


    #If statement for either plotting animation or static picture
    if Animation == "yes":
        print("Animate")
        #Adding the trailing line
        orbit = ax.plot(X[0], Y[0], Z[0], linewidth=0.9, c='k', linestyle='dashed')[0]
        orbit2 = ax.plot(X2[0], Y2[0], Z2[0], linewidth=0.9, c='k', linestyle='dashed')[0]
        #Adding the line between the station and sat
        Station_Sat_Line, = ax.plot([Station_Coords[0], Station_Coords[0]], [Station_Coords[1], Station_Coords[1]], [Station_Coords[2], Station_Coords[2]], c='r')
        ani = animation.FuncAnimation(fig, update,init_func=init ,frames=len(X), interval=20, blit=True, repeat=False)
        plt.show()
        print("saving animation!")
        ani.save(filename="Full3d.gif", writer='pillow')
        print("Done!")
        plt.close()
    else:

        # plot()
        orbit = ax.plot(X, Y, Z, linewidth=0.9, linestyle='dashed', c='k')
    
        # show
        plt.show()

    # finish
    return

def init():
    global scat
    global scat2
    # x2 = X[0]+1000*np.cos(u)*np.sin(v)
    # y2 = Y[0] +1000*np.sin(u)*np.sin(v)
    # z2 = Z[0]+1000*np.cos(v)
    # scat = ax.plot_surface(x2,y2,z2, alpha = 0.5)

    #Start point of satellite orbit
    #scat = ax.scatter(X[0], Y[0], Z[0], c='g', s=100)
    '''
    You NEED the commas here and in the return because it expects a sequence
    '''
    scat = ax.plot([], [], [], marker = 'o', ms=10, c='g')[0]
    scat2 = ax.plot([], [], [], marker = 'o', ms=10, c='g')[0]
    return scat, scat2,

def update(frame):
    if frame == len(X) - 1:
        scat.set_data([], [])
        scat.set_3d_properties([])
        scat2.set_data([], [])
        scat2.set_3d_properties([])
        orbit.set_data([], [])
        orbit.set_3d_properties([])
        orbit2.set_data([], [])
        orbit2.set_3d_properties([])
        Station_Sat_Line.set_data([], [])
        Station_Sat_Line.set_3d_properties([])
        #Earth.set_alpha(0)
        return scat, orbit, Station_Sat_Line, scat2, orbit2, Earth
    
    #print("update")
    XCoord, YCoord, ZCoord = X[frame], Y[frame], Z[frame]
    #scat, = ax.plot(XCoord, YCoord, ZCoord, marker = 'o', ms=10, c='g')
    scat.set_data([XCoord], [YCoord])
    scat.set_3d_properties([ZCoord])

    #updating position of second orbit
    XCoord2, YCoord2, ZCoord2 = X2[frame], Y2[frame], Z2[frame]
    #scat, = ax.plot(XCoord, YCoord, ZCoord, marker = 'o', ms=10, c='g')
    scat2.set_data([XCoord2], [YCoord2])
    scat2.set_3d_properties([ZCoord2])

    #Updating the trailing line 
    #orbit, = ax.plot(X[:frame], Y[:frame], Z[:frame], linewidth=0.9, linestyle='dashed', c='k')
    orbit.set_data(X[:frame], Y[:frame])
    orbit.set_3d_properties(Z[:frame])

    #Updating the trailing line for second sat
    #orbit, = ax.plot(X[:frame], Y[:frame], Z[:frame], linewidth=0.9, linestyle='dashed', c='k')
    orbit2.set_data(X2[:frame], Y2[:frame])
    orbit2.set_3d_properties(Z2[:frame])

    #Updating the station_sat line
    #Station_Sat_Line, = ax.plot([Station_Coords[0], XCoord], [Station_Coords[1], YCoord], [Station_Coords[2],ZCoord], c='r')
    Station_Sat_Line.set_data([Station_Coords[0], XCoord], [Station_Coords[1], YCoord])
    Station_Sat_Line.set_3d_properties([Station_Coords[2],ZCoord])
    return (scat, orbit, Station_Sat_Line, orbit2, scat2, Earth)

def CalcElevationAngle(StateVectors):
    #link budget calc
    link_budget_list = []
    #SatXYZ is lists of each coord for each time

    SatX, SatY, SatZ = StateVectors[0], StateVectors[1], StateVectors[2]
    print(SatX[0])
    StationCoords = np.array([Station_Coords[0], Station_Coords[1], Station_Coords[2]])
    ElevationAngle = []

    #calculating the vector between sat and ground station
    for j, value  in enumerate(SatX):
        Vector_Sat_Stn = np.array([SatX[j]-StationCoords[0], SatY[j]-StationCoords[1], SatZ[j]-StationCoords[2]])
        Dot_Prod = np.dot(Vector_Sat_Stn, StationCoords)
        Angle_between = math.acos(Dot_Prod/(np.linalg.norm(Vector_Sat_Stn)*np.linalg.norm(StationCoords)))
        #angle from 90 degree elevation angle
        
        Real_angle = math.pi/2-abs(Angle_between) #Elevation angle in radians
        Real_angle = math.degrees(Real_angle) #convert to degrees
        
        if Real_angle >= 0:
            #print('e')
            Magnitude_vector_sat_earth = np.linalg.norm([SatX[j], SatY[j], SatZ[j]])
            sat_altitude = Magnitude_vector_sat_earth-r
            #print(sat_altitude)
            laser_dB = link_budget(laser_power, wavelength, tx_aperture, rx_aperture, sat_altitude*1000, Real_angle, tx_transmission, rx_transmission, rx_obscuration, rx_smf_il)
        else:
            laser_dB = None
        link_budget_list.append(laser_dB)
        ElevationAngle.append(Real_angle)

    
    return ElevationAngle, link_budget_list

def plot_link_budget(link_list, time_arr):
    plt.plot(time_arr, link_list)
    plt.show()



def main():
    satellite, epoch = get_satellite(sat1_l1, sat1_l2)
    time_arr         = get_propagation_times(epoch)
    state_vectors    = get_state_vectors(satellite, time_arr)
    Elevation_Angle, link_bud = CalcElevationAngle(state_vectors)

    satellite2, epoch2 = get_satellite(sat2_l1, sat2_l2)
    time_arr2         = get_propagation_times(epoch2)
    state_vectors2    = get_state_vectors(satellite2, time_arr2)
    #Elevation_Angle2, link_bud2 = CalcElevationAngle(state_vectors2)
    
    print(Elevation_Angle)
    plot_xyz(state_vectors, r, "yes", state_vectors2)
    plot_link_budget(link_bud, time_arr)
    return

# main
if __name__ == '__main__':
    main()