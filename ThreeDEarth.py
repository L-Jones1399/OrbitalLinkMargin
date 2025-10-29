#code taken from https://gist.github.com/ravi4ram/a3b22fb69d966c6dfd0b2d6638bae7e2 

# Plot 3D orbit from TLE data

import numpy as np
import math
from sgp4 import api
from sgp4.conveniences import sat_epoch_datetime
import matplotlib.animation as animation
from BigLinkBudget import *
from matplotlib.figure import Figure
from astropy.time import Time
from astropy.coordinates import TEME, ITRS, CartesianRepresentation, CartesianDifferential
import astropy.units as u
from datetime import time, datetime
import pandas as pd
import matplotlib.dates as mdates


import matplotlib.pyplot as plt

dt = 10

#laser_power = convert_to_dBm(0.5) # Laser power in watts
laser_power = 1.0 # Laser power in watts
wavelength = 1.54E-06 # Wavelength in metres
tx_aperture = 0.4 # Transmitter aperture diameter in metres
rx_aperture = 0.4 # Receiver aperture diameter in metres
orbit_altitude = 450000 # Orbital altitude above the Earth in metres
elevation_angle = 45 # Elevation angle from the Earth in degrees
tx_transmission = 1 # Optical transmission ratio of the transmitter (unitless)
rx_transmission = 0.8 # Optical transmission ratio of the receiver (unitless)
rx_obscuration = 0.42 # Obscuration ratio of the receiver (unitless)
rx_smf_il = -4.0 # Receiver SMF IL (dB) (look up later)
rx_switch_il = -1.2 # Receiver switch IL (dB) (look up later)
smf_coupling = -5 # Quick and dirty dB loss figure for SMF coupling

# TLE
nm1 = '2023-132A     '          
#sat1_l1 = '1 25544U 98067A   25194.82082520  .00010000  00000+0  17963-3 0  9997'
#sat1_l2 = '2 25544  51.6335 175.0526 0002287   7.9305 352.1726 15.50527943519299'

sat1_l1 = '1 42719U 17027A   25287.72824329  .00000627  00000+0 -13342-1 0  9990' #cosmos 2518 high eccentricity - tundra/molinya
sat1_l2 = '2 42719  63.0691 279.0723 7105107 273.0854  14.9814  2.00754496 61472'

sat1_l1 = '1 29268U 06031A   25286.65870861  .00000515  00000+0  10580-3 0  9998' #Kompsat-2 sun synch
sat1_l2 = '2 29268  97.8323 118.0507 0015436 132.2214 317.8777 14.64413986 25392'

#Geostationary
sat1_l1 = '1 60133U 24119A   25251.86105140 -.00000251  00000+0  00000+0 0  9997'
sat1_l2 = '2 60133   0.0206  21.1929 0000552 233.9644 327.7985  1.00272146  4143'

#sat1_l1 = '1 52734U 22057C   22333.28670189  .00017291  00000-0  87416-3 0  9990'
#sat1_l2 = '2 52734  97.5289  85.9900 0007263 332.7546  27.3305 15.17097234 28371' #TBIRD Schlier sat

#sat1_l1 = '1 48859U 21054A   25289.21780061 -.00000094  00000+0  00000+0 0  9999'
#sat1_l2 = '2 48859  55.2757 342.2645 0024064 230.1083 138.8203  2.00556506 31868' #USA-319 Navstar 81 GPS

#sat1_l1 = "1 26677U 80050G   25226.43619803  .00000319  00000-0  00000+0 0  9994"
#sat1_l2 = "2 26677  72.4979 193.0344 6001883 168.2468   2.1949  1.93165565223152" #COSMOS 1188 DEB

#sat1_l1 = '1 46826U 20078A   25289.20312859 -.00000051  00000-0  00000+0 0  9995'
#sat1_l2 = "2 46826  54.0306 220.6753 0063294 203.5880 326.5084  2.00567133 36612" #Navstar80

#TLE for second satellite
nm2 = '2023-132A     '          
sat2_l1 = '1 25544U 98067A   25212.99608426  .00012187  00000+0  21918-3 0  9998'
sat2_l2 = '2 25544  51.6364  84.9763 0002287 136.8880 197.5582 15.50259437522105'
sat2_l1 = '1 52734U 22057C   22333.28670189  .00017291  00000-0  87416-3 0  9990'
sat2_l2 = '2 52734  97.5289  85.9900 0007263 332.7546  27.3305 15.17097234 28371' #TBIRD Schlier sat

sat2_l1 = '1 29268U 06031A   25286.65870861  .00000515  00000+0  10580-3 0  9998' #Kompsat-2 sun synch
sat2_l2 = '2 29268  97.8323 118.0507 0015436 132.2214 317.8777 14.64413986 25392'


sat2_l1 = '1 46826U 20078A   25289.20312859 -.00000051  00000-0  00000+0 0  9995'
sat2_l2 = "2 46826  54.0306 220.6753 0063294 203.5880 326.5084  2.00567133 36612" #Navstar80

sat2_l1 = '1 25544U 98067A   25194.82082520  .00010000  00000+0  17963-3 0  9997' #ISS
sat2_l2 = '2 25544  51.6335 175.0526 0002287   7.9305 352.1726 15.50527943519299'

sat2_l1 = '1 42719U 17027A   25287.72824329  .00000627  00000+0 -13342-1 0  9990' #cosmos 2518 high eccentricity - tundra/molinya
sat2_l2 = '2 42719  63.0691 279.0723 7105107 273.0854  14.9814  2.00754496 61472'

#Geostationary
#sat2_l1 = '1 60133U 24119A   25251.86105140 -.00000251  00000+0  00000+0 0  9997'
#sat2_l2 = '2 60133   0.0206  21.1929 0000552 233.9644 327.7985  1.00272146  4143'

# Station_Coords = [6378.135,0,0]
# Station_Coords = [3892.1565149716216, -5052.892345830759, -1.63787480995528]
# Station_Coords = [-4135.262609147723, 3995.877149802852, -2759.198259511845]
# Station_Coords = [4009.9342629525772, -1871.3894866971962, 4593.357668863071] #Middle of ISS orbit
# Station_Coords = [-544.911416761162, 677.8825887146957, 6318.556229236691] #Komsat2
# Station_Coords = [-2764.97993893081, -5688.2397108916075, 824.2699842713819] #Navstar 81
# Station_Coords = [-327.02619092, -6352.94086755,  -462.38758855] #Navstar81 rotated 15
# Station_Coords = [-1873.4920326059193, -6096.771504987369, -3.3012060208078227] #Cosmos 1188 deb
# Station_Coords = [5730.427334788963, -2541.6242313754915, -1175.9910311791712] #halfway through cosmos
# Station_Coords = [2817.3500753343774, -5647.9164424912915, 918.79512943582] #Navstar80
Station_Coords = [5634.03806084984, 2865.985704785897, -851.0858635197549] #halfwaynavstar80

#Station_Coords = [-4669.370227796119, -4344.834601643922, 0.19584229410561427] #Geostationary GOES satellite
#Station_Coords = [1629.5500629971712, -6166.455435422829, 0.18263326814885214] #Geostationary take 2
#Station_Coords = [2935.9828662193718, 411.5843391374828, 5647.230207745074] #halfway through Molniya orbit
#Station_Coords = [-2448866.296/1000,-4667978.099/1000, 3582785.595/1000] #Schiler2023 Table Mountain observatory

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
    #end = 3*60*60
    
    # time steps
    #dt = 10
    # create time array
    time_end = np.arange(0.0, end, dt)
    time_arr = start + time_end.astype('timedelta64[s]')
    # 
    return time_arr, time_end, t

# For 2 satellite systems to ensure the animation is the same length
def get_propagation_times2(epoch, t):
    # DeprecationWarning: parsing timezone aware datetimes is deprecated;
    epoch = epoch.replace(tzinfo=None)
    
    # Epoch (start time)
    start = np.datetime64(epoch)

    # Period (end time)
    #t = 2 * np.pi * (a**3/mu)**0.5
    # end time for full orbit    
    end = t
    #end = 3*60*60
    
    # time steps
    #dt = 10

    # create time array
    time_end = np.arange(0.0, end, dt)
    time_arr = start + time_end.astype('timedelta64[s]')
    # 
    return time_arr, time_end

# return state vectors for each datetime in the time array
def get_state_vectors(satellite, time_arr):
    time_arr2 = []
    position = []; velocity = [];
    for j in time_arr.tolist():
        jd, fr = api.jday(j.year, j.month, j.day, j.hour, j.minute, j.second)
        e, p, v = satellite.sgp4(jd, fr)
        t = Time(jd+fr, format='jd')
        start = time(8, 33)
        end = time(8,41)
        dt = t.to_datetime()
        target_day = datetime(2022, 11, 29).date()
        # if not (j.date() == target_day and  start <= j.time() <= end):
        #     #print("nay")
        #     continue
        # else:
        #     print("yay")
        #     time_arr2.append(j)

        #If e=/= that means there has been error
        if e==0:
            teme_p = CartesianRepresentation(p*u.km)
            teme_v = CartesianDifferential(v*u.km/u.s)
            teme = TEME(teme_p.with_differentials(teme_v), obstime=t)
            itrs_geo = teme.transform_to(ITRS(obstime=t))
            p = [itrs_geo.x.value, itrs_geo.y.value, itrs_geo.z.value]
            #position.append(itrs_geo.x.km)
            v = [itrs_geo.v_x.value, itrs_geo.v_y.value, itrs_geo.v_z.value]
            position.append(p); velocity.append(v)
        else:
            print(p); print(v)

    # slice into columns
    pos, vel   = list(zip(*position)), list(zip(*velocity))
    X, Y, Z    = np.array(pos[0]), np.array(pos[1]), np.array(pos[2])
    VX, VY, VZ = np.array(vel[0]), np.array(vel[1]), np.array(vel[2])

    #In reference frame TEME - True Equator Mean Equinox
    state_vectors = [X, Y, Z, VX, VY, VZ]
    state_vectors2 = [X, Y, Z]
    #
    print("Starting Coords: ", X[0],Y[0],Z[0])
    print("Middle Coords: ", X[int(len(X)/2)], Y[int(len(Y)/2)], Z[int(len(Z)/2)])
    #print(state_vectors2)
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


#Calculate the position of the ground station as the Earth rotates
def GroundRotation(time_end):
    global Ground_Station_Array
    #print("Ground rotation")
    #Angular velocity of Earth, 2 pi every 23hrs 56 min (sidereal day)
    omega = (2*math.pi)/((23*60*60)+(56*60)+4) #rad/sec

    x, y, z = Station_Coords

    Ground_Station_Array = []
    Ground_Station_Array.append(Station_Coords)

    #For every time step in the array
    for i, time in enumerate(time_end):
        angle_rad = omega*time_end[i] #Calculate the angle of rotation

        rotation_matrix = np.array([
            [np.cos(angle_rad), -np.sin(angle_rad), 0],
            [np.sin(angle_rad),  np.cos(angle_rad), 0],
            [0,                 0,                 1]
        ])

        rotated = rotation_matrix @ np.array([x,y,z])
        Ground_Station_Array.append(rotated.tolist())
        #print(rotated)
    
    #print(Ground_Station_Array)
    first = Station_Coords
    Ground_Station_Array = [first[:] for _ in Ground_Station_Array] 
    return Ground_Station_Array
    



# plot orbit
def plot_xyz(state_vectors, r, Animation, state_vectors2, GUI, GUIOneSat, GroundSat):
    global fig, ax, orbit, satellite, Earth
    global X, Y, Z
    global u2, v
    global Station_Sat_Line
    

    global X2, Y2, Z2
    global orbit2
    
    if GUIOneSat==True:
        X, Y, Z = state_vectors[0], state_vectors[1], state_vectors[2]
    else:
        X, Y, Z = state_vectors[0], state_vectors[1], state_vectors[2]
        X2, Y2, Z2 = state_vectors2[0], state_vectors2[1], state_vectors2[2]
    
    fig = plt.figure()
    #ax = plt3.Axes3D(fig)
    ax = plt.axes(projection='3d')
    #ax.set_facecolor('black')

    # set labels
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')

    # set correct aspect ratio
    ax.set_box_aspect([1,1,1])
    set_axes_equal_3d(ax)
    # set view
    #azm=29; ele=-36; roll = 180 #39 150,36
    azm = 150; ele = 36
    ax.view_init(elev=ele, azim=azm)
    # set limit
    size = 1.0
    #limit = max(max(X[0]), max(Y[0]), max(Z[0]))
    #limit = max(max(X), max(Y), max(Z))
    limit = max(max(X, key=abs), max(Y, key=abs), max(Z, key=abs), key=abs)

    

    ax.set_xlim(size*limit, -size*limit)
    ax.set_ylim(size*limit, -size*limit)
    ax.set_zlim(size*limit, -size*limit)

    # earth
    #ax.scatter(0, 0, 0, marker='o', c='deepskyblue', s=3*r, alpha=0.5)

    u2, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:50j]
    x = r*np.cos(u2)*np.sin(v)
    y = r*np.sin(u2)*np.sin(v)
    z = r*np.cos(v)
    Earth = ax.plot_surface(x, y, z, alpha = 0.3)

    # Station
    #ax.scatter(Station_Coords[0], Station_Coords[1], Station_Coords[2], c='r', s=100)

    #ax.scatter(-6354.282219769004, 551.0914877540154, -1.3122033207603634, c='r', s=100)


    #If statement for either plotting animation or static picture
    if Animation == "yes":
        print("Animate")
        #Adding the trailing line
        if GUIOneSat==True:
            orbit = ax.plot(X[0], Y[0], Z[0], linewidth=0.9, c='k', linestyle='dashed')[0]
            Station_Sat_Line, = ax.plot([Station_Coords[0], Station_Coords[0]], [Station_Coords[1], Station_Coords[1]], [Station_Coords[2], Station_Coords[2]], c='r')
            ani = animation.FuncAnimation(fig, update2,init_func=init ,frames=len(X), interval=20, blit=False, repeat=False)
        else:
            orbit = ax.plot(X[0], Y[0], Z[0], linewidth=0.9, c='k', linestyle='dashed')[0]
            orbit2 = ax.plot(X2[0], Y2[0], Z2[0], linewidth=0.9, c='k', linestyle='dashed')[0]
            Station_Sat_Line, = ax.plot([Station_Coords[0], Station_Coords[0]], [Station_Coords[1], Station_Coords[1]], [Station_Coords[2], Station_Coords[2]], c='r')
            if GroundSat == True:
                ani = animation.FuncAnimation(fig, update,init_func=init ,frames=len(X), interval=20, blit=False, repeat=False)
            else:
                ani = animation.FuncAnimation(fig, update3,init_func=init ,frames=len(X), interval=20, blit=False, repeat=False)
            
        #Adding the line between the station and sat
        #Station_Sat_Line, = ax.plot([Station_Coords[0], Station_Coords[0]], [Station_Coords[1], Station_Coords[1]], [Station_Coords[2], Station_Coords[2]], c='r')
        #ani = animation.FuncAnimation(fig, update,init_func=init ,frames=len(X), interval=20, blit=False, repeat=False)
        if GUI == True:
            return fig, ani
        elif GUI== False:
            plt.show()
            print("saving animation!")
            #ani.save(filename="Full3d.gif", writer='pillow')
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
    global station
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
    #Moving station
    station = ax.plot([], [], [], marker = 'o', ms=10, c='r')[0]

    return scat, scat2, station, 

def update(frame):
    if frame == len(X) - 1:
        scat.set_data([], [])
        scat.set_3d_properties([])
        scat2.set_data([], [])
        scat2.set_3d_properties([])
        station.set_data([], [])
        station.set_3d_properties([])
        # orbit.set_data([], [])
        # orbit.set_3d_properties([])
        # orbit2.set_data([], [])
        # orbit2.set_3d_properties([])
        # Station_Sat_Line.set_data([], [])
        # Station_Sat_Line.set_3d_properties([])
        #Earth.set_alpha(0)
        return scat, orbit, Station_Sat_Line, scat2, orbit2, Earth, station
    
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

    #Updating the ground station position
    NewStationX, NewStationY, NewStationZ = Ground_Station_Array[frame][0], Ground_Station_Array[frame][1], Ground_Station_Array[frame][2]
    station.set_data([NewStationX], [NewStationY])
    station.set_3d_properties([NewStationZ])

    '''
    Check if the condition for intersat link is True
    '''
    Station_Sat_Line.set_data([Ground_Station_Array[frame][0], XCoord], [Ground_Station_Array[frame][1], YCoord])
    Station_Sat_Line.set_3d_properties([Ground_Station_Array[frame][2],ZCoord])
    return (scat, orbit, Station_Sat_Line, orbit2, scat2, Earth, station)

#This one is for only one satellite
def update2(frame):
    if frame == len(X) - 1:
        # scat.set_data([], [])
        # scat.set_3d_properties([])
        # #scat2.set_data([], [])
        # #scat2.set_3d_properties([])
        # orbit.set_data([], [])
        # orbit.set_3d_properties([])
        # station.set_data([], [])
        # station.set_3d_properties([])
        # orbit2.set_data([], [])
        # orbit2.set_3d_properties([])
        # Station_Sat_Line.set_data([], [])
        # Station_Sat_Line.set_3d_properties([])
        #Earth.set_alpha(0)
        return scat, orbit, Station_Sat_Line, Earth, station
    
    #print("update")
    XCoord, YCoord, ZCoord = X[frame], Y[frame], Z[frame]
    #scat, = ax.plot(XCoord, YCoord, ZCoord, marker = 'o', ms=10, c='g')
    scat.set_data([XCoord], [YCoord])
    scat.set_3d_properties([ZCoord])

    # #updating position of second orbit
    # XCoord2, YCoord2, ZCoord2 = X2[frame], Y2[frame], Z2[frame]
    # #scat, = ax.plot(XCoord, YCoord, ZCoord, marker = 'o', ms=10, c='g')
    # scat2.set_data([XCoord2], [YCoord2])
    # scat2.set_3d_properties([ZCoord2])

    #Updating the trailing line 
    #orbit, = ax.plot(X[:frame], Y[:frame], Z[:frame], linewidth=0.9, linestyle='dashed', c='k')
    orbit.set_data(X[:frame], Y[:frame])
    orbit.set_3d_properties(Z[:frame])

    # #Updating the trailing line for second sat
    # #orbit, = ax.plot(X[:frame], Y[:frame], Z[:frame], linewidth=0.9, linestyle='dashed', c='k')
    # orbit2.set_data(X2[:frame], Y2[:frame])
    # orbit2.set_3d_properties(Z2[:frame])

    #Updating the station_sat line
    #Station_Sat_Line, = ax.plot([Station_Coords[0], XCoord], [Station_Coords[1], YCoord], [Station_Coords[2],ZCoord], c='r')
    Station_Sat_Line.set_data([Ground_Station_Array[frame][0], XCoord], [Ground_Station_Array[frame][1], YCoord])
    Station_Sat_Line.set_3d_properties([Ground_Station_Array[frame][2],ZCoord])

    NewStationX, NewStationY, NewStationZ = Ground_Station_Array[frame][0], Ground_Station_Array[frame][1], Ground_Station_Array[frame][2]
    station.set_data([NewStationX], [NewStationY])
    station.set_3d_properties([NewStationZ])
    return (scat, orbit, Station_Sat_Line, Earth, station)

def update3(frame):
    if frame == len(X) - 1:
        scat.set_data([], [])
        scat.set_3d_properties([])
        scat2.set_data([], [])
        scat2.set_3d_properties([])
        station.set_data([], [])
        station.set_3d_properties([])
        # orbit.set_data([], [])
        # orbit.set_3d_properties([])
        # orbit2.set_data([], [])
        # orbit2.set_3d_properties([])
        # Station_Sat_Line.set_data([], [])
        # Station_Sat_Line.set_3d_properties([])
        #Earth.set_alpha(0)
        return scat, orbit, Station_Sat_Line, scat2, orbit2, Earth, station
    
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

    NewStationX, NewStationY, NewStationZ = Ground_Station_Array[frame][0], Ground_Station_Array[frame][1], Ground_Station_Array[frame][2]
    station.set_data([NewStationX], [NewStationY])
    station.set_3d_properties([NewStationZ])

    '''
    Check if the condition for intersat link is True
    '''
    Station_Sat_Line.set_data([XCoord2, XCoord], [YCoord2, YCoord])
    Station_Sat_Line.set_3d_properties([ZCoord2,ZCoord])
    return (scat, orbit, Station_Sat_Line, orbit2, scat2, Earth, station)


def CalcElevationAngle(StateVectors,laser_power, wavelength, tx_aperture, rx_aperture, tx_transmission, rx_transmission, rx_obscuration, rx_smf_il, smf_coupling):
    #link budget calc
    link_budget_list = []
    #SatXYZ is lists of each coord for each time

    SatX, SatY, SatZ = StateVectors[0], StateVectors[1], StateVectors[2]
    #print(SatX[0])
    
    ElevationAngle = []
    Distance_array = []

    #calculating the vector between sat and ground station
    for j, value  in enumerate(SatX):
        Vector_Sat_Stn = np.array([SatX[j]-Ground_Station_Array[j][0], SatY[j]-Ground_Station_Array[j][1], SatZ[j]-Ground_Station_Array[j][2]])
        Distance_Sat_Stn = np.linalg.norm(Vector_Sat_Stn)
        Distance_array.append(Distance_Sat_Stn)
        StationCoords = np.array([Ground_Station_Array[j][0], Ground_Station_Array[j][1], Ground_Station_Array[j][2]])
        Dot_Prod = np.dot(Vector_Sat_Stn, StationCoords)
        #print(np.linalg.norm(Vector_Sat_Stn))
        #print(np.linalg.norm(StationCoords))
        #print(np.linalg.norm(Vector_Sat_Stn)*np.linalg.norm(StationCoords))
        try:
            Angle_between = math.acos(Dot_Prod/(np.linalg.norm(Vector_Sat_Stn)*np.linalg.norm(StationCoords)))

            #angle from 90 degree elevation angle
            Real_angle = math.pi/2-abs(Angle_between) #Elevation angle in radians
            Real_angle = math.degrees(Real_angle) #convert to degrees
        except:
            '''
            This is for the geostationary orbit case, when the coordinates of the satellite and station are so close that there is a floating point error or something
            Will basically never happen IRL but it may occur where you have defined a ground station orbit so that it is DIRECTLY under the satellite
            '''
            Angle_between = 0 
            Real_angle = 90 #Satellite is directly overhead so elevation angle = 90 degrees
        
        
        
        
        if Real_angle >= 5:
            #print('e')
            Magnitude_vector_sat_earth = np.linalg.norm([SatX[j], SatY[j], SatZ[j]])
            sat_altitude = Magnitude_vector_sat_earth-r
            #print(sat_altitude)
            print("TX: ", tx_aperture)
            laser_dB = link_budget(laser_power, wavelength, tx_aperture, rx_aperture, sat_altitude*1000, Real_angle, tx_transmission, rx_transmission, rx_obscuration, rx_smf_il, smf_coupling)
            #laser_dB = link_budget3(laser_power, wavelength, tx_aperture, rx_aperture, sat_altitude*1000, Real_angle, tx_transmission, rx_transmission, rx_obscuration, rx_smf_il, smf_coupling)
        else:
            laser_dB = None
        link_budget_list.append(laser_dB)
        ElevationAngle.append(Real_angle)

    #print(ElevationAngle)

    
    return ElevationAngle, link_budget_list, Distance_array

#def CalcLinkBudgetInterSat():
    #laser_dB = link_budget(laser_power, wavelength, tx_aperture, rx_aperture, sat_altitude*1000, Real_angle, tx_transmission, rx_transmission, rx_obscuration, rx_smf_il)

def plot_link_budget(link_list, time_arr, inGUI):
    #time_arr = pd.to_datetime(time_arr)
    print(type(time_arr[0]))
    if inGUI == True:
        plt.clf()
        fig2 = plt.figure()
        
        ax2 = fig2.add_subplot(111)
        ax2.plot(time_arr, link_list)
        ax2.set_xlabel("Time")
        ax2.set_ylabel("Link margin (dB)")
        #plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        fig2.show()
        #fig2.savefig("TestFIGURE.png")
        return fig2
    else:
        #print(link_list)
        plt.plot(time_arr, link_list)
        plt.xlabel("Time")
        plt.ylabel("Link margin (dB)")
        #plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        plt.show()    

def IsIntersecion(p1, p2):
    # Line direction
    d = p2 - p1
    
    # Quadratic coefficients: ||p1 + t*d||^2 = R^2
    a = np.dot(d, d)
    b = 2 * np.dot(p1, d)
    c = np.dot(p1, p1) - r**2
    
    # Discriminant
    discriminant = b**2 - 4*a*c
    
    if discriminant < 0:
        return False  # no intersection
    
    # Solve quadratic
    sqrt_disc = np.sqrt(discriminant)
    t1 = (-b - sqrt_disc) / (2*a)
    t2 = (-b + sqrt_disc) / (2*a)
    
    # Check if either intersection point is on the segment
    if (0 <= t1 <= 1) or (0 <= t2 <= 1):
        return True
    else:
        return False

def LinkBudgetInterSatellite(state_vector1, state_vector2, laser_power, wavelength, tx_aperture, rx_aperture, tx_transmission, rx_transmission, rx_obscuration, rx_smf_il, smf_coupling):
    print('intersatellite link budget')
    Sat1X, Sat1Y, Sat1Z = state_vector1[0], state_vector1[1], state_vector1[2]
    Sat2X, Sat2Y, Sat2Z = state_vector2[0], state_vector2[1], state_vector2[2]

    link_budget_list = []
    distance_list = []

    for i, value in enumerate(Sat1X):


        #Calculating the distance between the two satellites
        Vector_sat_sat = np.array([Sat1X[i]-Sat2X[i], Sat1Y[i]-Sat2Y[i], Sat1Z[i]-Sat2Z[i]])
        distance_between = np.linalg.norm(Vector_sat_sat)

        #Calculate whether the line between them intersects with the Earth
        p1 = np.array([Sat1X[i], Sat1Y[i], Sat1Z[i]])
        p2 = np.array([Sat2X[i], Sat2Y[i], Sat2Z[i]])

        #If is true, it does intersect the Earth. If False, does not.
        DoesCollide = IsIntersecion(p1, p2)

        if DoesCollide == False:
            #print("Does Not Intersect")
            laser_dB = link_budget2(laser_power, wavelength, tx_aperture, rx_aperture, distance_between*1000, tx_transmission, rx_transmission, rx_obscuration, rx_smf_il, smf_coupling)
            link_budget_list.append(laser_dB)
        else:
            laser_dB = None
            link_budget_list.append(laser_dB)

        distance_list.append(distance_between)
    return link_budget_list, distance_list

def main():
    satellite, epoch = get_satellite(sat1_l1, sat1_l2)
    time_arr, time_end, t1         = get_propagation_times(epoch)
    Ground_Station_Array = GroundRotation(time_end)

    satellite2, epoch2 = get_satellite(sat2_l1, sat2_l2)
    time_arr2 , time_end, t2        = get_propagation_times(epoch2)

    #if satellite 1 has larger period than 2
    if t1 > t2:
        time_arr2, time_end = get_propagation_times2(epoch2, t1)
    elif t2 > t1:
        time_arr, time_end = get_propagation_times2(epoch, t2)

    state_vectors    = get_state_vectors(satellite, time_arr)
    print(len(state_vectors[0]))
    print(state_vectors[0][int(len(state_vectors[0])/2)], state_vectors[1][int(len(state_vectors[0])/2)], state_vectors[2][int(len(state_vectors[0])/2)])
    #print(state_vectors[0][0], state_vectors[1][0], state_vectors[2][0])
    Elevation_Angle, link_bud, distance_array = CalcElevationAngle(state_vectors,laser_power, wavelength, tx_aperture, rx_aperture, tx_transmission, rx_transmission, rx_obscuration, rx_smf_il, smf_coupling)

    
    state_vectors2    = get_state_vectors(satellite2, time_arr2)
    Elevation_Angle2, link_bud2, distance_array2 = CalcElevationAngle(state_vectors2,laser_power, wavelength, tx_aperture, rx_aperture, tx_transmission, rx_transmission, rx_obscuration, rx_smf_il, smf_coupling)
    
    #print(Elevation_Angle)

    

    plot_xyz(state_vectors, r, "yes", state_vectors2, False, True, True)
    #print(type(time_arr))
    #print(type(time_arr[0]))
    fig2 = plot_link_budget(link_bud, time_arr, False)

    #print(link_bud)
    return

# main
if __name__ == '__main__':
    main()