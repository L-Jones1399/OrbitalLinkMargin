import tkinter
from ThreeDEarth import *
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from PIL import Image

# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2Tk)

#make the process DPI aware so to avoid rescaling issues when draging window to second monitor/without the main monitor
import ctypes
try:
    ctypes.windll.shcore.SetProcessDpiAwareness(1)  # system DPI aware
    # or use 2 for per-monitor DPI aware
except Exception:
    pass

root = tkinter.Tk()
#root.geometry("1000x600") #set the size of the GUI
root.wm_title("Link Margin Simulator") #Title


satellite, epoch = get_satellite(sat1_l1, sat1_l2)
time_arr, time_end, t1         = get_propagation_times(epoch)

satellite2, epoch2 = get_satellite(sat2_l1, sat2_l2)
time_arr2, time_end, t2         = get_propagation_times(epoch2)

#if satellite 1 has larger period than 2
if t1 > t2:
    time_arr2, time_end = get_propagation_times2(epoch2, t1)
elif t2 > t1:
    time_arr, time_end = get_propagation_times2(epoch, t2)

Ground_Station_Array = GroundRotation(time_end)
state_vectors    = get_state_vectors(satellite, time_arr)
print(len(state_vectors[0]))
#print(state_vectors[0][33], state_vectors[1][33], state_vectors[2][33])
Elevation_Angle, link_bud, distance_array = CalcElevationAngle(state_vectors,laser_power, wavelength, tx_aperture, rx_aperture, tx_transmission, rx_transmission, rx_obscuration, rx_smf_il, smf_coupling)



state_vectors2    = get_state_vectors(satellite2, time_arr2)
r = satellite.radiusearthkm

link_bud, distance_array = LinkBudgetInterSatellite(state_vectors, state_vectors2,laser_power, wavelength, tx_aperture, rx_aperture, tx_transmission, rx_transmission, rx_obscuration, rx_smf_il, smf_coupling)
fig, ani = plot_xyz(state_vectors, r, "yes", state_vectors2, True, False, False)
fig2 = plot_link_budget(link_bud, time_arr, True)
#fig2.savefig("TESTING2.png")

# fig = Figure(figsize=(5, 4), dpi=100)
# t = np.arange(0, 3, .01)
# ax = fig.add_subplot()
# line, = ax.plot(t, 2 * np.sin(2 * np.pi * t))
# ax.set_xlabel("time [s]")
# ax.set_ylabel("f(t)")

canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.
canvas.draw()
canvas.get_tk_widget().grid(row=0, column=0)
#canvas.get_tk_widget().pack(side=tkinter.TOP)


# pack_toolbar=False will make it easier to use a layout manager later on.
toolbar = NavigationToolbar2Tk(canvas, root, pack_toolbar=False)
toolbar.update()

canvas.mpl_connect(
    "key_press_event", lambda event: print(f"you pressed {event.key}"))
canvas.mpl_connect("key_press_event", key_press_handler)

#Function closes the Tkinter window and the matplotlib figure
def endProgram():
    root.destroy()
    plt.close("all")

def saveGif():
    print("Saving data to text file!")
    print("Filename is: ", filename_var.get())
    fname = r"C:\Users\zeus1\Documents\RMIT\Honours\Project\FSOC Code\DataFiles\\" + filename_var.get() +".txt"
    #DataFile = open(r"C:\Users\zeus1\Documents\RMIT\Honours\Project\FSOC Code\DataFiles\testfile.txt", "w")
    DataFile = open(fname, "w")
    DataFile.write("-----Start-----\n")
    text = TLE_1.get("1.0", "end-1c")
    lines = text.split('\n')
    DataFile.write(lines[0]+ '\n')
    DataFile.write(lines[1] +'\n')
    text = TLE_2.get("1.0", "end-1c")
    lines = text.split('\n')
    DataFile.write(lines[0]+ '\n')
    DataFile.write(lines[1] +'\n')
    #Determine the type of link
    if NoofSats.get() ==1:
        DataFile.write("One Satellite, Ground Link \n")
    elif NoofSats.get()==2 and TypeofLink.get()==1:
        DataFile.write("Two Satellites, Ground Link \n")
    elif NoofSats.get()==2 and TypeofLink.get()==2:
        DataFile.write("Two Satellites, Inter-Satellite Link \n")
    
    DataFile.write("Station XYZ: "+ ",".join(str(x) for x in Station_Coords))
    DataFile.write("\nLaser Power: "+ laserpower_var.get())
    DataFile.write("\nWavelength: " + str(wavelength))
    DataFile.write('\ntx_aperture: ' + txaperture_var.get())
    DataFile.write('\ntx_transmission: ' + txtransmission_var.get())
    DataFile.write('\nrx_aperture: ' + rxaperture_var.get())
    DataFile.write('\nrx_transmission: ' + rxtransmission_var.get())
    DataFile.write('\nrx_obscuration: ' + rxobscuration_var.get())
    DataFile.write('\nrx_smfil: ' + rxsmfil_var.get())
    DataFile.write('\nrx_switch: ' + rxswitchil_var.get())
    DataFile.write('\nsmf_coupling: ' + smfcoupling_var.get())
    DataFile.write('\n'+"-----End-----\n")

    for row in zip(time_arr, Elevation_Angle, distance_array, link_bud):
        DataFile.write("\t".join(map(str, row))+'\n')
    
    DataFile.close()
    print("Saving Link Budget Fig!")
    fname_graph = r"C:\Users\zeus1\Documents\RMIT\Honours\Project\FSOC Code\Graphs\\" + filename_var.get() +".png"
    fig2.savefig(fname_graph)
    if is_on == True:
        print("saving animation!")
        fname_gif = r"C:\Users\zeus1\Documents\RMIT\Honours\Project\FSOC Code\GIFs\\" + filename_var.get() +".gif"
        ani.save(filename=fname_gif, writer='pillow')
    print("Done!")




    
'''
Create some radiobuttons
'''
NoofSats = tkinter.IntVar(value=2)
rb1 = tkinter.Radiobutton(root, text="One satellite", variable=NoofSats, value=1)
rb2 = tkinter.Radiobutton(root, text="Two satellites", variable=NoofSats, value=2)
NoofSats.set(2)

TypeofLink = tkinter.IntVar(value=1)
rb3 = tkinter.Radiobutton(root, text="Ground-Sat", variable=TypeofLink, value=1)
rb4 = tkinter.Radiobutton(root, text="Intersat", variable=TypeofLink, value=2)
TypeofLink.set(2)


helv16Bold = tkinter.font.Font(family='Helvetica', size=16, weight=tkinter.font.BOLD)
helv16 = tkinter.font.Font(family='Helvetica', size=16)
button_quit = tkinter.Button(master=root, text="Quit",height=2, width=6, command=endProgram, font=helv16Bold)
button_save = tkinter.Button(master=root, text="Save", command=saveGif, font=helv16)

#Generate some labels to display information
#global sat_1_apogee_Var, sat_1_perigee_Var,sat_2_apogee_Var,sat_2_perigee_Var
sat_1_apogee_Var = tkinter.StringVar()
sat_1_apogee_Var.set("Satellite Apogee (km): " + str(round(satellite.alta*r,4)))
sat_1_perigee_Var = tkinter.StringVar()
sat_1_perigee_Var.set("Satellite Perigee (km): " + str(round(satellite.altp*r,4)))
Sat_1_apogee =  tkinter.Label(master=root, textvariable=sat_1_apogee_Var)
Sat_1_perigee =  tkinter.Label(master=root, textvariable=sat_1_perigee_Var)

sat_1_period_Var = tkinter.StringVar()
sat_1_period_Var.set("Orbital Period (min): "+str(round((2 * np.pi * ((satellite.a*r)**3/satellite.mu)**0.5)/60,4)))
sat_1_period = tkinter.Label(master=root, textvariable=sat_1_period_Var)

sat_2_apogee_Var = tkinter.StringVar()
sat_2_apogee_Var.set("Satellite Apogee (km): " + str(round(satellite2.alta*r,4)))
sat_2_perigee_Var = tkinter.StringVar()
sat_2_perigee_Var.set("Satellite Perigee (km): " + str(round(satellite2.altp*r,4)))

sat_2_period_Var = tkinter.StringVar()
sat_2_period_Var.set("Orbital Period (min): "+str(round((2 * np.pi * ((satellite2.a*r)**3/satellite.mu)**0.5)/60,4)))
sat_2_period = tkinter.Label(master=root, textvariable=sat_2_period_Var)

Sat_2_apogee =  tkinter.Label(master=root, textvariable=sat_2_apogee_Var)
Sat_2_perigee =  tkinter.Label(master=root, textvariable=sat_2_perigee_Var)

'''
Create some entry boxes for the satellite/ground infrastructure
'''
laserpower_var = tkinter.StringVar()
laserpower_lab = tkinter.Label(root, text="Laser Power (W)")
txaperture_var = tkinter.StringVar()
txaperture_lab = tkinter.Label(root, text="Transmitter Aperture (m)")
rxaperture_var = tkinter.StringVar()
rxaperture_lab = tkinter.Label(root, text="Receiver Aperture (m)")
txtransmission_var = tkinter.StringVar()
txtransmission_lab = tkinter.Label(root, text = "Transmitter Transmission (0<=x<=1)")
rxtransmission_var = tkinter.StringVar()
rxtransmission_lab = tkinter.Label(root, text = "Receiver Transmission (0<=x<=1)")
rxobscuration_var = tkinter.StringVar()
rxobscuration_lab = tkinter.Label(root, text="Receiver Obscuration (0<=x<=1)")
rxsmfil_var = tkinter.StringVar()
rxsmfil_lab = tkinter.Label(root, text = "Receiver Fibre Insertion Loss (dB)")
rxswitchil_var = tkinter.StringVar()
rxswitchil_lab = tkinter.Label(root, text="Receiver Switch Insertion Loss (dB)")
smfcoupling_var = tkinter.StringVar()
smfcoupling_lab = tkinter.Label(root, text="Single Mode Fibre Coupling Loss (dB)")

laserpower_entry = tkinter.Entry(root,textvariable=laserpower_var)
txaperture_entry = tkinter.Entry(root,textvariable=txaperture_var)
rxaperture_entry = tkinter.Entry(root,textvariable=rxaperture_var)
txtransmission_entry = tkinter.Entry(root,textvariable=txtransmission_var)
rxtransmission_entry = tkinter.Entry(root,textvariable=rxtransmission_var)
rxobscuration_entry = tkinter.Entry(root, textvariable=rxobscuration_var)
rxsmfil_entry =tkinter.Entry(root,textvariable=rxsmfil_var)
rxswitchil_entry = tkinter.Entry(root,textvariable=rxswitchil_var)
smfcoupling_entry = tkinter.Entry(root, textvariable=smfcoupling_var)

laserpower_entry.insert(0, laser_power)
txaperture_entry.insert(0, tx_aperture)
rxaperture_entry.insert(0, rx_aperture)
txtransmission_entry.insert(0, tx_transmission)
rxtransmission_entry.insert(0, rx_transmission)
rxobscuration_entry.insert(0, rx_obscuration)
rxsmfil_entry.insert(0, rx_smf_il)
rxswitchil_entry.insert(0, rx_switch_il)
smfcoupling_entry.insert(0, smf_coupling)

stationx_var = tkinter.StringVar()
stationx_entry = tkinter.Entry(root, textvariable=stationx_var)
stationy_var = tkinter.StringVar()
stationy_entry = tkinter.Entry(root, textvariable=stationy_var)
stationz_var = tkinter.StringVar()
stationz_entry = tkinter.Entry(root, textvariable=stationz_var)
station_lab = tkinter.Label(root, text="Station Coordinates (XYZ (km))")

stationx_entry.insert(0, Station_Coords[0])
stationy_entry.insert(0, Station_Coords[1])
stationz_entry.insert(0, Station_Coords[2])

filename_var = tkinter.StringVar()
filename_lab = tkinter.Label(root, text="Filename: ")
filename_entry = tkinter.Entry(root, textvariable=filename_var)

is_on = True
# Create Label
Switch_label = tkinter.Label(root, 
    text = "GIF will be saved", 
    fg = "green", 
    font = ("Helvetica", 11))
def switch():
    global is_on
    # Determine is on or off
    if is_on:
        savegif_slider.config(image = off)
        Switch_label.config(text = "GIF will not be saved", 
                        fg = "red")
        is_on = False
    else:
      
        savegif_slider.config(image = on)
        Switch_label.config(text = "GIF will be saved", fg = "green")
        is_on = True
# Define Our Images
#on_original = Image.open("on.png")
#on_resize = on_original.resize((on_original.width//5, on_original.height//5))
#on = tkinter.PhotoImage(on_resize)
on_original = tkinter.PhotoImage(file="on.png")
#on = on_original
on = on_original.subsample(7,7)

#off_original = Image.open("on.png")
#off_resize = on_original.resize((off_original.width//5, off_original.height//5))
#off = tkinter.PhotoImage(off_resize)
#off = tkinter.PhotoImage(file = "off.png")
off_original = tkinter.PhotoImage(file="off.png")
#off = off_original
off = off_original.subsample(7,7)

# Create A Button
savegif_slider = tkinter.Button(root, image = on, bd = 0,
                   command = switch)
#savegif_slider = tkinter.Scale(root, from_=0, to=1, orient="horizontal")
#savegif_slider.set(1)

TLE_1 = tkinter.Text(root, height=2, width = 70, wrap='none')
TLE_1.insert("1.1", sat1_l1+'\n')
TLE_1.insert("2.0", sat1_l2)
TLE_2 = tkinter.Text(root, height=2, width = 70, wrap='none')
TLE_2.insert("1.1", sat2_l1+'\n')
TLE_2.insert("2.0", sat2_l2)
    
MAX_LINES = 2
MAX_CHARS = 69

def validate_input(event):
    widget = event.widget
    # Ignore special control keys like arrows, delete, etc.
    if event.keysym in ("BackSpace", "Left", "Right", "Up", "Down", "Delete", "Home", "End"):
        return

    # Get current text
    text = widget.get("1.0", "end-1c")
    lines = text.split("\n")

    # Get cursor line & column
    index = widget.index("insert")
    line_no, col = map(int, index.split("."))

    # --- Block Enter if already at max lines ---
    if event.keysym == "Return" and len(lines) >= MAX_LINES:
        return "break"

    if event.char and event.char not in ("\n", "\r"):
        current_line = lines[line_no - 1] if line_no - 1 < len(lines) else ""
        if len(current_line) >= MAX_CHARS:
            return "break"

    # --- Prevent paste from exceeding limits ---
    if event.state & 0x4 and event.keysym.lower() == "v":  # Ctrl+V
        widget.after_idle(lambda w=widget: trim_text(w))

def trim_text(widget):
    """Trim excess lines and characters (for pasted input)."""
    text = widget.get("1.0", "end-1c")
    lines = text.split("\n")[:MAX_LINES]
    lines = [line[:MAX_CHARS] for line in lines]
    widget.delete("1.0", "end")
    widget.insert("1.0", "\n".join(lines))


#Enforce a 2 line and 69 character limit
TLE_1.bind("<Key>", validate_input)
TLE_2.bind("<Key>", validate_input)




'''
Update the graph when the button is pushed
'''
def updateGraph():
    global canvas
    global ani
    global Elevation_Angle
    global distance_array
    global link_bud
    global fig2
    print("updating graph")

    laser_power = float(laserpower_var.get())
    tx_aperture = float(txaperture_var.get())
    rx_aperture = float(rxaperture_var.get())
    tx_transmission = float(txtransmission_var.get())
    rx_transmission = float(rxtransmission_var.get())
    rx_obscuration = float(rxobscuration_var.get())
    rx_smf_il = float(rxsmfil_var.get())
    rx_switch_il = float(rxswitchil_var.get())
    smf_coupling = float(smfcoupling_var.get())

    Station_Coords[0] = float(stationx_var.get())
    Station_Coords[1] = float(stationy_var.get())
    Station_Coords[2] = float(stationz_var.get())

    #If one satellite is picked
    if NoofSats.get() == 1:
        print("One satellite")
        text = TLE_1.get("1.0", "end-1c")
        lines = text.split('\n')
        satellite, epoch = get_satellite(lines[0], lines[1])
        time_arr, time_end, t1         = get_propagation_times(epoch)
        Ground_Station_Array = GroundRotation(time_end)
        state_vectors    = get_state_vectors(satellite, time_arr)
        Elevation_Angle, link_bud,distance_array = CalcElevationAngle(state_vectors,laser_power, wavelength, tx_aperture, rx_aperture, tx_transmission, rx_transmission, rx_obscuration, rx_smf_il, smf_coupling)
        state_vectors2=None

        sat_1_apogee_Var.set("Satellite Apogee (km): " + str(round(satellite.alta*r,4)))
        sat_1_perigee_Var.set("Satellite Perigee (km): " + str(round(satellite.altp*r,4)))
        sat_1_period_Var.set("Orbital Period (min): "+str(round((2 * np.pi * ((satellite.a*r)**3/satellite.mu)**0.5)/60,4)))

        fig2 = plot_link_budget(link_bud, time_arr, True)

        fig, ani = plot_xyz(state_vectors, r, "yes", state_vectors2, True, True, False)
        
    elif NoofSats.get() == 2:
        
        
        text = TLE_1.get("1.0", "end-1c")
        lines = text.split('\n')

        text2 = TLE_2.get("1.0", "end-1c")
        lines2 = text2.split('\n')
        satellite, epoch = get_satellite(lines[0], lines[1])
        time_arr, time_end, t1         = get_propagation_times(epoch)

        satellite2, epoch2 = get_satellite(lines2[0], lines2[1])
        time_arr2, time_end, t2         = get_propagation_times(epoch2)

        #if satellite 1 has larger period than 2
        if t1 > t2:
            time_arr2, time_end = get_propagation_times2(epoch2, t1)
        elif t2 > t1:
            time_arr, time_end = get_propagation_times2(epoch, t2)

        Ground_Station_Array = GroundRotation(time_end)
        state_vectors    = get_state_vectors(satellite, time_arr)
        Elevation_Angle, link_bud,distance_array = CalcElevationAngle(state_vectors,laser_power, wavelength, tx_aperture, rx_aperture, tx_transmission, rx_transmission, rx_obscuration, rx_smf_il, smf_coupling)

        
        
        state_vectors2    = get_state_vectors(satellite2, time_arr2)

        sat_1_apogee_Var.set("Satellite Apogee (km): " + str(round(satellite.alta*r,4)))
        sat_1_perigee_Var.set("Satellite Perigee (km): " + str(round(satellite.altp*r,4)))
        sat_2_apogee_Var.set("Satellite Apogee (km): " + str(round(satellite2.alta*r,4)))
        sat_2_perigee_Var.set("Satellite Perigee (km): " + str(round(satellite2.altp*r,4)))
        sat_1_period_Var.set("Orbital Period (min): "+str(round((2 * np.pi * ((satellite.a*r)**3/satellite.mu)**0.5)/60,4)))
        sat_2_period_Var.set("Orbital Period (min): "+str(round((2 * np.pi * ((satellite2.a*r)**3/satellite.mu)**0.5)/60,4)))

        if TypeofLink.get() == 1:
            print("Two satellites, Ground-Sat link")
            fig2 = plot_link_budget(link_bud, time_arr, True)
            fig, ani = plot_xyz(state_vectors, r, "yes", state_vectors2, True, False, True)
        elif TypeofLink.get() == 2:
            print("Two satellites, Inter-sat link")
            Elevation_Angle = [0]*len(Elevation_Angle)
            link_bud, distance_array = LinkBudgetInterSatellite(state_vectors, state_vectors2,laser_power, wavelength, tx_aperture, rx_aperture, tx_transmission, rx_transmission, rx_obscuration, rx_smf_il, smf_coupling)
            fig2 = plot_link_budget(link_bud, time_arr, True)
            fig, ani = plot_xyz(state_vectors, r, "yes", state_vectors2, True, False, False)

        #If the link is intersatellite
        
            
    
    canvas.get_tk_widget().destroy()
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.get_tk_widget().grid(row=2,column=5, columnspan=2, rowspan=12, padx=10)
    canvas.draw()

button_submit = tkinter.Button(master=root, text="Submit", font=helv16, command=updateGraph)

# Packing order is important. Widgets are processed sequentially and if there
# is no space left, because the window is too small, they are not displayed.
# The canvas is rather flexible in its size, so we pack it last which makes
# sure the UI controls are displayed as long as possible.

button_quit.grid(row=1,column=6, sticky='e',padx=10,pady=10)
button_save.grid(row=14,column=6, sticky='e',padx=10, pady=10)
#slider_update.pack(side=tkinter.BOTTOM)
toolbar.grid(row=14,column=5, columnspan=2,sticky='n')

canvas.get_tk_widget().grid(row=2,column=5, columnspan=2, rowspan=12, padx=10)

rb1.grid(row=14, column=1)
rb2.grid(row=14,column=2)

rb3.grid(row=14, column=3, sticky='w')
rb4.grid(row=14, column=4, sticky='w')

button_submit.grid(row=12, column=3, sticky="s")

savegif_slider.grid(row=12, column=4)
Switch_label.grid(row=11, column=4)

TLE_1.grid(row=1, column=1, columnspan=4)
TLE_2.grid(row=5, column=1, columnspan=4)

Sat_1_apogee.grid(row=2, column=1, sticky="n")
Sat_1_perigee.grid(row=2, column=2, sticky="n")
sat_1_period.grid(row=2, column=3, sticky='n')

Sat_2_apogee.grid(row=6, column=1, sticky="n")
Sat_2_perigee.grid(row=6, column=2, sticky="n")
sat_2_period.grid(row=6, column=3, sticky='n')



txaperture_entry.grid(row=4, column=1, sticky='n')
txaperture_lab.grid(row=3, column=1, sticky='n')
txtransmission_entry.grid(row=4, column=2, sticky='n')
txtransmission_lab.grid(row=3, column=2, sticky='n')
rxaperture_entry.grid(row=8, column=1, sticky='n')
rxaperture_lab.grid(row=7, column=1, sticky='n')
rxtransmission_entry.grid(row=8,column=2, sticky='n')
rxtransmission_lab.grid(row=7,column=2, sticky='n')
rxobscuration_entry.grid(row=8,column=3, sticky='n')
rxobscuration_lab.grid(row=7,column=3, sticky='n')
rxsmfil_entry.grid(row=10,column=1, sticky='n')
rxsmfil_lab.grid(row=9,column=1, sticky='n')
rxswitchil_entry.grid(row=10,column=2, sticky='n')
rxswitchil_lab.grid(row=9,column=2, sticky='n')
smfcoupling_entry.grid(row=12,column=1, sticky='n')
smfcoupling_lab.grid(row=11,column=1, sticky='n')
laserpower_entry.grid(row=12, column=2, sticky='n')
laserpower_lab.grid(row=11, column=2, sticky='n')

station_lab.grid(row=13, column=1)
stationx_entry.grid(row=13, column=2, sticky="w")
stationy_entry.grid(row=13, column=3, sticky="w")
stationz_entry.grid(row=13, column=4, sticky="w")

filename_entry.grid(row=11, column=3)
filename_lab.grid(row=10, column=3)



# button_quit.pack(side=tkinter.BOTTOM)
# button_save.pack(side=tkinter.BOTTOM)
# #slider_update.pack(side=tkinter.BOTTOM)
# toolbar.pack(side=tkinter.BOTTOM, fill=tkinter.X)

# canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=True)

tkinter.mainloop()