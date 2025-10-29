import matplotlib.pyplot as plt
import pandas as pd
from io import StringIO
#import matplotlib.dates as mdates
from matplotlib.ticker import MultipleLocator

# Path to your file
file_path = r"C:\Users\zeus1\Documents\RMIT\Honours\Project\FSOC Code\DataFiles\INTERSAT_TAKE4.txt"

# Read and find where the '-----End-----' line occurs
with open(file_path, "r") as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    if "-----End-----" in line:
        start_index = i + 1
        break

# Read data after that line
data = pd.read_csv(
    StringIO("".join(lines[start_index:])),
    sep=r"\t",
    engine="python",
    header=None,
    names=["time", "angle", "distance", "link"]
)



# Convert time column to datetime
#data["time"] = pd.to_datetime(data["time"], format= '%Y-%m-%d %H:%M:%S.%f')
print(type(data['time'][0]))
data["time"] = pd.to_datetime(data["time"])
xyz = data['time'].tolist()
print(type(xyz[0]))
print(data['time'])
print(data['link'])

# Plot using true datetime axis
plt.figure(figsize=(12, 6))
plt.plot(data["time"], data["link"], label="Column 4", linewidth=1.5)

# Format x-axis to show only HH:MM
#plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
#plt.gca().yaxis.set_major_locator(MultipleLocator(5))
#plt.xlim(pd.Timestamp("2022-11-29 08:33"), pd.Timestamp("2022-11-29 08:41"))
#plt.ylim(-60, -45)

plt.xlabel("Time (HH:MM)")
plt.ylabel(r"dB W/${\text{m}^2}$")
plt.grid(True)
#plt.legend()
plt.tight_layout()
plt.savefig("Schieler_comp.png")
plt.show()
