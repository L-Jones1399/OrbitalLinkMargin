import os
import pandas as pd
from io import StringIO
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

folder_path = r"C:\Users\zeus1\Documents\RMIT\Honours\Project\FSOC Code\DataFiles\IntersatPlot"

DistanceList = []
AngleList = []
LinkList = []
for filename in os.listdir(folder_path):
    if filename.endswith(".txt"):
        file_path = os.path.join(folder_path, filename)
        with open(file_path, "r", encoding="utf-8") as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            if "-----End-----" in line:
                start_index = i + 1
                break
        
        data = pd.read_csv(
        StringIO("".join(lines[start_index:])),
        sep=r"\t",
        engine="python",
        header=None,
        names=["time", "angle", "distance", "link"]
        )

        filtered_df = data[data['link'] != None]
        filtered_df["time"] = pd.to_datetime(filtered_df["time"])
        plt.plot(filtered_df["time"], filtered_df['link'], label= filename.split('_', 1)[0])

plt.xlabel("Time (HH:MM)")
plt.ylabel("Link Margin (dB)")
plt.legend()
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
plt.savefig("IntersatPLots.pdf")


plt.show()
