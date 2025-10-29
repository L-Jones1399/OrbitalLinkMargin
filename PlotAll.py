import os
import pandas as pd
from io import StringIO
import matplotlib.pyplot as plt

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

        #Filter to keep angles >5
        filtered_df = data[data['angle'] >= 5]
        #Grab a specific column and turn it into a list
        distance = filtered_df["distance"].values.tolist()
        angle = filtered_df["angle"].values.tolist()
        link = filtered_df['link'].values.tolist()
        DistanceList += distance
        AngleList += angle
        LinkList += link


plt.scatter(DistanceList, AngleList, c = LinkList, cmap="viridis")
plt.xlabel("Distance (km)")
plt.ylabel("Elevation Angle (degrees)")
plt.colorbar(label='Link Margin (dB)')

plt.savefig("Colourplot.pdf")
plt.show()