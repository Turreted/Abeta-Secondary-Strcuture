import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import re
from matplotlib.pyplot import figure


figure(figsize=(8, 4.5), dpi=500)

df = pd.read_csv("/Users/gideon/Desktop/trajrmsd.dat", delimiter="\s+")
fig, ax = plt.subplots()

mol = "mol5"

# assumes 2000 fs per frame
fs_per_frame = 500

ax.plot(df[mol])
labels = [int(re.sub(u"\u2212", "-", i.get_text()))/fs_per_frame for i in ax.get_xticklabels()]
ax.set_xticklabels(labels)

plt.xlabel("Time (ns)")
plt.ylabel("RMSD (A)")
plt.title("RMSD Of Complex")
plt.savefig("RMSD.png", bbox_inches='tight', transparent=True, pad_inches=0, dpi=500) 
