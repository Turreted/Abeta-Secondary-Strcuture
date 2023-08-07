import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable

data_x = [0,1,2,3]
data_hight = [60,60,80,100]
data_color = [1000.,500.,1000.,900.]


data_color = [x / max(data_color) for x in data_color]
fig, ax = plt.subplots(figsize=(15, 4))

my_cmap = plt.cm._colormaps['GnBu']
colors = my_cmap(data_color)
rects = ax.bar(data_x, data_hight, color=colors)

sm = ScalarMappable(cmap=my_cmap, norm=plt.Normalize(0,max(data_color)))
sm.set_array([])

cbar = plt.colorbar(sm, ax=ax)
cbar.set_label('Color', rotation=270,labelpad=25)

plt.xticks(data_x)    
plt.ylabel("Y")

plt.show()