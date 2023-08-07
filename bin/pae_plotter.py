import matplotlib.pyplot as plt
import numpy as np
import json

pae_data = json.load(open("/Users/gideon/Desktop/pae_model_5_multimer_v3_pred_4.json", "r"))[0]

# Generate dictionary for predicted aligned error results from pkl file
pae_outputs = {'protein': (
    pae_data['predicted_aligned_error'],
    pae_data['max_predicted_aligned_error']
)}

print(pae_data['predicted_aligned_error'],)

# Output file_path for the plot
pae_output = 'pae.png'

# Plot predicted align error results for each aligned residue
pae, max_pae = list(pae_outputs.values())[0]

fig = plt.figure()  # generate figure
fig.set_facecolor('white')  # color background white
plt.imshow(pae, vmin=0., vmax=max_pae)  # plot pae
plt.colorbar(fraction=0.46, pad=0.04)  # create color bar

plt.title('Predicted Aligned Error')  # plot title
plt.xlabel('Scored residue')  # plot x-axis label
plt.ylabel('Aligned residue')  # plot y-axis label

plt.savefig(pae_output, dpi=1000, bbox_inches='tight', transparent=True, pad_inches=0)  # save plot to output directory

print('\n predicted aligned error plotted\n')