import matplotlib.pyplot as plt
import numpy as np
import json

plot = json.load(open("/Users/gideon/Desktop/confidence_model_5_multimer_v3_pred_4.json" , "r"))["confidenceScore"]
nchains = 8
nres = 40
total_len = nchains * nres

fig, ax = plt.subplots()
ax.plot(plot)
plt.ylim(0.0, 110.0)

xs = np.linspace(1, 21, 81)

for i in range(1, nchains):  
       line = plt.vlines(x=nres*i, ymin=0, ymax=len(xs), colors='green', ls=':', lw=2, label='vline_single - full height')

ax.set(xlabel='Residue Number', ylabel='pLDDT',
       title='Predicted LDDT Score per Residue')
ax.legend(['pLDDT', 'Chain Cutoff'])

fig.savefig("plddt.png", dpi=400, bbox_inches='tight', transparent=True, pad_inches=0)
