import matplotlib.pyplot as plt
import constants
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

newcolors = [np.array([1, a/255.0, a/255.0, 1]) for a in np.flip(np.arange(256))]
newcmp = ListedColormap(newcolors)

fig, ax= plt.subplots(1,1,figsize=(10,10))
pos=ax.imshow(np.array([[0,1.5,3.0], [0,1.5,3.0]]), cmap=newcmp, interpolation='none')

fig.colorbar(pos, ax=ax)

plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "color_legend.png"))