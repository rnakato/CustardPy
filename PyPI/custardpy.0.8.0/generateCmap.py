import numpy as np
from matplotlib.colors import LinearSegmentedColormap

def generate_cmap(colors):
    return LinearSegmentedColormap.from_list('custom_cmap', colors)
    
def generate_cmap_old(colors):
    values = range(len(colors))

    vmax = np.ceil(np.max(values))
    color_list = []
    for v, c in zip(values, colors):
        color_list.append((v/vmax, c))
    return LinearSegmentedColormap.from_list('custom_cmap', color_list)
