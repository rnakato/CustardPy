import numpy as np

def generate_cmap(colors):
    from matplotlib.colors import LinearSegmentedColormap
    values = range(len(colors))

    vmax = np.ceil(np.max(values))
    color_list = []
    for v, c in zip(values, colors):
        color_list.append((v/vmax, c))
    return LinearSegmentedColormap.from_list('custom_cmap', color_list)
