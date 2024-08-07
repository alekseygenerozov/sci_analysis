import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def annotate_multiple_ecdf(datasets, labels, x_offset=0.0, y_offset=0.0, ha='left', levels=None, linestyles=None, ax=None):
    """
    Annotate multiple ECDF plots with one annotation per line.

    Parameters:
    datasets (list of array-like): List of datasets for the ECDF plots.
    labels (list of str, optional): List of labels for each dataset.
    """
    if ax is None:
        ax = plt.gca()

    if levels is None:
        levels = [80] * len(labels)

    if linestyles is None:
        linestyles = ["-"] * len(labels)

    for i, data in enumerate(datasets):
        # Create the ECDF plot for each dataset
        l1 = sns.ecdfplot(data, linestyle=linestyles[i], ax=ax)

        # Get the ECDF data points
        ecdf_data = np.sort(data)
        ecdf_y = np.arange(1, len(data) + 1) / len(data)

        # Select a middle x-coordinate for annotation
        x_coord = np.percentile(data, levels[i])
        y_coord = np.interp(x_coord, ecdf_data, ecdf_y)

        # Add the annotation
        ax.annotate(f'{labels[i]}',
                     xy=(x_coord + x_offset, y_coord + y_offset), color=l1.lines[-1].get_color(), ha=ha)