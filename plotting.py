import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

from labelLine import labelLines


def annotate_multiple_ecdf(datasets, labels, x_offset=None, y_offset=0.0, ha=None, levels=None, colors=None, alphas=None, linestyles=None, ax=None, fontsize=None):
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

    if colors is None:
        colors = [None] * len(labels)

    if alphas is None:
        alphas = [None] * len(labels)

    if ha is None:
        ha = ["left"] * len(labels)

    if x_offset is None:
        x_offset = [0] * len(labels)

    for i, data in enumerate(datasets):
        # Create the ECDF plot for each dataset
        l1 = sns.ecdfplot(data, linestyle=linestyles[i], color=colors[i], alpha=alphas[i], ax=ax)

        # Get the ECDF data points
        ecdf_data = np.sort(data)
        ecdf_y = np.arange(1, len(data) + 1) / len(data)

        # Select a middle x-coordinate for annotation
        x_coord = np.percentile(data, levels[i])
        y_coord = np.interp(x_coord, ecdf_data, ecdf_y)

        # Add the annotation
        ax.annotate(f'{labels[i]}',
                     xy=(x_coord + x_offset[i], y_coord + y_offset), color=l1.lines[-1].get_color(), ha=ha[i], fontsize=fontsize)

def scaled_kde(data, desired_peak=1, npts=1000):
    # Sample data
    # Compute KDE using scipy
    kde = gaussian_kde(data)

    # Generate x values for the KDE plot
    x_values = np.linspace(min(data), max(data), npts)
    y_values = kde(x_values)

    # Find the peak of the KDE
    peak_value = np.max(y_values)
    # Desired peak value

    # Scale the y_values so the peak matches the desired value
    scaled_y_values = (y_values / peak_value) * desired_peak

    return x_values, scaled_y_values

def generate_ecdf(data):
    """
    Generate the x and y coordinates for an ECDF plot.
    """
    # Sort the data
    sorted_data = np.sort(data)
    # The y-axis values are the cumulative probabilities
    y_values = np.arange(1, len(sorted_data) + 1) / len(sorted_data)
    # The x-axis values are the sorted data points
    return sorted_data, y_values

def fill_between_ecdf(data1, data2, fig=None, fill_color=None):
    x1, y1 = generate_ecdf(data1)
    x2, y2 = generate_ecdf(data2)

    # 3. Use a common x-axis for interpolation
    # Create a comprehensive set of x-values from both datasets.
    # This is crucial for filling between step functions correctly.
    x_all = np.unique(np.concatenate([x1, x2]))

    # Use NumPy's interpolation to create common y-values for each ECDF
    y1_interp = np.interp(x_all, x1, y1)
    y2_interp = np.interp(x_all, x2, y2)

    if fig is None:
        # 4. Plot the ECDFs and fill the area between them
        fig = plt.figure()

    # Plot the individual ECDF curves as step functions
    # plt.step(x1, y1)
    # plt.step(x2, y2)

    # Fill the area between the two interpolated curves
    plt.fill_between(x_all, y1_interp, y2_interp, alpha=0.3, label='Area Between ECDFs', color=fill_color)

##Plot ecdf with error region...

def scatter_plus_regression(X, y, ax=None):
    """
    Make a scatter plot with a linear regression

    x, y--1D Arrays
    """
    ##Fitting linear regression...
    reg1 = LinearRegression()
    X = X.reshape(-1, 1)
    reg1.fit(X, y)
    ##Scoring
    y_true = y
    y_pred = reg1.predict(X)
    r2s = r2_score(y_true, y_pred)
    absc = np.linspace(X.min(), X.max())
    absc = absc.reshape(-1, 1)
    ords = reg1.predict(np.log10(absc))
    ##Plotting
    if ax is None:
        fig,ax = plt.subplots()
    l1, = ax.plot(absc, ords, "r--", label=rf"Slope=${reg1.coef_[0]:.2f}$"+ "\n" + r"$R^2$=" + f"{r2s:.2f}")
    labelLines([l1], align=False, y_offset=0.7, fontsize=16, xvals=[0.8])


