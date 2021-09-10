import matplotlib
import matplotlib.pyplot as plt

import numpy as np

# Make big figures:
gFIG_SIZE_in = [14, 10]

# Set plotting defaults:
gPLOT_PARAMS = {
    "legend.fontsize": "x-large",
    "figure.figsize": gFIG_SIZE_in,
    "axes.labelsize": "x-large",
    "axes.titlesize": "x-large",
    "xtick.labelsize": "x-large",
    "ytick.labelsize": "x-large"
}
matplotlib.rcParams.update(gPLOT_PARAMS)
plt.rcParams.update(gPLOT_PARAMS)

# Some single-place definitions of sizes for plots / figures:
gFONT_SIZE_UNITS = "pt"
gTITLE_FONT_SIZE = 36
gAXIS_LABEL_FONT_SIZE = 24
gTICK_LABEL_FONT_SIZE = 16
gTEXT_FONT_SIZE = 16

# To track global figure number creation:
gFIG_NUM = 0


def fix_plot_visuals(fig,
                     titlesize=gTITLE_FONT_SIZE,
                     labelsize=gAXIS_LABEL_FONT_SIZE,
                     ticklabelsize=gTICK_LABEL_FONT_SIZE,
                     textsize=gTEXT_FONT_SIZE,
                     tight_rect=None):
    """Fix the plot elements to be appropriate sizes for a slide / presentation."""

    if not textsize:
        textsize = ticklabelsize

    for ax in fig.get_axes():

        for ticklabel in (ax.get_xticklabels()):
            ticklabel.set_fontsize(ticklabelsize)
        for ticklabel in (ax.get_yticklabels()):
            ticklabel.set_fontsize(ticklabelsize)
        for c in ax.get_children():
            if c.__class__ == matplotlib.text.Text:
                c.set_fontsize(textsize)

        ax.xaxis.get_label().set_fontsize(labelsize)
        ax.yaxis.get_label().set_fontsize(labelsize)
        ax.title.set_fontsize(titlesize)

    for c in fig.get_children():
        if c.__class__ == matplotlib.legend.Legend:
            c.prop.set_size(ticklabelsize)
            c.get_title().set_size(ticklabelsize)

    if tight_rect:
        fig.tight_layout(rect=tight_rect)
    else:
        fig.tight_layout()


def save_figure(fig=None, name=None, prefix=None, suffix=None, fig_dir=None):
    """Save the current figure in the pyplot buffer to disk as both svg and png."""
    global gFIG_NUM

    if not fig_dir:
        fig_dir = ""
    elif not fig_dir.endswith("/"):
        fig_dir = f"{fig_dir}/"

    if not fig:
        fig = plt.gcf()

    if not name:
        # Get our name from the figure we're given:
        name = [ax.get_title() for ax in fig.axes if ax.get_title() and len(ax.get_title()) > 0][0]

    # Sanitize our name for writing to disk:
    name = name.replace(" ", "_").replace("\n", "_").replace("\t", "_").replace("=", "-").replace("(", "").replace(")", "")

    # Add prefix and suffix to the name if we have them defined:
    if prefix:
        name = f"{prefix}_{gFIG_NUM:02d}_{name}"
    if suffix:
        name = f"{name}_{suffix}"
    if not prefix:
        name = f"{gFIG_NUM:02d}_{name}"
    gFIG_NUM += 1

    # Save the figure:
    plt.savefig(f"{fig_dir}{name}.svg")
    plt.savefig(f"{fig_dir}{name}.png")


def signif(x, p=4):
    """Round the given value(s) in x to the number of significant figures given by p."""
    x_positive = np.where(np.isfinite(x) & (x != 0), np.abs(x), 10**(p-1))
    mags = 10 ** (p - 1 - np.floor(np.log10(x_positive)))
    return np.round(x * mags) / mags


def get_heat_cmap(cmap_name="jet", do_reverse=False):
    """Get the color map of the given name.  Reverse the colors if requested."""
    cmap = plt.get_cmap(cmap_name)
    if do_reverse:
        cmap = cmap.reversed()

    return cmap


def get_zero_white_cmap(cmap_name=None, base_cmap=None):
    """Add white to the bottom of a given colormap so that the minimum value in the map is white."""
    # create colormap
    # Taken from: https://stackoverflow.com/a/49367945
    # ---------------

    # create a colormap that consists of
    # - 1/5 : custom colormap, ranging from white to the first color of the colormap
    # - 4/5 : existing colormap

    # set upper part: 4 * 256/4 entries

    if cmap_name and base_cmap:
        raise RuntimeError("cmap_name and base_cmap are mutually exclusive.")

    if not base_cmap:
        base_cmap = plt.get_cmap(cmap_name)

    upper = base_cmap(np.arange(256))

    # set lower part: 1 * 256/4 entries
    # - initialize all entries to 1 to make sure that the alpha channel (4th column) is 1
    lower = np.ones((int(256 / 4), 4))
    # - modify the first three columns (RGB):
    #   range linearly between white (1,1,1) and the first color of the upper colormap
    for i in range(3):
        lower[:, i] = np.linspace(1, upper[0, i], lower.shape[0])

    # combine parts of colormap
    cmap = np.vstack((lower, upper))

    # convert to matplotlib colormap
    cmap = matplotlib.colors.ListedColormap(cmap, name=f"{cmap_name}ZeroWhite", N=cmap.shape[0])

    return cmap
