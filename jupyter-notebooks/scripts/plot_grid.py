import warnings
warnings.simplefilter('ignore')

import numpy
from matplotlib import pyplot
import seaborn
clear_bkgd = {'axes.facecolor':'none', 'figure.facecolor':'none'}
seaborn.set(style='ticks', context='notebook', rc=clear_bkgd)

import pygridgen

def plot_grid(grid, ax):
    ax.plot(grid.x.flatten(), grid.y.flatten(), 'k.', label='Grid nodes', zorder=5)
    ax.set_aspect('equal')
    ax.set_xlim([0, 4])
    ax.set_ylim([0, 4])
    ax.plot(grid.xbry, grid.ybry, '-', color='0.5', zorder=0)
    pos = numpy.nonzero(grid.beta == 1)
    neg = numpy.nonzero(grid.beta == -1)
    ax.plot(x[pos], y[pos], 'go', label='Positive', zorder=2, alpha=0.5)
    ax.plot(x[neg], y[neg], 'rs', label='Negative', zorder=2, alpha=0.5)
    ax.legend(numpoints=1)

