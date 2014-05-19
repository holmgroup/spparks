#!/usr/bin/env python

import numpy as np
import spktools.input as spkin
from scipy.ndimage.filters import generic_filter

def mode_filter(buffer):
  """ return the most common value in the neighborhood """
  return np.bincount(buffer.astype(np.int64)).argmax()

def diff_filter(buffer):
  return (buffer != buffer[4]).sum()

infile = 'dump1.potts'
nsites, x_dim, y_dim, z_dim, spins, pos = spkin.parse_potts(infile)
dimensionality = spkin.get_dimensionality(x_dim, y_dim, z_dim)

spins = spins.reshape(((x_dim[1]-x_dim[0]), (y_dim[1]-y_dim[0])))

cross = np.array([[0,1,0],[1,1,1],[0,1,0]])
X = np.array([[1,0,1],[0,1,0],[1,0,1]])
square = np.array([[1,1,1],[1,1,1],[1,1,1]])
filtered = generic_filter(spins, mode_filter, footprint=X, mode='wrap')

grain_boundary = generic_filter(filtered, diff_filter, footprint=square, mode='wrap')

import matplotlib.pyplot as plt
plt.imshow(spins, interpolation='None')
plt.colorbar()
plt.show()

plt.imshow(filtered, interpolation='None')
plt.colorbar()
plt.show()

microstructure = filtered
microstructure[grain_boundary > 2 ] = 0

plt.imshow(microstructure, interpolation='None')
plt.show()
