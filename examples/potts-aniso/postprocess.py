#!/usr/bin/env python

import numpy as np
import spktools.input as spkin

infile = 'dump1.potts'
nsites, x_dim, y_dim, z_dim, spins, pos = spkin.parse_potts(infile)
dimensionality = spkin.get_dimensionality(x_dim, y_dim, z_dim)


import matplotlib.pyplot as plt
plt.imshow(spins.reshape(((x_dim[1]-x_dim[0]), (y_dim[1]-y_dim[0]))), interpolation='None')
plt.colorbar()
plt.show()
