#!/usr/bin/env python

import numpy as np
import spktools.input as spkin
import spktools.texture.ori_gen as spktex
import spktools.special_grain as spkspec

infile = 'init1.potts'
outfile = 'input.potts'
nsites, x_dim, y_dim, z_dim, spins, pos = spkin.parse_potts(infile)
dimensionality = spkin.get_dimensionality(x_dim, y_dim, z_dim)
print('loaded data')

spkspec.create_special_grain(spins, pos,
                             x_dim[1]-x_dim[0],
                             y_dim[1]-y_dim[0],
                             z_dim[1]-z_dim[0],
                             size_advantage=6)
print('created special grain')

spin_map = spkspec.special_spin_map(spins)

print('created spin map')

high_fraction = 0.2
ori = spktex.random_two_tone(spin_map, high_fraction)
print('created texture map')

mtab_file = 'potts-test.mtab'
spkspec.write_mtab(ori, mtab_file, 1000)

with open(outfile, 'w') as of:
  spkin.print_header(of, dimensionality, nsites, x_dim, y_dim, z_dim)
  spkin.print_Sites(of, spins, spin_map)
print('wrote spin file')

# colors = np.copy(spins)
# for i in np.arange(spins.size):
#   if spins[i] == 1:
#     colors[i] = 0
#   else:
#     colors[i] = ori[spin_map[spins[i]]]

# print(len(spin_map))
# import matplotlib.pyplot as plt
# plt.imshow(colors.reshape(((x_dim[1]-x_dim[0]), (y_dim[1]-y_dim[0]))), interpolation='None')
# plt.colorbar()
# plt.show()
