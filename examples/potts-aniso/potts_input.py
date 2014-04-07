#!/usr/bin/env python
import numpy as np
import random

"""
  Process a SPPARKS potts dump file for use with potts-aniso.
  Reduce spins to the range 1 - nspins for use with lookup-tables
"""

def print_header(outfile, dim=3, N=1, x=[-.5,.5], y=[-.5,.5], z=[-.5,.5]):
  outfile.write('# SPPARKS sites file\n')
  outfile.write('{0} dimension\n'.format(dim))
  outfile.write('{0} sites\n'.format(N))
  outfile.write('{0} {1} xlo xhi\n'.format(x[0], x[1]))
  outfile.write('{0} {1} ylo yhi\n'.format(y[0], y[1]))
  outfile.write('{0} {1} zlo zhi\n'.format(z[0], z[1]))
  outfile.write('\n')

def print_Sites(outfile, spins, spin_map):
  outfile.write('Values\n')
  outfile.write('# segfault-preventing comment\n')
  for i, s in enumerate(spins):
    outfile.write('{0} {1}\n'.format(i+1, spin_map[s]))

def get_dim(dim_line):
  min_val, max_val = dim_line.split()
  return [float(min_val), float(max_val)]

def parse_potts(potts_file):
  nsites = 0
  spin = []
  x_dim, y_dim, z_dim = 0, 0, 0
  with open(potts_file) as pf:
    extract = False
    for line in pf:
      if 'NUMBER OF ATOMS' in line:
        nsites = int(next(pf))

      if 'BOX' in line:
        x_dim = get_dim(next(pf))
        y_dim = get_dim(next(pf))
        z_dim = get_dim(next(pf))
        # skip forward two lines:
        next(pf); line = next(pf)
        spin = np.zeros(nsites, dtype='int')
        extract = True

      if extract:
        # 0  1    2 3 4
        # id spin x y z
        # site_id i is offset by 1, hence spin.size == nsites+1
        i, s, x, y, z = map(int,line.split())
        spin[i-1] = s

  return nsites, x_dim, y_dim, z_dim, spin

def create_spin_map(spins):
  spin_set = set()
  for s in spins:
    spin_set.add(s)
  # if any spins in the original microstructure are 0
  # this is going to cause a bug -- 0 should be special.
  max_spin = len(spin_set) + 1
  
  spin_map = {}
  uniques = set([i for i in range(1, max_spin)])
  random.seed()
  for s in spin_set:
    if s not in spin_map:
      # hack-ish way of drawing from a set without replacement
      new_spin = random.sample(uniques, 1)[0]
      uniques.remove(new_spin)
      spin_map[s] = new_spin

  return spin_map

def get_dimensionality(x_dim, y_dim, z_dim):
  D = 0
  x_range = x_dim[1] - x_dim[0]
  y_range = y_dim[1] - y_dim[0]
  z_range = z_dim[1] - z_dim[0]
  for r in [x_range, y_range, z_range]:
    if r > 1.0:
      D += 1

  return D


if __name__ == '__main__':
  infile = 'init1.potts'
  outfile = 'input.potts'
  nsites, x_dim, y_dim, z_dim, spins = parse_potts(infile)
  dimensionality = get_dimensionality(x_dim, y_dim, z_dim)
                  
  spin_map = create_spin_map(spins)

  with open(outfile, 'w') as of:
    print_header(of, dimensionality, nsites, x_dim, y_dim, z_dim)
    print_Sites(of, spins, spin_map)
