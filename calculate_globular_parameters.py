#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 20:07:44 2024

@author: misha
"""

import MDAnalysis as mda
from MDAnalysis import transformations
from matplotlib import pyplot as plt
import argparse
import numpy as np

atom_types = ['1', '2']
block_name = { '1' : 'VCL', '2' : 'VI'}
line_style = {'1' : '-', '2' : '--'}

parser = argparse.ArgumentParser(
    description = 'Calculate globule parameters')

parser.add_argument('file', metavar = 'DATA', type = str, nargs = '+',
    help = 'simulation data')

parser.add_argument('--nbins', metavar = 'NBINS', type = int,
                    nargs = '?', default = 20, help = 'histogram bins')

args = parser.parse_args()

for file in args.file:

    u = mda.Universe(file)
    
    n_atom = u.atoms.n_atoms

    unwrap = transformations.unwrap(u.atoms)
    u.trajectory.add_transformations(unwrap)

    com = u.atoms.center_of_mass()
    rg = u.atoms.radius_of_gyration()

    print(f'N {n_atom} Rg = {rg:.3}')

    for atom_type in atom_types:
        atoms =u.select_atoms(f'type {atom_type}')
        positions = atoms.positions
        distances = np.linalg.norm(positions - com, axis = 1)
        max_dist = np.max(distances)
        print(f'N {n_atom} Max. {atom_type} = {max_dist:.3}')
        frequencies, bin_edges = np.histogram(distances, bins = args.nbins)
        bin_volumes_cumul = 4./3.*np.pi*np.power(bin_edges, 3)
        bin_volumes = bin_volumes_cumul[1:] - bin_volumes_cumul[:-1]
        bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
        densities = frequencies / bin_volumes
        plt.plot(bin_centers, densities / np.sum(densities), 
                 linestyle = line_style[atom_type],
                 label = f'N {n_atom} {atom_type}')
    
plt.legend()
plt.show()