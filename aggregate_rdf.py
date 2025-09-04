#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 19 20:17:45 2025

@author: misha
"""

import sys
import numpy as np
from copoly.misc import chain_ix, bring_chains_together
from mouse2.lib.neighbor import calculate_distances

import pdb

def aggregate_rdf(ag, ixs_per_chain, selection = None,
                  r_max = 50, n_bins = 25):

    u = ag.universe

    chain_indices = [chain_ix(ix,ixs_per_chain) 
                             for ix in ag.atoms.ix]

    unique_chains = list(set(chain_indices))
    unique_chains.sort()

    size_of_ag = len(unique_chains)

    bring_chains_together(ag, ixs_per_chain)

    rcm = ag.center_of_mass()
    
    if selection is not None:
        selected_ag = ag.select_atoms(selection)
        selection = selection.replace(" ", "_")
    else:
        selected_ag = ag
        selection = "all"

    x = selected_ag.positions[:,0]
    y = selected_ag.positions[:,1]
    z = selected_ag.positions[:,2]
    relative_distances = calculate_distances([x, y, z], rcm,
                                                     u.dimensions)
    max_distance = np.max(relative_distances)
    if max_distance > np.min(u.dimensions[:3])/2:
        #pdb.set_trace()
        sys.stderr.write(
                    "Distance from atom to c.m. exceeds cell/2")
        return { "n" : size_of_ag,
                f"rdf_{selection}_valid" : "False",
                f"rdf_{selection}_bins" : "size_exceeded",
                f"rdf_{selection}_densities" : "size_exceeded",
                }
    distr, bins = np.histogram(relative_distances, 
                                   range = (0, r_max),
                                   bins = n_bins)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    sph_volumes = 4/3*np.pi*bins**3
    bin_volumes = sph_volumes[1:] -sph_volumes[:-1]
    densities = np.divide(distr, bin_volumes)

    return { "n" : size_of_ag,
            f"rdf_{selection}_valid" : "True",
            f"rdf_{selection}_bins" : list(bin_centers),
            f"rdf_{selection}_densities" : list(densities),
            }
