#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 19 16:40:05 2025

@author: misha
"""

import numpy as np
from copoly.misc import chain_ix, bring_chains_together
from mouse2.lib.neighbor import calculate_distances

def rh(ag):
    u = ag.universe
    coordinates = ag.atoms.positions

    sum_inverse_distances = 0
    for j in range(ag.n_atoms - 1):
        first_atom = coordinates[j]
        second_atoms = coordinates[j+1:,:]
        coordinates_x = second_atoms[:,0]
        coordinates_y = second_atoms[:,1]
        coordinates_z = second_atoms[:,2]
        distances = calculate_distances([coordinates_x,
                                         coordinates_y,
                                         coordinates_z],
                                        first_atom,
                                        u.dimensions)
        inverse_distances = np.power(distances, -1)
        sum_inverse_distances += np.sum(inverse_distances)

    rh = (sum_inverse_distances/(ag.n_atoms**2))**(-1)
    return rh


def gyration_tensor(ag):
    u = ag.universe
    s = np.zeros((3,3))
    coords = ag.atoms.positions
    for i in range(3):
        for j in range(3):
            for n in range(ag.n_atoms - 1):
                box1 = u.dimensions[i]
                box2 = u.dimensions[j]
                ref_coords = coords[n]
                other_coords = coords[n+1:, :]
                ref_coord1 = ref_coords[i]
                ref_coord2 = ref_coords[j]
                other_coords1 = other_coords[:,i]
                other_coords2 = other_coords[:,j]
                distances1 = calculate_distances([other_coords1,
                                            np.zeros_like(other_coords1),
                                            np.zeros_like(other_coords1)],
                                                 [ref_coord1, 0., 0],
                                                 [box1, box1, box1])
                distances2 = calculate_distances([other_coords2,
                                            np.zeros_like(other_coords2),
                                            np.zeros_like(other_coords2)],
                                                 [ref_coord2, 0., 0],
                                                 [box2, box2, box2])
                s[i][j] += np.sum(distances1 * distances2)
    s = s / ag.n_atoms**2
    return s


def aggregate_rg_rh(ag, ixs_per_chain):

    chain_indices = [chain_ix(ix,ixs_per_chain) 
                             for ix in ag.atoms.ix]

    unique_chains = list(set(chain_indices))
    unique_chains.sort()

    size_of_ag = len(unique_chains)

    bring_chains_together(ag, ixs_per_chain)
    
    hydrodynamic_radius = rh(ag)

    gyr_tens = gyration_tensor(ag)
    
    eival = np.linalg.eigvalsh(gyr_tens)
    
    rg_sq = np.sum(eival)
    
    rg = np.sqrt(rg_sq)

    normalized_asphericity = (eival[2] - 0.5 * (eival[0] + eival[1]))/rg_sq

    anisotropy = 1.5 * np.sum(eival**2) / (np.sum(eival))**2 - 0.5


    radius_of_gyration = ag.radius_of_gyration()
    aggregate = {
                     "n" : size_of_ag,
                     "rg" : radius_of_gyration,
                     "rg_pairwise" : rg,
                     "rh" : hydrodynamic_radius,
                     "chains_ixs" : unique_chains,
                     "n_atoms" : ag.n_atoms,
                     "rg/rh" : rg / hydrodynamic_radius,
                     "asphericity" : normalized_asphericity,
                     "anisotropy" : anisotropy }
    return aggregate