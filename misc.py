#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 18 15:30:44 2025

@author: misha
"""
import numpy as np
from mouse2.lib.aggregation import determine_aggregates
import pdb

def ixs_to_query(ixs):
    query = "index " + " ".join(map(str, ixs))
    return query


def chains_list(u):
    aggregation_data = determine_aggregates(u, r_neigh = 0,
                                            bonded_as_neighbors = True)
    aggregates_list = list(aggregation_data['data'].values())[0]
    return aggregates_list


def chain_ix(atom_ix, chains_atoms_ixs):
    """
    Takes the "chains_atoms_ixs" list of lists.
    Returns the index of the list in the list of lists which contains
    atom_ix

    """
    ixs = [x[0] for x in enumerate(chains_atoms_ixs) if atom_ix in x[1]]
    if len(ixs) == 1:
        return ixs[0]
    elif len(ixs) == 0:
        raise NameError(f"Atom with index {atom_ix} not found in the chains list")
    elif len(ixs) > 1:
        raise NameError(f"Atom with index {atom_ix} listed in the chains {ixs}")
        
def bring_chains_together(ag, ixs_per_chain):
    
    box_x = ag.universe.dimensions[0]
    box_y = ag.universe.dimensions[1]
    box_z = ag.universe.dimensions[2]

    chain_indices = list(set([chain_ix(ix,ixs_per_chain) 
                         for ix in ag.atoms.ix]))
    first_chain = ag.select_atoms(ixs_to_query(
        ixs_per_chain[chain_indices[0]]))
    #pdb.set_trace()
    first_chain_com = first_chain.center_of_mass()
    for another_chain_ix in chain_indices[1:]:
        #pdb.set_trace()
        another_chain = ag.select_atoms(ixs_to_query(
            ixs_per_chain[another_chain_ix]))
        another_chain_com = another_chain.center_of_mass()
        dx = first_chain_com[0] - another_chain_com[0]
        jumps_x = box_x * round(dx / box_x)
        dy = first_chain_com[1] - another_chain_com[1]
        jumps_y = box_y * round(dy / box_y)
        dz = first_chain_com[2] - another_chain_com[2]
        jumps_z = box_z * round(dz / box_z)
        jump =  np.array([jumps_x, jumps_y, jumps_z])
        another_chain.positions -= jump