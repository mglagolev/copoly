#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 19 17:33:12 2025

@author: misha
"""
import numpy as np
from copoly.asa import calculate_asa
from copoly.misc import chain_ix

def aggregate_asa(ag, ixs_per_chain, r_probe, n_points, radii):
        u = ag.universe
        chain_indices = [chain_ix(ix,ixs_per_chain) 
                             for ix in ag.atoms.ix]
        size_of_ag = len(set(chain_indices))

        asas = calculate_asa(ag, r_probe, n_points,
                                 u.dimensions, radii)

        return { "n" : size_of_ag,
                "s1_mean" : np.mean(asas["1"]),
                "s1_sum" : np.sum(asas["1"]),
                "s2_mean" : np.mean(asas["2"]),
                "s2_sum" : np.sum(asas["2"]) }