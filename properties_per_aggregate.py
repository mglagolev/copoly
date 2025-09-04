import sys
import json
import argparse
import MDAnalysis as mda
from MDAnalysis import transformations
import pandas as pd
from mouse2.lib.aggregation import determine_aggregates
from aggregate_asa import aggregate_asa
from rg_rh import aggregate_rg_rh
from aggregate_rdf import aggregate_rdf
from copoly.misc import chains_list

r_probe = 0.

radii = {
 '1': 1,
 '2': 1,
}

def properties_per_aggregate(u, ixs_per_chain = None,
                                 r_probe = None, n_points = None,
                                 radii = None,
                                 filename = None):
    """
    Calculate the properties independently for each aggregate,
    then return the data in a list of dictionaries containing also
    the sizes of the aggregates.
    """
    
    unwrap = transformations.unwrap(u.atoms)
    u.trajectory.add_transformations(unwrap)

    aggregate_properties = []

    result = determine_aggregates(u, r_neigh = 1.2)
    aggregates_list = list(result["data"].values())[0]

    for aggregate in aggregates_list:
        values = {"filename" : filename}
        atom_ixs_str = list(map(str, aggregate))
        index_query = " ".join(atom_ixs_str)
        ag = u.select_atoms(f"index {index_query}")
        values |= aggregate_asa(ag, ixs_per_chain,
                                         r_probe, n_points, radii)
        values |= aggregate_rg_rh(ag, ixs_per_chain)
        values |= aggregate_rdf(ag, ixs_per_chain)
        values |= aggregate_rdf(ag, ixs_per_chain, selection = "type 1")
        values |= aggregate_rdf(ag, ixs_per_chain, selection = "type 2")
        aggregate_properties.append(values)
    return aggregate_properties


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='Aggregate sizes in chains')

    parser.add_argument('--first', metavar='FIRST_STEP', type=int, nargs=1,
                        help='first step')

    parser.add_argument('--last', metavar='LAST_STEP', type=int, nargs=1,
                        help='last step')

    parser.add_argument('--step', metavar='STEP', type=int, nargs='?',
                        default = 1,
                        help='interval between the considered steps')

    parser.add_argument('--output', metavar='XLSX', type=str, nargs=1,
                        help='Excel file for output data')

    parser.add_argument('--sasa-points', metavar='n', type=int, nargs='?',
                        default = 960,
                        help='Number of points for SASA assessment')

    args = parser.parse_args()

    all_data = []

    for i in range(args.first[0], args.last[0] + 1, args.step):
        sys.stderr.write(f"{i}\n")
        filename = f"{i}.data"
        u = mda.Universe(filename)

        parameters = {
            "ixs_per_chain" : chains_list(u),
            "r_probe" : r_probe,
            "n_points" : args.sasa_points,
            "radii" : radii,
            "filename" : filename
                     }

        this_universe_data = properties_per_aggregate(u, **parameters)

        all_data += this_universe_data

    df = pd.DataFrame(all_data)
    try:
        df['chains_ixs'] = df['chains_ixs'].apply(json.dumps)
    except KeyError:
        pass
    try:
        df['rdf_type_1_bins'] = df['rdf_type_1_bins'].apply(json.dumps)
        df['rdf_type_1_densities']= df['rdf_type_1_densities'].apply(json.dumps)
    except KeyError:
        sys.stderr.write("type1 rdfs not found")
    try:
        df['rdf_type_2_bins'] = df['rdf_type_2_bins'].apply(json.dumps)
        df['rdf_type_2_densities'] = df['rdf_type_2_densities'].apply(json.dumps)
    except KeyError:
        sys.stderr.write("type2 rdfs not found")
    try:
        df['rdf_all_bins'] = df['rdf_all_bins'].apply(json.dumps)
        df['rdf_all_densities'] = df['rdf_all_densities'].apply(json.dumps)
    except KeyError:
        sys.stderr.write("all rdfs not found")
        
    df["comment"] = None
    df.loc[0,"comment"] = (f"first: {args.first[0]}, "
                                   + f"last: {args.last[0]}, "
                                   + f"step: {args.step}, "
                                   + f"n_points: {args.sasa_points}, "
                                   + f"r_probe: {r_probe}, "
                                   + f"radii: {radii}")
    df.to_excel(args.output[0], index = False)
    pickle_filename = f"{args.output[0]}.pkl"
    df.to_pickle(pickle_filename)