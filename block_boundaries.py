#!/usr/bin/python3
import argparse
import pandas as pd
import sys

parser = argparse.ArgumentParser(
    description = 'Convert pandas in pickle chain sequences to strings')

parser.add_argument('filename', metavar = 'FILE', type = str,
    help = 'pickle filename')

parser.add_argument('--nmol', metavar = 'n', type = int, nargs = '?',
                    const = 1, help = "number of molecules to read")
                    
parser.add_argument('--nmol-offset', metavar = 'offset', type = int,
                    nargs = '?', const = 0,
                    help = "skip this number of sequences")

args = parser.parse_args()


data = pd.read_pickle(args.filename)

for imol in range(args.nmol_offset,args.nmol + args.nmol_offset):
    print(f'Seq nr. {imol}:')
    seq = data[imol]
    seq_nonzero = [atype for atype in seq if atype != 0]
    indexed_seq_nonzero = list(zip(range(len(seq_nonzero)), seq_nonzero))
    atypes = list(set(seq_nonzero))
    atypes.sort()
    print(atypes)
    for atype in atypes:
        maximum = [i[0] for i in indexed_seq_nonzero if i[1] == atype][-1]
        print(f'Type {atype} max {maximum}')
