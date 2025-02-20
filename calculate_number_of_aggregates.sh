#!/bin/bash
AGGREGATES="${HOME}/Documents/Codes/git/copoly/aggregates-tails.py"
N_AGGR="${HOME}/Documents/Codes/git/copoly/number_of_aggregates.py"
for i in $(seq $1 $2)
do
	echo -n "${i} "; ${AGGREGATES} ${i}.data | python3 ${N_AGGR}
done #2>/dev/null
