#!/bin/bash
AGGREGATES="${HOME}/Documents/Codes/git/mouse2/mouse2/aggregates.py"
#N_AGGR="${HOME}/Documents/Codes/git/copoly/average_aggregation_number.py"
N_AGGR="${HOME}/Documents/Codes/git/copoly/number_of_aggregates.py"
for i in $(seq $1 $2)
do
	echo -n "${i} "; ${AGGREGATES} --selection "type 1" ${i}.data | python3 ${N_AGGR}
done 2>/dev/null
