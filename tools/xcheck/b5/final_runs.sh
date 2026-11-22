#!/bin/sh
# Rerun the chain with all computed B5 values (b5_values.json): own mode m = 134..137, compat 134..136,
# then the full own-mode sweep m = 92..136 (adding a rule can only lower the chain, so m <= 133 must stay closed).
cd "$(dirname "$0")"
python3 run_chain.py own 134 135 136 137 > logs/chain_own_134_137.log 2>&1
python3 run_chain.py compat 134 135 136 > logs/chain_compat_134_136.log 2>&1
python3 run_chain.py own $(seq 92 136) > logs/sweep_own_92_136.log 2>&1
