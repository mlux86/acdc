#!/bin/bash

fasta="$1"

if [[ -z "$fasta" ]]; then
 	echo "No fasta file given!"
 	exit 1
fi

krakenFile=$(mktemp)
krakenOut=$(mktemp)

kraken --db /home/mlux/minikraken_20141208/ --output "$krakenFile" "$fasta"

kraken-translate --db /home/mlux/minikraken_20141208/ "$krakenFile" > "$krakenOut"

rm "$krakenFile"

echo "$krakenOut"
