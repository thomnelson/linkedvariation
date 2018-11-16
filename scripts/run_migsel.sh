#!/bin/bash

scriptdir="/Users/thom/Dropbox/1.Oregon/Data/active/simulations/forward/SLiM/scriptsEidos/"
outdir="/Users/thom/Dropbox/1.Oregon/Data/active/simulations/forward/SLiM/simulations/migrationSelectionBalance/outputs/"
eidos="MigrationSelectionBalance_cmdLine.eidos"

nreps=250
chrlen=100000
binsize=1000
stepsize=100

reps=$(seq $nreps)
for rep in $reps
do echo "SIMULATION REPLICATE "${rep}
  slim -define mu=5e-7 -define N2=500 -define r=5e-5 -define s=0.5 -define repl=$rep ${scriptdir}${eidos} 
  ./ms2tab.py ${outdir}"replicate"${rep}"_gen10000.ms" > ${outdir}"replicate"${rep}"_gen10000.tsv"
Rscript tab2stats.R ${outdir}"replicate"${rep}"_gen10000.tsv" $chrlen $binsize $stepsize $rep
done

# SEE SLIM MANUAL P 274!