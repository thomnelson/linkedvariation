#!/bin/bash

maindir=$(pwd)
tmpdir=${maindir}"/tmp/"
outdir=${maindir}"/output/"
eidos=${maindir}"/selection_in_structured_pops_bigcore.slim"
s=0.2
r=5e-5
mu=5e-7
N1=1000
N2=$N1
mig_red_factor=1
kpops=1
burnin=5000
runGens=10000

slim -define s=$s -define r=$r -define mu=$mu \
     -define N1=$N1 -define N2=$N2 \
     -define mig_red_factor=$mig_red_factor \
     -define kpops=$kpops \
     -define burnin=$burnin -define runGens=$runGens \
     -define simID=1 \
     $eidos
