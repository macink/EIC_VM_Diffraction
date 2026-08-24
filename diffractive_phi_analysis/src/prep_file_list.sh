#!/bin/bash

#incoherent phi
#prod="BeAGLE1.03.02-1.1_phi_eAu_10x100_q2_1to10000_hiAcc_run"
#did="epic:/RECO/26.04.1/epic_craterlake/EXCLUSIVE/DIFFRACTIVE_PHI_ABCONV/BeAGLE1.03.02-1.1/eAu/10x100/q2_1to10000"

#coherent phi
prod="sartre1.39-1.0_coherent_phi_eAu_bsat_10x100_ab"
did="epic:/RECO/26.04.1/epic_craterlake/EXCLUSIVE/DIFFRACTIVE_PHI_ABCONV/sartre1.39-1.0/eAu/coherent/bsat/10x100"

#coherent rho
#prod="sartre1.39-1.1_coherent_rho_eAu_bsat_10x100_q2_1to20_hiAcc"
#did="epic:/RECO/26.04.1/epic_craterlake/EXCLUSIVE/DIFFRACTIVE_RHO_ABCONV/sartre1.39-1.1/eAu/coherent/bsat/10x100/q2_1to20"

folder="${prod}"

filelistall=file.all.list
> $filelistall

echo "[i] Querying Rucio..."
rucio replica list file $did --pfns \
    | grep "root://" \
    | grep "$prod" \
    | grep -v "sca2302.jlab.org" \
    | grep -v "/mss/" \
    >> $filelistall

echo "[i] Found $(wc -l < $filelistall) files"

split -l 20 --numeric-suffixes --suffix-length=3 \
    $filelistall --additional-suffix=.list subList_

mkdir -pv "$folder"
rm -rf "$folder"/*
mv $filelistall "$folder"
mv subList_* "$folder"

echo "[i] Done. Lists written to $folder/"
