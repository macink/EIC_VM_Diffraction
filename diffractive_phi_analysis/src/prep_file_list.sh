#!/bin/bash

#datadir=/volatile/eic/EPIC/RECO/25.10.2/epic_craterlake/EXCLUSIVE/DIFFRACTIVE_PHI_ABCONV/sartre1.39-1.0/eAu/coherent/bsat/10x100
#datadir=/volatile/eic/EPIC/RECO/25.10.3/epic_craterlake/EXCLUSIVE/DIFFRACTIVE_PHI_ABCONV/BeAGLE1.03.02-1.1/eAu/10x100/q2_1to10000
#datadir=/volatile/eic/EPIC/RECO/25.10.3/epic_craterlake/EXCLUSIVE/DIFFRACTIVE_RHO_ABCONV/sartre1.39-1.1/eAu/coherent/bsat/10x100/q2_1to20
datadir=/volatile/eic/EPIC/RECO/25.10.2/epic_craterlake_without_zdc/DIS/BeAGLE1.03.02-1.0/eAu/10x100/q2_1to10

#prod=sartre1.39-1.0_coherent_phi_eAu_bsat_10x100_ab
#prod=BeAGLE1.03.02-1.1_phi_eAu_10x100_q2_1to10000_hiAcc_run
#prod=sartre1.39-1.1_coherent_rho_eAu_bsat_10x100_q2_1to20_hiAcc
prod=BeAGLE1.03.02-1.0_DIS_eAu_10x100_q2_1to10_ab_run


filelistall=file.all.list
cp blank $filelistall

#files=$(find /gpfs02/eic/macink/BeAGLE_100x10/input_files -type f -name "*BeAGLE_eAu_phi_hiacc_10x100_ab_*.root" ! -name "*.hist.root")
#> $filelistall
#for file in $files; do
#    echo $file >> $filelistall
#done

files=`xrdfs root://dtn-eic.jlab.org ls $datadir`

for file in $files; do
    if [[ $file == *"BeAGLE1.03.02-1.0_DIS_eAu_10x100_q2_1to10_ab_run"*".eicrecon.edm4eic.root" ]]; then
    #if [[ $file == *"sartre1.39-1.0_coherent_phi_eAu_bsat_10x100_ab."*".eicrecon.edm4eic.root" ]]; then
    #if [[ $file == *"BeAGLE1.03.02-1.1_phi_eAu_10x100_q2_1to10000_hiAcc_run"*".eicrecon.edm4eic.root" ]]; then
    #if [[ $file == *"sartre1.39-1.1_coherent_rho_eAu_bsat_10x100_q2_1to20_hiAcc."*".eicrecon.edm4eic.root" ]]; then
        echo root://dtn-eic.jlab.org/$file >> $filelistall
    fi
done


split -l 20 --numeric-suffixes --suffix-length=3 $filelistall --additional-suffix=.list subList_
#split -l 1 --numeric-suffixes --suffix-length=3 $filelistall --additional-suffix=.list subList_

if [ ! -d $prod ]; then
    mkdir -pv $prod
fi
rm -rf $prod/*
mv $filelistall $prod
mv subList* $prod

