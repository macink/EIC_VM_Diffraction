#!/bin/bash

datadir=/volatile/eic/EPIC/RECO/25.06.1/epic_craterlake/EXCLUSIVE/DIFFRACTIVE_PHI_ABCONV/Sartre/Coherent

#prod=sartre_bnonsat_Au_phi_ab_eAu_q2_15
prod=sartre_bnonsat_Au_phi_ab_eAu
#prod=BeAGLE1.03.02-1.0_phi_eCa_18x137.5_q2_1to10_ab
#prod=sartre1.39-1.0_coherent_phi_eCa_bsat_18x137.5_q2_1to10_ab

filelistall=file.all.list
cp blank $filelistall
#files=`xrdfs root://dtn-eic.jlab.org  ls $datadir`
#for file in $files; do
#    echo root://dtn-eic.jlab.org/$file >> $filelistall
#done

files=`xrdfs root://dtn-eic.jlab.org ls $datadir`
#echo "$files"
for file in $files; do
    #if [[ $file == *"sartre_bnonsat_Au_phi_ab_eAu_q2_15_1."*".eicrecon.edm4eic.root" ]]; then
    if [[ $file == *"sartre_bnonsat_Au_phi_ab_eAu_1."*".eicrecon.edm4eic.root" ]]; then
    #if [[ $file == *"BeAGLE1.03.02-1.0_phi_eCa_18x137.5_q2_1to10_ab."*".eicrecon.edm4eic.root" ]]; then
    #if [[ $file == *"sartre1.39-1.0_coherent_phi_eCa_bsat_18x137.5_q2_1to10_ab."*".eicrecon.edm4eic.root" ]]; then
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

