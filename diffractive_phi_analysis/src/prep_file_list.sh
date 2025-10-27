#!/bin/bash

datadir=/volatile/eic/EPIC/RECO/25.06.1/epic_craterlake/EXCLUSIVE/DIFFRACTIVE_PHI_ABCONV/sartre1.39-1.0/eAu/coherent/bsat/10x100
prod=sartre1.39-1.0_coherent_phi_eAu_bsat_10x100_ab

filelistall=file.all.list
cp blank $filelistall

files=`xrdfs root://dtn-eic.jlab.org ls $datadir`

for file in $files; do
    if [[ $file == *"sartre1.39-1.0_coherent_phi_eAu_bsat_10x100_ab."*".eicrecon.edm4eic.root" ]]; then
        echo root://dtn-eic.jlab.org/$file >> $filelistall
    fi
done


split -l 20 --numeric-suffixes --suffix-length=3 $filelistall --additional-suffix=.list subList_

if [ ! -d $prod ]; then
    mkdir -pv $prod
fi
rm -rf $prod/*
mv $filelistall $prod
mv subList* $prod


