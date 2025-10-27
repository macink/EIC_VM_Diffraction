#!/bin/bash

inputs=("sartre_bnonsat_Au_phi_ab_eAu_q2_15_1.7944.eicrecon.edm4eic.root_output.root") 
output_dir="/eic/u/macink/EICreconOutputReader/output_sartre_eAu_diffractive_phi_10x100"

mkdir -p "${output_dir}" 

for input in "${inputs[@]}"; do
    output="${output_dir}/$(basename ${input%.*})"

    echo "Generating plot for file: ${input}"
    echo "Output at [ ${output}.root ]"


    echo "Generating plots for file: ${output}.root"

    root -b -q macros/generate_all_plots_Lnorm.C\(\"${output}.root\"\)
    root -b -q macros/generate_all_plots_noNorm.C\(\"${output}.root\"\)
    root -b -q macros/plot_QA.C\(\"${output}.root\"\)
    
    echo "Finished processing file: ${input}"
    echo
done

echo "All files processed!"