#!/bin/bash

inputs=("25.04.1_eAu_merged.root") 
output_dir="/home/macilla_vanilla/eic/EICreconOutputReader"

mkdir -p "${output_dir}" 

for input in "${inputs[@]}"; do
    output="${output_dir}/$(basename ${input%.*})"

    echo "Generating plot for file: ${input}"
    echo "Output at [ ${output}.root ]"


    echo "Generating plots for file: ${output}.root"
    root -b -q macros/generate_all_plots.C\(\"${output}.root\"\)
    root -b -q macros/plot_diffractive_event_kinematics.C\(\"${output}.root\"\)
    root -b -q macros/plot_diffractive_vm_resolution.C\(\"${output}.root\"\)

    echo "Finished processing file: ${input}"
    echo
done

echo "All files processed!"
