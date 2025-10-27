#!/bin/bash

input=$1
output_dir="/home/macilla_vanilla/EIC_VM_Diffraction/diffractive_phi_analysis"
#"/eic/u/macink/EICreconOutputReader/output_sartre_eAu_diffractive_phi_10x100"

mkdir -p "${output_dir}"

output="${output_dir}/$(basename ${input})"
echo "Running analysis for file: ${input}"
echo "Output at [ ${output}_output.root ]"
    
# Run the diffractive VM analysis
root -b -q diffractive_vm_full_analysis.cxx\(\"${input}\",\"${output}\"\)
echo "Finished processing file: ${input}"
echo
    

    

