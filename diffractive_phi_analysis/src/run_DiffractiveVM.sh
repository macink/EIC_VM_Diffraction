#!/bin/bash

input=$1
output_dir="/eic/u/macink/EICreconOutputReader/output_DIS_10x100"

mkdir -p "${output_dir}"

output="${output_dir}/$(basename ${input})"
echo "Running analysis for file: ${input}"
echo "Output at [ ${output}_output.root ]"
    
# Run the diffractive VM analysis
root -b -q diffractive_vm_full_analysis.cxx\(\"${input}\",\"${output}\"\)
echo "Finished processing file: ${input}"
echo
    
