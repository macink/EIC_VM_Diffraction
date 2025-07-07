#!/bin/bash

input=$1
output_dir="/eic/u/macink/EICreconOutputReader/output_25.06.1"

mkdir -p "${output_dir}"

output="${output_dir}/$(basename ${input})"
echo "Running analysis for file: ${input}"
echo "Output at [ ${output}_output.root ]"
    
# Run the diffractive VM analysis
root -b -q diffractive_vm_analysis.cxx\(\"${input}\",\"${output}\"\)
# Use this if you want to process just one file (supports wildcard)
#root -b -q diffractive_vm_analysis_1file.cxx\(\"${input}\",\"${output}\"\)
echo "Finished processing file: ${input}"
echo
    

