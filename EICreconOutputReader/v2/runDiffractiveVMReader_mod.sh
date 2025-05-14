#!/bin/bash

input=$1
output_dir="/eic/u/macink/EICreconOutputReader/25.04.1_eAu_output"

mkdir -p "${output_dir}"

output="${output_dir}/$(basename ${input})"
echo "Running analysis for file: ${input}"
echo "Input edm4eic data is at [ ${input} ]"
echo "Output at [ ${output}_output.root ]"
    
root -b -q diffractive_vm_simple_analysis_mod.cxx\(\"${input}\",\"${output}\"\)

echo "Finished processing file: ${input}"
echo


