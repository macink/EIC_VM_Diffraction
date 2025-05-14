#!/bin/bash

cd ${1}

/eic/u/macink/eic-shell  << EOF
./diffractive_vm_simple_analysis_mod ${2} ${3} 
exit
EOF