#!/bin/bash

cd ${1}

/eic/u/macink/eic-shell  << EOF
./diffractive_vm_full_analysis ${2} ${3} 
exit
EOF