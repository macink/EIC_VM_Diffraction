#!/bin/bash

cd ${1}

/eic/u/macink/eic-shell --version 25.12.0-stable << EOF
./diffractive_vm_full_analysis ${2} ${3} 
exit
EOF
