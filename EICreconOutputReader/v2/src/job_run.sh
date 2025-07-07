#!/bin/bash

cd ${1}

/eic/u/macink/eic-shell  << EOF
./run_DiffractiveVM.sh ${2} ${3} 
exit
EOF