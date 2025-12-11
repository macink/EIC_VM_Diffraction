#!/bin/bash

submit=$1

###############################################################
### Batch production for real data using file list          ###
###############################################################

if [ $submit -eq 2 ]; then
    #configs=(sartre1.39-1.0_coherent_phi_eAu_bsat_10x100_ab)
    #configs=(BeAGLE1.03.02-1.1_phi_eAu_10x100_q2_1to10000_hiAcc_run)
    #configs=(sartre1.39-1.1_coherent_rho_eAu_bsat_10x100_q2_1to20_hiAcc)
    configs=(BeAGLE1.03.02-1.0_DIS_eAu_10x100_q2_1to10_ab_run)
    pwd=$PWD
    for config in "${configs[@]}"; do
        echo "[i] Running config = $config"

        odir=/eic/u/macink/EICreconOutputReader/analysis/$config
        logdir=$odir/log
        if [ ! -d $odir ]; then
            mkdir -pv $odir
        fi
        rm -rf $odir/*
        mkdir $logdir
        echo "Checking files: " $(ls /eic/u/macink/EICreconOutputReader/$config/subList*)

        executable=job_run.sh
        cp -v ${executable} $odir/.
        cp -v analysis $odir/.
        cp -v pleaseIncludeMe.h $odir/.
        cp -v diffractive_vm_full_analysis.cxx $odir/. 
        cp -v run_DiffractiveVM.sh $odir/.

        echo $odir/$executable
        # Initializing Condor File
        condor_file=CondorFile_$config
        echo "" > ${condor_file}
        echo "Universe    = vanilla" >> ${condor_file}
        echo "Executable  = ${odir}/${executable}" >> ${condor_file}
        echo "GetEnv  =  True" >> ${condor_file}

        files=`ls /eic/u/macink/EICreconOutputReader/$config/subList*`
        for file in $files; do
            listNum=`basename ${file} | sed "s/.list//g" | cut -f 2 -d _`
            LogFile=${logdir}/log_${listNum}.out
            ErrFile=${logdir}/log_${listNum}.err
            OutFile=${odir}/output_${listNum}.root

            echo "" >> ${condor_file}
            echo "Output    = ${LogFile}" >> ${condor_file}
            echo "Error     = ${ErrFile}" >> ${condor_file}
            echo "Arguments = ${odir} ${file} ${OutFile}" >> ${condor_file}
            echo "Queue" >> ${condor_file}
        done
        mv ${condor_file} $odir/.
        cd $odir
        condor_submit ${condor_file}
        cd $pwd
    done
fi
