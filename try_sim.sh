#!/bin/bash
start=`date +%s`
for i in {1..10}
do
    echo "start sim number $i"
    projectPath="/home/berdnikovich/try_mHiC/result_simulation_simmple/40_sim_const/sim_6_5_4"
    simnum="sim$i"
    resdir="${projectPath}/${simnum}"
    #creation of file-system
    mkdir ${resdir}
    mkdir ${resdir}/s4
    mkdir ${resdir}/s6
    mkdir ${resdir}/up_prior
    mkdir ${resdir}/down_prior
    mkdir ${resdir}/visualisation

    ##starts simulation of data
    echo "Simulation" 
    python3 Stage_1.py "${resdir}/s4" "${resdir}/visualisation" #?''
    python3 Stage2_generate_multi.py "${resdir}/s4" "${resdir}/visualisation" #?''
    python3 Stage3_bindata.py "${resdir}/s4" "${resdir}/visualisation"

    ##split data 
    echo "Data-splitting"
    python3 separator_uni.py "${resdir}/s4"
    
    ##build prior
    echo "Building prior"
    BL=500000
    validI1="${resdir}/s4/Uni_up_${BL}.Red_c"
    validI2="${resdir}/s4/Uni_down_${BL}.Red_c"
    updir="${resdir}/up_prior"
    downdir="${resdir}/down_prior"
    splineBin=50
    resolution=500000

    python3 s5_prior_RedC.py -f "${resdir}/s4/fragments_$BL" -i $validI1 -o ${updir}  -b $splineBin -r $resolution 
    python3 s5_prior_RedC.py -f "${resdir}/s4/fragments_$BL" -i $validI2 -o ${downdir}  -b $splineBin -r $resolution 

    ##conversion of two prior distribution to one input file for s6
    echo "Conversion of two prior distribution"
    python3 convert_prior_distrib.py "${resdir}/up_prior/s5_prior.mhic" "${resdir}/down_prior/s5_prior.mhic" "${resdir}/s5_prior12" 

    ##start step6
    echo "Step 6"
    multi="${resdir}/s4/MULTI_contacts_ALL_${BL}.Red_c"
    multiKeys="${resdir}/s4/MULTI_contacts_ALL_${BL}.Red_c.multiKeys"
    prior="${resdir}/s5_prior12"
    uni="${resdir}/s4/UNI_contacts_ALL_${BL}.Red_c"
    threshold=0.5
    filename="outFile.multi"

    awk -v OFS="_" '{print $2, $3, $4, $5}' $multi | sort -u >$multiKeys
    python3 s6_em.py -p $prior -u $uni -m $multi -mk $multiKeys -t $threshold -o "${resdir}/s6" -f $filename
    python3 parse_output.py $resdir

done
end=`date +%s`
runtime=$((end-start))
echo "Runtime is $runtime"