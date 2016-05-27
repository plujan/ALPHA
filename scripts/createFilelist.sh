#!/bin/bash

# step 0: execute script with PhedEx request number
#python2.6 filelists/get_ds_file_info.py -n 669099 > query.txt
#python2.6 scripts/get_ds_file_info.py -d "/*/*RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2*/MINIAOD*" > query.txt
python2.6 scripts/get_ds_file_info.py -d "/*/Run2016B*/MINIAOD" > query.txt

# step 1: get name of the temporary file
tmpname=$(cat query.txt | awk '{print $4}' | sed -e 's/to//g')
echo File to be read: $tmpname

## step 2: get list of the samples names (with postfix)
cat query.txt | grep MINIAOD | awk '{print $1}' | sed -e 's/\/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2//g' | sed -e 's/\/MINIAODSIM//g' | sed -e 's/\/MINIAOD//g' | sed -e 's/\///g' > samplelist.txt

# step 3: get filelist form the dump, filter it, and dump it into appropriate files
cat samplelist.txt | while read sample
do
    trimname=${sample%_v*}
    cat $tmpname | grep store | grep $trimname > filelists/$sample.txt
    sed -i -e 's/^/dcap:\/\/t2-srm-02.lnl.infn.it\/pnfs\/lnl.infn.it\/data\/cms\//' filelists/$sample.txt
    echo Created filelist $sample.txt
done

# final step: clean up
rm query.txt
rm samplelist.txt

# source scripts/createFilelist.sh
