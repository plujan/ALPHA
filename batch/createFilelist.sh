#!/bin/bash

if [[ "$1" != "MC" ]] && [[ "$1" != "DATA" ]]
then
    echo specify DATA or MC
    return
fi

############ from here

mccampaign="PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6"

# step 0: execute script with PhedEx request number
if [[ "$1" == "MC" ]]
then
    dataset=/*/*/MINIAODSIM 
    python2.6 batch/get_ds_file_info.py -d "/*/*RunIISummer16MiniAODv2*$mccampaign*/MINIAODSIM" > query.txt #append Summer16
else
    if [[ "$2" == "B" ]] || [[ "$2" == "C" ]] || [[ "$2" == "D" ]] || [[ "$2" == "E" ]] || [[ "$2" == "F" ]] || [[ "$2" == "G" ]] || [[ "$2" == "H" ]]
    then
        python2.6 batch/get_ds_file_info.py -d "/*/Run2016$2*/MINIAOD" > query.txt   # for Data
    else
        echo "for DATA also specify the era, e.g.:"
        echo "sh batch/createFilelist.sh DATA D"
        exit
    fi
fi

# step 1: get name of the temporary file
tmpname=$(cat query.txt | awk '{print $4}' | sed -e 's/to//g')
echo File to be read: $tmpname

# # step 2: get list of the samples names (with postfix)
cat query.txt | grep MINIAOD | awk '{print $1}' | sed -e "s/\/RunIISummer16MiniAODv2-$mccampaign//g" | sed -e 's/\/MINIAODSIM//g' | sed -e 's/\/MINIAOD//g' | sed -e 's/\///g' > samplelist.txt

############ to here

# step 3: get filelist form the dump, filter it, and dump it into appropriate files
cat samplelist.txt | while read sample
do
    if [[ "$1" == "MC" ]]
    then
        checkthis=${mccampaign}
        
        if [[ $sample == *"_ext"* ]]
        then
            trimname=$(echo $sample | sed -e 's/_ext[0-9]-v[0-9]//g')
            versionname=""${sample: -8}""
        else
            trimname=$(echo $sample | sed -e 's/-v[0-9]//g')
            versionname=""${sample: -3}""
        fi

        dirname="Summer16" # for MC
        recoversion=$trimname
    else
        checkthis=""

        dirname="Run2016" # for Data
        trimname=${sample%Run2016*} # for Data
        recoversion=$(echo $sample | cut -d "-" -f2)
        versionname=${sample: -2} # for Data
    fi

    cat $tmpname | grep store | grep $trimname | grep $checkthis$versionname | grep $recoversion | sort > filelists/$dirname/$sample.txt
    sed -i -e 's/^/dcap:\/\/t2-srm-02.lnl.infn.it\/pnfs\/lnl.infn.it\/data\/cms\//' filelists/$dirname/$sample.txt
    echo Created filelist filelists/$dirname/$sample.txt

done

# final step: clean up
rm query.txt
rm samplelist.txt

# # removed unused lists
# if [[ "$1" == "MC" ]]
# then
    # rm filelists/Spring16/QCD*GenJets5*
    # #rm filelists/Spring16/QCD*BGenFilter*
    # #rm filelists/Spring16/DYB*
    # rm filelists/Spring16/*Contin*
    # rm filelists/Spring16/GJets*
    # rm filelists/Spring16/GGJets*
    # rm filelists/Spring16/WWTo4Q_13TeV-powheg_v0-v1.txt
    # rm filelists/Spring16/TT_TuneCUETP8M1_alphaS01273_13TeV-powheg-scaleup-pythia8_v0-v1.txt
    # rm filelists/Spring16/GluGluSpin0ToZG*
    # rm filelists/Spring16/WpWp*
    # rm filelists/Spring16/WprimeToMuNu*
# fi

# source batch/createFilelist.sh
