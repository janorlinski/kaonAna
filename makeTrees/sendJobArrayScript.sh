#!/bin/bash

# Submission script for GridEngine (GE). Each job will 
# be executed via the jobScript.sh
# This jobScript supports up to 7 parameters. Edit 
# the user specific part of the script according to 
# your program.
#
# Input to the script is a filelist with 1 file per line.
# For each file a job is started. With the parameter 
# nFilesPerJob a comma separated filelist will be 
# generated and handed to the job script. This feature
# is usefull when running many small jobs. Each
# job has its own logfile. All needed directories for the
# logfiles will be created if non existing.
#
# IMPORTANT: the hera/prometheus cluster jobs will only
# see the /hera file system. All needed scripts, programs
# and parameters have to be located on /hera or the job
# will crash. This script syncs your working dir to the submission
# dir on /hera . Make sure your scripts use the submission dir!
# Software should be taken from /cvmfs/hades.gsi.de/install/
#
# job log files will be named like inputfiles. If nFilesPerJob > 1
# the log files will contain the partnumber.
#
######################################################################
#   CONFIGURATION

user=$(whoami)
currentdir=$(pwd | xargs -i basename {})                           
currentDir=$(pwd)

day=${1}
dataset=ag158ag_3200A_${2}
particle=${3}

submmissionbase=/lustre/hades/user/${user}/sub/flow_exp

submissiondir=${submmissionbase}/${currentdir}
    outputdir=/lustre/hades/user/${user}/flow/output/makeTrees/${dataset}/${day}   # outputdir for files AND logFiles
pathoutputlog=${outputdir}/logs     # protocol from batch farm for each file
#    outputdir=/lustre/hades/user/${user}/flow/output/${generation}/${dataset}   # outputdir for files AND logFiles
#pathoutputlog=/lustre/hades/user/${user}/flow/output/${generation}/${dataset}/logs     # protocol from batch farm for each file

#nFilesPerJob=1                                                 # number of files to be analyzed by 1 job (default==1)
#nFilesPerJob=10                                                 # number of files to be analyzed by 1 job (default==1)
#nFilesPerJob=15                                                 # number of files to be analyzed by 1 job (default==1)
#nFilesPerJob=20                                                 # number of files to be analyzed by 1 job (default==1)
nFilesPerJob=40                                                 # number of files to be analyzed by 1 job (default==1)
#nFilesPerJob=50                                                 # number of files to be analyzed by 1 job (default==1)
#nFilesPerJob=100                                                 # number of files to be analyzed by 1 job (default==1)
#nFilesPerJob=200                                                 # number of files to be analyzed by 1 job (default==1)

jobscript=${submissiondir}/jobArrayScript.sh     # exec script (full path, call without dot, set it executable!)
#wrap=${submissiondir}/wrap.sh     # exec script (full path, call without dot, set it executable!)


#fillNtuple -> histograms Theta_Dphi
#filename=${particle}_May21_v${4}${5}${6}                          # filename of log file if nFilesPerJob > 1 (partnumber will be appended)
#par1=/lustre/hades/user/lchlad/flow/scripts/apr12/behruz_lib/defall.sh            # optional par1 : environment script

#proton
#filename=${particle}_May21_Pid${4}_dummy${5}_OCpid${6}                          # filename of log file if nFilesPerJob > 1 (partnumber will be appended)
#par1=/lustre/hades/user/lchlad/flow/scripts/apr12/behruz_lib/defall.sh            # optional par1 : environment script

#kCharged     
filename=${particle}_v300_PidOpt${4}_TrackCuts${5}_OCpid${6}                          # filename of log file if nFilesPerJob > 1 (partnumber will be appended)
par1=/lustre/hades/user/lchlad/flow/scripts/apr12/behruz_lib/defall.sh            # optional par1 : environment script
#par1=/lustre/hades/user/lchlad/flow/scripts/apr12/behruz_lib/defall_debian10.sh            # optional par1 : environment script

#kZero
#filename=${particle}_v204_PidOpt${4}_Precuts${5}_MomCorr${6}                         # filename of log file if nFilesPerJob > 1 (partnumber will be appended)
#par1=/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/behruz_root6_lib/defall.sh     # optional par1 : environment script

par2=${submissiondir}/loopDST/build/analysis                               # optional par2 : executable
par3=""                                                         # optional par3 : input file list
par4=${outputdir}                                               # optional par4 : outputfile (part number will be appended (_num.root))
par5=${4}                                                 # optional par5 : number of events
par6=${5}                                                       # optional par6
par7=${6}                                                       # optional par7

#resources="--partition=long --mem=4000 --time=0-12:00:00"                    
resources="--mem=2000 --time=0-4:00:00"                          
#resources="--partition=debug --mem=4000 --time=0-0:10:00"                          

jobarrayFile="jobarray_${day}_${dataset}_${4}_${5}_${6}.dat"
                                                              
filelist=./filelists/day_${day}_${dataset}.list   # file list in local dir! not in submissiondir!!!
######################################################################

nFiles=$( cat $filelist | wc -l)

#---------------------------------------------------------------------
# create sub dirs
if [ ! -d $submmissionbase ]
then
    echo "===> CREATE SUBMISSIONBASEDIR : $submmissionbase"
    mkdir -p $submmissionbase
else
    echo "===> USE SUBMISSIONBASEDIR : $submmissionbase"
fi

#---------------------------------------------------------------------
# output dirs

if [ ! -d $outputdir ]
then
   echo "===> CREATE OUTPUTDIR : $outputdir"
   mkdir -p $outputdir
else
   echo "===> USE OUTPUTDIR : $outputdir"
fi

if [ ! -d $pathoutputlog ]
then
   echo "===> CREATE LOGDIR : $pathoutputlog"
   mkdir -p $pathoutputlog
else
   echo "===> USE LOGDIR : $pathoutputlog"
fi
#---------------------------------------------------------------------


ctF=0          # counter for file number
ctJ=0          # counter for job number
partNumber=0   # counter for part number

#---------------------------------------------------------------------
# read the files list into an job array
if [ -f $jobarrayFile ]
then
  rm -f $jobarrayFile
fi

echo "===> CREATING JOB ARRAY FILE"

declare -a jobarray
ct1=0
for file in $(cat $filelist)
do
   jobarray[$ct1]=$file
   ((ct1+=1))
done
#---------------------------------------------------------------------


#---------------------------------------------------------------------
# loop over the job array and submit parts with
# nFilesPerJob to SLURM

while ((ctF<$nFiles))
do

     #---------------------------------------------------------------------
     # build comma separated file list
     # per job
     if [ $nFilesPerJob -gt 1 ]
     then
        infileList=${jobarray[${ctF}]}
        ((ctF+=1))
        for (( ctList=1;ctList<$nFilesPerJob; ctList++ ))
        do   	
            if [ $ctF -lt ${nFiles} ]
            then
               infileList="${infileList},${jobarray[${ctF}]}"
               ((ctF+=1))
            fi
        done
     else 
        infileList=${jobarray[${ctF}]}
        ((ctF+=1))
     fi

     ######################################################################
     #  SEND NEW JOB (USER SPECIFIC)
     
     par3=${infileList}

           #defall.sh prog filelist outdir  nev
     echo "${par1} ${par2} ${par3} ${par4} ${par5} ${par6} ${par7}" >>  $jobarrayFile
     

     ######################################################################
     
done
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# sync the local modified stuff 
# to the submission dir
echo "===> SYNC CURENTDIR TO SUBMISSIONDIR : rsync  -vHa $currentDir ${submmissionbase}"
rsync  -vHa $currentDir ${submmissionbase}/

syncStat=$?

if [ ! $syncStat -eq 0 ]
then
     echo "===> ERROR : SYNCHRONIZATION ENCOUNTERED PROBLEMS"
else

  echo "-------------------------------------------------"


  nFiles=$( cat $jobarrayFile | wc -l)
  ctsend=0
  block=500
  while ((${ctsend} * ${block} < ${nFiles}))
  do
     ((start=${ctsend}*${block}+1))
     ((stop= ${start}+${block}-1))
     ((rest=${nFiles}-${start}))
     if [ $rest -le $block ]
     then
        ((stop=$start+$rest))
     fi

     command="--array=${start}-${stop} ${resources} -D ${submissiondir}  --output=${pathoutputlog}/slurm-%A_%a.out -- ${jobscript} ${submissiondir}/${jobarrayFile} ${pathoutputlog} ${filename}"
     echo $command
     sbatch $command

     ((ctsend+=1))
  done

  echo "${nFiles} jobs are submitted"
fi
