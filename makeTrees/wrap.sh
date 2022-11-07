#!/bin/bash

case "$#" in

3)  
jobscript=$1
jobarrayFile=$2
pathoutputlog=$3

singularity exec \
-B /cvmfs/hadessoft.gsi.de/install/debian8/install:/cvmfs/hades.gsi.de/install \
-B /cvmfs/hadessoft.gsi.de/param:/cvmfs/hades.gsi.de/param \
-B /cvmfs/hadessoft.gsi.de/install/debian8/oracle:/cvmfs/it.gsi.de/oracle \
-B /lustre \
/cvmfs/vae.gsi.de/debian8/containers/user_container-production.sif  ${jobscript} ${jobarrayFile} ${pathoutputlog}
    ;;
4) 
jobscript=$1
jobarrayFile=$2
pathoutputlog=$3
arrayoffset=$4

singularity exec \
-B /cvmfs/hadessoft.gsi.de/install/debian8/install:/cvmfs/hades.gsi.de/install \
-B /cvmfs/hadessoft.gsi.de/param:/cvmfs/hades.gsi.de/param \
-B /cvmfs/hadessoft.gsi.de/install/debian8/oracle:/cvmfs/it.gsi.de/oracle \
-B /lustre \
/cvmfs/vae.gsi.de/debian8/containers/user_container-production.sif  ${jobscript} ${jobarrayFile} ${pathoutputlog} ${arrayoffset}
    ;;
*) echo "Unsupported number of arguments" 
   echo "Usage : wrap.sh jobScript.sh jobarrayfile pathoutputlog ?par4?"
   exit
   ;;
esac
