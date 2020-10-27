#!/bin/bash  

#$ -P P_lsst

#$ -l os=cl7 
#$ -l sps=1
#$ -l s_fsize=4G
#$ -l s_cpu=20:00:00
#$ -l s_rss=10G 

#$ -M arcelin@apc.in2p3.fr
#$ -m be   ## envoie un email quand le job commence et termine 

#$ -o /sps/lsst/users/barcelin/job_outputs/dc2/ 
#$ -e /sps/lsst/users/barcelin/job_outputs/dc2/

source /pbs/home/b/barcelin/pbs_throng_link/lsst_stack/loadLSST.bash
setup lsst_distrib

cd /pbs/home/b/barcelin/pbs_throng_link/dc2_img_generation/scripts/

##python generate_dc2_img.py 4637 validation 10000 ## Vadidation
##python generate_dc2_img.py 4855 test 10000 ## Test

python generate_dc2_img.py 5074 training 10000 ## Training
