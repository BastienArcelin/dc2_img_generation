#!/bin/bash  

#$ -P P_lsst

#$ -l os=cl7 
#$ -l sps=1
#$ -l s_fsize=2G
#$ -l s_cpu=2:00:00
#$ -l s_rss=1G 

#$ -M arcelin@apc.in2p3.fr
#$ -m be   ## envoie un email quand le job commence et termine 

#$ -o /sps/lsst/users/barcelin/job_outputs/dc2/ 
#$ -e /sps/lsst/users/barcelin/job_outputs/dc2/

source /pbs/home/b/barcelin/pbs_throng_link/lsst_stack/loadLSST.bash
setup lsst_distrib

cd /pbs/home/b/barcelin/pbs_throng_link/dc2_img_generation/script/

##python generate_dc2_img.py 4637 validation 10000 ## Vadidation
##python generate_dc2_img.py 4855 test 10000 ## Test
python generate_dc2_img.py 4232 /training_24.5_v2 100 24.5 ## Training
