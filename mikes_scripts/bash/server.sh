export LSF_DOCKER_VOLUMES='/storage1/fs1/kutluay/Active/:/storage1/fs1/kutluay/Active /scratch1/fs1/kutluay/:/scratch1/fs1/kutluay' # to access storage and scratch in docker

LSF_DOCKER_NETWORK=host bsub -G compute-kutluay -q general-interactive -Is -a 'docker(ubuntu)' /bin/bash

export LSF_DOCKER_PRESERVE_ENVIRONMENT=false # to run own docker (if set, ubuntu jobs cannot do bsub anymore)

bsub -G compute-kutluay -q general-interactive -Is -a 'docker(yiqingw/kutluay_docker:latest)' Rscript testR.R


##################
# command to run interactive job with own docker #
# "PATH=..." is to run own docker, better than "LSF_DOCKER_PRESERVE_ENVIRONMENT=false"
PATH="/opt/conda/bin:$PATH" bsub -G compute-kutluay -q general-interactive -Is -a 'docker(yiqingw/kutluay_docker_test)' /bin/bash
##################


export LSF_DOCKER_NETWORK=host # to run docker in job



# to run non-interactive job with own docker example #
bsub -n 4 -G compute-kutluay -M 50GB -q general -J 'rseq_alignment_test_cellular' -u yiqingw@wustl.edu -R 'rusage[mem=50000] span[hosts=1]' -a 'docker(yiqingw/kutluay_docker)' bash /storage1/fs1/kutluay/Active/rnaseq_alignment/align_hgenome_star_ris.sh $index $i $baseoutput'/'${array[6]} $gtf



## Build and push docker ##
docker build -t your_dockerhub_username/image_name:tag .
docker push your_dockerhub_username/image_name:optional_tag

