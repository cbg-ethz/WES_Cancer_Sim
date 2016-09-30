#!/bin/bash

# always set
pathToSnakeFile=`pwd`
snakemake=/usr/local/beerenwinkel/pipeline/bin/snakemake

# fill in
snakeFile=$pathToSnakeFile/rank_combination.snake
logFile=$pathToSnakeFile/log.log
queue=mpi04-ht.q


dryRun="-n" # yes
dryRun="" # no

command="${snakemake} -s ${snakeFile} --latency-wait 120 --cluster \"qsub -q ${queue} -N {params.jobname} -o {params.lsfoutfile} -e {params.lsferrfile}  \" -j 50 -p ${dryRun} 2>&1 | tee -a $logFile"
echo $command
eval $command

##/PATH/TO/SNAKEMAKE/bin/snakemake -s /PATH/TO/NGS-PIPE/snake/example.snake --configfile /PATH/TO/NGS-PIPE/snake/config.json -p --latency-wait 120 --cluster "<qsub/bsub>" -j 50 2>&1 | tee -a /PATH/TO/PROJECT/log.log

