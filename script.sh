tmux a -t RNAseq

WORK=/mnt/datagenetique/ANALYSIS/BIT/PROJECTS/BL/Projet/SKILLS2
DATA=$WORK/data
RSLT=$WORK/result
SCRIPT=$WORK/script
NF=$SCRIPT/nextflow
cd $WORK/script

$NF -c $WORK/script/RNAseq.config run $WORK/script/RNAseq.nf -profile singularity
