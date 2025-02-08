nextflow run Lipinski-B/RNAseq




WORK=/home/bobo/Bureau/Git ; cd $WORK
DATA=$WORK/data
RSLT=$WORK/result
SCRIPT=$WORK/script
nextflow run $SCRIPT/RNAseq/RNAseq.nf --input $DATA --output $RSLT -profile docker -resume


docker build -t mon-app:1.0 .


docker login

docker tag rnaseq:1.1 lipinskiboris/rnaseq:1.1
docker push lipinskiboris/rnaseq:1.1

docker tag lipinskiboris/rnaseq:1.1 rnaseq:1.1 



sudo apt-get install libssl-dev libfontconfig1-dev libxml2-dev

sudo apt-get install libharfbuzz-dev libfribidi-dev 

sudo apt-get install libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev