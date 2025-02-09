nextflow run Lipinski-B/RNAseq




WORK=/home/bobo/Bureau/Git ; cd $WORK
DATA=$WORK/data
RSLT=$WORK/result
SCRIPT=$WORK/script
nextflow run $SCRIPT/RNAseq/RNAseq.nf --input $DATA --output $RSLT --protocol "PE" -profile docker -resume




WORK=/home/bobo/Bureau/Git ; 
DATA=$WORK/data
RSLT=$WORK/result_git
SCRIPT=$WORK/script
nextflow run Lipinski-B/RNAseq --input $DATA --output $RSLT --protocol "PE" -profile docker -resume





docker build -t mon-app:1.0 .



VERSION=1.2
docker build -t rnaseq:$VERSION .
docker tag rnaseq:$VERSION lipinskiboris/rnaseq:$VERSION
docker login
docker push lipinskiboris/rnaseq:$VERSION




docker tag lipinskiboris/rnaseq:$VERSION rnaseq:$VERSION 



sudo apt-get install libssl-dev libfontconfig1-dev libxml2-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev