### Here is the code I used to try to gapfill the assembly after Chromonomer with our PacBio HiFi data

## I had to start off by setting up my directories and locally installing some different packages
BASE=/project/conservationgenomicsfraik/CCT
SOURCE=/project/rainbow_trout
DATA=${SOURCE}/CoastalCT_hifiasm/with_nanopore_l100k
INDIR=${BASE}/Input_Files
HIFI=${DATA}/Geoff/cctrout
TMPDIR=/90daydata/conservationgenomicssparks/CCT_Assembly
##SRC=/home/guangtu.gao/rainbow_trout/Geoff/cctrout
#PROJ=m64310e_220430_052108.hifi_reads
#bam2fastq -o $PROJ -c 9 "$SRC"/"$PROJ".bam

## Load required modeules
module load pbsuite/15.8.24
module load busco5/5.7.1
module load seqkit/2.4.0
# pyenv install 2.7
pyenv global 2.7
eval "$(pyenv init -)"
# https://github.com/pyenv/pyenv
python -V
# pip install networkx

