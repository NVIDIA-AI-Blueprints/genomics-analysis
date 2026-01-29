#!/bin/bash 

echo "Downloading pangenome data..." 

# Make data directory (if it doesn't exist)
DATA_DIR=$(pwd)/data
mkdir -p ${DATA_DIR}

# Donwload index files
cd ${DATA_DIR} && \
    echo -e "dist\nmin\ngbz" | xargs -I {} \
    wget "https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.{}"

# Extract the list of paths corresponding to GRCh38
docker run --rm \
    --volume ${DATA_DIR}:/workdir \
    --workdir /workdir \
    quay.io/vgteam/vg:v1.59.0 \
    vg paths -x hprc-v1.1-mc-grch38.gbz -L -Q GRCh38 > ${DATA_DIR}/hprc-v1.1-mc-grch38.paths

# Filter paths list
cd ${DATA_DIR} && \
    grep -v _decoy hprc-v1.1-mc-grch38.paths \
    | grep -v _random \
    | grep -v chrUn_ \
    | grep -v chrEBV \
    | grep -v chrM \
    | grep -v chain_ > hprc-v1.1-mc-grch38.paths.sub