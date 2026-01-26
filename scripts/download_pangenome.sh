#!/bin/bash 

echo "Downloading pangenome data..." 

# Donwload index files
echo -e "dist\nmin\ngbz\nhapl" | xargs -I {} \
    wget "https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.{}"

# Extract the list of paths corresponding to GRCh38
docker run --rm --volume $(pwd):/workdir \
    --workdir /workdir \
    quay.io/vgteam/vg:v1.59.0 \
    vg paths -x hprc-v1.1-mc-grch38.gbz -L -Q GRCh38 > hprc-v1.1-mc-grch38.paths

# Filter paths list
grep -v _decoy hprc-v1.1-mc-grch38.paths \
    | grep -v _random \
    | grep -v chrUn_ \
    | grep -v chrEBV \
    | grep -v chrM \
    | grep -v chain_ > hprc-v1.1-mc-grch38.paths.sub