#!/bin/bash
set -euo pipefail

# SPDX-FileCopyrightText: Copyright (c) 2024-2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: Apache-2.0

# Make directories for data and references files
DATA_DIR=$(pwd)/../data
mkdir -p $DATA_DIR 

# Download the data 
# Note: The tarball downloads unnecessary files and can be streamlined in the future 
cd $DATA_DIR && \
    wget -nc --progress=bar:force 2>&1 https://s3.amazonaws.com/parabricks.sample/parabricks_sample.tar.gz 

cd $DATA_DIR/ref && \
    printf "%s\n" dist min gbz | xargs -I {} \
    wget -nc --progress=bar:force 2>&1 "https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.{}"

# Pre-process the data 
cd $DATA_DIR && \
    pigz -dc parabricks_sample.tar.gz | tar xvf - && \
    mv parabricks_sample/Data/sample_* . 
