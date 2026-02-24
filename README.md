# Genomics Analysis Blueprint

- [Overview](#overview)
- [Notebooks](#notebooks)
  - [Germline WES](germline_wes.ipynb)
  - [Pangenome](pangenome.ipynb)
  - [Variant Effect Prediction](variant_effect_prediction.ipynb)
- [System Requirements](#system-requirements)
- [Terms of Use](#terms-of-use)

## Overview

This repository contains notebooks for GPU-accelerated genomic analysis using [Parabricks](https://docs.nvidia.com/clara/parabricks/latest/index.html), NVIDIA's software suite for secondary genomic analysis. Parabricks accelerates standard bioinformatics tools by 100x or more compared to CPU-only pipelines — reducing runtimes from hours to minutes — without changing outputs or requiring workflow modifications.

The goal of this repository is to help users quickly explore three workflows: germline variant calling on whole exome data, pangenome-based alignment and variant detection, and variant effect prediction using NVIDIA's CodonFM RNA foundation model. Workflows can be run on any CUDA-capable GPU system or through the quick deploy capability of [Brev.dev](https://developer.nvidia.com/brev) Launchables. Sequencing data is publicly available from the [Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) Consortium.

These notebooks are designed for any bioinformatics scientist or developer who wants to experience GPU-accelerated genomics firsthand — and because Parabricks scales efficiently across larger GPUs and multi-GPU systems, the same workflows demonstrated here can be applied directly to cohort-scale or production datasets without re-engineering. For more information, see the [latest Parabricks documentation and release information](https://docs.nvidia.com/clara/parabricks/latest/index.html).

## Notebooks

This repository contains three end-to-end notebooks that demonstrate GPU-accelerated genomics workflows using NVIDIA Parabricks. Each notebook is self-contained and walks through downloading data, running analysis tools, and interpreting results.

- **[germline_wes.ipynb](germline_wes.ipynb)** — Runs a standard germline variant calling workflow on whole exome sequencing (WES) data. Downloads the NA12878 sample from the Genome in a Bottle consortium, aligns reads to the GRCh38 reference using GPU-accelerated BWA-MEM via Parabricks `fq2bam`, and calls variants with GPU-accelerated DeepVariant, producing a final `.vcf` file.

- **[pangenome.ipynb](pangenome.ipynb)** — Demonstrates a pangenome analysis workflow as an alternative to single-reference alignment. Downloads the HPRC v1.1 pangenome graph, aligns short-read FASTQ samples using GPU-accelerated Giraffe, and calls variants with Pangenome-Aware DeepVariant — a variant of DeepVariant that uses the pangenome graph to improve alignment accuracy and variant detection across diverse populations.

- **[variant_effect_prediction.ipynb](variant_effect_prediction.ipynb)** — Runs a full variant effect prediction pipeline starting from raw FASTQ files. Uses Parabricks to align reads and call variants, processes GENCODE gene annotations to extract protein-coding sequences, maps detected variants onto transcripts, and uses CodonFM (NVIDIA's RNA foundation model) to predict the functional impact of each variant via log likelihood ratios.

## System Requirements

| Requirement | Notes |
| -------- | ------- |
| GPU  | We recommend L40S to balance cost and performance.  <br> - Higher performance: A100 <br> - Better cost: L4 and T4 <br>|
| GPU Memory | 48 GB is recommended. <br> - All tools require at least 16 GB of available GPU memory. <br> - For GPUs with 16-48 GB memory, the --low-memory flag is required. |
| System RAM | At least 100 GB. |
| CPU | At least 24 CPU threads. |
| Driver | NVIDIA Driver version 525.60.13 or greater. See documentation about [forward compatibility](https://docs.nvidia.com/deploy/cuda-compatibility/#forward-compatibility). |
| OS | Any Linux OS that supports nvidia-docker2 and Docker version 20.10 or higher. |

Users may have to wait 5-10 minutes for the instance to start depending on cloud availability. 

## Terms of use
**Governing Terms**: The Parabricks container is governed by the [NVIDIA Software License Agreement](https://www.nvidia.com/en-us/agreements/enterprise-software/nvidia-software-license-agreement/) and the [Product-Specific Terms for NVIDIA AI Products](https://www.nvidia.com/en-us/agreements/enterprise-software/product-specific-terms-for-ai-products/). This Genomics Analysis github repository is provided under [Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0).
