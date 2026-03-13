# Genomics Analysis Developer Example

Easily run essential genomics workflows to save time leveraging Parabricks and CodonFM. 

- [Overview](#overview)
- [Experience Workflow](#experience-workflow)
  - [Architecture Diagram](#architecture-diagram)
  - [Notebook Outline](#notebook-outline)
- [How to Run](#how-to-run)
  - [Hardware Requirements](#hardware-requirements)
  - [Software Requirements](#software-requirements)
  - [Pre-configured Instances](#pre-configured-instances)
  - [Manual Installation](#manual-installation)
- [References](#references)
- [Terms of Use](#terms-of-use)
- [Ethical Considerations](#ethical-considerations)

## Overview

This developer example enables bioinformaticians to run GPU-accelerated genomics workflows in minutes on any cloud through Brev.dev. [NVIDIA® Parabricks®](https://docs.nvidia.com/clara/parabricks/latest/index.html) powers both linear and graph-based read alignment along with variant calling via DeepVariant. [CodonFM](https://github.com/NVIDIA-Digital-Bio/CodonFM), NVIDIA's RNA foundation model, can then be used to predict the functional impact of each detected variant on specific genes.

## Experience Workflow

This developer example shows how to use GPU accelerated tools for alignment (linear and graph), variant calling, and variant effect prediction. 

### Architecture Diagram

The exact steps to run this workflow are outlined below: 

![](images/genomics-analysis-arch-diagram.svg)

### Notebook Outline 

All the code can be found in Jupyter notebooks in the [`notebooks`](notebooks) directory of the Github repo. 

#### `germline_wes.ipynb`
Runs a standard germline variant calling workflow on whole exome sequencing (WES) data. Downloads the NA12878 sample from the Genome in a Bottle consortium, aligns reads to the GRCh38 reference using GPU-accelerated BWA-MEM via Parabricks fq2bam, and calls variants with GPU-accelerated DeepVariant, producing a final .vcf file.

#### `pangenome.ipynb`
Demonstrates a pangenome analysis workflow as an alternative to single-reference alignment. Downloads the HPRC v1.1 pangenome graph, aligns short-read FASTQ samples using GPU-accelerated Giraffe, and calls variants with Pangenome-Aware DeepVariant — a variant of DeepVariant that uses the pangenome graph to improve alignment accuracy and variant detection across diverse populations.

#### `variant_effect_prediction.ipynb`
Runs a full variant effect prediction pipeline starting from raw FASTQ files. Uses Parabricks to align reads and call variants, processes GENCODE gene annotations to extract protein-coding sequences, maps detected variants onto transcripts, and uses CodonFM (NVIDIA's RNA foundation model) to predict the functional impact of each variant via log likelihood ratios.

## How to Run 

### Hardware Requirements 

The L40s with at least 48GB of GPU memory is recommended for the best combination of cost and performance. Users can also try L4 or T4 (better cost) or RTX Pro 6000 (better performance). 

NVIDIA Parabricks can be run on any NVIDIA GPU that supports CUDA® architecture 75, 80, 86, 89, 90, 100, or 120 and has at least 16GB of GPU RAM. 

Parabricks has been tested specifically on the following NVIDIA GPUs:  

* T4  
* A10, A30, A40, A100, A6000  
* L4, L40  
* H100, H200  
* GH200
* B200, B300
* GB200, GB300
* RTX PRO 6000 Blackwell Server Edition
* RTX PRO 4500
* DGX Spark
* DGX Station

The minimum amount of CPU RAM and CPU threads depends on the number of GPUs. Please refer to the table below: 

| GPUs | Minimum CPU RAM (GB) | Minimum CPU Threads |
|------|----------------|---------------------|
| 2 | 100 | 24 |
| 4 | 196 | 32 |
| 8 | 392 | 48 |

### Software Requirements 

* Any NVIDIA driver that is compatible with CUDA 12.9 (535, 550, 570, 575, or similar). Please check [here](https://docs.nvidia.com/deploy/cuda-compatibility/#forward-compatibility) for more details on forward compatibility.  
* Any Linux operating system that supports Docker version 20.10 (or higher) with the NVIDIA GPU runtime.

### Pre-configured Instances 

These notebooks are available as a launchable on [Brev](https://login.brev.nvidia.com/signin). This is a one-click method, that automatically installs dependencies, provisions hardware, and loads this repository. 

 [![ Click here to deploy.](https://brev-assets.s3.us-west-1.amazonaws.com/nv-lb-dark.svg)](https://brev.nvidia.com/launchable/deploy?launchableID=env-3AjR5pVHTtMm2ToM3KTy6wtQWzE)

### Manual installation 

For users who prefer to run on their own hardware, installation instructions are provided below: 

Prerequisites: `Python3`

```
# Create Python virtual environment and activate it
python3 -m venv .venv 
source .venv/bin/activate

# Run the setup script 
./scripts/local_setup.sh

# Start Jupyter lab 
jupyter lab
```

## References

* [Parabricks Documentation](https://docs.nvidia.com/clara/parabricks/latest/index.html) 
* [CodonFM Blog](https://developer.nvidia.com/blog/introducing-the-codonfm-open-model-for-rna-design-and-analysis/)
* [Parabricks Pangenome Alignemnt Blog](https://developer.nvidia.com/blog/discover-new-biological-insights-with-accelerated-pangenome-alignment-in-nvidia-parabricks/)

## Terms of Use

**Governing Terms**: The Parabricks container is governed by the [NVIDIA Software License Agreement](https://www.nvidia.com/en-us/agreements/enterprise-software/nvidia-software-license-agreement/) and the [Product-Specific Terms for NVIDIA AI Products](https://www.nvidia.com/en-us/agreements/enterprise-software/product-specific-terms-for-ai-products/). This Genomics Analysis Blueprint github repository is provided under [Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0).

## Ethical Considerations

NVIDIA believes Trustworthy AI is a shared responsibility, and we have established policies and practices to enable development for a wide array of AI applications. When downloaded or used in accordance with our terms of service, developers should work with their supporting model team to ensure the models meet requirements for the relevant industry and use case and addresses unforeseen product misuse. For more detailed information on ethical considerations for the models, please see the Model Card++ Explainability, Bias, Safety & Security, and Privacy Subcards. Please report security vulnerabilities or NVIDIA AI Concerns [here](https://www.nvidia.com/en-us/support/submit-security-vulnerability/). 