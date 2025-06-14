# DeepHap
![_1e45051c-3812-47c2-a6c4-4c38260dfe83](https://user-images.githubusercontent.com/91246296/227034145-859a89ea-3c8d-4311-93f5-9e9f38c4f77e.jpg)

A Deep Learning Approach for Haplarithmisis data analysis and anomaly detection



### Abstract

Haplarithmisis is a powerful computational pipeline that detects haplotypes, copy numbers, and identifies the parental origin of chromosomal abnormalities using continuous B allele frequency (BAF) values. However, interpreting haplarithm plots can be a challenge due to the complexity of the data, the presence of noise, and difficulty in identifying and interpreting breakpoints and anomalies. To address these challenges, we developed DeepHap, a deep learning-based platform for haplarithmisis data analysis and anomaly detection.

### Features

    Breakpoint detection module using a two-layer long short-term memory (LSTM) model
    Robust mode-based denoising module
    Rule-based anomaly assignment module


### Getting Started


### Pre-requisites

    Python 3.5
    matplotlib==3.2.1
    scipy==1.4.1

### Set up
Install dependencies within your virtual environment

```bash
    pip install -r source/requirements.txt
```
    
### Usage

To use DeepHap, users can follow the instructions provided in the user manual, which is included in the repository. The repository contains the following files:

    main.py: the main script that implements DeepHap.
    user_manual.pdf: the user manual that includes instructions on how to use DeepHap.
    sample_data/: a directory that contains sample data for testing DeepHap.

### Update V2

```bash
    python.exe .\hadi_cli.py -w 300  input_dir  output_dir
```

### Contributing

Contributions to DeepHap are welcome. Please see the CONTRIBUTING.md file for more information.

### License

DeepHap is released under the MIT License. Please see the LICENSE file for more information.

### Authors

 Mohsen Yazdani (msilico1@email.com)
    
### Acknowledgments

We would like to thank the Laboratory of Bioinformatics and Drug Design (LBD), University of Tehran for funding this project. We also acknowledge the contributions of our colleagues at Maastricht University who provided valuable feedback and suggestions during the development of DeepHap.
