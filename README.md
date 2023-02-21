# DeepHap
A Deep Learning Approach for Haplarithmisis data analysis and anomaly detection



#### Introduction 
Haplarithmisis is a computational pipeline that uses continuous B allele frequency (BAF) values to detect haplotypes and copy numbers, as well as identify the parental origin of chromosomal abnormalities. However, the interpretation of haplarithm plots can be challenging due to the complexity of the data, the presence of noise, and the difficulty of identifying and interpreting breakpoints and anomalies. In this study, we developed DeepHap, a platform for haplarithmisis data analysis and anomaly detection

##### What is DeepHap?
DeepHap consists of three main components: a breakpoint detection module using a two-layer long short-term memory (LSTM) model, a robust mode-based denoising module, and a rule-based anomaly assignment module.  


![deephap_small](https://user-images.githubusercontent.com/91246296/220306055-78eb8480-1a0d-4390-a258-faf94d8d3198.JPG)



#### Pre-requisites

Python 3.5
matplotlib==3.2.1,
scipy==1.4.1.

#### Set up
Install dependencies within your virtual environment

 ```bash
    pip install -r source/requirements.txt
    ```
