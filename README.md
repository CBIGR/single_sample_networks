# Evaluation of sample-specific network inference methods 
Code accompanying the paper "Evaluation of single-sample network inference methods for precision oncology" 

## Folder resources 
- Driver genes per cancer subtype, downloaded from [IntOGen](https://www.intogen.org/search) and [COSMIC](https://cancer.sanger.ac.uk/cosmic) databases, HumanNet-XN network, downloaded from [HumanNetv2](http://www.inetbio.org/humannet/), and the preprocessed expression data for lung and brain.

## Folder scripts
Scripts used to create the results in this paper. 
- Preparation: The cell line selection and preprocessing file to prepare the DepMap expression data for network inference. Expression data downloaded from [DepMap](https://depmap.org/portal/) version 20Q4. 
- Network_inference: All scripts used to infer the sample-specific networks (SSN, LIONESS, iENA, SSPGI, CSN) and the make_one_Edgelist file to create one file combining all samples from separate files per sample, used for SSN and iENA. 
- Postprocessing: Select top n edges and scale network files for postprocessing of the sample-specific networks. 
- Results: All scripts to create the results and figures in the paper. 
