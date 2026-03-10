# Utilities for manipulating Cas9 ground truth datasets

These tools accompany our released Cas9 datasets with encoded ground-truth heterogeneity. More information about the datasets can be found [here](). 

# Installing 
The resampling scripts provided here have very few dependencies and could be run in most conda environments. The required packages are numpy and pandas. 


# Using 
The two scripts provided here can be used to subsample from any provided .star file based on dataset of origin. Users just need to pass the input .star file (below, input.star), along with the desired number of particles to be sampled from each origin dataset and an output file name.   

For example, the following generates a new file (```split_10_80_10.star```) with 10% 2 bp particles, 80% 8 bp particles, and 10% 14 bp particles.  
```
python generate_distributions.py --infile input.star --two 10000 --eight 80000 --fourteen 10000 --outfile split_10_80_10.star
```  

The ```generate_distributions.py``` script only allows 2, 8, and 14 bp particles to be selected, whereas the ```generate_distributions_fullstack.py``` script allows particles to be selected from all datasets. 
