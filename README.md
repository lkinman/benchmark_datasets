# Utilities for manipulating Cas9 ground truth datasets

These tools accompany our released Cas9 datasets with encoded ground-truth heterogeneity. More information about the datasets can be found [here](). 

# Installing 
The resampling scripts provided here have very few dependencies and could be run in most conda environments. The required packages are numpy and pandas. 


# Using 
The two scripts provided here can be used to subsample from any provided .star file based on dataset of origin. Users just need to pass the input .star file (below, input.star), along with the desired number of particles to be sampled from each origin dataset and an output file name.  

The ```--num_particles``` flag is used to indicate the desired number of particles from each dataset. Users should provide a list of 13 integer values, with the desired number of particles from the 2 bp dataset listed first and the desired number of particles from the 14 bp dataset listed last.  
  
For example, the following generates a new file (```split_10_80_10.star```) with 10000 2 bp particles, 80000 8 bp particles, and 10000 14 bp particles:  
```
python generate_distributions.py --infile input.star --num_particles 10000 0 0 0 0 0 80000 0 0 0 0 0 10000 --outfile split_10_80_10.star
```  
  
Users can also select subsets based on the conformational label from deep classification with either ```--conformation_list``` or ```--conformation-range```. The former allows users to provide a discrete list of labels they want to select, whereas the latter selects all particles with label within a given range.  

For example, the following generates a new file (```split_10_80_10_labels_0_20_40.star```) with 10000 2 bp particles, 80000 8 bp particles, and 10000 14 bp particles, where all particles have label 0, 20, or 40:  
```
python generate_distributions.py --infile input.star --conformation_list 0 20 40 --num_particles 10000 0 0 0 0 0 80000 0 0 0 0 0 10000 --outfile split_10_80_10_labels_0_20_40.star
```

whereas the following creates a similar distribution of files but with particle labels from 0 to 20 (inclusive):  
```
python generate_distributions.py --infile input.star --conformation_range 0 20 --num_particles 10000 0 0 0 0 0 80000 0 0 0 0 0 10000 --outfile split_10_80_10_labels_0to20.star
```
