
=========
voxppdec
=========

author : steeve laquitaine  
mail : steeve@stanford.edu  
keywords : machine learning, Bayesian statistics, generative modeling of the data, reconstruction  

Required packages : sldataMunging

This small package uses a matrix of Ni response instances (observations) x Nv voxels (data features) and a vector of Ni stimulus features (motion stimulus direction in deg), models the voxel responses with a generative probabilistic population code (van Bergen 2015, Nature Neuroscience) that assumes that voxels' responses are the linear combination of the responses 
of 8 noisy half-rectified cosine channels raised to the power of 5, characterized by different phases (stimulus feature selectivities) to an observed motion direction.
Three sources of noises in the voxel responses are modelled : multivariate Gaussian noise shared by channels with same tunings, multivariate Gaussian noise shared 
globally among voxels irrespective of their tunings, multivariate noise specific to individual voxels.

usage : 

```matlab
>> run slsimvoxppdec.m
```
