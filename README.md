
=========
voxppdec
=========

author : steeve laquitaine  
mail : steeve@stanford.edu  
keywords : machine learning, Bayesian statistics, generative modeling of the data, reconstruction

This small package uses a matrix of Ni response instance x Nv voxels and a vector of Ni stimulus feature (motion stimulus direction in deg) 
models the voxel responses with a probabilistic population code (van Bergen 2015, Nature Neuroscience) that assumes that voxels' responses are the linear combination of the responses 
of 8 noisy half-rectified cosine channels raised to the power of 5 with different phases (feature preference) to an observed motion direction.
Three sources of noises are modelled : multivariate Gaussian noise shared by channels with same tunings, multivariate Gaussian noise shared 
globally among voxels irrespective of their tunings, multivariate noise specific to individual voxels.

usage : 

```matlab
>> run slsimvoxppdec.m
```
