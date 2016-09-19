


%analysis00.m
%
%
%
%author : steeve laquitaine
%  date : 160919
%purpose: model data sample of actual brain voxel responses with a 5 fold 
%         cross-validated forward model that assumes voxel responses are 
%         the linear sum of 5 sinewave "channel" function with 5 different 
%         direction preference (phase) raised to power 4.
%         Data are concatenated over scans from diffferent prior
%         conditions (prior mean are 135 and 225 deg)



%Each matrix stacks data from 3 sessions from the same subject.
%We collect both matrices and stack them together
%prior 225
load('data/prior225/V1/instanceMatrix.mat')
instances = d.instances;
svec = d.myRandomDir;
%%get most motion-responsive voxels
%%note that it seems to hurt the representation to 
%%remove voxels
%load('data_sample/r2_V1.mat')
%instances = instances(:,r2 >= 0.05);

%% cross validated likelihood decoding
pp.phi_k = unique(svec);
pp.exponent = 4;
[LLH_f,pp] = slvoxppKFoldCVdec(instances,svec,5,pp);


