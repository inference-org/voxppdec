


%slvoxppmodelTrain
%
%
%author : steeve laquitaine
%purpose: train the parameters of a probabilistic population model of 
%         voxel responses
%

function [W_tr,rho_tr,tau_tr,sigma_tr,pp] = slvoxppmodelTrain(b_train,svec_train)

%% average channel responses to vector of stimulus directions
%Ni instances x 8 stimulus feature - tuned channels
pp = slsimPPchannels(0);
C = pp.f_k_s(:,svec_train)';

%% 1st step : train the model's weights (OLS)
%note that the more different stimulus directions we use to train the model
%and the more the fit is constrained and the closer the weights "Wtrained"
%are to veridical (W). This is also true with the number of repeats of 
%the directions: the more stimulus directions are repeated and the closer 
%the weights to veridical.
W_tr = sltrainPPmodel(b_train,C);

%% 2nd step : train the model's noise parameters (fminsearch Nelder mead)
%The more rectangular b the instance x voxel matrix (Ni >> Nv) 
%the better the training of rho.
%The quality of sigma depends on the quality of the trained weights
%We get eta using the conjugate gradient method. Eta might be inaccurate
%when W'*W is ill-conditioned
[rho_tr,tau_tr,sigma_tr] = sltrainPPmodelStep2(b_train,svec_train,pp.f_k_s,W_tr);
