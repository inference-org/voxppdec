

%simulate voxel ppc decoding
%
%
%author : steeve laquitaine
%
%
%status : 
%       - The weights are very poorly trained
%       - sigma: when the weights are veridical, sigma trains well
%                when weigths are bad, sigma tends to 0
%       -  rho : sometimes train well 
%       -  tau : 
%
%simulate training instance matrix "b_wres" : 
%200 instances x 100 voxels responses matrix to 200 trials with stimulus 
%5 different stimulus directions
clear;
Ni = 626;
Nv = 320;
svec = randsample([15 85 155 225 295],Ni,'true');
W = rand(Nv,8); %for random weights
sigma = 0.3;
rho = 0.1;
tau = rand(Nv,1);
save('tau','tau')
[b_train,C,~,pp] = slsimPPCinstances(svec,W,sigma,tau,rho);

%% 1st step : train the model's weights (OLS)
%note that the more different stimulus directions we use to train the model
%and the more the fit is constrained and the closer the weights "Wtrained"
%are to veridical (W). This is also true with the number of repeats of 
%the directions: the more stimulus directions are repeated and the closer 
%the weights to veridical.
W_tr = sltrainfmmodel(b_train,C);
slppPlotWeights(W_tr,W,pp.phi_k);

%% 2nd step : train the model's noise parameters (fminsearch Nelder mead)
%The more rectangular b the instance x voxel matrix (Ni >> Nv) 
%the better the training of rho.
%The quality of sigma depends on the quality of the trained weights
%We get eta using the conjugate gradient method. Eta might be inaccurate
%when W'*W is ill-conditioned
[rho_tr,tau_tr,sigma_tr,nglogl] = sltrainPPmodelStep2(b_train,svec,pp.f_k_s,W_tr);

%% calculate posterior
%W_tr, rho_tr, tau_tr, sigma_tr, 
%we first calculate the them from each trial BOLD pattern
%We first try to infer the posteriors out of the training set for sanity
%check
b_test = b_train;
%Nv x Nv mu MV gaussian mean 
mu = W_tr*pp.f_k_s;
%Nv x Nv Omega global noise matrix
Om = rho_tr*(tau_tr*tau_tr') + (1-rho_tr)*times(eye(Nv,Nv),tau_tr*tau_tr')+(sigma_tr^2)*(W_tr*W_tr');
[~,e] = cholcov(Om);
if e~=0
    fprintf('%s \n','(slsimvoxppdec) Covariance matrix Omega is not symmetric, positive definite')
    dbstack
    keyboard    
end
if det(Om)==0
    fprintf('%s \n','(slsimvoxppdec) Covariance matrix Omega singular')
    dbstack 
    keyboard   
end

%likelihood probability distribution
for si = 1 : 360
    pb_si(:,si) = mvnpdf(b_test,mu(:,si)',Om);
end
pbgivs = bsxfun(@rdivide,pb_si,sum(pb_si,2));

%average likelihood by displayed directions
%calculate
s_disp = unique(svec);
for i = 1 : length(s_disp)
    pb_si_pers(i,:) = nanmean(pbgivs(svec == s_disp(i),:),1);
end
%plot
figure('color','w');
cl = linspecer(length(s_disp));
for i = 1 : length(s_disp)
    plot(pb_si_pers(i,:),'color',cl(i,:))
    hold on; plot([s_disp(i) s_disp(i)],...
        [0 max(pb_si_pers(:))],'color',cl(i,:),'linestyle',':')
end
box off






