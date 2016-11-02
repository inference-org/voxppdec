

%slvoxppmodelTest
%
%
%author : steeve laquitaine
%  date : 161101
%purpose: decode likelihood of stimulus feature (e.g., motion directions) 
%         from Ni instances by Nv voxels responses with model parameters
%         W_tr, rho_tr,tau_tr,sigma_tr
%
% usage :
%
%       LLHs = slvoxppmodelTest(b_test)


function LLHs = slvoxppmodelTest(b_test,W_tr,rho_tr,tau_tr,sigma_tr,pp)

%number of data feature (voxels)
Nv = size(W_tr,1);

%calculate likelihoods
%W_tr, rho_tr, tau_tr, sigma_tr, 
%we first calculate them from each trial BOLD pattern
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
LLHs = bsxfun(@rdivide,pb_si,sum(pb_si,2));