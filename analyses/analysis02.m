
%slvoxppDataAnalysis_cv.m
%
%
% author: steeve laquitaine
%purpose: train and test probabilitic population code model (WJM 2005, NN) 
%         on one prior dataset (225) and compare with train on one prior (225), 
%         test on another (135). Are the likelihoods systematically biased 
%         toward 225 deg ? If yes the patterns of instances contain a
%         representation of the prior mean.
%         

% set train (uniform prior) and test (prior 125) datasets
%load uniform prior dataset
cd ~/Dropbox/voxppRes/data/V1
load ClassifStckSessmyRandomDir_V1
cd ~/Dropbox/voxppRes/

% ------------------------  TRAIN ------------------------------

%get Ni instances x Nv voxels matrix from which to reconstruct likelihood
instances = cell2mat(o.instances');
Ni = size(instances,1);
Nv = size(instances,2);
%get trial vector of displayed directions
s_u = nan(1,length(s.variable));
svec = [];
for i = 1 : length(s.variable)
    %get unique displayed stimulus directions
    s_u = str2double(s.variable{i}(find(s.variable{i}=='=')+1:end));    
    %trial vector of displayed directions
    nI = size(o.instances{i},1);
    svec = [svec; s_u(ones(nI,1))];
end

%balance instances by stimulus feature (motion direction)
%to prevent biasing training (particularly fminsearch) toward 
%most frequent instances
%sort instances in structure by stimulus feature
[i_struc,v_u] = slmat2structByVar(instances,svec);

%balance instances between simulus features
instancesBalstruc = setNumInstancesByClass(i_struc,'balanceByRemovI');

%get vectors of associated stimuli features
v_u_all = v_u(:,ones(1,size(instancesBalstruc{1},1)))';

%reshape as a training matrix and vector
b_train = cell2mat(instancesBalstruc');
svec = v_u_all(:);

%hypothetical sinewave channels (no noise)
%K channels x 360 hypothetical directions space
pp = slsimPPchannels(0);

%average channel responses to vector of stimulus directions
%Ni instances x 8 channels
C = pp. f_k_s(:,svec)';
W_tr = sltrainPPmodel(b_train,C);
figure('Color','w'); imagesc(W_tr); xlabel('channels by direction preference (deg)')
set(gca,'xtick',1:length(pp.phi_k),'xticklabel',pp.phi_k)

%2nd step : train the model's noise parameters
[rho_tr,tau_tr,sigma_tr,nglogl] = sltrainPPmodelStep2(b_train,svec,pp.f_k_s,W_tr);


% ------------------------  TEST ------------------------------

%% test 
cd ~/Dropbox/voxppRes/data/prior135/V1
load instanceMatrix.mat
cd ~/Dropbox/voxppRes/

%% calculate likelihoods in prior 135 test dataset
%we know W_tr, rho_tr, tau_tr, sigma_tr, from training on prior 225 dataset
%load prior 225 dataset for training

%set test dataset (no need to balance here
%training is over)
b_test = d.instances;
svec_test = d.myRandomDir;

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

%% plot average decoded llhs by directions
s_disp = unique(svec_test);
figure('color','w');
cl = linspecer(length(s_disp));
for i = 1 : length(s_disp)
    %posterior with sem over direction instance repeats
    errorarea(pbgivs(svec_test == s_disp(i),:)','k',cl(i,:))
    %displayed directions
    nanmean(pbgivs(svec_test == s_disp(i),:))
    hold on; plot([s_disp(i) s_disp(i)],...
        [min(nanmean(pbgivs(svec_test == s_disp(i),:))) max(nanmean(pbgivs(svec_test == s_disp(i),:)))],'color',cl(i,:),...
        'linestyle',':','linewidth',2)
end
box off
xlabel('Hypothetical motion directions (deg)')
ylabel('Likelihood (probability)')
title({'Average LLHs decoded by ppc from V1 bold patterns by directions (colors)',...
    '(area is sem over directions repeats)',...
    'Trained on p225 (balanced by dir) and tested on p135',...
    ['rho:' num2str(rho_tr) ' - sigma:' num2str(sigma_tr) ' - mean(tau):' num2str(mean(tau_tr))]})
xlim([0 360])


