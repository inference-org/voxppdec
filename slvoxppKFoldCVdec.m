%slvoxppKFoldCVdec.m
%
%
%
% author: steeve laquitaine
%purpose: models voxel responses to a stimulus feature (motion directions)
%         with a probabilistic population code and use the model to decode
%         the likelihood of the responses given hyothetical feature values
%         (direction space). Decoding is 5 fold cross-validated. Training 
%         on 4/5 of the data and testing on the remaining fold. The number of 
%         Voxel responses instances is balanced between feature values 
%         (directions) to prevent biased training toward the more frequent
%         feature value.
%
%usage :
%
%             slfmriInitAnalysisTaskDotDirfMRI05 
%             d = slfmriGetDBoverSessions(o,varargin); 
%             b = d.instances; svec = d.myRandomDir;
%             [LLH_f,pp] = slvoxppKFoldCVdec(b,svec,5);


function [LLH_f,pp] = slvoxppKFoldCVdec(instances,svec,Nf,pp)

%divide the dataset into 5 folds
%5 folds per stimulus feature (motion direction)
%to make sure data for all stimulus feature are present
%during training. The number of trials is equalized by 
%condition to prevent biasing training toward certain
%features.
%sort instances in structure by stimulus feature
[i_struc,v_u] = slmat2structByVar(instances,svec);

%balance instances between simulus features
instancesBal = setNumInstancesByClass(i_struc,'balanceByRemovI');

%get vectors of associated stimuli features
v_u_all = v_u(:,ones(1,size(instancesBal{1},1)))';

%get row starting and ending indices for each fold 
%by stimulus feature (col)
[f_st,f_end,~,Nivi] = slgetkFoldrows(instancesBal,Nf);

%decode from each train/test fold combinationW_tr_f(:,:,f) = W_tr;
Nv = size(instances,2);
Ni = size(instances,1);
Nk = 8;
rho_tr_f = nan(1,Nf);
tau_tr_f = nan(Nf,Nv);
sigma_tr_f = nan(1,Nf);
W_tr_f = nan(Nv,Nk,Nf);
LLH_f = nan(Nivi,360,Nf);
svecTest_f = nan(Nivi,Nf);
    
%Train and test each fold combination
for f = 1 : Nf
    
    %this fold test and train instances
    %picked up from 1/Nf of each stimulus feature
    testInst = []; trainInst = []; svecTrain = []; svecTest = [];
    for s_u = 1 : length(unique(svec))        
        
        %this fold test instances
        testInst = [testInst; instancesBal{s_u}(f_st(f) : f_end(f),:)];        
        %this fold's train instance positions
        trainIx = setdiff(1:Nivi,f_st(f) : f_end(f));
        %this fold's train instances
        trainInst = [trainInst; instancesBal{s_u}(trainIx ,:)];        
        %this fold's train stimulus feature
        svecTrain = [svecTrain; v_u_all(trainIx,s_u)];
        %this fold's test stimulus feature
        svecTest = [svecTest; v_u_all(f_st(f):f_end(f),s_u)];
    end
    
    %decode this test fold likelihoods (train then test)
    [W_tr,rho_tr,tau_tr,sigma_tr,pp] = slvoxppmodelTrain(trainInst,svecTrain,pp);        
    LLH = slvoxppmodelTest(testInst,W_tr,rho_tr,tau_tr,sigma_tr,pp);    
    
    %backup this fold
    W_tr_f(:,:,f) = W_tr;
    rho_tr_f(f) = rho_tr;
    tau_tr_f(f,:) = tau_tr';    
    sigma_tr_f(f) = sigma_tr;
    LLH_f(:,:,f) = LLH;
    svecTest_f(:,f) = svecTest;
end

%backup each fold data
pp.W_tr_f = W_tr_f;
pp.rho_tr_f = rho_tr_f;
pp.tau_tr_f = tau_tr_f;
pp.sigma_tr_f = sigma_tr_f;

%backup means over folds
pp.meanOverFolds_W_tr_f = mean(W_tr_f,3);
pp.meanOverFolds_rho_tr_f = mean(rho_tr_f);
pp.meanOverFolds_tau_tr_f = mean(tau_tr_f,1);
pp.meanOverFolds_sigma_tr_f = mean(sigma_tr_f);

%backup stds over folds
pp.stdOverFolds_W_tr_f = std(W_tr_f,0,3);
pp.stdOverFolds_rho_tr_f = std(rho_tr_f);
pp.stdOverFolds_tau_tr_f = std(tau_tr_f,1);
pp.stdOverFolds_sigma_tr_f = std(sigma_tr_f);

%calculate average LLH by stimulus feature (direction) over folds and trials
figure('color','w');
cl = linspecer(length(v_u));
LLH_s_f = cell(length(v_u));
meanLLH_all = nan(length(v_u),360);
stdLLH_all = nan(length(v_u),360);

%loop over displayed directions
for i = 1 : length(v_u)
    
    %get all LLH for this displayed direction    
    LLH_s_f_tmp = [];
    for f = 1 : Nf
        LLH_s_f_tmp = [LLH_s_f_tmp; LLH_f(svecTest_f(:,f)==v_u(i),:,f)];
    end
    LLH_s_f{i} = LLH_s_f_tmp;
    
    %LLH mean and std
    meanLLH_all(i,:) = mean(LLH_s_f_tmp,1);
    stdLLH_all(i,:) = std(LLH_s_f_tmp);    
    
    %likelihood with sem over direction repeats
    errorarea(LLH_s_f_tmp','k',cl(i,:))    
    
    %stimulus feature (directions)
    hold on; plot([v_u(i) v_u(i)],...
        [min(meanLLH_all(:)) max(meanLLH_all(:))],'color',cl(i,:),...
        'linestyle',':','linewidth',2)    
end

box off
xlabel('Hypothetical motion directions (deg)')
ylabel('Likelihood (probability)')
title({'Average LLHs decoded by cross-validated ppc from V1 bold patterns by directions (colors)',...
    '(area is sem over directions repeats)',...
    ['mean rho:' num2str(mean(rho_tr_f)) ' - mean sigma:' num2str(mean(sigma_tr_f)) ' - mean(tau):' num2str(mean(mean(tau_tr_f)))]})
xlim([0 360])