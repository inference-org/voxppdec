

%slvoxppmodelTest
%
%
%author : steeve laquitaine
%purpose: train the parameters of a probabilistic population model of 
%         voxel responses
%
% usage :
%
%       pbgivs = slvoxppmodelTest(b_test)


function pbgivs = slvoxppmodelTest(b_test,W_tr,rho_tr,tau_tr,sigma_tr,pp)

%number of data feature (voxels)
Nv = size(W_tr,1);

%calculate likelihoods
%W_tr, rho_tr, tau_tr, sigma_tr, 
%we first calculate the them from each trial BOLD pattern
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


% %plot average decoded llhs by directions
% figure('color','w');
% cl = linspecer(length(s_disp));
% for i = 1 : length(s_disp)
%     
%     %posterior with sem over direction instance repeats
%     errorarea(pbgivs(svec == s_disp(i),:)','k',cl(i,:))
%     
%     %displayed directions
%     hold on; plot([s_disp(i) s_disp(i)],...
%         [min(pb_si_pers(:)) max(pb_si_pers(:))],'color',cl(i,:),...
%         'linestyle',':','linewidth',2)
% end
% box off
% xlabel('Hypothetical motion directions (deg)')
% ylabel('Likelihood (probability)')
% title({'Average LLHs decoded by ppc from V1 bold patterns by directions (colors)',...
%     '(area is sem over directions repeats)',...
%     ['rho:' num2str(rho_tr) ' - sigma:' num2str(sigma_tr) ' - mean(tau):' num2str(mean(tau_tr))]})
% xlim([0 360])