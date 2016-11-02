
%slppmodelPlotDecodedLLHs.m
%
%
%author : steeve laquitaine
%purpose: plot decoded likelihood of stimulus feature (e.g., motion directions) 
%         from Ni instances by Nv voxels responses with model parameters
%         contained in model.
%usage :
%
%       slppmodelPlotDecodedLLHs(LLHs,directions_test,model,'V1','priorUnif')

function slppmodelPlotDecodedLLHs(LLHs,directions_test,model,roi,prior)

%plot decoded likelihoods averaged by directions
s_disp = unique(directions_test);
figure('color','w');
cl = linspecer(length(s_disp));

%each direction
for i = 1 : length(s_disp)
    
    %likelihood with sem over direction instances
    errorarea(LLHs(directions_test == s_disp(i),:)','k',cl(i,:))
    
    %average
    nanmean(LLHs(directions_test == s_disp(i),:))
    hold on; plot([s_disp(i) s_disp(i)],...
        [min(nanmean(LLHs(directions_test == s_disp(i),:))) max(nanmean(LLHs(directions_test == s_disp(i),:)))],'color',cl(i,:),...
        'linestyle',':','linewidth',2)
end
box off
xlabel('Hypothetical motion directions (deg)')
ylabel('Likelihood (probability)')
title({['Average LLHs decoded by ppc from' roi 'responses by directions (colors)'],...
    '(area:sem over directions)',...
    ['Trained on uniform prior - tested on ' prior],...
    ['rho:' num2str(model.rho_tr) ' - sigma:' num2str(model.sigma_tr) ' - mean(tau):' num2str(mean(model.tau_tr))]})
xlim([0 360])

% %dynamically plot each trial likelihood
% figure('color','w')
% for i = 1 : size(LLHs,1)
%     plot(LLHs(i,:))
%     drawnow
%     pause(0.5)
% end