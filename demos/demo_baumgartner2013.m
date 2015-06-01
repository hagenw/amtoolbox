%DEMO_BAUMGARTNER2013 Demo for sagittal-plane localization model from Baumgartner et al. (2013)
%
%   `demo_baumgartner2013(flag)` demonstrates how to compute and visualize 
%   the baseline prediction (localizing broadband sounds with own ears) 
%   for a listener of the listener pool and the median SP using the 
%   sagittal-plane localization model from Baumgartner et al. (2013).
%
%   .. figure::
%
%      Baseline prediction
% 
%      This demo computes the baseline prediction (localizing broadband 
%      sounds with own ears) for an exemplary listener (NH58).
%
%      Predicted polar response angle probability of subject NH58 as a  
%      function of the polar target angle with probabilities encoded by
%      brigthness.
%
%   See also: baumgartner2013

% AUTHOR : Robert Baumgartner

%% Settings

flag='NH58';  % subject ID
lateral = 0;  % lateral target angle in degrees
   

%% Search subject

s = data_baumgartner2013('pool');
for ids = 1:length(s)
    if strcmp(flag,s(ids).id)
        break
    end
end


%% Run model

[targets,tang] = extractsp(lateral,s(ids).Obj);
[p,rang] = baumgartner2013(targets,s(ids).Obj,'u',s(ids).u,'lat',lateral);

figure;
plotbaumgartner2013(p,tang,rang);
title(['Baseline prediction for ' s(ids).id]);

[qe,pe,pb] = baumgartner2013pmv2ppp(p,tang,rang,'print');

