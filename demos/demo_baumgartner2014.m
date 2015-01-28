%DEMO_BAUMGARTNER2014 Demo for sagittal-plane localization model from Baumgartner et al. (2014)
%
%   `demo_baumgartner2014(flag)` demonstrates how to compute and visualize 
%   the baseline prediction (localizing broadband sounds with own ears) 
%   for a listener of the listener pool and the median plane using the 
%   sagittal-plane localization model from Baumgartner et al. (2014).
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
%   See also: baumgartner2014 exp_baumgartner2014 baumgartner2014virtualexp
%   localizationerror

% AUTHOR : Robert Baumgartner

%% Settings

subID = 'NH58';   % subject ID of exemplary listener
lat = 0;          % lateral target angle in degrees
runs = 3;         % # of virtual experimental runs
   

%% Get listener's data

s = data_baumgartner2014('pool');   % load data of listener pool
ids = find(ismember({s.id},subID));  % index of exemplary listener


%% Run model with individual sensitivity S

[p,rang,tang] = baumgartner2014(s(ids).Obj,s(ids).Obj,'S',s(ids).S,'lat',lat);


%% Run virtual experiment

m = baumgartner2014virtualexp(p,tang,rang,'runs',2);


%% Calcualte performance measures 

amtdisp('Performance Predictions:')
amtdisp('------------------------')

% via expectancy values:
[qe,pe] = baumgartner2014pmv2ppp(p,tang,rang,'print'); 

% and/or via responses drawn from virtual experiments
[f,r] = localizationerror(m,'sirpMacpherson2000');
perMacpherson2003 = localizationerror(m,f,r,'perMacpherson2003');
amtdisp(['Local polar error rate (%)        ' num2str(perMacpherson2003,'%4.1f')])


%% Plot results

figure;
plotbaumgartner2014(p,tang,rang,m(:,6),m(:,8));
title(['Baseline prediction for ' s(ids).id]);
