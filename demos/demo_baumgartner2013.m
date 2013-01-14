%DEMO_BAUMGARTNER2013 Demo for sagittal-plane localization model from Baumgartner et al. (2013)
%   Usage: p = demo_baumgartner2013( flag )
%
%   `demo_baumgartner2013(flag)` demonstrates how to compute and visualize 
%   the baseline prediction (localizing broadband sounds with own ears) 
%   for a listener of the listener pool and the median SP using the 
%   sagittal-plane localization model from Baumgartner et al. (2013).
%
%   .. figure::
%
%     This is the headline of the figure, one line at most
%
%     This is the description of the figure. It can be
%     several lines.
%
%     and span several paragraphs.
%
%   See also: baumgartner2013

% AUTHOR : Robert Baumgartner

flag='NH58';
    
s = data_baumgartner2013('pool');
for ids = 1:length(s)
    if strcmp(flag,s(ids).id)
        break
    end
end

[spdtfs,tang] = extractsp(0,s(ids).dtfs,s(ids).pos);
[p,rang] = baumgartner2013(spdtfs,spdtfs,'u',s(ids).u);

figure;
plotbaumgartner2013(p,tang,rang);
title(['Baseline prediction for ' s(ids).id]);

[qe,pe,pb] = pmv2ppp(p,'print');

