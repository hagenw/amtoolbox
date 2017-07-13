function sigma = baumgartner2016_comparisonprocess(tar,tem)
%baumgartner2016_comparisonprocess - Comparison with direction-specific templates
%   Usage:     sigma = baumgartner2016_comparisonprocess(tar,tem)
%
%   Input parameters:
%     tar     : discharge rate profiles of target sounds (fields: tar.m for
%               magnitude and tar.sd for standard deviation)
%     tem     : discharge rate profiles of templates (fields as for tar)
%
%   Output parameters:
%     sigma   : internal distance metric
%
%   `baumgartner2016_comparisonprocess(...)` compares discharge rate profiles
%   on the basis of the quotient between rate difference and auditory nerve
%   variance (May and Huang, 1997; Reiss et al., 2011)
%
%   References: may1997 reiss2011

% AUTHOR: Robert Baumgartner

%% Comparison process, Eq.(4)

% Unbiased statistical index of rate discrimination acc. to Eq. 1 of
% may1997
sigma=zeros(size(tem.m,2),size(tar,2),size(tem.m,3),size(tem.m,4),size(tar,5)); % init
for itime = 1:size(tar.m,5)
  for itang = 1:size(tar.m,2)
    isd = repmat(tar.m(:,itang,:,:,itime),[1,size(tem.m,2),1,1]) - tem.m;
    sd = sqrt( repmat(tar.sd(:,itang,:,:,itime),[1,size(tem.sd,2),1,1]).^2 + tem.sd.^2);
    dprime = isd./sd;
    sigma(:,itang,:,:,itime) = nanmean(abs(dprime),1);
  end
end

end