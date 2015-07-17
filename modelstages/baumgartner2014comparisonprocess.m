function sigma = baumgartner2014comparisonprocess(tar,tem)
%BAUMGARTNER2014COMPARISONPROCESS - Comparison with direction-specific templates
%   Usage:     sigma = baumgartner2014comparisonprocess(tar,tem)
%
%   Input parameters:
%     tar     : internal representation(s) of target profile(s)
%     tem     : internal templates
%
%   Output parameters:
%     sigma   : internal distance metric
%
%   `baumgartner2014comparisonprocess(...)` performs spectro-spatial
%   comparison process between internal representation of incoming sound
%   (target) and the templates of the sagittal plane
%
%   References: baumgartner2014modeling

% AUTHOR: Robert Baumgartner

%% Comparison process, Eq.(4)
sigma=zeros(size(tem,2),size(tar,2),size(tem,3),size(tem,4)); % init
for it = 1:size(tar,2)
  isd = repmat(tar(:,it,:,:),[1,size(tem,2),1,1]) - tem;
  sigma(:,it,:,:) = mean(abs(isd),1);
end

end