function sigma = baumgartner2014_comparisonprocess(tar,tem)
%baumgartner2014_comparisonprocess - Comparison with direction-specific templates
%   Usage:     sigma = baumgartner2014_comparisonprocess(tar,tem)
%
%   Input parameters:
%     tar     : internal representation(s) of target profile(s)
%     tem     : internal templates
%
%   Output parameters:
%     sigma   : internal distance metric
%
%   `baumgartner2014_comparisonprocess(...)` performs spectro-spatial
%   comparison process between internal representation of incoming sound
%   (target) and the templates of the sagittal plane
%
%   References: baumgartner2014modeling

% AUTHOR: Robert Baumgartner

%% Comparison process, Eq.(4)
sigma=zeros(size(tem,2),size(tar,2),size(tem,3),size(tem,4),size(tar,5)); % init
for itime = 1:size(tar,5)
  for itang = 1:size(tar,2)
    isd = repmat(tar(:,itang,:,:,itime),[1,size(tem,2),1,1]) - tem;
    sigma(:,itang,:,:,itime) = mean(abs(isd),1);
  end
end

end