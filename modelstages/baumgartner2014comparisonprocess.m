function sigma = baumgartner2014comparisonprocess(tar,tem,kv,flags)
%BAUMGARTNER2014COMPARISONPROCESS - Comparison with direction-specific templates
%   Usage:     sigma = baumgartner2014comparisonprocess(tar,tem)
%              sigma = baumgartner2014comparisonprocess(tar,tem,flags,kv)
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
%   `baumgartner2014comparisonprocess` accepts the following optional parameters:
%
%     'flags',flags  Transfer flags. If not defaults of baumgartner2014 are used.
%
%     'kv',kv        Transfer key-value pairs. If not defaults of baumgartner2014 are used.
%
%   References: baumgartner2014modeling

% AUTHOR: Robert Baumgartner

if not(exist('flags','var') || exist('kv','var'))
  definput.import={'baumgartner2014'};
  [flags,kv]=ltfatarghelper({},definput,{});
end

%% Comparison process, Eq.(4)
sigma=zeros(size(tem,2),size(tar,2),size(tem,3),size(tem,4)); % init
for it = 1:size(tar,2)
  isd = repmat(tar(:,it,:,:),[1,size(tem,2),1,1]) - tem;
  sigma(:,it,:,:) = mean(abs(isd),1);
end

end