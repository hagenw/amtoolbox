function si = baumgartner2014similarityestimation(sigma,kv,flags)
%BAUMGARTNER2014SIMILARITYESTIMATION - Similarity estimation with listener-specific sensitivity
%   Usage:     si = baumgartner2014similarityestimation(sigma)
%              si = baumgartner2014similarityestimation(sigma,flags,kv)
%
%   Input parameters:
%     sigma   : internal distance metrics
%
%   Output parameters:
%     si      : similarity indices
%
%   `baumgartner2014similarityestimation(...)` maps internal distance
%   metrics to similarity indices according to listener-specific
%   sensitivity
%
%   `baumgartner2014similarityestimation` accepts the following optional parameters:
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

%% Similarity estimation, Eq.(5)

si=zeros(size(sigma)); % init
for ch = 1:size(sigma,3)
  for it = 1:size(sigma,2)
    si(:,it,ch) = 1+eps - (1+exp(-kv.gamma*(sigma(:,it,ch)-kv.S))).^-1;
  end
end

end