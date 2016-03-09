function si = baumgartner2014similarityestimation(sigma,varargin)
%BAUMGARTNER2014SIMILARITYESTIMATION - Similarity estimation with listener-specific sensitivity
%   Usage:     si = baumgartner2014similarityestimation(sigma)
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
%     'S',S          Set the listener-specific sensitivity threshold 
%                    (threshold of the sigmoid link function representing 
%                    the psychometric link between transformation from the
%                    distance metric and similarity index) to *S*. 
%                    Default value is 1.
%
%     'gamma',G      Set the degree of selectivity 
%                    (slope of the sigmoid link function representing 
%                    the psychometric link between transformation from the
%                    distance metric and similarity index) to *G*. 
%                    Default value is 6.
%
%   References: baumgartner2014modeling

% AUTHOR: Robert Baumgartner

definput.import={'baumgartner2014'};
[flags,kv]=ltfatarghelper({},definput,varargin);

%% Similarity estimation, Eq.(5)

si=zeros(size(sigma)); % init
for ch = 1:size(sigma,3)
  for it = 1:size(sigma,2)
    si(:,it,ch) = 1+eps - (1+exp(-kv.gamma*(sigma(:,it,ch)-kv.S))).^-1;
  end
end

end