function sibin = baumgartner2014binauralweighting(simon,kv,flags)
%BAUMGARTNER2014BINAURALWEIGHTING - Binaural combination of monaural similarity estimates
%   Usage:     sibin = baumgartner2014binauralweighting(simon)
%              sibin = baumgartner2014binauralweighting(simon,kv,flags)
%
%   Input parameters:
%     simon   : monaural similarity indices
%
%   Output parameters:
%     sibin   : monaural similarity indices
%
%   `baumgartner2014binauralweighting(...)` combines the monaural
%   similarity indices to binaural similartiy indices while accounting for
%   ipsilateral predominance.
%
%   `baumgartner2014binauralweighting` accepts the following optional parameters:
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

%% Binaural weighting, Eq.(6)

if size(simon,3) == 2
    binw = 1./(1+exp(-kv.lat/kv.bwcoef)); % weight of left ear signal with 0 <= binw <= 1
    sibin = binw * simon(:,:,1) + (1-binw) * simon(:,:,2);
end

end