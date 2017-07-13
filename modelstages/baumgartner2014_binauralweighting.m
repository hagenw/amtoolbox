function sibin = baumgartner2014_binauralweighting(simon,varargin)
%baumgartner2014_binauralweighting - Binaural combination of monaural similarity estimates
%   Usage:     sibin = baumgartner2014_binauralweighting(simon)
%
%   Input parameters:
%     simon   : monaural similarity indices
%
%   Output parameters:
%     sibin   : monaural similarity indices
%
%   `baumgartner2014_binauralweighting(...)` combines the monaural
%   similarity indices to binaural similartiy indices while accounting for
%   ipsilateral predominance.
%
%   `baumgartner2014_binauralweighting` accepts the following optional parameters:
%
%     'bwcoef',bwc   Set the binaural weighting coefficient *bwc*.
%                    Default value is 13 degrees.
%
%     'lat',lat      Set the apparent lateral angle of the target sound to
%                    *lat*. Default value is 0 degree (median SP).
%
%   References: baumgartner2014modeling

% AUTHOR: Robert Baumgartner

definput.import={'baumgartner2014'};
[flags,kv]=ltfatarghelper({},definput,varargin);

%% Binaural weighting, Eq.(6)

if size(simon,3) == 2
    binw = 1./(1+exp(-kv.lat/kv.bwcoef)); % weight of left ear signal with 0 <= binw <= 1
    sibin = binw * simon(:,:,1) + (1-binw) * simon(:,:,2);
end

end