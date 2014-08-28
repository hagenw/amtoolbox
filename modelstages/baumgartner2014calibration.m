function scalib = baumgartner2014calibration(s,kv,TolX)
%baumgartner2014calibration  Calibration of listener-specific sensitivity 
% thresholds to experimental performance
%   Usage: scalib = baumgartner2014calibration(s)
%
%   Input parameter:
%     s       : strucure containing subject's data. It must include the 
%               fields *Obj*, *pe_exp*, and *qe_exp*, representing the
%               listener's HRTF as SOFA object, the baseline local
%               polar RMS error, and the baseline quadrant error rate,
%               respectively.
%
%   Output parameter:
%     scalib  : strucure containing subject's data with calibrated u
%
%   `baumgartner2014calibration` returns listener data with
%   listener-specific sensitivity thresholds calibrated by joint
%   optimization of PE and QE to minimize mismatch between experimental
%   and predicted results.
%
%   `baumgartner2014calibration` accepts the following optional parameters:
%
%     kv        Key-value pairs according to `baumgartner2014`
%
%     TolX      Optimization tolerance of listener-specific sensitivity.
%               Default is 1e-3.
%
%   See also: baumgartner2014

% AUTHOR : Robert Baumgartner

if not(exist('kv','var'))
  definput.import={'baumgartner2014'};
  [~,kv]=ltfatarghelper({},definput,varargin);
end
kv.latseg = [-20,0,20];

if not(exist('TolX','var'))
  TolX = 0.001;
end
 
scalib = s;
for ss = 1:length(s)
  
  scalib(ss).S = fminsearch(@(S) evaldist(s(ss),S,kv),s(ss).S,...
    optimset('MaxIter',50,'TolX',TolX)...
    );
%   [~,scalib(ss).qe_pred,scalib(ss).pe_pred] = evaldist(s(ss),S,kv);
  disp([num2str(ss,'%2.0u') ' of ' num2str(length(s),'%2.0u') ' calibrated.'])

end


end


function [distmetric,qeM,peM] = evaldist(s,S,kv)

if S <= 0
  distmetric = Inf;
  return
end

%% LocaMo
qeM = zeros(length(s),1);
peM = zeros(length(s),1);
for ll = 1:length(s)

  for ii = 1:length(kv.latseg)

    s(ll).sphrtfs{ii} = 0;     % init
    s(ll).p{ii} = 0;        % init

    [s(ll).sphrtfs{ii},polang] = extractsp( kv.latseg(ii),s(ll).Obj );
    [s(ll).p{ii},respangs] = baumgartner2014(...
        s(ll).sphrtfs{ii},s(ll).sphrtfs{ii},s(ll).fs,...
        'S',S,'lat',kv.latseg(ii),'polsamp',polang,...
        'mrsmsp',kv.mrsmsp,'gamma',kv.gamma,'do',kv.do);

    [ qe(ii),pe(ii) ] = pmv2ppp( ...
        s(ll).p{ii} , polang , respangs , s(ll).target{ii});

    qeM(ll) = qeM(ll) + qe(ii)*s(ll).Ntargets{ii}/sum([s(ll).Ntargets{:}]);
    peM(ll) = peM(ll) + pe(ii)*s(ll).Ntargets{ii}/sum([s(ll).Ntargets{:}]);

  end

  dQE(ll) = s(ll).qe_exp - qeM(ll);
  dPE(ll) = s(ll).pe_exp - peM(ll);

end

[qe_chance,pe_chance] = pmv2ppp(ones(49,44));
distmetric =  (dQE/qe_chance).^2 + (dPE/pe_chance).^2; % Joint distance metric of QE and PE (standardized scatter)

end