function scalib = baumgartner2014calibration(s,kv)
%baumgartner2014calibration  Calibration of listener-specific sensitivity 
% thresholds to experimental performance
%   Usage: scalib = baumgartner2014calibration(s)
%          scalib = baumgartner2014calibration(s,latseg,dlat)
%
%   Input parameter:
%     s       : strucure containing subject's data
%
%   Output parameter:
%     scalib  : strucure containing subject's data with calibrated u
%
%   `baumgartner2014calibration` returns listener data with
%   listener-specific sensitivity thresholds calibrated by joint
%   optimization of PE and QE to minimize mismatch between experimental
%   and predicted results.

% AUTHOR : Robert Baumgartner

kv.latseg = 0;

scalib = s;
for ss = 1:length(s)
  
  scalib(ss).S = fminsearch(@(S) evaldist(s(ss),S,kv),s(ss).S,...
    optimset('MaxIter',50,'TolX',0.001)...
    );
  fprintf('%2.0u of %2.0u calibrated.\n',ss,length(s))

end


end


function distmetric = evaldist(s,S,kv)

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
        'mrsmsp',kv.mrsmsp,'gamma',kv.gamma,'do',kv.do,'print');

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