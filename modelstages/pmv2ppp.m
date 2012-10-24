function [ varargout ] = pmv2ppp( p,varargin )
%PMV2PPP Retrieve psychoacoustic performance parameters (PPPs) from
% predicted probability mass vectors (PMVs) of polar response angles.
%
%   Usage:  [ qe,pe,pb ] = pmv2ppp( p,tang,rang )
%           [ qe,pe,pb ] = pmv2ppp( p,tang,rang,exptang )
%
%   `pmv2ppp(...)` retrieves commonly used PPPs for sagittal-plane (SP)
%   localization like quadrant error (QEs), local polar RMS error (PE), 
%   and local polar bias (PB) from response PMVs predicted by a localization 
%   model. PPPs are retreived either for a specific polar target angle or
%   as an average across all available target angles. The latter is the
%   default.
%
%   Input arguments:
%     p          : prediction matrix (response PMVs)
%     tang       : possible polar target angles. As default, ARI's MSP 
%                  polar angles in the median SP is used.
%     rang       : polar angles of possible response angles.
%                  As default regular 5°-sampling is used (-30:5:210).
%
%   `pmv2ppp` needs the following optional parameter in order to retrieve
%   the PPPs for a specific (set of) target angles:
%
%     'exptang', exptang   experimental polar target angles
%
%  Output arguments:
%     qe         : quadrant error rate
%     pe         : local polar RMS error in degrees
%     pb         : polar bias in degrees
%
%   `pmv2ppp` accepts the following flag:
%
%     'print'      Display the outcomes.
%
%  References: baumgartner2013assessment baumgartner2012modelling

% AUTHOR : Robert Baumgartner

definput.flags.print = {'noprint','print'};
definput.keyvals.rang=-30:5:210;
definput.keyvals.tang=[-30:5:70,80,100,110:5:210];
definput.keyvals.exptang=[];
[flags,kv]=ltfatarghelper({'tang','rang','exptang'},definput,varargin);

    
p = p./repmat(sum(p),length(kv.rang),1);  % ensure probability mass vectors
tang = kv.tang(:);
rang = kv.rang(:);
nt = length(tang);

qet = zeros(nt,1); % QE for each target angle
pet = zeros(nt,1); % PE for each target angle
pbt = zeros(nt,1); % PB for each target angle
for ii = 1:nt % for all target positions
    d = abs(wrapTo180(tang(ii)-rang)); % unwraped angular distance between tang & rang
    qet(ii) = sum( p(d>=90,ii) );
    pc = p(d<90,ii);                   % pmv for conditional probability excluding QEs
    pc = pc/sum(pc);                   % normalization to sum=1
    pet(ii) = sqrt( sum( pc .* (d(d<90)).^2 )); % RMS of expected difference
    pbt(ii) = sum( pc .* rang(d<90) ) - tang(ii); % expectancy value of rang - tang
end

if ~isempty(kv.exptang)
    extang = [-90; tang(:); 270]; % extended tang for targets outside
    qetb = (qet(1)+qet(end))/2;  % boundaries for extang
    etb = (pet(1)+pet(end))/2;
    pbtb = (pbt(1)+pbt(end))/2;
    qet = interp1(extang,[qetb; qet(:); qetb],kv.exptang);
    pet = interp1(extang,[etb; pet(:); etb],kv.exptang);
    pbt = interp1(extang,[pbtb; pbt(:); pbtb],kv.exptang);
end

qe = mean(qet)*100;
pe = mean(pet);
pb = mean(pbt);

varargout{1} = qe;
varargout{2} = pe;
varargout{3} = pb;


if flags.do_print
    fprintf('Quadrant errors (%%) \t\t %4.1f \n',qe)
    if nargout > 1
      fprintf('Local polar RMS error (deg) \t %4.1f \n',pe)
    end
    if nargout > 2
      fprintf('Local polar bias (deg) \t\t %4.1f \n',pb)
    end
end

end