function inoutsig = adaptloop(inoutsig,fs,varargin);
%ADAPTLOOP   Adaptation loops.
%   Usage: outsig = adaptloop(insig,fs,limit,minlvl,tau);
%          outsig = adaptloop(insig,fs,limit,minlvl);
%          outsig = adaptloop(insig,fs,limit);
%          outsig = adaptloop(insig,fs);
%
%   ADAPTLOOP(insig,fs,limit,minlvl,tau) applies non-linear adaptation to an
%   input signal insig sampled at a sampling frequency of fs Hz. limit is used
%   to limit the overshoot of the output, minlvl determines the lowest
%   audible threshhold of the signal and tau are time constants involved
%   in the adaptation loops. The number of adaptation loops is determined
%   by the length of tau.
%
%   ADAPTLOOP(insig,fs,limit,minlvl) does as above, but uses the values for
%   tau determined in Dau. 1996.
%
%   ADAPTLOOP(insig,fs,limit) does as above with a minimum threshhold minlvl
%   equal to 1e-5.
%
%   ADAPTLOOP(insig,fs) does as above with an overshoot limit of limit=10.
%
%   ADAPTLOOP takes the following flags at the end of the line of input
%   arguments:
%
%      dau - Choose the parameters as in the Dau 1996 and 1997 models. This
%           consists of 5 adaptation loops with an overshoot limiting of 10
%           and a minimum level of 1e-5. This is a correction in regard to
%           the published version of Dau 96, which did not use overshoot
%           limiting. The adaptation loops have an exponential spacing.This
%           flag is the default.
%
%      breebart - Choose the parameters as in the Breebart 2001 model. This
%           consists of 5 adaptation loops without overshoot limiting and a
%           minimum level of XXX. The adapation loops have a linear spacing.
%
%R  dau1996qmeI puschel1988pza breebaart2001binaural

% Copyright (c) 1999 - 2004 Stephan Ewert. All rights reserved.

%   AUTHOR : Stephan Ewert, Morten L. Jepsen, Peter L. Soendergaard

% ------ Checking of input parameters and default parameters ---------

if nargin<2
  error('Too few input parameters.');
end;

defnopos.keyvals.limit=10;
defnopos.keyvals.minlvl=1e-5;
defnopos.keyvals.tau=[0.005 0.050 0.129 0.253 0.500];

defnopos.groups.dau     ={'limit',10,'minlvl',1e-5, ...
                    'tau',[0.005 0.050 0.129 0.253 0.500]};
defnopos.groups.breebart={'limit',10,'minlvl',1e-5,...
                    'tau',linspace(0.005,0.5,5)};

[flags,keyvals,limit,minlvl,tau]  = ltfatarghelper({'limit','minlvl','tau'},defnopos,varargin);

if ~isnumeric(tau) || ~isvector(tau) || any(tau<=0)
  error('%s: tau must be a vector with positive values.',upper(mfilename));
end;

if ~isnumeric(minlvl) || ~isscalar(minlvl) || minlvl<=0
  error('%s: minlvl must be a positive scalar.',upper(mfilename));
end;

if ~isnumeric(limit) || ~isscalar(limit) 
  error('%s: "limit" must be a scalar.',upper(mfilename));
end;  

% -------- Computation ------------------

inoutsig=comp_adaptloop(inoutsig,fs,limit,minlvl,tau);
