function inoutsig = adaptloop(inoutsig,fs,limit,minlvl,tau);
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
%M
%R  dau1996qmeI puschel1988pza

% Copyright (c) 1999 - 2004 Stephan Ewert. All rights reserved.

%   AUTHOR : Stephan Ewert, Morten L. Jepsen, Peter L. Soendergaard

% ------ Checking of input parameters and default parameters ---------

error(nargchk(2,5,nargin));
  
% Default parameters for tau measured in seconds.
if nargin<5
  tau=[0.005 0.050 0.129 0.253 0.500];
else
  if ~isnumeric(tau) || ~isvector(tau) || tau<=0
    error('%s: tau must be a vector with positive values.',upper(mfilename));
  end;
end;

if nargin<4
  minlvl =1e-5;
else
  if ~isnumeric(minlvl) || ~isscalar(minlvl) || minlvl<=0
    error('%s: minlvl must be a positive scalar.',upper(mfilename));
  end;
end;

if nargin<3
  limit = 10;
else
  if ~isnumeric(limit) || ~isscalar(limit) 
    error('%s: "limit" must be a scalar.',upper(mfilename));
  end;  
end;

% -------- Computation ------------------

inoutsig=comp_adaptloop(inoutsig,fs,limit,minlvl);
