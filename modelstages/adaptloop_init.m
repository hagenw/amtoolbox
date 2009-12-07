function s = adaptloop_init(nsigs,fs,limit,minlvl,tau);
%ADAPTLOOP   Adaptation loops.
%   Usage: s = adaptloop(nsigs,fs,limit,minlvl,tau);
%          s = adaptloop(nsigs,fs,limit,minlvl);
%          s = adaptloop(nsigs,fs,limit);
%          s = adaptloop(nsigs,fs);
%
%   ADAPTLOOP(nsigs,fs,limit,minlvl,tau) applies non-linear adaptation to an
%   input signal insig sampled at a sampling frequency of fs Hz. limit is used
%   to limit the overshoot of the output, minlvl determines the lowest
%   audible threshhold of the signal and tau are time constants involved
%   in the adaptation loops. The number of adaptation loops is determined
%   by the length of tau.
%
%   ADAPTLOOP(nsigs,fs,limit,minlvl) does as above, but uses the values for
%   tau determined in Dau. 1996.
%
%   ADAPTLOOP(nsigs,fs,limit) does as above with a minimum threshhold minlvl
%   equal to 1e-5.
%
%   ADAPTLOOP(nsigs,fs) does as above with an overshoot limit of limit=10.
%M
%R  dau1996qmeI puschel1988pza

%   AUTHOR: Peter L. Soendergaard

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

if nargin<4
  tau=[0.005 0.050 0.129 0.253 0.500];
end;

% -------- Computation ------------------

s.loops=length(tau);

% Calculate filter coefficients for the loops.

% b0 from RC-lowpass recursion relation y(n)=b0*x(n)+a1*y(n-1)
% a1 coefficient of the upper IIR-filter
a1=exp(-1./(tau*fs));

% To get a range from 0 to 100 model units
s.corr = minlvl^(1/32);
s.mult = 100/(1-s.corr); 

% Determine steady-state levels. The values are repeated to fit the
% number of input signals.
state=repmat(minlvl.^(1./(2.^((1:s.loops).'))),nsigs);
  

% Overshoot Limitation.
  
% Max. possible output value
maxvalue = (1 - state.^2) * limit - 1;

% Factor in formula to speed it up 
factor = maxvalue * 2; 			

% Exponential factor in output limiting function
expfac = -2./maxvalue;
offset = maxvalue - 1;

% Load it all up

s.nsigs  = nsigs;
s.state  = state;
s.a1     = a1;
s.factor = factor;
s.expfac = expfac;
s.offset = offset;
s.minlvl = minlvl;
s.limit  = limit;