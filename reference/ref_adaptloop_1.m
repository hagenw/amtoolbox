function bigout = ref_adaptloop_1(bigin,fs,limit,minlvl,tau);
%ADAPTLOOP   Adaptation loops.
%   Usage: out = adaptloop(in,fs,limit,minlvl,tau);
%          out = adaptloop(in,fs,limit,minlvl);
%          out = adaptloop(in,fs,limit);
%
%   This implementation of adaptloop turned out to be slow.
  
% Copyright (c) 1999 - 2004 Stephan Ewert. All rights reserved.

%   AUTHOR : Stephan Ewert, Morten L. Jepsen, Peter L. Søndergaard
%   Last changed on $Date: 2009-01-14 17:46:50 +0100 (ons, 14 jan 2009) $
%   Last change occured in $Rev: 1 $

% ------ Checking of input parameters and default parameters ---------

error(nargchk(2,5,nargin));
  
% Default parameters for tau measured in seconds.
if nargin<5
  tau=[0.005 0.050 0.129 0.253 0.500];
else
  if ~isnumeric(tau) || ~isvector(tau) || any(tau<=0)
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
  if ~isnumeric(limit) || ~isscalar(limit) || limit<=0
    error('%s: "limit" must be a positive scalar.',upper(mfilename));
  end;  
end;

% -------- Computation ------------------

% Determine the signal length and width
len=size(bigin,1);
W=size(bigin,2);

% Allocate output array.
bigout=zeros(len,W);

% Determine the number of adaptation loops
nloops=length(tau);

% Calculate filter coefficients for the loops.

% b0 from RC-lowpass recursion relation y(n)=b0*x(n)+a1*y(n-1)
% a1 coefficient of the upper IIR-filter
b0=1./(tau*fs);
a1=exp(-b0);
b0=1-a1;

% To get a range from 0 to 100 model units
corr = minlvl^(1/(2^nloops));
mult = 100/(1-corr); 

% Apply minimum level to the input
bigin = max(bigin,minlvl);

% Determine steady-state levels. The values are repeated to fit the
% number of input signals.
tmp2=repmat(minlvl.^(1./(2.^((1:nloops).'))),W);



if limit <=1 
  % No overshoot limitation

  for ii=1:len
    tmp1=bigin(ii,:);
    
    % Compute the adaptation loops.
    for jj=1:nloops
      tmp1=tmp1./tmp2(jj,:);
      tmp2(jj,:) = a1(jj)*tmp2(jj,:) + b0(jj)*tmp1;         
    end;    
    
    % store the result.
    bigout(ii,:)=tmp1;
  end;  
else 
  
  % Overshoot Limitation.
  
  % Max. possible output value
  maxvalue = (1 - tmp2.^2) * limit - 1;
  
  % Factor in formula to speed it up 
  factor = maxvalue * 2; 			
  
  % Exponential factor in output limiting function
  expfac = -2./maxvalue;
  offset = maxvalue - 1;
  
  for ii=1:len
      tmp1=bigin(ii,:);
      
      for jj=1:nloops
        
        tmp1=tmp1./tmp2(jj,:);

        for w=1:W          
          if ( tmp1(w) > 1 )
            tmp1(w) = factor(jj,w)/(1+exp(expfac(jj,w)*(tmp1(w)-1)))-offset(jj,w);
          end
          
        end;
        tmp2(jj,:) = a1(jj)*tmp2(jj,:) + b0(jj)*tmp1;
        
      end;
      
      % store the result.
      bigout(ii,:)=tmp1;    
      
  end;
end


% Scale to model units
bigout = (bigout-corr)*mult;



