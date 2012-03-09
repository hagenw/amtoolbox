function outsig = ref_adaptloop_2(insig,fs,limit,minlvl,tau);
%ADAPTLOOP   Adaptation loops.
%   Usage: out = adaptloop(in,fs,limit,minlvl,tau);
%          out = adaptloop(in,fs,limit,minlvl);
%          out = adaptloop(in,fs,limit);
%
%   This is an implementation in the style of difference equations.
  
%   AUTHOR : Stephan Ewert, Morten L. Jepsen, Peter L. Soendergaard
%   Last changed on $Date: 2009-01-14 17:46:50 +0100 (ons, 14 jan 2009) $
%   Last change occured in $Rev: 1 $

  
% ------ Checking of input parameters and default parameters ---------

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

% Determine the signal length and width
siglen=size(insig,1);
nsigs=size(insig,2);

% Allocate output array.
outsig=zeros(siglen,nsigs);

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
insig = max(insig,minlvl);

% Determine steady-state levels. The values are repeated to fit the
% number of input signals.
state=minlvl.^(1./(2.^((1:nloops))));

s=zeros(siglen+1,nloops);
s(1,:)=state;


% Back up the value, because state is overwritten
stateinit=state;

if limit <=1 
  % No overshoot limitation

  for w=1:nsigs
    state=stateinit;
    
    for ii=1:siglen      
      % Compute the adaptation loops.
      %for jj=1:nloops
      %  insig(ii,w)=insig(ii,w)/state(jj);
      %  state(jj) = a1(jj)*state(jj) + b0(jj)*insig(ii,w);
      %end;    

      % One formulation of the equations.
      %s(ii+1,1) = a1(1)*s(ii,1) + b0(1)*insig(ii,w)/s(ii,1);
      %s(ii+1,2) = a1(2)*s(ii,2) + b0(2)*insig(ii,w)/(s(ii,1)*s(ii,2));
      %s(ii+1,3) = a1(3)*s(ii,3) + b0(3)*insig(ii,w)/(s(ii,1)*s(ii,2)*s(ii,3));
      %s(ii+1,4) = a1(4)*s(ii,4) + b0(4)*insig(ii,w)/(s(ii,1)*s(ii,2)*s(ii,3)*s(ii,4));
      %s(ii+1,5) = a1(5)*s(ii,5) + b0(5)*insig(ii,w)/(s(ii,1)*s(ii,2)*s(ii,3)*s(ii,4)*s(ii,5));
      %outsig(ii,w)=insig(ii,w)/(s(ii,1)*s(ii,2)*s(ii,3)*s(ii,4)*s(ii,5));
    
      
      
      s(ii+1,1) = a1(1)*s(ii,1) + (1-a1(1))*insig(ii,w)/s(ii,1);
      s(ii+1,2) = a1(2)*s(ii,2) + (1-a1(2))/(1-a1(1))*(s(ii+1,1)-a1(1)*s(ii,1))/s(ii,2);      
      s(ii+1,3) = a1(3)*s(ii,3) + b0(3)/b0(2)*(s(ii+1,2)-a1(2)*s(ii,2))/s(ii,3);
      s(ii+1,4) = a1(4)*s(ii,4) + b0(4)/b0(3)*(s(ii+1,3)-a1(3)*s(ii,3))/s(ii,4);
      s(ii+1,5) = a1(5)*s(ii,5) + b0(5)/b0(4)*(s(ii+1,4)-a1(4)*s(ii,4))/s(ii,5);

      outsig(ii,w) = (s(ii+1,5)-a1(5)*s(ii,5)) / b0(5);
            
    end;
  
  end;
  
else 
  
  % Overshoot Limitation.
  
  % Max. possible output value
  maxvalue = (1 - state.^2) * limit - 1;
  
  % Factor in formula to speed it up 
  factor = maxvalue * 2; 			
  
  % Exponential factor in output limiting function
  expfac = -2./maxvalue;
  offset = maxvalue - 1;

  for w=1:nsigs
    
    % Determine steady-state levels. The values are repeated to fit the
    % number of input signals.
    state=stateinit;
  
    for ii=1:siglen

      tmp1=insig(ii,w);
      
      for jj=1:nloops
        
        tmp1=tmp1/state(jj);
        
        if ( tmp1 > 1 )
          tmp1 = factor(jj)/(1+exp(expfac(jj)*(tmp1-1)))-offset(jj);
        end
        
        state(jj) = a1(jj)*state(jj) + b0(jj)*tmp1;
      
      end;
      
      % store the result.
      outsig(ii,w)=tmp1;    
      
    end;
  end
end;

% Scale to model units
outsig = (outsig-corr)*mult;



