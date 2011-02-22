function [bout,aout] = convolveba(bin,ain,nrep);
%
%   CONVOLVEBA(bin,ain,nrep) returns the coefficients of the filter given by
%   bin and ain repeater nrep times.

 bout=1;
 aout=1;
  
 for ii=1:nrep
   bout=conv(bin,bout);
   aout=conv(ain,aout);    
 end;