function outsig=rampdown(siglen)
%RAMPDOWN  Falling ramp function
%   Usage: outsig=rampdown(siglen);
%
%   RAMPDOWN(siglen) will return a falling ramp function of length siglen. The
%   ramp is a sinsosoidal starting from one and ending at zero. The ramp
%   is centered such that the first element is always one and the last
%   element is not quite zero, such that the ramp fits with following zeros.
%
%   See also: rampup
   
outsig=cos((0:siglen-1)/siglen*pi/2).';
  