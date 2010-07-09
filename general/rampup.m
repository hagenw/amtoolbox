function outsig=rampup(siglen)
%RAMPUP  Rising ramp function
%   Usage: outsig=rampup(siglen);
%
%   RAMPUP(siglen) will return a rising ramp function of length siglen. The
%   ramp is a sinsosoidal starting from zero and ending at one. The ramp
%   is centered such that the first element is always 0 and the last
%   element is not quite 1, such that the ramp fits with following ones.
%
%   See also: rampdown
  
outsig=sin((0:siglen-1)/siglen*pi/2).';
  