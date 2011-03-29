function y = dbspl(insig,varargin)
%DBSPL RMS value of signal (in dB)
%   Usage: y = dbspl(insig);
%          y = dbspl(insig,'ac');
%
%   DBSPL(insig) computes the SPL (sound pressure level) of the input signal
%   measured in dB, using the convention that a pure tone at 100 dB SPL has
%   an RMS value of 1.
%  
%   DBSPL takes the following flags at the end of the line of input
%   parameters:
%
%-     'ac'     : Consider only the AC component of the signal (i.e. the mean is
%                 removed).
%
%-     'dim',d  : Work along specified dimension. The default value of []
%                 means to work along the first non-singleton one.
%
%   See also: setdbspl
%
%R  moore2003introduction

%   AUTHOR : Hagen Wierstorf
  
% ------ Computation --------------------------

% The level of a signal in dB SPL is given by the following formula:
% level = 20*log10(p/p_0)
% To get to the standard used in the toolbox.
y = 20*log10( rms(insig,varargin{:}) )+100;
