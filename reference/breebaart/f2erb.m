function y=f2erb(x)

% ==========================================================
% f2erb - function to comute ERBs corresponding to CF
%
% Written by Jeroen Breebaart - jeroen.breebaart@philips.com
% (C) 2001 Philips Research Labs, Eindhoven.
%
% Usage: y=f2erb(x)
%
% x	: center frequency [Hz]
% y	: erbrate
%
% The centerfrequency (in Hz) of filter y is given by:
%
%		  x=(exp(0.11*y)-1)/.00437.
%
% Consequently, the erbrate of a frequency cf is given by
%
%	y=log(0.00437x+1)/0.11;
%


y=log(0.00437*x+1)/.11;


