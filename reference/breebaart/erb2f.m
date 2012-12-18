function y=erb2f(x)


% ==========================================================
% erb2f - function to compute frequencies corresponding to ERBs
%
% Written by Jeroen Breebaart - jeroen.breebaart@philips.com
% (C) 2001 Philips Research Labs, Eindhoven.
%
% Usage: y=erb2f(x)
%
% x	: erbrate [ERB]
% y	: center frequency [Hz]
%
% The centerfrequency (in Hz) of filter x is given by:
%
%		  cf=(exp(0.11*x)-1)/.00437.
%


y=(exp(0.11*x)-1)/.00437;



