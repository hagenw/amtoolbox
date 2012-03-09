% ==========================================================
% gtfbank - 3th-order gammatone filterbank with phase and
% delay-aligned filters
%
% Written by Jeroen Breebaart - jeroen.breebaart@philips.com
% (C) 2003 Philips Research Labs, Eindhoven.
%
% Usage: y=gtfbank(x,fs,numfilt,spacing)
%
% x	        : input signal, this must be a [n by 1] matrix
% fs        : sampling rate of input signals
% numfilt   : total number of filters
% spacing   : filter spacing expressed in ERBs
%
% y	        : complex output signal (n x numfilt matrix)
%
%
% The output consists of a set of complex bandpass filtered signals
% according to the ERB scale (Glasberg & Moore).
% The centerfrequency (in Hz) of filter m is given by:
%
%		  cf=(exp(0.11*m)-1)/.00437.
%
% If the spacing parameter equals 1, the center frequencies
% are as given above. For a spacing parameter of 2,
% the center frequencies are given by ERB rates 1, 3, 5, etc.
%
%
% If the center frequency of some filters is beyond the nyquist
% frequency defined by the sampling rate, an error message is issued.
% as a rule of thumb, the full bandwidth can be covered
% at a sampling rate of 44.1 kHz if numfilt*spacing equals 40.
%
% The ERB bandwidth of the filter equals the ERB according
% to Glasberg & Moore:
%
%		  bw=spacing*24.7*(cf*0.00437+1)=24.7*exp(0.11*m)
%
% The output is given as complex numbers. For waveform output,
% just take the real part of the complex output. The ABS of
% each filter output results in the Hilbert envelope. The ANGLE
% of each output results in the instananeous phase.
%
% ==========================================================
% History
%
% 20/08/2001: 	finished c-code for optimized processing
%
% 07/01/2003:   Included variable spacing, number of filters,
%               complex output.
%               Delay is also removed.
% 05/09/2003:	Fixed a bug in bandwidth parameter computation.
%		Cos/sin operations replaced by tables.



% This file does not include source code but is intended as
% help file only.



