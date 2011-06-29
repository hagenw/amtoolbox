% ==========================================================
% Fadapt - adaptation loops according to Pueschel. The 
% input signal is fed through 5 consequtive adaptation loops
% with time constants linearly spaced between 5 and 500 ms.
%
% Written by Jeroen Breebaart - jeroen.breebaart@philips.com
%
% Usage: y=fadapt(x,fs)
%
% x	: input signal, this must be a [n by 1] matrix
% y	: output signal (equal size of input signal)
% fs	: sample rate of input/output signal
%
% Convention: rms=1 ==> 0 dB SPL
%
% ==========================================================
% History
%
% 02/06/2001: 	matlab code replaced by linked c-code to
% 					improve speed.
%
% 08/01/2002:   Bug fix: corrected initialization of states
%

% This file does not include source code but is intended as
% help file only.

