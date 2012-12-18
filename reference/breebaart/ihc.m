% ==========================================================
% ihc - simple inner haircell model 
%
% Written by Jeroen Breebaart - jeroen.breebaart@philips.com
%
% Usage: ihc(x,fs)
%
% x	: input signal, this must be a [n by 1] matrix
% y	: output signal (equal size of input signal)
% fs	: sample rate of input/output signal
%
% ==========================================================
% History
%
% 02/06/2001: 	matlab code replaced by linked c-code to
% 					improve speed.
%
% COMMENT BY PETER: The code in ihc.c produces /almost/ the same result
% as ihcenvelope(x,fs,'breebaart');


% This file does not include source code but is intended as
% help file only.



