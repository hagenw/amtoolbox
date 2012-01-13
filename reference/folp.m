% folp.m - digital first-order lowpass-filter transfer function
%          (bilinear transform)
%
% Usage: [b,a] = folp(f0,fs)
%
% f0 = cutoff frequency of the lowpass filter in Hz
% fs = sampling rate in Hz
%
% b = [b0,b1, ... ,bN] = numerator coefficients
% a = [1,a1, ... ,aN] = denominator coefficients
%
% This function computes [b,a] = butter(1,f0/(fs/2));
function [b,a] = folp(f0,fs)

w0 = 2*pi*f0/fs;	% frequency conversion
W0 = tan(w0/2);		% prewarping

b = [W0, W0]/(1 + W0);
a = [1,(W0 - 1)/(1 + W0)];

%OLDFORMAT
