% This script compares the implentation of gammatone filter coefficiens
% supplied for the DRNL with the standard gammatone implementation.
%
% The DRNL filters seems to have must higher order, and they are not
% all-pole.


fc=1000;
fs=16000;
n=4;
bw=100;

[b_ref,a_ref]=coefGtDRNL(1000,100,4,16000);

[b,a]=gammatone(fc,fs,n,bw,'complex');
b=real(b);
a=real(a);

