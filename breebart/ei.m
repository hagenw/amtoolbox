function z = ei( l, r , fs , tau , IID )
%EI  XXX
%   Usage: y = ei( l , r, fs, tau, IID)
%
%   Input parameters:
%        l,r	    : input signals, these must be a [n by 1] matrix
%        fs        : sampling rate of input signals
%        tau       : characteristic delay in seconds (positive: left is leading)
%        IID       : characteristic IID in dB (positive: left is louder)
%
%   Output parameters:
%        y	   : EI-type cell output as a function of time
%
%

% Written by Jeroen Breebaart - jeroen.breebaart@philips.com
% (C) 2003 Philips Research Labs, Eindhoven.


% parameters:
tc          = 30e-3;            % Temporal smoothing constant
a           = 0.1;              % non-linear I/O parameter 'a' 
b           = 0.00002;          % non-linear I/O parameter 'b'
ptau        = 2.2e-3;           % time constant for p(tau) function

% apply characteristic delay:
n = round( abs(tau) * fs );

if tau > 0,
    l = [zeros(n,1) ; l(1:end-n)];
else
    r = [zeros(n,1) ; r(1:end-n)];
end

% Old code, harder to read.
%c1 = 10^( +IID/40 );
%c2 = 10^( -IID/40 );
%l = l * c1;
%r = r * c2;

% apply characteristic IID:
l=gaindb(l, IID/2);
r=gaindb(l,-IID/2);

% compute instanteneous EI output:
x = (l - r).^2;

% temporal smoothing:
A=[1 -exp(-1/(fs*tc))];
B=[1-exp(-1/(fs*tc)) ];
y= filtfilt(B,A,x);% / ( (1-exp(-1/(fs*tc)))/2 );

% compressive I/O: Scale signal by 200. This approximately
% results in JNDs of 1 in the output
z = exp(-tau/ptau) * a * log( b * y + 1);


