% example.m : example script to call binaural model.

% input signal parameters:
fs      = 44100;            % Sampling rate
T       = 0.3;              % Duration
Spl1    = 75;               % SPL of input signal 1
Spl2    = 75;               % SPL of input signal 2
rho     = 0.0;              % normalized correlation of signals

% model parameters:
ERBspce = 1.0;              % ERB spacing of filterbank
ERBnum  = 40;               % Number of ERB filters


% generate signals:
t       = [0:1/fs:T];
n1      = randn( length(t) , 1);
n2      = randn( length(t) , 1);
x1      = 10^(Spl1/20) * ( n1*sqrt((1+rho)/2) + n2*sqrt((1-rho)/2) );
x2      = 10^(Spl2/20) * ( n1*sqrt((1+rho)/2) - n2*sqrt((1-rho)/2) );

%x1      = sinusoid(1000*T, 1000, Spl1, fs, 0);
%x2      = sinusoid(1000*T, 1000, Spl1, fs, 0);

% feed signals through gammatone filterbank:
[b,a]   = gammatone(erbspacebw(0,erbtofreq(ERBnum),ERBspce),fs,'complex')
x1g = 2*real(filterbankz(b,a,x1));
x2g = 2*real(filterbankz(b,a,x2));
%x1g     = real( gtfbank(x1, fs, ERBnum, ERBspce) );
%x2g     = real( gtfbank(x2, fs, ERBnum, ERBspce) );

% feed signals through inner haircell model:

for idx = 1 : ERBnum,
  x1h_ref(:,idx) = ihc( x1g(:,idx) , fs);
  x2h_ref(:,idx) = ihc( x2g(:,idx) , fs);
end

x1h = ihcenvelope(x1g,fs,'breebaart');
x2h = ihcenvelope(x2g,fs,'breebaart');

% feed signals through adaptation loops:
for idx=1:ERBnum
    x1a_ref(:,idx) = fadapt( x1h_ref(:,idx) , fs);
    x2a_ref(:,idx) = fadapt( x2h_ref(:,idx) , fs);
end

% feed signals through adaptation loops:
for idx=1:ERBnum
    x1a(:,idx) = fadapt( x1h(:,idx) , fs);
    x2a(:,idx) = fadapt( x2h(:,idx) , fs);
end

% compute EI representation for a single cell as a function of time:
ei_single = ei( x1a(:,10) , x2a(:,10) , fs , 0 , 0 );
ei_single_ref = ei( x1a_ref(:,10) , x2a_ref(:,10) , fs , 0 , 0 );



%OLDFORMAT
