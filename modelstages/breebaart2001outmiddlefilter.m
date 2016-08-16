function out = breebaart2001outmiddlefilter(in,fs)
%breebaart2001outmiddlefilter simulates the outer- and middle ear transfer
%function from Breebaart et. al. 2001
% 
%   Usage: out = breebaart2001outmiddlefilter(in,fs)
%
%   Input parameters:
%        in  : input acoustic signal.
%        fs     : sampling rate.
% 
%   `breebaart2001outmiddlefilter(in,fs)` filters the input signal *in* with the
%   transfer function of the outer and middle ear  sampled with a frequency
%   of *fs* Hz as described in Breebaart (2001).
%
%   References: breebaart2001a

%   AUTHOR: Martina Kreuzbichler


% coefficients
q = 2 - cos(2*pi*4000/fs) - sqrt((cos(2*pi*4000/fs)-2).^2-1);
r = 2 - cos(2*pi*1000/fs) - sqrt((cos(2*pi*1000/fs)-2).^2-1);

% initialize vectors
len=length(in)+2; 
% get the length + overhead
y=zeros(len,1);    
% define y and x
x=zeros(len,1);     
x(3:len)=in;     

% main loop
for n=3:len
    y(n)=(1-q)*r*x(n) - (1-q)*r*x(n-1) + (q+r)*y(n-1) - q*r*y(n-2);
end

% set the output
out=y(3:end);

% to see bode
% imp=[1; zeros(1024,1)];
% out= outmiddleartransfunct(imp,48000);
% spect(fft(out),48000);



