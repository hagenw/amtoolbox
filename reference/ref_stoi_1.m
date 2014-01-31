function d = ref_stoi_1(x, y, fs_signal)
%   d = stoi(x, y, fs_signal) returns the output of the short-time
%   objective intelligibility (STOI) measure described in [1, 2], where x 
%   and y denote the clean and processed speech, respectively, with sample
%   rate fs_signal in Hz. The output d is expected to have a monotonic 
%   relation with the subjective speech-intelligibility, where a higher d 
%   denotes better intelligible speech. See [1, 2] for more details.
%
%   This is a cleaned up version of ref_stoi, producing the exact same
%   results.

  
if length(x)~=length(y)
    error('x and y should have the same length');
end

% initialization
% clean speech column vector
x = x(:);

% processed speech column vector
y = y(:);                             

% sample rate of proposed intelligibility measure
fs = 10000;

% window support
N_frame	= 256;                          

% FFT size
K = 512;

% Number of 1/3 octave bands
J = 15; 

% Center frequency of first 1/3 octave band in Hz.
mn = 150;          

% Get 1/3 octave band matrix
H = thirdoct(fs, K, J, mn);

% Number of frames for intermediate intelligibility measure (Length analysis window)
N = 30;           

% lower SDR-bound
Beta = -15;

% speech dynamic range
dyn_range = 40;

% resample signals if other samplerate is used than fs
if fs_signal ~= fs
    x = resample(x, fs, fs_signal);
    y = resample(y, fs, fs_signal);
end

% remove silent frames
[x y] = removeSilentFrames(x, y, dyn_range, N_frame, N_frame/2);

% apply 1/3 octave band TF-decomposition
g=fftshift(hanning(N_frame));

x_hat = dgtreal(x,g,N_frame/2,K); 	% apply short-time DFT to clean speech
y_hat = dgtreal(y,g,N_frame/2,K); 	% apply short-time DFT to processed speech

x_hat = x_hat(:, 2:end-2); % Remove border columns
y_hat = y_hat(:, 2:end-2);        	

% Collect the frequency bands into the 1/3 octave band TF-representation
% This is done through multiplication with a sparse matrix consisting of
% 0 and 1's
X = sqrt(H*abs(x_hat).^2);
Y = sqrt(H*abs(y_hat).^2);

% loop al segments of length N and obtain intermediate intelligibility measure for all TF-regions

% init memory for intermediate intelligibility measure
d_interm  	= zeros(J, length(N:size(X, 2)));

% constant for clipping procedure
c           = 10^(-Beta/20);                          

for m = N:size(X, 2)
  % region with length N of clean TF-units for all j
  X_seg = X(:, (m-N+1):m);
  
  % region with length N of processed TF-units for all j
  Y_seg = Y(:, (m-N+1):m);
  
  % obtain scale factor for normalizing processed TF-region for all j
  alpha   = sqrt(sum(X_seg.^2, 2)./sum(Y_seg.^2, 2));

  % obtain \alpha*Y_j(n) from Eq.(2) [1]
  aY_seg 	= Y_seg.*repmat(alpha, [1 N]);
  for j = 1:J
    % apply clipping from Eq.(3)   	
    Y_prime             = min(aY_seg(j, :), X_seg(j, :)+X_seg(j, :)*c);

    % obtain correlation coeffecient from Eq.(4) [1]
    d_interm(j, m-N+1)  = taa_corr(X_seg(j, :).', Y_prime(:));
  end
end
        
% combine all intermediate intelligibility measures as in Eq.(4) [1]
d = mean(d_interm(:));                              

%%
function  [A cf] = thirdoct(fs, N_fft, numBands, mn)
%   [A CF] = THIRDOCT(FS, N_FFT, NUMBANDS, MN) returns 1/3 octave band matrix
%   inputs:
%       FS:         samplerate 
%       N_FFT:      FFT size
%       NUMBANDS:   number of bands
%       MN:         center frequency of first 1/3 octave band
%   outputs:
%       A:          octave band matrix
%       CF:         center frequencies

f  = linspace(0, fs, N_fft+1);
f  = f(1:(N_fft/2+1));
k  = 0:(numBands-1); 
cf = 2.^(k/3)*mn;
fl = sqrt((2.^(k/3)*mn).*2.^((k-1)/3)*mn);
fr = sqrt((2.^(k/3)*mn).*2.^((k+1)/3)*mn);
A  = sparse(numBands, length(f));

for i = 1:(length(cf))
  [a b] = min((f-fl(i)).^2);
  fl(i) = f(b);
  fl_ii = b;
  
  [a b] = min((f-fr(i)).^2);
  fr(i) = f(b);
  fr_ii = b;
  A(i,fl_ii:(fr_ii-1))	= 1;
end

rnk         = sum(A, 2);
numBands = find((rnk(2:end)>=rnk(1:(end-1))) & (rnk(2:end)~=0)~=0, 1, 'last' )+1;
A           = A(1:numBands, :);
cf          = cf(1:numBands);

%%
function [x_sil y_sil] = removeSilentFrames(x, y, range, N, K)
%   [X_SIL Y_SIL] = REMOVESILENTFRAMES(X, Y, RANGE, N, K) X and Y
%   are segmented with frame-length N and overlap K, where the maximum energy
%   of all frames of X is determined, say X_MAX. X_SIL and Y_SIL are the
%   reconstructed signals, excluding the frames, where the energy of a frame
%   of X is smaller than X_MAX-RANGE

x       = x(:);
y       = y(:);

frames  = 1:K:(length(x)-N);
w       = hanning(N);
msk     = zeros(size(frames));

for j = 1:length(frames)
    jj      = frames(j):(frames(j)+N-1);
    msk(j) 	= 20*log10(norm(x(jj).*w)./sqrt(N));
end

msk     = (msk-max(msk)+range)>0;
count   = 1;

x_sil   = zeros(size(x));
y_sil   = zeros(size(y));

for j = 1:length(frames)
    if msk(j)
        jj_i            = frames(j):(frames(j)+N-1);
        jj_o            = frames(count):(frames(count)+N-1);
        x_sil(jj_o)     = x_sil(jj_o) + x(jj_i).*w;
        y_sil(jj_o)  	= y_sil(jj_o) + y(jj_i).*w;
        count           = count+1;
    end
end

x_sil = x_sil(1:jj_o(end));
y_sil = y_sil(1:jj_o(end));


%%
function rho = taa_corr(x, y)
%   RHO = TAA_CORR(X, Y) Returns correlation coeffecient between column
%   vectors x and y. Gives same results as 'corr' from statistics toolbox.
xn    	= x-mean(x);
yn   	= y-mean(y);
rho   	= dot(xn/norm(xn),yn/norm(yn));