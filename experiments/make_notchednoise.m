function output = make_notchednoise(fs,fc,dur,L,bw,delta)
%MAKE_NOTCHEDNOISE  Generates a notched-noise-type masker
%   Usage: output = make_notchednoise(fs,fc,dur,L,bw,delta);
%          output = make_notchednoise(fs,fc,dur,L,bw,[deltaL deltaR]);
% 
%   Generates a notched-noise-type masker as used in psychoacoustical
%   studies investigating the auditory filters' shape (original method
%   described in Patterson, 1974). The noise is composed of two noise bands
%   of width `bw` and a stopband centered at `fc` with a deviation from `fc`
%   given by `delta`.
%
%   `output = make_notchednoise(fs,fc,dur,L,bw,delta)` generates a
%   notched-noise masker with duration dur (in sec) and overall level `L`
%   (in dB SPL) with a sampling rate of `fs` Hz. The deviation from center
%   frequency `fc` is symmetric and is given by `delta` such that the
%   stopband is `[fc-delta*fc fc+delta*fc]`. The left and right noise bands
%   have a bandwidth of `bw*fc`, in Hz. If `delta=0` then a broadband noise
%   is returned.
%
%   `output = make_notchednoise(fs,fc,dur,L,bw,[deltaL deltaR])` generates a
%   notched-noise masker with an asymmetric configuration. `deltaL` and
%   `deltaR` denote the left and right deviations from `fc`, respectively.
%   In this case the stopband is `[fc-deltaL*fc fc+deltaR*fc]`.
% 
%   References: patterson1974auditory moore1990auditory

if nargin<6
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isnumeric(delta) || ~isempty(find(delta<0,1))
  error('%s: delta must be a positive scalar or array.',upper(mfilename));
end;

%% Generate broadband Gaussian noise
n = round(dur*fs);
if mod(n,2) ~= 0
    n = n+1;
end

noise = randn(1,n);
noise = noise./max(noise);% Normalize to 1

n_ramp = round(10E-3*fs);% Default = 10-ms Hanning on/off ramps

noise = rampsignal(noise,n_ramp);

%ramp = hann(n_ramp);
% Apply temporal windowing
%noise(1:n_ramp/2) = noise(1:n_ramp/2).*ramp(1:end/2)';% Onset
%noise(end-n_ramp/2+1:end) = noise(end-n_ramp/2+1:end).*ramp(end/2+1:end)';% Offset
% Zero padding to account for FIR delay (l = IR length)

l = 1024;% Length of filter impulse response
noiseZP = [zeros(1,l) noise zeros(1,l)];

% Set overall level in dB SPL
% noiseZP = setdbspl(noiseZP,L);
noiseZP = gaindb(noiseZP,L);

%% If delta = 0 then filter is not required
if isscalar(delta) && delta == 0
    output = noiseZP;
end
%% Multiband filter design (FIR filter)
if isscalar(delta) && delta > 0
%     Symmetric notch
    b1l = ((fc*(1-delta-bw))*2)/fs;% Low edge of left noise band
    if b1l < 0
        b1l = 0;
    end
    b1h = ((fc*(1-delta))*2)/fs;% High edge of left noise band
    b2l = ((fc*(1+delta))*2)/fs;% Low edge of right noise band
    b2h = ((fc*(1+delta+bw))*2)/fs;% High edge of right noise band
    if b2h > 1
        b2h = 1;
    end
    % Compute and analyze filter
    f = [0 b1l b1l b1h b1h b2l b2l b2h b2h 1];
    m = [0 0 1 1 0 0 1 1 0 0];
    b = firls(l,f,m);
    % Plot for verification (uncomment if needed)
    % [h,w] = freqz(b,1,1024);
    % figure
    % plot(f,m,w/pi,abs(h),'--'),grid on
    % xlabel('Normalized frequency'), ylabel('Magnitude'), title('Filter response')
    % Apply filter:
    output = filter(b,1,noiseZP);
elseif ~isscalar(delta) 
    if length(delta) > 2
        error('%s: Stopband is not correctly specified.',upper(mfilename));
    end
%     Asymmetric notch
    b1l = ((fc*(1-delta(1)-bw))*2)/fs;% Low edge of left noise band
    if b1l < 0
        b1l = 0;
    end
    b1h = ((fc*(1-delta(1)))*2)/fs;% High edge of left noise band
    b2l = ((fc*(1+delta(2)))*2)/fs;% Low edge of right noise band
    b2h = ((fc*(1+delta(2)+bw))*2)/fs;% High edge of right noise band
    if b2h > 1
        b2h = 1;
    end
    % Compute and analyze filter
    f = [0 b1l b1l b1h b1h b2l b2l b2h b2h 1];
    m = [0 0 1 1 0 0 1 1 0 0];
    b = firls(l,f,m);
    % Plot for verification (uncomment if needed)
    % [h,w] = freqz(b,1,1024);
    % figure
    % plot(f,m,w/pi,abs(h),'--'),grid on
    % xlabel('Normalized frequency'), ylabel('Magnitude'), title('Filter response')
    % Apply filter:
    output = filter(b,1,noiseZP);
end

%% Plot results (uncomment if needed)
% fft_noise1 = fft(noiseZP)./(n+2*l);
% fft_noise2 = fft(output)./(n+2*l);
% figure
% % Time domain
% subplot(2,1,1)
% plot(noiseZP), hold on
% plot(output,'--r'), hold off
% legend('Input','Output')
% xlabel('Samples'),ylabel('Amplitude'),title('Time domain')
% % Frequency domain
% subplot(2,1,2)
% plot(linspace(0,fs,n+2*l),20*log10(abs(fft_noise1))), hold on
% plot(linspace(0,fs,n+2*l),20*log10(abs(fft_noise2)),'--r'), hold off
% legend('Broadband noise','Filtered noise')
% xlabel('Frequency (Hz)'),ylabel('Squared modulus (dB SPL)'),title('Frequency domain')

% eof
