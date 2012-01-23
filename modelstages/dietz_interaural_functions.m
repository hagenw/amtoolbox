function [outp] = dietz_interaural_functions(s1, s2, tau, fc, signal_level_dB_SPL, compr, coh_cycles, fs)
% Calculates interaural parameters of two complex-valued
% Hilbert-transformed signals, supplied by the GFB Analyzer()
%
%   Input parameters
%     s1, s2 : input signals
%     tau    : lowpass filter parameter, such that a = exp(-1./(fs*tau)), 
%              lowpass = filter([1-a], [1, -a], x)
%  
%              Note: If tau is a scalar, lowpass is done on every channel
%              with the value of tau. If the filter parameter tau is a
%              vector, it has to have as many elements as the number of 
%              channels of s1 and s2; each value tau(i) is applied to
%              the signals in s1(i,:) and s2(i,:).
%     fc : central frequencies
%     fs     : sampling frequencies
%
%   Output parameters:
%     itf : transfer function
%     itf_equal : transfer function without amplitude
%     ipd : phase difference in rad
%     ipd_lp : based on lowpass-filtered itf, phase difference in rad
%     ild : level difference in dB
%     itd, itd_C, itd_lp, itd_C_lp - time difference based on instantaneous
%                  and central frequencies, with and without low-passed itf
%     f_inst_1 : instantaneous frequencies in the channels of the filtered s1
%     f_inst_2 : instantaneous frequencies in the channels of the filtered s2

a = exp( -1./(fs*tau) );

outp.itf = s2 .* conj(s1);

outp.ipd_lp = zeros(size(s1));
f_inst_1 = zeros(size(s1));
f_inst_2 = zeros(size(s1));

% interaural phase difference
outp.ipd = angle(outp.itf);

% interaural coherence
outp.ic = new_ic(outp.itf, coh_cycles./fc, fs);

% interaural level difference for the case of tau - scalar

s1_low = lowpass(abs(s1),a);
s2_low = lowpass(abs(s2),a);    % take envelope at higher frequencies

s1_low( find(abs(s1_low) < eps) ) = eps;    % avoid log(0)
s2_low( find(abs(s2_low) < eps) ) = eps;    % avoid division by zero

outp.ild = 20*log10(s1_low./s2_low)./compr;
if length(a) == 1
  a(1:length(fc)) = a;
  tau(1:length(fc)) = tau;
end
for k = 1:length(fc);
  outp.ipd_lp(:,k) = angle(lowpass(outp.itf(:,k),a(k)))';
end

% interaural time difference, based on central and instantaneous frequencies
for k = 1:length(fc)
  f_inst_1(:,k) = calc_f_inst(s1(:,k),fs,tau(k),0);
  f_inst_2(:,k) = calc_f_inst(s2(:,k),fs,tau(k),0);
  outp.itd_C(:,k) = 1/(2*pi)*outp.ipd(:,k)/fc(k);
  outp.itd_C_lp(:,k) = 1/(2*pi)*outp.ipd_lp(:,k)/fc(k);
end
outp.f_inst = max(eps,0.5*(f_inst_1 + f_inst_2));

% to avoid division by zero

% based on instantaneous frequencies
outp.itd = 1/(2*pi)*outp.ipd./outp.f_inst;    
outp.itd_lp = 1/(2*pi)*outp.ipd_lp./outp.f_inst;

%level = signal_level_dB_SPL*compr + 20*log10(abs(s1)+abs(s2));



% weighting of channels for cumulative ixd determination
% sqrt(2) is due to half-wave rectification (included 28th Sep 07)
outp.rms = signal_level_dB_SPL*compr + 20*log10(sqrt(2)*min(rms(abs(s1)),rms(abs(s2))));
outp.rms = max(outp.rms,0); % avoid negative weights


end;

%% lowpass
function y = lowpass(x, a)
% y = lowpass(x, a)
% This is a simple low-pass filter y(n) = (1-a)*x(n) - a*y(n-1)
% Meaning of parameter a:
%   a - damping coefficient (0 - no filtering, 1 - flat output)
%   tau = 1/(2*pi*f_c)      where f_c is the cutoff frequency of the filter
%   a = exp(-1/(f_s*tau))   where fs - sampling frequency
%
% Input
%   x - input signal
%   1] The signal x may be either row or column vector, the output y is ALWAYS a column vector
%   2] If x is a matrix, time is considered along the LONGER dimension of x
%      Then, the output y is always a column vector or a matrix with time
%      along the rows
%
% Output
%  y - filtered signal
%
% Example see ...\examples\example_lowpass_tester.m

[rows, columns] = size(x);
if rows < columns
  x = x.';
  y = zeros(columns,rows);
else
  y = zeros(rows, columns);
end
[rows, columns] = size(y);

y = filter([1-a], [1, -a], x);
end



%% f_inst - instantaneous frequency
function f_inst = calc_f_inst(sig,fs,tau,norm)
%
% function f_inst = calc_f_inst(sig,fs,tau,norm);
%
% Calculates instantaneous frequency from a complex (analytical) signal
% using first order differences
%
% input parameters:
%   sig  : complex (analytical) input signal
%   fs   : sampling frequency of sig
%   tau  : exponential decay time of temporal averaging filter
%   norm : exponent for amplitude weighting for temporal averaging filter
%          process(0: no level weigthing)
%
% output values:
%   f_inst:   vector of estimated inst. frequency values (w temporal
%            averaging
%
% copyright: Universitaet Oldenburg
% author   : vh
% date     : 12/04
%

sig = sig';

alpha = exp(-1/tau/fs);
b = [1-alpha];
a = [1 -alpha];

f_inst = sig./(abs(sig)+eps);
f_inst = abs(sig).^norm.*[0 f_inst(2:end).*conj(f_inst(1:end-1))];
f_inst = filter(b,a,f_inst);
f_inst = angle(f_inst')/2/pi*fs;

end



function ic = new_ic(itf,tau_coherence,fs)

%tau_coherence = 15e-3; % good value for ipd_fine

c_coh = exp(-1./(fs.*tau_coherence));
if length(tau_coherence)==1
    ic = abs(filter(1-c_coh,[1 -c_coh],itf))./abs(filter(1-c_coh,[1 -c_coh],abs(itf)));
elseif length(tau_coherence)==size(itf,2)
    ic = zeros(size(itf));
    for n = 1:length(tau_coherence)
        ic(:,n) = abs(filter(1-c_coh(n),[1 -c_coh(n)],itf(:,n)))./ ...
            abs(filter(1-c_coh(n),[1 -c_coh(n)],abs(itf(:,n))));
    end
else
    error('wrong number of tau_coherence values')
end
end


