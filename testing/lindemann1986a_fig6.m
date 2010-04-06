%LINDEMANN1986a_FIG6 Reproduces fig. 6 from lindemann1986a
%
%   Reproduces fig.6: cross-correlation for different inhibition factors c_s
%
%R lindemann1986a
%

%   AUTHOR: Hagen Wierstorf

% Sampling rate
fs = 44100;
% Frequency of the sinusoid
f = 500;
% Plot results
plotting = true;

% Calculate crosscorrelations for 21 ITD points between 0~ms and 1~ms
ii = 1;
cc_a = zeros(21,45);
cc_b = cc_a;
cc_c = cc_a;
for itd = linspace(0,1,21); 
    % Generate ITD shifted sinusoid
    sig = itdsin(f,itd,fs)(1:441,:);
    % Calculate crosscorrelation pattern for
    % fig. 6 (a)
    cc_a(ii,:) = mean(lindemann(sig,0,0,fs),3);
    % fig. 6 (b)
    cc_b(ii,:) = mean(lindemann(sig,0,0.3,fs),3);
    % fig. 6 (c)
    cc_c(ii,:) = mean(lindemann(sig,0,1,fs),3);
    
    ii = ii+1;
end

if(plotting)
    % Generate time axis
    t = linspace(1,0,21);
    tau = linspace(-1,1,45);
    %
    figure;
    mesh(tau,t,cc_a);
    view(0,57);
    xlabel('correlation-time tau (ms)');
    ylabel('interaural time difference (ms)');
    title('c_s = 0');
    %
    figure;
    mesh(tau,t,cc_b);
    view(0,57);
    xlabel('correlation-time tau (ms)');
    ylabel('interaural time difference (ms)');
    title('c_s = 0.3');
    %
    figure;
    mesh(tau,t,cc_c);
    view(0,57);
    xlabel('correlation-time tau (ms)');
    ylabel('interaural time difference (ms)');
    title('c_s = 1');
end

