% demonstration of the Braasch (2013) model
% using a lead/lag pair with ITDs of +-0.3ms
% at different lag/lead ratios and interstimulus 
% intervals
% uses auditoryfilterbank from AMT toolbox

Fs=48000;          % Sampling Frequency [Hz]
ISI=[0 1 3];       % Inter-stimulus Interval [ms]
Ratio=[-60 0 4 8]; % lead/lag level ratio [dB]
StimDuration=200;  % stimulus duration [ms]

ITD=zeros(length(Ratio),length(ISI));
ILD=zeros(length(Ratio),length(ISI));
for m=1:length(Ratio)
    for n=1:length(ISI)
            y = stimulus(Fs, 2, 200, 500, 800, 0.3, -0.3, ISI(n), 20, 20, 0,db2amp(Ratio(m))); 
            [ITD(m,n),ILD(m,n)]=acmod(y,Fs);
    end % of for 
end
    
for m=1:length(Ratio)
    for n=1:length(ISI)
            y = stimulus(Fs, 2, 200, 500, 800, -0.3, 0.3, ISI(n), 20, 20, 0,db2amp(Ratio(m))); 
            [ITD2(m,n),ILD2(m,n)]=acmod(y,Fs);
    end % of for 
end

% display output
figure;
subplot(2,1,1)
plot(ISI,ITD(1,:),'k+-');
hold on
plot(ISI,ITD(2,:),'kx--');
plot(ISI,ITD(3,:),'k*:');
plot(ISI,ITD(4,:),'k^-.');
hold off
legend('no lag','0 dB','4 dB', '8 dB')
axis([-0.2 3.2 -0.3 0.3])
xlabel('Interstimulus interval [ms]')
ylabel('ITD [ms]')
title('lead left, lag right')

subplot(2,1,2)
plot(ISI,ITD2(1,:),'k+-');
hold on
plot(ISI,ITD2(2,:),'kx--');
plot(ISI,ITD2(3,:),'k*:');
plot(ISI,ITD2(4,:),'k^-.');
hold off
legend('no lag','0 dB','4 dB', '8 dB')
axis([-0.2 3.2 -0.3 0.3])
xlabel('Interstimulus interval [ms]')
ylabel('ITD [ms]')
title('lead right, lag left')

