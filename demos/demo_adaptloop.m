%DEMO_ADAPTLOOP  Show the effect of adaptation.
%
%   This script demonstrates the effect of adaptation applied to a test
%   signal with and without noise.
%
%   The test signal is made of a sinosoidal ramp up and down between 0
%   and 1.
%
%   FIGURE 1 Clean test signal
%
%     This figure shows the effect of adaptation on the clean test signal with and
%     without overshoot limiting.
%
%   FIGURE 2 Noisy test signal
%
%     This figure shows the effect of adaptation on the noisy test signal with and
%     without overshoot limiting.
%
%   See also: adaptloop


siglen=10000;
fs=10000;

% This is the default minimum level (0 dB) of the adaptation loops. The
% loops assume that a signal is never silent, and sets all values below
% minlvl equal to minlvl. For plotting purposes, we do the same explicitly.
minlvl=setdbspl(0);

part=siglen/10;

insig=[zeros(2*part,1);
       rampup(part);
       ones(2*part,1);
       rampdown(part);
       zeros(4*part,1)];

insig=max(insig,minlvl);

figure(1);
subplot(3,1,1);
plot(20*log10(insig));
title('Input signal');
ylabel('level / Db');

subplot(3,1,2);
plot(adaptloop(insig,fs,0));
title('Adaptation.');
ylabel('level / model units');

subplot(3,1,3);
plot(adaptloop(insig,fs));
title('Adaptation w. limiting.');
ylabel('level / model units');


% Add a low level of noise
insig=abs(insig+0.001*randn(siglen,1));
insig=max(insig,minlvl);

figure(2);

subplot(3,1,1);
plot(20*log10(insig));
title('Input signal with added Gaussian noise.');
ylabel('level / Db');

subplot(3,1,2);
plot(adaptloop(insig,fs,0));
title('Adaptation.');
ylabel('level / model units');

subplot(3,1,3);
plot(adaptloop(insig,fs));
title('Adaptation w. limiting.');
ylabel('level / model units');


%OLDFORMAT
