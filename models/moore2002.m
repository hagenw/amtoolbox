% function xxx = moore1997(yyy)
% params:
% - stationary / non stationary signal
%
% todo: two-channel inSig, andere fs (-> window sizes), inSig Normalisierung, fft
% Normalisierung nach Fenstern durch psd okay? (mehr Überlappen -> mehr Energie)

%% init
fs = 32000; % todo: für andere fs müssen Fenstergrößen/fftLen angepasst werden! (max(hannLenSmp) = fftLen)
fVec = 20:fs/2;
inSig = randn(fs,1);
t = linspace(0,1,fs);
inSig = sin(2*pi*1000*t).';

%% model
% todo: normalize inSig to SPL (pascalize??)

data = data_moore1997(fVec,0);

% filter order as in moore2002
order = 4096;

% create FIR filter
tfLinear = db2mag(data.tfOuterMiddle);
outerMiddleFilter = fir2(order, linspace(0, 1, length(fVec)), tfLinear);

% h = freqz(outerMiddleFilter,1,fVec,fs);
% figure
% semilogx(fVec,mag2db(abs(h)))
% grid on

earSig = filter(outerMiddleFilter,1,inSig);

fftLen = 2048; % according to moore2002

% todo: signal length < 1ms
% compute windows
hannLenMs = [2,4,8,16,32,64]; % hanning window size (moore2002) in ms
hannLenSmp = round(hannLenMs./1000 * fs); % windows size in samples

% calculate hann windows
hannWin2 = hann(hannLenSmp(1));
hannWin4 = hann(hannLenSmp(2));
hannWin8 = hann(hannLenSmp(3));
hannWin16 = hann(hannLenSmp(4));
hannWin32 = hann(hannLenSmp(5));
hannWin64 = hann(hannLenSmp(6));

timeStep = 0.001; % 1ms steps as in moore2002
updateRate = round(timeStep*fs);
nOverlap = round(hannLenSmp - updateRate); % number of overlapping samples

numBlocks = ceil(length(earSig)./updateRate);

earSigPad = earSig;
earSigPad(end+1:end+hannLenSmp(6)) = zeros(hannLenSmp(6),1);  % zero padding
% compute ffts in 1ms intervals
% for ii=1:numBlocks
%     % fft for each window
%     lower = (ii-1)*updateRate+1;
%     upper = (ii-1)*updateRate+hannLenSmp(1);
%     spect1(ii,:) = fft(earSigPad(lower:upper) .* hannWin2, fftLen);
%     
%     upper = (ii-1)*updateRate+hannLenSmp(2);
%     spect2(ii,:) = fft(earSigPad(lower:upper) .* hannWin4, fftLen);
%     
%     upper = (ii-1)*updateRate+hannLenSmp(3);
%     spect3(ii,:) = fft(earSigPad(lower:upper) .* hannWin8, fftLen);
% 
%     upper = (ii-1)*updateRate+hannLenSmp(4);
%     spect4(ii,:) = fft(earSigPad(lower:upper) .* hannWin16, fftLen);
% 
%     upper = (ii-1)*updateRate+hannLenSmp(5);
%     spect5(ii,:) = fft(earSigPad(lower:upper) .* hannWin32, fftLen);
% 
%     upper = (ii-1)*updateRate+hannLenSmp(6);
%     spect6(ii,:) = fft(earSigPad(lower:upper) .* hannWin64, fftLen);
% 
% end

binWidth = fs/(fftLen+2); %  bandwidth in Hz represented by 1 fft frequency bin
oneHz = (fftLen+2)/fs;  % number of frequency bins representing 1Hz
% 
% % truncate ffts to match the frequency ranges specified in moore2002
% % and put PSD for each window and each time interval in matrix -> window
% % normalization
% spect = zeros(numBlocks,fftLen/2+1);
% spect(:,round(4050*oneHz)+1:fftLen/2+1) = abs(spect1(:,round(4050*oneHz)+1:fftLen/2+1)).^2/sum(hannWin2.^2); % 4050-fs/2
% spect(:,round(2540*oneHz)+1:round(4050*oneHz)) = abs(spect2(:,round(2540*oneHz)+1:round(4050*oneHz))).^2/sum(hannWin4.^2); % 2540-4050Hz
% spect(:,round(1250*oneHz)+1:round(2540*oneHz)) = abs(spect3(:,round(1250*oneHz)+1:round(2540*oneHz))).^2/sum(hannWin8.^2); % 1250-2540Hz
% spect(:,round(500*oneHz)+1:round(1250*oneHz)) = abs(spect4(:,round(500*oneHz)+1:round(1250*oneHz))).^2/sum(hannWin16.^2); % 500-1250Hz
% spect(:,round(80*oneHz)+1:round(500*oneHz)) = abs(spect5(:,round(80*oneHz)+1:round(500*oneHz))).^2/sum(hannWin32.^2); % 80-500Hz
% spect(:,1:round(80*oneHz)) = abs(spect6(:,1:round(80*oneHz))).^2/sum(hannWin64.^2); % 0-80Hz
% compInt = 2*spect./fs;
% compdB = 10*log10(compInt); % intensity level of each frequency component in dB
% compFq = linspace(0,fs/2,fftLen/2+1);

% plot first fft time segment
% figure
% semilogx(linspace(0,fs/2,fftLen/2+1), compdB(1,:));

%% using spectrogram (stft) instead of fft
%   -> same result (not zero padded! -> different lengths)
% [S,F,T,P,Fc,Tc] = spectrogram(X,WINDOW,NOVERLAP,F,Fs,SPECTRUMTYPE,FREQRANGE)
[s1,compFq,~,p1] = spectrogram(earSigPad,hann(hannLenSmp(1)),nOverlap(1),fftLen,fs,'psd','onesided');
[s2,~,~,p2] = spectrogram(earSigPad,hann(hannLenSmp(2)),nOverlap(2),fftLen,fs,'psd','onesided');
[s3,~,~,p3] = spectrogram(earSigPad,hann(hannLenSmp(3)),nOverlap(3),fftLen,fs,'psd','onesided');
[s4,~,~,p4] = spectrogram(earSigPad,hann(hannLenSmp(4)),nOverlap(4),fftLen,fs,'psd','onesided');
[s5,~,~,p5] = spectrogram(earSigPad,hann(hannLenSmp(5)),nOverlap(5),fftLen,fs,'psd','onesided');
[s6,~,~,p6] = spectrogram(earSigPad,hann(hannLenSmp(6)),nOverlap(6),fftLen,fs,'psd','onesided');
compFq = compFq.';
compInt = zeros(numBlocks,fftLen/2+1);  % calculate component intensity (psd)
compInt(:,round(4050*oneHz)+1:fftLen/2+1) = p1(round(4050*oneHz)+1:fftLen/2+1,1:numBlocks).'; % 4050-fs/2
compInt(:,round(2540*oneHz)+1:round(4050*oneHz)) = p2(round(2540*oneHz)+1:round(4050*oneHz),1:numBlocks).'; % 2540-4050Hz
compInt(:,round(1250*oneHz)+1:round(2540*oneHz)) = p3(round(1250*oneHz)+1:round(2540*oneHz),1:numBlocks).'; % 1250-2540Hz
compInt(:,round(500*oneHz)+1:round(1250*oneHz)) = p4(round(500*oneHz)+1:round(1250*oneHz),1:numBlocks).'; % 500-1250Hz
compInt(:,round(80*oneHz)+1:round(500*oneHz)) = p5(round(80*oneHz)+1:round(500*oneHz),1:numBlocks).'; % 80-500Hz
compInt(:,1:round(80*oneHz)) = p6(1:round(80*oneHz),1:numBlocks).'; % 0-80Hz
compdB = 10*log10(compInt); % intensity level of each frequency component in dB
%% calculating excitation patterns
% calculate ERB numbers corresponding to ERB mid frequencies
erbStep = 0.25;    % according to moore1997, moore2002
erbFcMin = 50;
erbFcMax = 15000;
erbNMin = fc2erbN(erbFcMin);
erbNMax = fc2erbN(erbFcMax);
erbN = erbNMin:erbStep:erbNMax;    % numbers of erb bands
erbFc = erbN2fc(erbN);               % center frequency of erb bands

erbLoFreq = erbN2fc(erbN-0.5); % lower limit of each ERB filter
erbHiFreq = erbN2fc(erbN+0.5); % upper limit of each ERB filter

% calculate intensity for each ERB (dB/ERB) and each time step
for ii=1:length(erbFc)
    range = round(erbLoFreq(ii)*oneHz):round(erbHiFreq(ii)*oneHz);
    erbInt(:,ii) = sum(compInt(:,range),2);   % intensity sum in each erb
end
erbdB = 10*log10(erbInt);   % intensity level in each erb

% p determines bandwidth and slope of the erb filter and is generally
% asymmetrical
% p is roughly symmetrical for an excitation level of 51dB per ERB
p51 = 4*erbFc./f2erb(erbFc);    % p for erb center frequencies and a level of 51dB
p511 = 4*1000/f2erb(1000);    % p for fc=1kHz and a level of 51dB (at 1kHz filters are symmetrical)

pU = p51;   % pU for all erbFc and all time steps
g = abs(repmat(compFq,149,1) - repmat(erbFc.',1,length(compFq)))./repmat(erbFc.',1,length(compFq));    % normalized deviation of each f to erbFc for each erb band
for ii = 1:numBlocks % time steps
    erbdB2f(ii,:) = interp1([0 erbFc fs/2], [min(erbdB(ii,:)) erbdB(ii,:) min(erbdB(ii,:))], compFq);   % map erbFc to compFq
    pL(ii,:,:) = repmat(p51,length(compFq),1) - 0.35.*(repmat(p51,length(compFq),1)./p511).*(repmat(erbdB2f(ii,:).',1,length(erbN))-51); 
    p(ii,:,:) = pL(ii,:,:);
    for jj = 1:length(erbN)
        p(ii,round(erbFc(jj)*oneHz)+1:end,jj) = pU(jj); % p(f>erbFc) = pU
    end
    w(ii,:,:) = (1+squeeze(p(ii,:,:)).*g.').*exp(-squeeze(p(ii,:,:)).*g.');    % calculate weighting function
    e(ii,:,:) = squeeze(w(ii,:,:)).*compInt(ii,:).';
    eL(ii,:) = sum(squeeze(e(ii,:,:)),1);   % sum excitation level in each erb
end
eLdB = 10*log10(eL);

figure
for ii=10:10:140
    plot(compFq,w(500,:,ii))
    hold on
end
figure
for ii=1:10:length(e(1,1,:))
    plot(fc2erbN(compFq),e(500,:,ii))
    hold on
end
figure
plot(erbFc,eLdB(1,:))

%% calculating specific loudness 

dataFc = data_moore1997(erbFc);
tQdB = dataFc.tQ;
tQ = 10.^(tQdB./10);
data500 = data_moore1997(500); % internal excitation level at threshold for freq. 500Hz and above
tQdB500 = data500.tQ;
gdB = tQdB500-tQdB;    % low level gain in cochlea amplifier
g = 10.^((tQdB500-tQdB)/10);
dataA = data_moore1997([],[],gdB);
a = dataA.a;    % parameter for linearization around absolute threshold
alpha = dataA.alpha;    % compressive exponent
c = 0.046871; % constant to get loudness scale to sone

specLoud = zeros(size(eL));

specLoud1 = c*(2*eL./(eL+repmat(tQ,numBlocks,1))).^1.5 .*((repmat(g,numBlocks,1).* ...
    eL+repmat(a,numBlocks,1)).^alpha-repmat(a,numBlocks,1).^alpha);
specLoud2 = c * ((repmat(g,numBlocks,1) .*eL+repmat(a,numBlocks,1)).^alpha - ...
    repmat(a,numBlocks,1).^alpha);
specLoud3 = c*(eL./1.04e6).^0.5;
specLoud(eL<repmat(tQ,numBlocks,1)) = specLoud1(eL<repmat(tQ,numBlocks,1));
specLoud(eL<=10^10 & eL>repmat(tQ,numBlocks,1)) = specLoud2(eL<=10^10 & eL>repmat(tQ,numBlocks,1));
specLoud(eL>10^10) = specLoud3(eL>10^10);

%% instant loudness, short term loudness (STL), long term loudness (LTL)
instLoudness = sum(specLoud,2) * erbStep;     % integrate over the erbs

% STL and LTL:
aSTL = 0.045;
rSTL = 0.02;
STL = zeros(length(instLoudness),1);
aLTL = 0.01;
rLTL = 0.0005;
LTL = zeros(length(instLoudness),1);
for ii = 2:length(instLoudness)
    if instLoudness(ii)>STL(ii-1)
        STL(ii) = aSTL*instLoudness(ii)+(1-aSTL)*STL(ii-1);
    else
        STL(ii) = rSTL*instLoudness(ii)+(1-rSTL)*STL(ii-1);
    end
    if STL(ii)>LTL(ii-1)
        LTL(ii) = aLTL*STL(ii)+(1-aLTL)*LTL(ii-1);
    else
        LTL(ii) = rLTL*STL(ii)+(1-rLTL)*LTL(ii-1);
    end
end
figure
plot(STL);
figure
plot(LTL);
