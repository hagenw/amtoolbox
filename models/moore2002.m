function [results] = moore2002(inSig,fs)
%
%
% todo: two-channel inSig, different fs (-> window sizes), inSig normalization

    %% init -> for testing, todo: remove
%     fs = 32000; % todo: allow different fs -> adapt window sizes (max(hannLenSmp) = fftLen)
%     t = linspace(0,1,fs);
%     sig = sin(2*pi*1000*t).';
%     inSig = setdbspl(sig,100);  % set dbspl, needs AMToolbox!

    %% model
    fVec = 20:fs/2;
    data = data_moore2002('tfOuterMiddle1997','fieldType','free','fVec',fVec);

    % filter order as in moore2002
    order = 4096;
    % create FIR filter
    tfLinear = 10.^(data.tfOuterMiddle/10);
    outerMiddleFilter = fir2(order, linspace(0, 1, length(fVec)), tfLinear);
    earSig = filtfilt(outerMiddleFilter,1,inSig);   % why does filter(..) not work?

    fftLen = 2048; % according to moore2002

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

    numBlocks = ceil(length(earSig)./updateRate);

    earSigPad = earSig;
    earSigPad(end+1:end+hannLenSmp(6)) = zeros(hannLenSmp(6),1);  % zero padding
%     compute ffts in 1ms intervals
%     todo: preallocation for speed...
    for ii=1:numBlocks
        % fft for each window
        lower = (ii-1)*updateRate+1;
        upper = (ii-1)*updateRate+hannLenSmp(1);
        spect1(ii,:) = fft(earSigPad(lower:upper) .* hannWin2, fftLen);
        
        upper = (ii-1)*updateRate+hannLenSmp(2);
        spect2(ii,:) = fft(earSigPad(lower:upper) .* hannWin4, fftLen);
        
        upper = (ii-1)*updateRate+hannLenSmp(3);
        spect3(ii,:) = fft(earSigPad(lower:upper) .* hannWin8, fftLen);
    
        upper = (ii-1)*updateRate+hannLenSmp(4);
        spect4(ii,:) = fft(earSigPad(lower:upper) .* hannWin16, fftLen);
    
        upper = (ii-1)*updateRate+hannLenSmp(5);
        spect5(ii,:) = fft(earSigPad(lower:upper) .* hannWin32, fftLen);
    
        upper = (ii-1)*updateRate+hannLenSmp(6);
        spect6(ii,:) = fft(earSigPad(lower:upper) .* hannWin64, fftLen);
    end

    binWidth = fs/(fftLen+2); %  bandwidth in Hz represented by 1 fft frequency bin
    oneHz = (fftLen+2)/fs;  % number of frequency bins representing 1Hz
    
    % truncate ffts to match the frequency ranges specified in moore2002
    % and put PSD for each window and each time interval in matrix -> window
    % normalization
    spect = zeros(numBlocks,fftLen/2+1);
    spect(:,round(4050*oneHz)+1:fftLen/2+1) = abs(spect1(:,round(4050*oneHz)+1:fftLen/2+1)).^2/sum(hannWin2.^2); % 4050-fs/2
    spect(:,round(2540*oneHz)+1:round(4050*oneHz)) = abs(spect2(:,round(2540*oneHz)+1:round(4050*oneHz))).^2/sum(hannWin4.^2); % 2540-4050Hz
    spect(:,round(1250*oneHz)+1:round(2540*oneHz)) = abs(spect3(:,round(1250*oneHz)+1:round(2540*oneHz))).^2/sum(hannWin8.^2); % 1250-2540Hz
    spect(:,round(500*oneHz)+1:round(1250*oneHz)) = abs(spect4(:,round(500*oneHz)+1:round(1250*oneHz))).^2/sum(hannWin16.^2); % 500-1250Hz
    spect(:,round(80*oneHz)+1:round(500*oneHz)) = abs(spect5(:,round(80*oneHz)+1:round(500*oneHz))).^2/sum(hannWin32.^2); % 80-500Hz
    spect(:,1:round(80*oneHz)) = abs(spect6(:,1:round(80*oneHz))).^2/sum(hannWin64.^2); % 0-80Hz
    compInt = 2*spect./fs;  %   psd
    compdB = 10*log10(compInt./(20e-6)^2); % intensity level of each frequency component in dB
    compFq = linspace(0,fs/2,fftLen/2+1);

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
    erbdB = 10*log10(erbInt./(20e-6)^2);   % intensity level in each erb using reference SPL of 20 uPa

    % p determines bandwidth and slope of the erb filter and is generally
    % asymmetrical
    % p is roughly symmetrical for an excitation level of 51dB per ERB
    p51 = 4*erbFc./f2erb(erbFc);    % p for erb center frequencies and a level of 51dB
    p511 = 4*1000/f2erb(1000);    % p for fc=1kHz and a level of 51dB (at 1kHz filters are symmetrical)

    pU = p51;   % pU for all erbFc and all time steps
    g = abs(repmat(compFq,149,1) - repmat(erbFc.',1,length(compFq)))./repmat(erbFc.',1,length(compFq));    % normalized deviation of each f to erbFc for each erb band
    % todo: preallocation for speed...
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
    results.eLdB = 10*log10(eL./(20e-6)^2);
    results.erbN = erbN;
    % plot excitation pattern at some time
    % figure
    % plot(erbN,results.eLdB(50,:))
    % grid on

    %% calculating specific loudness 

    dataSL = data_moore2002('specLoud','fVec',erbFc);
    tQdB = dataSL.tQ;
    tQ = 10.^(tQdB./10);
    tQdB500 = dataSL.tQ500;
    gdB = dataSL.g;    % low level gain in cochlea amplifier
    g = 10.^((tQdB500-tQdB)/10);
    a = dataSL.a;    % parameter for linearization around absolute threshold
    alpha = dataSL.alpha;    % compressive exponent
    c = dataSL.c; % constant to get loudness scale to sone

    specLoud = zeros(size(eL));

    specLoud1 = c*(2*eL./(eL+repmat(tQ,numBlocks,1))).^1.5 .*((repmat(g,numBlocks,1).* ...
        eL+repmat(a,numBlocks,1)).^alpha-repmat(a,numBlocks,1).^alpha);
    specLoud2 = c * ((repmat(g,numBlocks,1) .*eL+repmat(a,numBlocks,1)).^alpha - ...
        repmat(a,numBlocks,1).^alpha);
    specLoud3 = c*(eL./1.04e6).^0.5;
    specLoud(eL<repmat(tQ,numBlocks,1)) = specLoud1(eL<repmat(tQ,numBlocks,1));
    specLoud(eL<=10^10 & eL>repmat(tQ,numBlocks,1)) = specLoud2(eL<=10^10 & eL>repmat(tQ,numBlocks,1));
    specLoud(eL>10^10) = specLoud3(eL>10^10);

    %% monaural/binaural loudness (= instantaneous loudness), short term loudness (STL), long term loudness (LTL)
    results.monauralLoudness = sum(specLoud,2) * erbStep;     % integrate over the erbs
    results.binauralLoudness = 2*results.monauralLoudness;  % integrate moore2007 (Modeling binaural loudness) for better results

    % STL and LTL:
    aSTL = 0.045;
    rSTL = 0.02;
    STL = zeros(length(results.binauralLoudness),1);
    aLTL = 0.01;
    rLTL = 0.0005;
    LTL = zeros(length(results.binauralLoudness),1);
    for ii = 2:length(results.binauralLoudness)
        if results.binauralLoudness(ii)>STL(ii-1)
            STL(ii) = aSTL*results.binauralLoudness(ii)+(1-aSTL)*STL(ii-1);
        else
            STL(ii) = rSTL*results.binauralLoudness(ii)+(1-rSTL)*STL(ii-1);
        end
        if STL(ii)>LTL(ii-1)
            LTL(ii) = aLTL*STL(ii)+(1-aLTL)*LTL(ii-1);
        else
            LTL(ii) = rLTL*STL(ii)+(1-rLTL)*