function [results] = moore1997(inSig,fs)
    % moore1990 model for stationary signals

    %% init -> for testing, remove
%     fs = 32000; % todo: allow different fs -> adapt window sizes (max(hannLenSmp) = fftLen)
%     t = linspace(0,1,fs);
%     sig = sin(2*pi*1000*t).';
%     inSig = setdbspl(sig,100);  % setdbspl needs AMToolbox!

    %% model
    fVec = 20:fs/2;
    data = data_moore2002('tfOuterMiddle1997','fieldType','free','fVec',fVec);

    % filter order as in moore2002
    order = 4096;
    % create FIR filter
    tfLinear = 10.^(data.tfOuterMiddle/10);
    outerMiddleFilter = fir2(order, linspace(0, 1, length(fVec)), tfLinear);
    earSig = filtfilt(outerMiddleFilter,1,inSig);   % why does filter(..) not work?

    % compute fft
    spect = fft(earSig);
    fftLen = length(spect);

    binWidth = fs/(fftLen+2); %  bandwidth in Hz represented by 1 fft frequency bin
    oneHz = (fftLen+2)/fs;  % number of frequency bins representing 1Hz

    numBins = round(fftLen/2+1);
    compInt = 2*abs(spect(1:numBins)).^2/(numBins*fs);  % psd
    compdB = 10*log10(compInt./(20e-6)^2); % intensity level in dBSPL
    compFq = linspace(0,fs/2,numBins);
    nPoints = length(compFq);
    compErb = fc2erbN(compFq);

    % calculate ERB numbers corresponding to ERB mid frequencies
    erbStep = 0.25;
    erbFcMin = 50;
    erbFcMax = 15000;
    erbNMin = fc2erbN(erbFcMin);
    erbNMax = fc2erbN(erbFcMax);
    erbN = erbNMin:erbStep:erbNMax;    % numbers of erb bands
    erbFc = erbN2fc(erbN);               % center frequency of erb bands

    erbLoFreq = erbN2fc(erbN-0.5); % lower limit of each ERB filter
    erbHiFreq = erbN2fc(erbN+0.5); % upper limit of each ERB filter

    %calculate intensity for each ERB (dB/ERB)
    for ii=1:length(erbFc)
        range = round(erbLoFreq(ii)*oneHz):round(erbHiFreq(ii)*oneHz);
        erbInt(ii) = sum(compInt(range));   % intensity sum in each erb
    end

    erbdB = 10*log10(erbInt./(20e-6)^2);   % intensity level in each erb using reference SPL of 20 uPa
    p511 = 4*1000/f2erb(1000);    % p for fc=1kHz and a level of 51dB (at 1kHz filters are symmetrical)
    erbdB2F = interp1([0 erbFc fs/2], [min(erbdB) erbdB min(erbdB)], compFq);   % map erbFc to compFq

    for e = 1:length(erbN)
        erb = f2erb(erbFc(e));
        p51 = 4*erbFc(e)/erb;
        intensity = 0;
        for comp = 1:nPoints
            g = (compFq(comp)-erbFc(e))/erbFc(e);
            if g<0
                p = p51 - 0.35*(p51/p511) * (erbdB2F(comp)-51);
            else
                p = p51;
            end
            g = abs(g);
            w = (1+p*g)*exp(-p*g);
            intensity = intensity+w*compInt(comp);  %intensity per erb
        end
        eL(e) = intensity;
    end
    results.eLdB = 10*log10(eL./(20e-6)^2); % get dB SPL (20uPa reference)
    results.erbN = erbN;
    
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
    specLoud1 = c*(2*eL./(eL+tQ)).^1.5 .*((g.* eL + a).^alpha-a.^alpha);
    specLoud2 = c * ((g .*eL + a).^alpha - a.^alpha);
    specLoud3 = c*(eL./1.04e6).^0.5;
    specLoud(eL<tQ) = specLoud1(eL<tQ);
    specLoud(eL<=10^10 & eL>tQ) = specLoud2(eL<=10^10 & eL>tQ);
    specLoud(eL>10^10) = specLoud3(eL>10^10);

    %% monaural/binaural loudness (= instantaneous loudness), short term loudness (STL), long term loudness (LTL)
    results.monauralLoudness = sum(specLoud) * erbStep;     % integrate over the erbs
    results.binauralLoudness = 2*results.monauralLoudness;  % integrate moore2007 model (Modeling binaural loudness) for better results

end