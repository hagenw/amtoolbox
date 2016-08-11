function signalout = breebaart2001siggen(inttype,fc,sl,sdur,sphase,nbw,nl,ndur,nphase,hannramp,fs)
%BREEBAART2001SIGGEN computes the signals for the breebaart2001 experiment

% generate noise
noise(:,1) = bandpassnoisefreq(fc,fs,ndur,nl,nbw);

% set diotic phase difference for noise masker
if nphase == 0
    noise(:,2) = noise(:,1);
elseif nphase == pi;
    noise(:,2) = -noise(:,1);
else % set correlation
    noise(:,2) = bandpassnoisefreq(fc,fs,ndur,nl,nbw);
    noiseL = 0.5*sqrt(2)*sqrt(1+nphase)*noise(:,1) + 0.5*sqrt(2)*...
        sqrt(1-nphase)*noise(:,2);
    noiseR = 0.5*sqrt(2)*sqrt(1+nphase)*noise(:,1) - 0.5*sqrt(2)*...
        sqrt(1-nphase)*noise(:,2);
    noise = [noiseL, noiseR];
end

% apply hanning ramps to noise
n_ramp = round(hannramp*fs); 
noise= rampsignal(noise,n_ramp);

% set noiselevel
noise = setdbspl(noise,nl);

if strcmp(inttype,'target')

    signal(:,1) = sin(2*pi*(0:(sdur*fs)-1)'*fc/fs);

    % set diotic phase difference for signal
    if sphase == 0
        signal(:,2) = signal(:,1);
    else
        signal(:,2) = sin(2*pi*(0:(sdur*fs)-1)'*fc/fs-sphase);
    end

    % apply hanning ramps to signal  
    signal = rampsignal(signal,n_ramp);

    % set signal level
    signal = setdbspl(signal,sl);

    if sdur < ndur
        zerolength = round((ndur-sdur)*fs/2);
        signal = [zeros(zerolength,2); signal; zeros(zerolength,2)];
    end
    
    signalout = signal + noise;

elseif strcmp(inttype,'reference')
    signalout = noise;
end
    

