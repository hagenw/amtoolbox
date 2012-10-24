function varargout = baumgartner2013( target,template,varargin )
%BAUMGARTNER2013 Model for localization in saggital planes
%   Usage:    [p,respang] = baumgartner2013( target,template,varargin )
%
%   Input parameters:
%     target  : binaural impulse response(s) referring to the directional 
%               transfer function(s) (DFTs) of the target sound(s).
%     template: binaural impulse responses of all available
%               listener-specific DTFs of the sagittal plane referring to
%               the perceived lateral angle of the target sound
%
%   Output parameters:
%     p       : predicted probability mass vectors for response angles 
%               with respect to target positions
%               1st dim: polar response angle
%               2nd dim: polar taget angle
%     respang : polar response angles (after regularization of angular 
%               sampling)
%
%   `baumgartner2013(...)` is a model for sound-source localization in
%   sagittal planes (SPs). It bases on the comparison of internal sound 
%   representation with a template and results in a probabilistic
%   prediction of polar angle response.
%
%   `baumgartner2013` accepts the following optional parameters:
%
%     'fs',fs        Define the sampling rate of the impulse responses. 
%                    Default value is 48000 Hz.
%
%     'stim',stim    Define the stimulus (source signal without directional
%                    features). As default an impulse is used.
%
%     'fsstim',fss   Define the sampling rate of the stimulus. 
%                    Default value is 48000 Hz.
%
%     'flow',flow    Set the lowest frequency in the filterbank to
%                    flow. Default value is 700 Hz.
%
%     'fhigh',fhigh  Set the highest frequency in the filterbank to
%                    fhigh. Default value is 18000 Hz.
%
%     'lat',lat      Set the perceived lateral angle of the target sound to
%                    lat. Default value is 0° (median SP).
%
%     'u',u          Set the listener-specific uncertainty (standard
%                    deviation of the Gaussian transformation from the
%                    distance metric of the comparison process to the
%                    similarity index) to u. Default value is 2 dB.
%
%     'space',s      Set spacing of auditory filter bands to s numbers of
%                    equivalent rectangular bandwidths (ERBs). 
%                    Default value is 1 ERB.
%
%     'bwsteep',bws  Set the steepness factor bws of the sigmoid function 
%                    applied for binaural weighting of monaural similarity 
%                    indices. Default value is 13°.
%
%     'polsamp',ps   Define the the polar angular sampling of the current
%                    SP. As default the sampling of ARI's HRTF format at
%                    the median SP is used, i.e.,
%                    ps = [-30:5:70,80,100,110:5:210] .
%
%   `baumgartner2013` accepts the following flags:
%
%     'gt'           Use the Gammatone filterbank for peripheral processing. 
%                    This is the default.
%
%     'cqdft'        Use a filterbank approximation based on DFT with 
%                    constant relative bandwidth for peripheral processing. 
%                    This was used by Langendijk and Bronkhorst (2002).
%
%     'ihc'          Incorporate the transduction model of inner hair 
%                    cells used by Dau et al. (1996). This is the default.
%
%     'noihc'        Do not incorporate the IHC stage.
%
%     'regular'      Apply spline interpolation in order to regularize the 
%                    angular sampling of the polar response angle. 
%                    This is the default.
%
%     'noregular'    Disable regularization of angular sampling.
%
%
%   Example:
%   ---------
%
%   To compute and visualize the baseline prediction (localizing broadband  
%   sounds with own ears) for listener NH58 and the median SP use :::
%
%     p = demo_baumgartner2013('NH58');
%
%   See also: plotbaumgartner2013, data_baumgartner2013,
%   exp_baumgartner2013
%
%   References: baumgartner2013assessment baumgartner2012modelling langendijk2002contribution patterson1988efficient dau1996qmeI
%
% AUTHOR: Robert Baumgartner, Acoustics Research Institute, Vienna, Austria

%% Check input options 

definput.flags.fbank = {'gt','cqdft','drnl'};
definput.flags.headphonefilter = {'','headphone'};
definput.flags.middleearfilter = {'','middleear'};
definput.flags.ihc = {'ihc','noihc'};
definput.flags.regularization = {'regular','noregular'};

% CP-Falgs:
definput.flags.cp={'std','xcorr'};

definput.keyvals.fs=48000;      % Hz
definput.keyvals.space=1;       % No. of ERBs (Cams)
definput.keyvals.do=0;
definput.keyvals.u=2;           % listener-specific uncertainty in dB
definput.keyvals.lat=0;         % deg
definput.keyvals.flow=1000;     % Hz
definput.keyvals.fhigh=18000;   % Hz
definput.keyvals.lvlstim = 40; 	% dBSPL
definput.keyvals.lvltem = 40;  	% dBSPL
definput.keyvals.stim=[];
definput.keyvals.fsstim=[];
definput.keyvals.bwsteep=13;    % steepness in degrees of binaural weighting function
definput.keyvals.polsamp=[-30:5:70 80 100 110:5:210];  % polar sampling (for regular)

[flags,kv]=ltfatarghelper(...
  {'fs','space','do','u','lat','flow','fhigh',...
  'lvlstim','lvltem','stim','fsstim','bwsteep','polsamp'},definput,varargin);


%% Stimulus 
if isempty(kv.stim) 
    kv.stim = [1;0];%[1;zeros(size(target,1),1)];    % impulse
    kv.fsstim = kv.fs;
end

if flags.do_headphone% || flags.do_drnl
    hpfilt = headphonefilter(kv.fs);
    kv.stim = convolve(kv.stim,hpfilt(:));
end

if flags.do_middleear% || flags.do_drnl
    miearfilt = middleearfilter(kv.fs);
    kv.stim = convolve(kv.stim,miearfilt(:));
end


%% DTF filtering
if ~isequal(kv.fs,kv.fsstim)
    disp('Sorry, sampling rate of stimulus and HRIRs must be equal!')
    return
end

tmp = convolve(target,kv.stim);
target = reshape(tmp,[size(tmp,1),size(target,2),size(target,3)]);


%% Cochlear filter bank -> internal representations
if flags.do_cqdft
    
    bpo = kv.space*6; % bands per octave (1 oct. approx. as 6 ERBs)
    ireptem = cqdft(template,kv.fs,kv.flow,kv.fhigh,bpo);
    ireptar = cqdft(target,kv.fs,kv.flow,kv.fhigh,bpo);

elseif flags.do_gt

    % Determine filterbank
    fc = audspacebw(kv.flow,kv.fhigh,kv.space,'erb');
    Nfc = length(fc);   % # of bands
    [bgt,agt] = gammatone(fc,kv.fs,'classic');
    
    % Filtering
    ireptar = ufilterbankz(bgt,agt,target(:,:)); % channel (3rd) dimension resolved!
    ireptem = ufilterbankz(bgt,agt,template(:,:)); % resolve 3rd dim
    
    % IHC transduction
    if flags.do_ihc
        ireptar = ihcenvelope(ireptar,kv.fs,'ihc_dau');
        ireptem = ihcenvelope(ireptem,kv.fs,'ihc_dau');
    end
    
    % Set back the channel dimension
    ireptar = reshape(ireptar,[size(target,1),... 
        Nfc,size(target,2),size(target,3)]);
    ireptem = reshape(ireptem,[size(template,1),Nfc,size(template,2),size(template,3)]);
    
    % Averaging over time (RMS)
    ireptar = db(squeeze(rms(ireptar)));      % in dB !!!!!!
    ireptem = db(squeeze(rms(ireptem)));
        
elseif flags.do_drnl   
    
    % Set level
    for ch = 1:size(template,3)
        target(:,:,ch) = setdbspl(target(:,:,ch),kv.lvlstim);
        template(:,:,ch) = setdbspl(template(:,:,ch),kv.lvltem);
    end
    
    % Filtering
    [ireptar,fc] = drnl(target(:,:),kv.fs,'flow',kv.flow,'fhigh',kv.fhigh);  % includes middle ear
    ireptem = drnl(template(:,:),kv.fs,'flow',kv.flow,'fhigh',kv.fhigh);
       
    % IHC transduction
    if flags.do_ihc 
        ireptar = ihcenvelope(ireptar,kv.fs,'ihc_dau');
        ireptem = ihcenvelope(ireptem,kv.fs,'ihc_dau');
    end
    
    % Set back the channel dimension
    ireptar = reshape(ireptar,[size(ireptar,1),... 
        length(fc),size(target,2),size(target,3)]);
    ireptem = reshape(ireptem,[size(ireptem,1),length(fc),size(template,2),size(template,3)]);
    
    % Averaging over time (RMS)
    ireptar = db(squeeze(rms(ireptar)));
    ireptem = db(squeeze(rms(ireptem)));

end

if size(ireptar,2) ~= size(target,2) % retreive polar dimension if squeezed out
    ireptar = reshape(ireptar,[size(ireptar,1),size(target,2),size(target,3)]);
end


%% Comparison process -> monaural similarity indices (SIs)
si=zeros(size(template,2),size(target,2),size(template,3)); % initialisation
for it = 1:size(target,2)
	si(:,it,:) = langendijkcomp(ireptar(:,it,:),ireptem,'s',kv.u, ...
    'argimport',flags,kv);
end


%% Binaural weighting -> binaural SIs
if size(si,3) == 2
    binw = 1./(1+exp(-kv.lat/kv.bwsteep)); % weight of left ear signal with 0 <= binw <= 1
    si = binw * si(:,:,1) + (1-binw) * si(:,:,2);
end


%% Interpolation (regularize polar angular sampling)
if flags.do_regular
    respang0 = ceil(min(kv.polsamp)*0.2)*5;    % ceil to 5°
    respangs = respang0:5:max(kv.polsamp);
    siint = zeros(length(respangs),size(si,2));
    for tt = 1:size(si,2)
        siint(:,tt) = interp1(kv.polsamp,si(:,tt),respangs,'spline');
    end
    si = siint;
    si(si<0) = 0; % SIs must be positive (necessary due to spline interp)
else
    respangs = kv.polsamp;
end


%% Normalization to PMV
p = si ./ repmat(sum(si),size(si,1),1);


%% Output
varargout{1} = p;
if nargout == 2
    varargout{2} = respangs;
end
  
end