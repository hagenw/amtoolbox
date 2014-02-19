function varargout = baumgartner2014( target,template,varargin )
%BAUMGARTNER2014 Model for localization in saggital planes
%   Usage:    [p,respang] = baumgartner2014( target,template )
%             [p,respang,tang] = baumgartner2014( target,template )
%             [p,respang,tang] = baumgartner2014( target,template,fs,S,lat,stim,fsstim )
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
%               1st dim: response angle
%               2nd dim: target angle
%     respang : polar response angles (after regularization of angular 
%               sampling)
%     tang    : polar target angles (usefull if sagittal-plane HRTFs are
%               extracted directly from SOFA object)
%
%   `baumgartner2014(...)` is a model for sound-source localization
%   in sagittal planes (SPs). It bases on the comparison of internal sound 
%   representation with a template and results in a probabilistic
%   prediction of polar angle response.
%
%   `baumgartner2014` accepts the following optional parameters:
%
%     'fs',fs        Define the sampling rate of the impulse responses. 
%                    Default value is 48000 Hz.
%
%     'S',S          Set the listener-specific sensitivity threshold 
%                    (threshold of the sigmoid link function representing 
%                    the psychometric link between transformation from the
%                    distance metric and similarity index) to *S*. 
%                    Default value is 1.
%
%     'lat',lat      Set the apparent lateral angle of the target sound to
%                    *lat*. Default value is 0° (median SP).
%
%     'stim',stim    Define the stimulus (source signal without directional
%                    features). As default an impulse is used.
%
%     'fsstim',fss   Define the sampling rate of the stimulus. 
%                    Default value is 48000 Hz.
%
%     'flow',flow    Set the lowest frequency in the filterbank to
%                    *flow*. Default value is 700 Hz.
%
%     'fhigh',fhigh  Set the highest frequency in the filterbank to
%                    *fhigh*. Default value is 18000 Hz.
%
%     'space',sp     Set spacing of auditory filter bands (i.e., distance 
%                    between neighbouring bands) to *sp* in number of
%                    equivalent rectangular bandwidths (ERBs). 
%                    Default value is 1 ERB.
%
%     'conGain',cg   Set the contralateral gain *cg* of the sigmoid function 
%                    applied for binaural weighting of monaural similarity 
%                    indices. Default value is 13 degrees.
%
%     'polsamp',ps   Define the the polar angular sampling of the current
%                    SP. As default the sampling of ARI's HRTF format at
%                    the median SP is used, i.e.,
%                    ps = [-30:5:70,80,100,110:5:210] degrees.
%
%     'mrsmsp',mrs   Set the motoric response scatter mrs within the median 
%                    sagittal plane. Default value is 17° in accordance
%                    with scatter of unimodal response distribution
%                    proposed in Langendijk and Bronkhorst (2002).
%
%   `baumgartner2014` accepts the following flags:
%
%     'gammatone'    Use the Gammatone filterbank for peripheral processing. 
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
%   See also: plotbaumgartner2013, data_baumgartner2014
%
%   References: baumgartner2013assessment baumgartner2012modelling langendijk2002contribution patterson1988efficient dau1996qmeI

    
% AUTHOR: Robert Baumgartner, Acoustics Research Institute, Vienna, Austria

%% Check input options 

definput.flags.fbank = {'gammatone','cqdft','drnl','zilany2007humanized','zilany5'};   % disp('zilany2007humanized used!')
definput.flags.headphonefilter = {'','headphone'};
definput.flags.middleearfilter = {'','middleear'};
definput.flags.ihc = {'noihc','ihc'};    
definput.flags.Ifw = {'nointensityweighting','intensityweighting'};
definput.flags.regularization = {'regular','noregular'};
definput.flags.motoricresponsescatter = {'mrs','nomrs'};
definput.flags.sensitivitymapping = {'sigmoidmapping','normstdmapping'};
definput.flags.settings = {'notprint','print'};

% CP-Falgs:
definput.flags.cp={'fwstd','std','xcorr'};

definput.keyvals.fs=48000;      % Hz
definput.keyvals.S=0.5;         % listener-specific sensitivity parameter
definput.keyvals.lat=0;         % deg
definput.keyvals.stim=[];
definput.keyvals.fsstim=[];
definput.keyvals.space=1;       % No. of ERBs (Cams) 
definput.keyvals.do=1;
definput.keyvals.flow=700;      % Hz
definput.keyvals.fhigh=18000;   % Hz
definput.keyvals.lvltar = 50; 	% dBSPL
definput.keyvals.lvltem = 60;  	% dBSPL
definput.keyvals.SL = [];       % db/ERB; spectral density of target sound re absolut detection threshold
definput.keyvals.conGain=13;    % steepness in degrees of binaural weighting function
definput.keyvals.polsamp=[-30:5:70 80 100 110:5:210];  % polar sampling (for regular)
definput.keyvals.mrsmsp=20;     % degrees
definput.keyvals.gamma=6;       % slope of psychometric function

[flags,kv]=ltfatarghelper(...
  {'fs','S','lat','stim','fsstim','space','do','flow','fhigh',...
  'lvltar','lvltem','SL','conGain','polsamp','mrsmsp','gamma'},definput,varargin);

% Print Settings
if flags.do_print 
  if flags.do_mrs
    fprintf('Settings: DCN = %1.0f; Gamma = %1.0u; Epsilon = %1.0f deg \n',kv.do,kv.gamma,kv.mrsmsp)
  else
    fprintf('Settings: DCN = %1.0f; Gamma = %1.0u; Epsilon = 0 deg \n',kv.do,kv.gamma)
  end
end

% HRTF format
if isstruct(target) % Targets given in SOFA format
  kv.fs = target.Data.SamplingRate;
  [target,tang] = extractsp( kv.lat,target );
end

if isstruct(template) % Template given in SOFA format
  [template,kv.polsamp] = extractsp( kv.lat,template );
end


%% Error handling
if size(template,2) ~= length(kv.polsamp)
  fprintf('\n Error: Second dimension of template and length of polsamp need to be of the same size! \n')
  return
end

if kv.S <= 0
  fprintf('\n Error: Listener-specific uncertainty has to be larger than zero! \n')
  return
end


%% Stimulus 
if isempty(kv.stim) 
    kv.stim = [1;0];%[1;zeros(size(target,1),1)];    % impulse
    kv.fsstim = kv.fs;
elseif isempty(kv.fsstim) 
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


%% Set level
idnztar = target~=0;    % to ignore pausings
idnztem = template~=0;  % to ignore pausings
% aht = setdbspl([1;0],kv.lvltar-kv.SL);  % absolut hearing threshold
for ch = 1:size(template,3)
    target(idnztar(:,:,ch)) = setdbspl(target(idnztar(:,:,ch)),kv.lvltar);
% aht(idnztar(:,:,ch)) = setdbspl(target(idnztar(:,:,ch)),kv.lvltar-kv.SL);
    template(idnztem(:,:,ch)) = setdbspl(template(idnztem(:,:,ch)),kv.lvltem);
end
    
    
%% Cochlear filter bank -> internal representations

if kv.space == 1
  [ireptar,fc] = auditoryfilterbank(target(:,:),kv.fs,...
      'flow',kv.flow,'fhigh',kv.fhigh);
  ireptem = auditoryfilterbank(template(:,:),kv.fs,...
      'flow',kv.flow,'fhigh',kv.fhigh);
else
  fc = audspacebw(kv.flow,kv.fhigh,kv.space,'erb');
  [bgt,agt] = gammatone(fc,kv.fs,'complex');
  ireptar = 2*real(ufilterbankz(bgt,agt,target(:,:)));  % channel (3rd) dimension resolved!
  ireptem = 2*real(ufilterbankz(bgt,agt,template(:,:)));
end
Nfc = length(fc);   % # of bands

% Set back the channel dimension
ireptar = reshape(ireptar,[size(target,1),Nfc,size(target,2),size(target,3)]);
ireptem = reshape(ireptem,[size(template,1),Nfc,size(template,2),size(template,3)]);

% Averaging over time (RMS)
ireptar = 20*log10(squeeze(rms(ireptar)));      % in dB!
ireptem = 20*log10(squeeze(rms(ireptem)));

if size(ireptar,2) ~= size(target,2) % retreive polar dimension if squeezed out
    ireptar = reshape(ireptar,[size(ireptar,1),size(target,2),size(target,3)]);
end


%% Comparison process -> monaural similarity indices (SIs)

si=zeros(size(ireptem,2),size(ireptar,2),size(ireptem,3)); % initialisation
for ch = 1:size(ireptar,3)

  if kv.do == 1 % DCN model
      nrep.tem = dcn(ireptem(:,:,ch),kv);
      nrep.tar = dcn(ireptar(:,:,ch),kv);
  elseif kv.do == 2 
      nrep.tem = diff(ireptem(:,:,ch),kv.do);
      nrep.tar = diff(ireptar(:,:,ch),kv.do);
  else
      nrep.tem = ireptem(:,:,ch);
      nrep.tar = ireptar(:,:,ch);
  end

  for it = 1:size(ireptar,2)

    % Distance Metric
    isd = repmat(nrep.tar(:,it),[1,size(nrep.tem,2),1]) - nrep.tem; 
    if kv.do == 0
      sigma = sqrt(squeeze(var(isd)));
    else
%         pnorm = 1;
%         sigma = sum( abs(isd).^pnorm .* repmat(fw,1,length(kv.polsamp)) ).^(1/pnorm);
      sigma = mean(abs(isd));
    end

    % Similarity Percept
    if not(exist('flags','var')) || flags.do_normstdmapping
      si(:,it,ch) = normpdf(sigma,0,kv.S);
    else % flags.do_sigmoidmapping
      si(:,it,ch) = 1+eps - (1+exp(-kv.gamma*(sigma-kv.S))).^-1;
    end
  end
end


%% Binaural weighting -> binaural SIs
if size(si,3) == 2
    binw = 1./(1+exp(-kv.lat/kv.conGain)); % weight of left ear signal with 0 <= binw <= 1
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


%% Motoric response scatter
if flags.do_mrs && flags.do_regular && kv.mrsmsp > 0
  
    angbelow = -90:5:min(respangs)-5;
    angabove = max(respangs)+5:5:265;
    respangs = [angbelow,respangs,angabove];
    si = [zeros(length(angbelow),size(si,2)) ; si ; zeros(length(angabove),size(si,2))];
    
    mrs = kv.mrsmsp/cos(deg2rad(kv.lat)); % direction dependent scatter (derivation: const. length rel. to the circumferences of circles considered as cross sections of a unit sphere)
    
    x = 0:2*pi/72:2*pi-2*pi/72;
    kappa = 1/deg2rad(mrs)^2; % concentration parameter (~1/sigma^2 of normpdf)
    mrspdf = exp(kappa*cos(x)) / (2*pi*besseli(0,kappa)); % von Mises PDF 
    for tt = 1:size(si,2)
      %si(:,tt) = circonv(si(:,tt),mrspdf,360/5);
      si(:,tt) = pconv(si(:,tt),mrspdf(:));
    end
    
end


%% Normalization to PMV
p = si ./ repmat(sum(si),size(si,1),1);


%% Output
varargout{1} = p;
if nargout >= 2
    varargout{2} = respangs;
    if nargout >= 3
        varargout{3} = tang;
    end
end
  
end

function t4 = dcn(an,kv)
%DCN Phenomenological model of dorsal cochlear nucleus (DCN)
%   Usage:      out = dcn(in)
%
%   Input parameters:
%     an      : spectral profile in dB
%
%   Output parameters:
%     t4      : activity of type IV unit

%% Parameter Settings
c2 = 1; % inhibitory coupling between type II and type IV neurons
c4 = 1; % coupling between an and type IV neuron
dilatation = 1; % of tonotopical 1-ERB-spacing between type IV and II neurons

%% Calculations
Nb = size(an,1); % # auditory bands
dt4t2 = round(dilatation/kv.space); % tonotopical distance between type IV and II neurons
t4 = zeros(Nb-dt4t2,size(an,2),size(an,3)); % type IV output
for b = 1:Nb-dt4t2
  t4(b,:,:) = c4 * an(b+dt4t2,:,:) - c2 * an(b,:,:);
end

t4 = max(t4,0); %disp('only rising edges')
end