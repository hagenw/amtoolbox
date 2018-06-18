function [E,varargout] = baumgartner2017( target,template,varargin )
%BAUMGARTNER2017 Model for sound externalization
%   Usage:    [E,cues] = baumgartner2017( target,template )
%
%   Input parameters:
%     target  : binaural impulse response(s) referring to the directional 
%               transfer function(s) (DFTs) of the target sound(s).
%               Option 1: given in SOFA format -> sagittal plane DTFs will 
%               be extracted internally. 
%               Option 2: binaural impulse responses of all available
%               listener-specific DTFs of the sagittal plane formated 
%               according to the following matrix dimensions: 
%               time x direction x channel/ear
%     template: binaural impulse responses of all available
%               listener-specific DTFs of the sagittal plane referring to
%               the perceived lateral angle of the target sound.
%               Options 1 & 2 equivalent to *target*.
%
%   Output parameters:
%     E       : predicted degree of externalization
%     cues    : outcomes of individual cues 
%
%   `baumgartner2017(...)` is a model for sound externalization.
%   It bases on the comparison of the intra-aural internal representation
%   of the incoming sound with a template and results in a probabilistic
%   prediction of polar angle response.
%
%   `baumgartner2017` accepts the following optional parameters:
%
%     'cueWeights',cW  Set the weights of individual cues to determine the 
%                      final externalization score. Cue-specific weights 
%                      (entered as a vector) are ordered as follows: 
%                      1. monaural spectral similarity (c.f., Baumgartner et
%                         al., 2014). This is the default.
%                      2. interaural spectral similarity of ILDs (c.f.,
%                         Hassager et al., 2016)
%                      3. spectral standard deviation of ILDs (c.f.,
%                         Georganti et al., 2013)
%                      4. temporal standard deviation of ILDs (c.f., Catic
%                         et al., 2015)
%                      5. interaural broadband time-intensity coherence
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
%                    *lat*. Default value is 0 degree (median SP).
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
%     'do',do        Set the differential order of the spectral gradient 
%                    extraction to *do*. Default value is 1 and includes  
%                    restriction to positive gradients inspired by cat DCN
%                    functionality.
%
%     'bwcoef',bwc   Set the binaural weighting coefficient *bwc*.
%                    Default value is 13 degrees.
%
%     'range',c1     Set the range factor of the externalization scores to *c1*.
%                    Default value is 3.78 from Hassager et al. (2016).
%
%     'offset',c2    Set the offset of the externalization score to *c2*.
%                    Default value is 1 from Hassager et al. (2016).
%
%     'JND',JND      Set the just noticeable difference *JND* from the 
%                    internal template. Default value is 1.5 (dB).
%                    
%
%   Requirements: 
%   -------------
%
%   1) SOFA API from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
% 
%   2) Data in hrtf/baumgartner2017
%
%   3) Circular Statistics Toolbox from http://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox--directional-statistics-
%
%
%   See also: baumgartner2016_spectralanalysis,
%   baumgartner2016_gradientextraction, baumgartner2014_binauralweighting

% AUTHOR: Robert Baumgartner, Acoustics Research Institute, Vienna, Austria

%% Check input

definput.import={'baumgartner2017','baumgartner2014','baumgartner2014_pmv2ppp','localizationerror','amt_cache'};

[flags,kv]=ltfatarghelper(...
  {'fs','S','lat','stim','space','do','flow','fhigh',... %'fsstim'
  'bwcoef','polsamp','rangsamp','mrsmsp','gamma'},definput,varargin);

if not(isstruct(target)) && ismatrix(target)
  target = permute(target,[1,3,2]);
%   warning(['Matrix dimensions of target should be: time x direction x channel/ear.' ...
%    'Since 3rd dimension was empty, 2nd dimension was used as channel dimension.'])
end

if not(isstruct(template)) && ismatrix(template)
  template = permute(template,[1,3,2]);
%   warning(['Matrix dimensions of template should be: time x direction x channel/ear.' ... 
%     'Since 3rd dimension was empty, 2nd dimension was used as channel dimension.'])
end

%% Print Settings

if flags.do_print 
  if flags.do_nomrs
    kv.mrsmsp = 0;
  end
  amt_disp(['Settings: PSGE = ' num2str(kv.do,'%1.0f') '; Gamma = ' ...
    num2str(kv.gamma,'%1.0u') '; Epsilon = ' num2str(kv.mrsmsp,'%1.0f') ' deg'])
end


%% Determine lateral angle and extract HRTFs of sagittal plane
 
if isstruct(target) % Targets given in SOFA format
  kv.fs = target.Data.SamplingRate;
  [target,tang] = extractsp( kv.lat,target );
% else
%   fncache = ['latLookup_',template.GLOBAL_ListenerShortName];
%   latLookup = amt_cache('get',fncache,flags.cachemode);
%   if isempty(latLookup)
%     latLookup = itd2angle_lookuptable(template,template.Data.SamplingRate,'dietz2011');
%     amt_cache('set',fncache,latLookup)
%   end
%   tarSig = squeeze(target);
%   kv.lat = wierstorf2013_estimateazimuth(tarSig,latLookup,'fs',kv.fs,'dietz2011','rms_weighting');
%   disp(kv.lat)
end

if isstruct(template) % Template given in SOFA format
  [template,rang] = extractsp( kv.lat,template );
end

% Error handling
% if size(template,2) ~= length(rang)
%   fprintf('\n Error: Second dimension of template and length of polsamp need to be of the same size! \n')
%   return
% end

%% Optional: Middle ear filtering
if flags.do_middleEarFilter
  b=middleearfilter(kv.fs);
  target = filter(b,1,target);
  template = filter(b,1,template);
end

%% Optional: HRTF filtering

dimtar = size(target); % for lconv dim check

if not(isempty(kv.stim))
  target = lconv(target,kv.stim);
end

% check that lconv preserved matrix dimensions (earlier bug in lconv)
if size(target,2) ~= dimtar(2)
  target = reshape(target,[size(target,1),dimtar(2:end)]);
end

% frameLength = round((kv.tempWin*kv.fs));

%% ITD
tem.itd = itdestimator(shiftdim(template,1),'fs',kv.fs,'MaxIACCe');
tar.itd = itdestimator(shiftdim(target,1),'fs',kv.fs,'MaxIACCe');

%% Filterbank
[tem.mp,fc] = baumgartner2016_spectralanalysis(template,70,'tiwin',kv.tempWin,'GT_minSPL',-100,'gammatone','redo');
[tar.mp,fc] = baumgartner2016_spectralanalysis(target,70,'tiwin',kv.tempWin,'GT_minSPL',-100,'gammatone','redo');

%% Spectral cues
tem.psg = baumgartner2016_gradientextraction(tem.mp,fc);
tem.ild = -diff(tem.mp,1,3);
% ref.ild = mean(tem.ild,5);
tar.psg = baumgartner2016_gradientextraction(tar.mp,fc);
tar.ild = -diff(tar.mp,1,3);

%% spectral SD of ISLD (Georganti et al., 2013) -> dprime possible
% tar.sdild = std(tar.ild,0,1);
% tem.sdild = std(tem.ild,0,1);
% ref.sdild = std(ref.ild,0,1);
% ISLDspecSD = mean(tar.sdild./tem.sdild,5);

%% temporal SD of ISLD (Catic et al., 2015)
ISLDtemSD = 1 - mean(std(tar.ild,0,5)./std(tem.ild,0,5));

%% Interaural broadband time-intenstiy coherence (IBTIC) -> dprime possible
IBTIC = abs(tar.itd/tem.itd - mean(tar.ild(:))/mean(tem.ild(:)));

%% Spectral comparison

for iSC = 1:3 % first monaural then interaural
  if iSC == 1 % monaural spectral gradients
    tem.nrep = tem.psg.m;
    tar.nrep = tar.psg.m;
  elseif iSC == 2 % interaural spectral differences
    tem.nrep = tem.ild;
    tar.nrep = tar.ild;
  elseif iSC == 3 % spectral SD of interaural differences
    tem.nrep = std(tem.ild,0,1);
    tar.nrep = std(tar.ild,0,1);
  end
  
  % comparison with time average of spectral template
  nrep = {tar.nrep};
  if flags.do_dprime
    nrep{2} = tem.nrep;
  end
  tem.nrep = mean(tem.nrep,5);
  sigma = cell(length(nrep),1);
  for inrep = 1:length(nrep)  
    tmp = repmat(tem.nrep,[1,1,1,1,size(nrep{inrep},5)]);
    delta = abs(tmp-repmat(nrep{inrep},[1,size(tmp,2),1,1,1]));
    delta(delta < kv.JND) = 0; % limit minimum ILD difference according to JND
    delta = delta./(eps+abs(tmp)); % normalization (Weber fraction)
    sigma{inrep} = mean(delta); % average across frequency bands
    if iSC == 1 % do_intraaural
      sigma{inrep} = baumgartner2014_binauralweighting(sigma{inrep},'argimport',flags,kv);
    end
  end

  % temporal integration
  if length(sigma{1}) == 1
    si = exp(-kv.S*sigma{1});
  elseif flags.do_dprime % signal detection theory applied to time histograms
  %   figure; histogram(sigma{1}); hold on ; histogram(sigma{2}); legend('target','reference')
    allsigma = [sigma{1}(:);sigma{2}(:)];
    msigma = mean(allsigma);
    sdsigma = std(allsigma);
    mzsigma(1) = mean((sigma{1}-msigma) ./ sdsigma);
    mzsigma(2) = mean((sigma{2}-msigma) ./ sdsigma);
    dprime = max(mzsigma(1)-mzsigma(2),0);
    distmetric = dprime;
  else % temporal weighting according to amount of sensory information available
%     si = exp(-kv.S*sigma{1});
    % figure; plot(squeeze(bsi))
    tweight = mean(mean(abs(tar.nrep)),3); % temporal weighting
    tweight = tweight-min(tweight,[],5); % bounded between 0
    tweight = 2*tweight./max(tweight,[],5); % and 2
    distmetric = sigma{1}.*tweight;
    distmetric = mean(distmetric,5);
  end
  
  si = distmetric;%exp(-kv.S*distmetric);
  
  if iSC == 1
    MSS = si; % monaural spectral similarity
  elseif iSC == 2
    ISS = si; % interaural spectral similarity
  elseif iSC == 3
    ISLDspecSD = si;
  end
  
end

%% Cue integration/weighting
if flags.do_intraaural
  kv.cueWeights = 1;
elseif flags.do_interaural
  kv.cueWeights = [0,1];
end
cues = [MSS; ISS; ISLDspecSD; ISLDtemSD; IBTIC];
kv.S = postpad(kv.S(:),length(cues));
kv.cueWeights = postpad(kv.cueWeights(:),length(cues))/sum(kv.cueWeights);
si = exp(-kv.S.*cues);
bsi = nansum(kv.cueWeights .* si);

E = kv.range*bsi +kv.offset;%max(bsi);%min(1,max(bsi));%geomean(bsi);
if nargout >= 2
  varargout{1} = cues;
end

  
end