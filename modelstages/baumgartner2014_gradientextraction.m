function [gp,gfc] = baumgartner2014_gradientextraction(mp,fc,varargin)
%baumgartner2014_gradientextraction - Extraction of positive spectral gradients
%   Usage:      [gp,gfc] = baumgartner2014_gradientextraction(mp,fc)
%
%   Input parameters:
%     mp      : spectral magnitude profile in dB
%     fc      : center frequencies
%
%   Output parameters:
%     gp      : positive spectral gradient profile
%     gfc     : center frequencies of gradient profile
%
%   `baumgartner2014_gradientextraction(...)` is a spectral cue extractor
%    inspired by functionality of dorsal cochlear nucleus in cats.
%
%   `baumgartner2014_gradientextraction` accepts the following optional parameters:
%
%     'c2',c2   Inhibitory coupling between type II mpd type IV neurons. 
%               Default is 1.
%     'c4',c4   Inhibitory coupling between AN and type IV neuron. 
%               Default is 1.
%     'spacing',sp  Tonotopical spacing between type IV mpd II neurons in ERBs. 
%                   Default is 1 ERB.
%
%   `baumgartner2014_gradientextraction` accepts the following flags:
%
%     'positive'  Perform positive spectral gradient extraction (default).
%     'negative'  Perform negative spectral gradient extraction.
%     'both'      Perform spectral gradient extraction.
%
%   References: baumgartner2014modeling

% AUTHOR: Robert Baumgartner

%% Parameter Settings

definput.keyvals.c2 = 1;
definput.keyvals.c4 = 1;
definput.keyvals.spacing = 1;
definput.flags.gradient = {'positive','negative','both'};

[flags,kv]=ltfatarghelper({'c2','c4','spacing'},definput,varargin);

% if not(exist('c2','var'))
%   c2 = 1; % inhibitory coupling between type II mpd type IV neurons
% end
% c4 = 1; % coupling between AN and type IV neuron
% dilatation = 1; % of tonotopical 1-ERB-spacing between type IV mpd II neurons


%% Calculations
Nb = size(mp,1); % # auditory bands
erb = audfiltbw(fc);
dgpt2 = round(mean(erb(2:end)./diff(fc))*kv.spacing); % tonotopical distance between type IV mpd II neurons
gp = zeros(Nb-dgpt2,size(mp,2),size(mp,3),size(mp,4),size(mp,5)); % type IV output
for b = 1:Nb-dgpt2
  gp(b,:,:,:,:) = kv.c4 * mp(b+dgpt2,:,:,:,:) - kv.c2 * mp(b,:,:,:,:);
end

if flags.do_positive
  gp(gp<0) = 0; % gp = (gp + c2*abs(gp))/2;
elseif flags.do_negative
  gp(gp>0) = 0;
end

% gfc = fc(dgpt2+1:end); % same as AN
gfc = sqrt(fc(1:Nb-dgpt2).*fc(dgpt2+1:end)); % use geometric mean

end