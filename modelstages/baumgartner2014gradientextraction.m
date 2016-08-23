function [gp,gfc] = baumgartner2014gradientextraction(mp,fc,c2)
%BAUMGARTNER2014GRADIENTEXTRACTION - Extraction of positive spectral gradients
%   Usage:      [gp,gfc] = baumgartner2014gradientextraction(mp,fc)
%
%   Input parameters:
%     mp      : spectral magnitude profile in dB
%     fc      : center frequencies
%
%   Output parameters:
%     gp      : positive spectral gradient profile
%     gfc     : center frequencies of gradient profile
%
%   `baumgartner2014gradientextraction(...)` is a spectral cue extractor
%    inspired by functionality of dorsal cochlear nucleus in cats.
%
%   References: baumgartner2014modeling

% AUTHOR: Robert Baumgartner

%% Parameter Settings
if not(exist('c2','var'))
  c2 = 1; % inhibitory coupling between type II mpd type IV neurons
end
c4 = 1; % coupling between AN and type IV neuron
dilatation = 1; % of tonotopical 1-ERB-spacing between type IV mpd II neurons

erb = audfiltbw(fc);

%% Calculations
Nb = size(mp,1); % # auditory bands
dgpt2 = round(mean(erb(2:end)./diff(fc))*dilatation); % tonotopical distance between type IV mpd II neurons
gp = zeros(Nb-dgpt2,size(mp,2),size(mp,3),size(mp,4),size(mp,5)); % type IV output
for b = 1:Nb-dgpt2
  gp(b,:,:,:,:) = c4 * mp(b+dgpt2,:,:,:,:) - c2 * mp(b,:,:,:,:);
end

gp = (gp + c2*abs(gp))/2; %disp('only rising edges')
% gp(gp<0) = nmp;

gfc = fc(dgpt2+1:end);

% if nargout > 1

end