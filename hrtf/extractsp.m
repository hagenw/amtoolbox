function varargout = extractsp(lat,hM,pos)
% EXTRACTSP extracts sagittal plane (SP) HRTFs from measurement data stored
% in ARI's HRTF format
%
% Usage:    [sphrtfs,polangs] = extractsp( lat,hM,pos )
%
% Input parameters:
%     lat     : lateral angle of the SP
%     hM      : matrix containing head-related impulse responses.
%               Dimensions: time,position,channel 
%               (for more details see doc: HRTF format description)
%     pos     : source-position matrix referring to 2nd dimension of hM and 
%               formated acc. to meta.pos (ARI format).
%               6th col: lateral angle. 7th col: polar angle
%
% Output parameters:
%     sphrtfs : all available HRTFs in the current SP, sorted acc. to
%               ascending polar angle
%     polangs : corresponding polar angles
%
%   `extractsp(...)` extracts all HRTFs available for a specific SP or
%   lateral angle. In order to result in a representative HRTF template,
%   demands are made on:
%   1) lateral tolerance to be as small as possible, but min. 2 and max. 5,
%   2) polar angles to range from max. -30° to min. 210°,
%   3) gaps of polar angles to be max. 30° large.
%
% AUTHOR: Robert Baumgartner, Acoustics Research Institute, Vienna, Austria

dlat = 2;      % initial lateral tolerance (+/-) in deg
pol = [0,0];   % initial polar angles

while (min(pol) >= -30 || max(pol) <= 210 ... % ensure that important polar range is included
        || max(diff(pol))>30)...    % and gaps are <= 30°
        && dlat <= 5                % but keep tolerance smaller than 5 deg

  idx=find(pos(:,6)>=-(dlat+0.01)/2+lat & pos(:,6)<=(dlat+0.01)/2+lat);
  [pol,polidx]=unique(real(pos(idx,7)));   % sorted polar angles
  sphrtfs=double(hM(:,idx,:));  % unsorted DTFs of SP
  sphrtfs=sphrtfs(:,polidx,:);          % sorted DTFs of SP

  dlat = dlat + 1;  % increase dlat

end

varargout{1} = sphrtfs;
if nargout == 2
  varargout{2} = pol;
end

end