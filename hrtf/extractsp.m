function varargout = extractsp(lat,varargin)
%EXTRACTSP Sagittal plane (SP) HRTFs from measurement data
%   Usage:    [sphrtfs,polangs] = extractsp( lat,hM,pos )
%             [sphrtfs,polangs] = extractsp( lat,Obj )
%
%   Input parameters:
%     lat     : lateral angle of the SP
%     Obj     : HRIR Data in SOFA format
%     hM      : matrix containing head-related impulse responses in ARI
%               Format
%               Dimensions: time,position,channel 
%               (for more details see doc: HRTF format description)
%     pos     : source-position matrix referring to 2nd dimension of hM and 
%               formated acc. to meta.pos (ARI format).
%               6th col: lateral angle. 7th col: polar angle
%
%   Output parameters:
%     sphrtfs : all available HRTFs in the current SP, sorted acc. to
%               ascending polar angle
%     polangs : corresponding polar angles
%
%   `extractsp(...)` extracts all HRTFs available for a specific SP or
%   lateral angle. In order to result in a representative HRTF template,
%   demands are made on:
%
%     1) lateral tolerance to be as small as possible, but min. 2 and
%        max. 5.
% 
%     2) polar angles to range from max. -30 deg. to min. 210 deg.,
%
%     3) gaps of polar angles to be max. 30 deg. large.
    
% AUTHOR: Robert Baumgartner, Acoustics Research Institute, Vienna, Austria

%% Check input options 

if nargin == 2 || isstruct(varargin{1}) % SOFA input
  
  Obj = varargin{1};
	hM = permute(double(Obj.Data.IR),[3 1 2]);

  pos(:,1)=bsxfun(@times,Obj.SourcePosition(:,1),ones(Obj.API.M,1));
  pos(:,2)=bsxfun(@times,Obj.SourcePosition(:,2),ones(Obj.API.M,1));
  [pos(:,6), pos(:,7)]=sph2hor(pos(:,1),pos(:,2));
  
elseif nargin == 3 % aRI Format
  
  hM = varargin{1};
  pos = varargin{2};
  
else return
end
  

%% Extract SP

dlat = 2;      % initial lateral tolerance (+/-) in deg
pol = [0,0];   % initial polar angles

while (min(pol) > -30 || max(pol) < 210 ... % ensure that important polar range is included
        || max(diff(pol))>30)...            % and gaps are <= 30°
        && dlat <= 5                        % but keep tolerance smaller than 5 deg

  idx=find(pos(:,6)>=-(dlat+0.01)/2+lat & pos(:,6)<=(dlat+0.01)/2+lat);
  latActual = pos(idx,6); % actual lateral angles
  [pol,polidx]=unique(real(pos(idx,7)));   % sorted polar angles
  latActual = latActual(polidx);
  sphrtfs=double(hM(:,idx,:));  % unsorted DTFs of SP
  sphrtfs=sphrtfs(:,polidx,:);          % sorted DTFs of SP

  dlat = dlat + 1;  % increase dlat

end

varargout{1} = sphrtfs;
if nargout >= 2
  varargout{2} = pol;
  if nargout == 3
    varargout{3} = latActual;
  end
end

end