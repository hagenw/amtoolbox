function [ri,rang] = baumgartner2014sensorimotormapping(si,varargin)
%BAUMGARTNER2014SENSORIMOTORMAPPING - Response scatter induced by
%localization task
%   Usage:     [ri,rang] = baumgartner2014sensorimotormapping(si)
%              [ri,rang] = baumgartner2014sensorimotormapping(si,flags,kv)
%
%   Input parameters:
%     ri      : response index
%
%   Output parameters:
%     si      : similarity index
%     rang    : response polar angles
%
%   `baumgartner2014sensorimotormapping(...)` performs polar-angle
%   interpolation and emulates task-induced response scatter.
%
%   `baumgartner2014sensorimotormapping` accepts the following optional parameters:
%
%     'flags',flags  Transfer flags. If not defaults of baumgartner2014 are used.
%
%     'kv',kv        Transfer key-value pairs. If not defaults of baumgartner2014 are used.
%
%   References: baumgartner2014modeling

% AUTHOR: Robert Baumgartner

definput.import={'baumgartner2014'};

if length(varargin) == 1 && isstruct(varargin{1})
  kv = varargin{1};
  flags=ltfatarghelper({},definput,varargin);
elseif length(varargin) == 2 && isstruct(varargin{1}) && isstruct(varargin{2})
  kv = varargin{1};
  flags = varargin{2};
else
  [flags,kv]=ltfatarghelper({},definput,varargin);
end

%% Interpolation (regularize polar angular sampling)
if flags.do_regular
    rang0 = ceil(min(kv.polsamp)*1/kv.rangsamp)*kv.rangsamp;    % ceil to kv.rangsamp deg
    rang = rang0:kv.rangsamp:max(kv.polsamp);
    siint = zeros(length(rang),size(si,2));
    for tt = 1:size(si,2)
        siint(:,tt) = interp1(kv.polsamp,si(:,tt),rang,'spline');
    end
    si = siint;
    si(si<0) = 0; % SIs must be positive (necessary due to spline interp)
else
    rang = kv.polsamp;
end


%% Sensorimotor mapping, Eq.(7)
if flags.do_mrs && flags.do_regular && kv.mrsmsp > 0
  
  angbelow = -90:kv.rangsamp:min(rang)-kv.rangsamp;
  angabove = max(rang)+kv.rangsamp:kv.rangsamp:265;
  rang = [angbelow,rang,angabove];
  ri = [zeros(length(angbelow),size(si,2)) ; si ; zeros(length(angabove),size(si,2))];

  mrs = kv.mrsmsp/cos(deg2rad(kv.lat)); % direction dependent scatter (derivation: const. length rel. to the circumferences of circles considered as cross sections of a unit sphere)

  x = 0:2*pi/length(rang):2*pi-2*pi/length(rang);
  kappa = 1/deg2rad(mrs)^2; % concentration parameter (~1/sigma^2 of normpdf)
  mrspdf = exp(kappa*cos(x)) / (2*pi*besseli(0,kappa)); % von Mises PDF 
  for tt = 1:size(ri,2)
    ri(:,tt) = pconv(ri(:,tt),mrspdf(:));
  end
    
else
  ri = si;
end
  

end