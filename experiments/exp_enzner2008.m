function exp_enzner2008(varargin) 
%EXP_ENZNER2008  Creates figures like [Enzner2008, Fig. 2], [Enzner2009, Fig. 4]
%   Usage: exp_enzner2008(flag)
%
%   The following flags can be specified:
%
%     'fig2_2008'   Plot Fig. 2 from Enzner et al. (2008)
%
%     'fig4_2009'   Plot Fig. 4 from Enzner et al. (2009)
%
%   Examples:
%   ---------
%
%   To display Figure 2 from the 2008 paper use :::
%
%     exp_enzner2008('fig2_2008');
%
%   To display Figure 4 from the 2009 paper use :::
%
%     exp_enzner2008('fig4_2009');
%
%   See also: enzner2008

  
%  Authors: Michael Weinert (Michael.Weinert@rub.de), Gerald Enzner (Gerald.Enzner@rub.de)
%  Date: 21-01-2013

%addpath(fullfile(amtbasepath,'hrtf','continuous-azimuth HRIR'))

enzner2008(1,1,varargin); % enzner2008(mu, delta_phi, varargin)


