function plotlindemann(crosscorr,t,varargin)
%PLOTLINDEMANN Plots the binaural output pattern of the lindemann model
%   Usage: plotlindemann(crosscorr,t,f,tstr);
%          plotlindemann(crosscorr,t,f);
%          plotlindemann(crosscorr,t,tstr);
%          plotlindemann(crosscorr,t);
%
%   Input parameters:
%       crosscorr   - cross-correlation matrix, output from the lindemann
%                     function
%       t           - time vector of the analysed stimuli (used for t axis)
%
%   PLOTLINDEMANN(crosscorr,t) plots the cross-correlation output from the
%   lindemann function as a so called binaural activity map. This means the
%   correlation value is plotted dependending on time of the stimulus and
%   the correlation-time delay. t is the time axis of the plot. f determines
%   the frequency channel to plot by using the channel in which the
%   frequency f belongs.
%
%   If crosscorr has more than one time step a 3D activity map is plotted, else
%   a 2D plot of the cross-correlation is done.
%
%   The function takes the following flags at the end of the line of
%   input arguments:
%
%-      'fc',fc - plot only the frequency channel with its center frequency
%                   is nearest to the frequency f. The default value of []
%                   means to plot the mean about all frequency channels
%
%-      'title',t - display t as the title overriding the default.
%
%   You may also supply the parameters in the input arguments in the
%   following order: PLOTLINDEMANN(crosscorr,t,fc);
%  
%   See also: lindemann, lindemannbincorr
%

%   AUTHOR: Hagen Wierstorf


% ------ Checking of input  parameters -----------------------------------

if nargin<2
  error('%s: Too few input arguments.',upper(mfilename));
end;

if ~isnumeric(crosscorr)
    error('%s: crosscorr has to be numeric!',upper(mfilename));
end

if ( ~isnumeric(t) || ~isvector(t) )
    error('%s: t has to be a vector!',upper(mfilename));
end

definput.keyvals.title=[];
definput.keyvals.fc=[];

[flags,keyvals]  = ltfatarghelper({'fc','title'},definput,varargin);

if isempty(keyvals.fc)
  binpattern = mean(crosscorr,3);
else
  % Minimum and maximum frequency in the lindemann model (see lindemann.m)
  flow = erbtofreq(5);
  fhigh = erbtofreq(40);
  if ~isscalar(keyvals.fc)
    error('%s: fc has to be a scalar!',upper(mfilename));
  elseif keyvals.fc<flow || keyvals.fc>fhigh
    error('%s: fc has to be between %.0f Hz and %.0f Hz.',...
          upper(mfilename),flow,fhigh);
  end  

  % Calculate the frequency channel to plot
  % NOTE: it starts with the fifth channel in the lindemann model, so we have
  % to subtract 4 to index the binpattern correctly.
  fc = round(freqtoerb(keyvals.fc));
  binpattern = crosscorr(:,:,fc-4);

end;

% ------ Computation -----------------------------------------------------
    
% Calculate tau (delay line time) axes
tau = linspace(-1,1,size(crosscorr,2));

% ------ Plotting --------------------------------------------------------
if size(crosscorr,1)==1
    % If we have only one time step (stationary case) plot 2D
    plot(tau,binpattern);
else
    mesh(tau,t,binpattern);
    ylabel('t (s)');
end

xlabel('correlation-time delay (ms)');
% Create title, if fc is given but not tstr
if isempty(keyvals.title) && ~isempty(keyvals.fc)
    keyvals.title = sprintf('fc = %i',fc);
end

% Plot title
if ~isempty(keyvals.title)  
    title(keyvals.title);
end
