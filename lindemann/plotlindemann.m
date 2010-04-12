function plotlindemann(crosscorr,t,varargin)
%PLOTLINDEMANN Plots the binaural output pattern of the lindemann model
%   Usage: plotlindemann(crosscorr,t,f,tstr)
%          plotlindemann(crosscorr,t,f)
%          plotlindemann(crosscorr,t,tstr)
%          plotlindemann(crosscorr,t)
%
%   Input parameters:
%       crosscorr   - cross-correlation matrix, output from the lindemann
%                     function
%       t           - time of the analysed stimuli (used for t axis)
%       f           - plot only the frequency channel with its center frequency
%                     is nearest to the frequency f. Default: mean about all
%                     frequency channels
%       tstr        - title string for the plot. Default: ''
%
%   PLOTLINDEMANN(crosscorr,t,f,tstr) plots the cross-correlation output 
%   from the lindemann function as a so called binaural activity map. This
%   means the correlation value is plotted dependend on time of the stimulus
%   and the correlation-time delay. t is the max value of the time axis of the
%   plot and tstr the title of the plot. f determines the frequency channel to
%   plot by using the channel in which the frequency f belongs.
%
%   see also: lindemann, bincorr
%

%   AUTHOR: Hagen Wierstorf


% ------ Checking of input  parameters -----------------------------------

error(nargchk(2,4,nargin));

if ~isnumeric(crosscorr)
    error('%s: crosscorr has to be numeric!',upper(mfilename));
end

if ~isnumeric(t) || ~isscalar(t) || t<=0
    error('%s: t has to be a positive scalar!',upper(mfilename));
end

% Parsing of varargin (fc, tstr)
% FIXME: can this be done nicer? Can I also use amtarghelper?
for ii = 1:nargin-2
    if isnumeric(varargin{ii}) && ~exist('f')
        f = varargin{ii};
        % Minimum and maximum frequency in the lindemann model (see lindemann.m)
        flow = erbtofreq(5);
        fhigh = erbtofreq(40);
        if ~isscalar(f)
            error('%s: f has to be a scalar!',upper(mfilename));
        elseif f<flow || f>fhigh
            error('%s: f has to be between %.0f Hz and %.0f Hz!',...
                upper(mfilename),flow,fhigh);
        end
    elseif ischar(varargin{ii}) && ~exist('tstr')
        tstr = varargin{ii};
    else
        error(['%s: the optional parameters have to be numeric (f) or/and ',...
            'a string (tstr)!'],upper(mfilename));
    end
end


% ------ Computation -----------------------------------------------------
    
% Calculate tau (delay line time) axes
tau = linspace(-1,1,size(crosscorr,2));
% Calculate t axes
t = linspace(0,t,size(crosscorr,1));
% Calculate mean binaural activation pattern
if ~exist('f')
    binpattern = mean(crosscorr,3);
else
    % Calculate the frequency channel to plot
    % NOTE: it starts with the fifth channel in the lindemann model, so we have
    % to subtract 4 to index the binpattern correctly.
    fc = round(freqtoerb(f));
    binpattern = crosscorr(:,:,fc-4);
end

% ------ Plotting --------------------------------------------------------
figure;
mesh(tau,t,binpattern);
xlabel('correlation-time delay (ms)');
ylabel('t (ms)');
% Create title, if fc is given but not tstr
if ~exist('tstr') && exist('fc')
    tstr = sprintf('fc = %i',fc);
end
% Plot title
if exist('tstr')
    title(tstr);
end
