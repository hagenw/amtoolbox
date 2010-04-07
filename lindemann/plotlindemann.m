function plotlindemann(crosscorr,t,fc)
%PLOTLINDEMANN Plots the binaural output pattern of the lindemann model
%   Usage: plotlindemann(crosscorr,t,fc)
%
%   Input parameters:
%       crosscorr   - cross-correlation matrix, output from the lindemann
%                     function
%       t           - time of the analysed stimuli (used for t axis)
%       fc          - frequency channel to plot. Default: mean about all
%                     channels
%
%   PLOTLINDEMANN(crosscorr,t,fc) plots the cross-correlation output from the
%   lindemann function as a so called binaural activity map. This mean the
%   correlation value is plotted dependend on time of the stimulus and the
%   correlation-time delay.
%
%   see also: lindemann, bincorr
%

%   AUTHOR: Hagen Wierstorf


% ------ Checking of input  parameters -----------------------------------

error(nargchk(2,3,nargin));

if ~isnumeric(crosscorr)
    error('%s: crosscorr has to be numeric!',upper(mfilename));
end

if ~isnumeric(t) || ~isscalar(t) || t<=0
    error('%s: t has to be a positive scalar!',upper(mfilename));
end

if nargin>2 && ( ~isnumeric(fc) || ~isscalar(fc) || fc<=0 )
    error('%s: fc has to be a positive scalar');
end

% ------ Computation -----------------------------------------------------
    
% Calculate tau (delay line time) axes
tau = linspace(-1,1,size(crosscorr,2));
% Calculate t axes
t = linspace(0,t,size(crosscorr,1));
% Calculate mean binaural activation pattern
if nargin==2
    binpattern = mean(crosscorr,3);
else
    binpattern = crosscorr(:,:,fc);
end

% plot result
figure;
mesh(tau,t,binpattern);
xlabel('correlation-time delay (ms)');
ylabel('t (ms)');
