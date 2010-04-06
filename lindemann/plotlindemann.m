function plotlindemann(crosscorr,t,fc)

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
