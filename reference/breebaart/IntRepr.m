function IR = IntRepr(x, s)
% ==========================================================
% IntRepr - Derive internal representation of input signal x
%
% function IR = IntRepr(x, s);
%
%
% Version 1.5 : 23 Mars 2006 : Auditory Channel now selected here. no need
% to process the unwanted ones.
% Version 1.6 : 13 Avril 2006 : Inlcude the lenght check before entering
% gammatone filterbank. Throwing away the extra length afterwards does not
% seem to be optimal though.
% Version 1.7 : 28 Septembre 2006: Include the middle ear filter.
% Version 1.8 : 21 Novembre 2006: Include delay compensation of the fftfilt
%
% Nicolas Le Goff - n.legoff@tm.tue.nl

%% Get the left and right signal to feed to the model for this interval.
xleft  = x(:,1);
xright = x(:,2);

%% send signals to middle ear
if s.ModelParameters.Middleear == 1;
    % require initialisastion of middle ear.
    % and compensation for the delay
    xleftpad  = [xleft; zeros(1024,1)];
    xrightpad = [xright; zeros(1024,1)];
    xleftpadout  = fftfilt(s.ModelParameters.middlearcoefs,xleftpad);
    xrightpadout = fftfilt(s.ModelParameters.middlearcoefs,xrightpad);
    xleftout  = xleftpadout(1025:end);
    xrightout = xrightpadout(1025:end);
else
    xleftout =xleft;
    xrightout =xright;
end

%% Feed signals through gammatone filterbank:
xleftgt  = real( gtfbank( xleftout, s.fs, s.ModelParameters.ERBnumber, s.ModelParameters.ERBspace));
xrightgt = real( gtfbank(xrightout, s.fs, s.ModelParameters.ERBnumber, s.ModelParameters.ERBspace));   


%% Feed signals through inner haircell model:
xlefth  = zeros(round(s.maskerlength*s.fs/1000), length(s.ModelParameters.ChannelToLookAt));
xrighth = zeros(round(s.maskerlength*s.fs/1000), length(s.ModelParameters.ChannelToLookAt));
for idx = 1 : length(s.ModelParameters.ChannelToLookAt),
    xlefth(:,idx)  = ihc( xleftgt (:,s.ModelParameters.ChannelToLookAt(idx)), s.fs);
    xrighth(:,idx) = ihc( xrightgt(:,s.ModelParameters.ChannelToLookAt(idx)), s.fs);
end

%% Feed signals through adaptation loops:
xleftb = zeros(round(s.maskerlength*s.fs/1000), length(s.ModelParameters.ChannelToLookAt));
xrightb = zeros(round(s.maskerlength*s.fs/1000), length(s.ModelParameters.ChannelToLookAt));
for idx = 1 : length(s.ModelParameters.ChannelToLookAt),
    xleftb(:,idx) = fadapt( xlefth(:,idx), s.fs); % 5,128.75,252.5,376.25,500
    xrightb(:,idx) = fadapt( xrighth(:,idx), s.fs);

end


%% remove the undershoot on offset of the adaptation loops

if s.ModelParamaters.AllAbovezeroSwitch ==1
    xlefta  = zeros(round(s.maskerlength*s.fs/1000), length(s.ModelParameters.ChannelToLookAt));
    xrighta = zeros(round(s.maskerlength*s.fs/1000), length(s.ModelParameters.ChannelToLookAt));
    for idx = 1 : length(s.ModelParameters.ChannelToLookAt),
        xlefta(:,idx)  = AllAboveZero(xleftb(:,idx));
        xrighta(:,idx) = AllAboveZero(xrightb(:,idx));
    end    
else
    xlefta=xleftb;
    xrighta=xrightb;
end

%keyboard
%% Apply dynamic-range compression
if s.ModelParameters.CompressorSwitch == 1;
    %fprintf('toto')
    xleftcomp  = zeros(round(s.maskerlength*s.fs/1000), length(s.ModelParameters.ChannelToLookAt));
    xrightcomp = zeros(round(s.maskerlength*s.fs/1000), length(s.ModelParameters.ChannelToLookAt));
    for idx = 1 : length(s.ModelParameters.ChannelToLookAt),
        xleftcomp(:,idx)  = compressor( xlefta(:,idx), s.ModelParameters.CompressionThreshold, s.ModelParameters.CompressionFactor);
        xrightcomp(:,idx) = compressor(xrighta(:,idx), s.ModelParameters.CompressionThreshold, s.ModelParameters.CompressionFactor);
    end
else
    xleftcomp = xlefta;
    xrightcomp = xrighta;
end


% Return internal representation:
IR(1,:,:) = xleftcomp;
IR(2,:,:) = xrightcomp;


%OLDFORMAT
