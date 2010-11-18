% Validation of Langendijk et al. (2002)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Robert Baumgartner, OEAW Acoustical Research Institute
% latest update: 2010-08-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               SETTINGS                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
listener='P6';  % ID of listener (P3 or P6)
fs = 48000;     % sampling frequency
bw = 06;        % bandwidth of averagingfb as partial of an octave
do = 0;         % differential order
cp = 'std';     % comparison process; 'std' (default) or 'xcorr' 
s  = 2;         % standard deviation of transforming Gaussian function
bal= 1;         % balance of left to right channel;   default: 1
fstart =2000;   % start frequency; minimum: 0,5kHz;   default: 2kHz
fend   =16000;  % end frequency;                      default: 16kHz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(['langendijk2002-' listener]); 
% loads hM data for all conditions which were calculated as follows:
% temp=data_langendijk2002([listener '-dtf']);
% pol=temp(1,:);
% med=temp(2:end,:);
% temp=data_langendijk2002([listener '-b']);
% targetb=temp(1,:); responseb=temp(2,:);
% medir=gr2ir(med,'b',fs);
% temp=data_langendijk2002([listener '-2o']);
% target2o=temp(1,:); response2o=temp(2,:);
% medir2o=gr2ir(med,'2o',fs);
% temp=data_langendijk2002([listener '-1ol']);
% target1ol=temp(1,:); response1ol=temp(2,:);
% medir1ol=gr2ir(med,'1ol',fs);
% temp=data_langendijk2002([listener '-1om']);
% target1om=temp(1,:); response1om=temp(2,:);
% medir1om=gr2ir(med,'1om',fs);
% temp=data_langendijk2002([listener '-1oh']);
% target1oh=temp(1,:); response1oh=temp(2,:);
% medir1oh=gr2ir(med,'1oh',fs);

% pdf calcualtion
h = waitbar(0,'Please wait...');
pb  = langendijk( medir   ,medir,bw,do,cp,s,bal,fstart,fend); % baseline
waitbar(1/5)
p2o = langendijk( medir2o ,medir,bw,do,cp,s,bal,fstart,fend); % 2-oct (4-16kHz)
waitbar(2/5)
p1ol= langendijk( medir1ol,medir,bw,do,cp,s,bal,fstart,fend); % 1-oct (low:4-8kHz)
waitbar(3/5)
p1om= langendijk( medir1om,medir,bw,do,cp,s,bal,fstart,fend); % 1-oct (middle:5.7-11.3kHz)
waitbar(4/5)
p1oh= langendijk( medir1oh,medir,bw,do,cp,s,bal,fstart,fend); % 1-oct (high:8-16kHz)
waitbar(5/5)

% likelihood estimations
la=zeros(5,1);le=zeros(5,1);ci=zeros(5,2);
idb=1:2:length(targetb); % in order to get comparable likelihoods
[la(1),le(1),ci(1,:)] = likelilangendijk( pb,pol,pol,targetb(idb),responseb(idb) );
[la(2),le(2),ci(2,:)] = likelilangendijk( p2o,pol,pol,targetc,response2o );
[la(3),le(3),ci(3,:)] = likelilangendijk( p1ol,pol,pol,targetc,response1ol );
[la(4),le(4),ci(4,:)] = likelilangendijk( p1om,pol,pol,targetc,response1om );
[la(5),le(5),ci(5,:)] = likelilangendijk( p1oh,pol,pol,targetc,response1oh );
close(h)

% pdf plots with actual responses
plotlangendijk(pb,pol,pol,[listener '; ' 'baseline']);
hold on; h=plot( targetb, responseb, 'ko'); set(h,'MarkerFaceColor','w')
plotlangendijk(p2o,pol,pol,[listener '; ' '2-oct (4-16kHz)']);
hold on; h=plot( targetc, response2o, 'ko'); set(h,'MarkerFaceColor','w')
plotlangendijk(p1ol,pol,pol,[listener '; ' '1-oct (low: 4-8kHz)']);
hold on; h=plot( targetc, response1ol, 'ko'); set(h,'MarkerFaceColor','w')
plotlangendijk(p1om,pol,pol,[listener '; ' '1-oct (middle: 5.7-11.3kHz)']);
hold on; h=plot( targetc, response1om, 'ko'); set(h,'MarkerFaceColor','w')
plotlangendijk(p1oh,pol,pol,[listener '; ' '1-oct (high: 8-16kHz)']);
hold on; h=plot( targetc, response1oh, 'ko'); set(h,'MarkerFaceColor','w')

% likelihood statistic
plotlikelilangendijk(la,le,ci)