% Implementation of Langendijk's localization model for spatstrat
% experiments
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Robert Baumgartner, OEAW Acoustical Research Institute
% latest update: 2010-08-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

listener={'CTrofeit';'IKöcher';'JStubner';'TStimmer';'EEnzinger';'SRamforth';'AFlorianz';'DMayer';'FHasenhindl';'APrey';'GNowak';'BElwischger';'VHolczmann';'WSchitter'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               SETTINGS                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global settings
plotflag=1;     % switch for plots
name=char(listener(5));
lat=00;         % lateral angle of sagital plane
delta=0;        % lateral variability in degree
deltat=20;      % lateral variability of target angle in degree
% model settings
bw = 06;        % bandwidth of averagingfb as partial of an octave
do = 0;         % differential order
cp = 'std';   % comparison process; 'std' (default) or 'xcorr' 
s  = 2;         % standard deviation of transforming Gaussian function
bal= 1;         % balance of left to right channel; default: 1
fstart =1000;   % start frequency; minimum: 0,5kHz;   default: 2kHz
% fend   =16000;   % end frequency;                      default: 16kHz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% good results for: cp='xcorr' & do=0 OR cp='std' & do=1

ad=load(['P:\Projects\HRTF_measurement\Database\Measurements\' name '\hrtf_M_dtf 512_dummy.mat']);
aw=load(['P:\Projects\HRTF_measurement\Database\Measurements\' name '\hrtf_M_dtf 512_warped.mat']);
b=load(['P:\Projects\HRTF_measurement\Database\Measurements\' name '\hrtf_M_dtf 256.mat']);
idx=find(b.posLIN(:,6)>=-(delta+0.01)/2+lat & b.posLIN(:,6)<=(delta+0.01)/2+lat); % median plane with lateral delta-variation; +0.01 because of rounding errors in posLIN
[pol,polidx]=unique(b.posLIN(idx,7));   % sorted polar angles
sagir2=double(b.hM(:,idx,:));           % unsorted impulse responses on sagital plane
ir2=sagir2(:,polidx,:);                 % sorted

sagird=double(ad.hM(:,idx,:)); % unsorted
ird=sagird(:,polidx,:);        % sorted
sagirw=double(aw.hM(:,idx,:)); 
irw=sagirw(:,polidx,:);        

% ir2(:,7,:)=zeros(size(ir1,1),1,2);
% ir2(1,7,:)=1;
% ir2(:,38,:)=zeros(size(ir1,1),1,2);
% ir2(1,38,:)=1;
% ir1(:,12,:)=zeros(size(ir1,1),1,2);
% ir1(1,12,:)=1;


% Results of localization tests
testb=GetData(name, 'Learn M 60', 'Loca\Spatstrat\', 2:6); %Ausgabe(spaltenweise): Target Az, Target El, Resp Az, Resp El, Target Lat, Target Pol, Resp Lat, Resp Pol
% testb=GetData(name, 'PreTest M 60', 'Loca\Spatstrat\',1);
testd=GetData(name, 'PreTest M 60_dummy', 'Loca\Spatstrat\');
testw=GetData(name, 'PreTest M 60_warped', 'Loca\Spatstrat\');
% data=wavread(['P:\Resources\Experimental_Data\Main Experiments\NH\' name '\Loca\Spatstrat\training.wav']);
datad=wavread('P:\Projects\CI_HRTF\Loca SpatStrat NH\Audio\dummy.wav');  % dummy simulus
dataw=wavread('P:\Projects\CI_HRTF\Loca SpatStrat NH\Audio\warped.wav'); % warped simulus

% aus AverageData03FrontBack.m (in P:\Projects\CI_HRTF\Loca BtE CI
% HRTF\Results) fuer GetData.m-Aufruf
% f={'mod', 'project', 'cond', 'name', 'ID', 'trials'};
% data={ ...
%        'NH' 'Loca\Spatstrat\' 'Learn M' 'EEnzinger' 'NH12' 38;  ...
%        'NH' 'Loca\Spatstrat\' 'Learn M' 'CTrofeit'  'NH15' 24; ...
%        'NH' 'Loca\Spatstrat\' 'Learn M' 'JMalik' 'NH16' 42; ...
%        'NH' 'Loca\Spatstrat\' 'Learn M' 'FWippel'  'NH17' 16; ...
%        'NH' 'Loca\Spatstrat\' 'Learn M' 'KAntonicek' 'NH18' 15;  ...
%      };       
% s=cell2struct(data,f,2);
% 
% clear bi bii nn qe chr
% for ii=1:length(s)
%   test=GetData(s(ii).name, s(ii).cond, s(ii).project, 5:12, s(ii).mod);

%   m1=SelectData(test,range(1,1),range(1,2),range(1,3),range(1,4),rangedim,0);  
%   idx=find(m1(:,7)>=range(1,1) & m1(:,7)<=range(1,2)); % lateral range
%   m1=m1(idx,:);
%   idx=find(abs(90-m1(:,6))>30); % +/-30° um top (90°) herum weg...
%   m1=m1(idx,:);
%   if rnd, m1(:,8)=rand(size(m1,1),1)*240-30; end
%   idx=find(abs(90-m1(:,8))>30);
%   m1=m1(idx,:);
%   idx=find(sign(90-m1(:,8))~=sign(90-m1(:,6)));
%   bi(ii)=size(idx,1);
%   nn(ii)=size(m1,1);
% %   disp([s(ii).ID ':' num2str(bi(ii)) ' / ' num2str(nn(ii)) ' = ' num2str(bi(ii)/nn(ii))]);
% end
% disp([s(1).mod ': ' num2str(sum(bi)) ' / ' num2str(sum(nn)) ' = ' num2str(sum(bi)/sum(nn))]);

% baseline test
idxtb=find(testb(:,5)>=-(deltat+0.01)/2+lat & testb(:,5)<=(deltat+0.01)/2+lat);
[targetb,polidxtb]=sort(testb(idxtb,6)); % sorted target angles of loc. test
responseb=testb(idxtb,8); % unsorted response angles
responseb=responseb(polidxtb);
% dummy test
idxtd=find(testd(:,5)>=-(deltat+0.01)/2+lat & testd(:,5)<=(deltat+0.01)/2+lat);
[targetd,polidxtd]=sort(testd(idxtd,6)); % sorted target angles of loc. test
responsed=testd(idxtd,8); % unsorted response angles
responsed=responsed(polidxtd);
% warped test
idxtw=find(testw(:,5)>=-(deltat+0.01)/2+lat & testw(:,5)<=(deltat+0.01)/2+lat);
[targetw,polidxtw]=sort(testw(idxtw,6)); % sorted target angles of loc. test
responsew=testw(idxtw,8); % unsorted response angles
responsew=responsew(polidxtw);


pb = langendijk( ir2,ir2,bw,do,cp,s,bal,fstart,16000);        % baseline
pd = langendijk( ird,ir2,bw,do,cp,s,bal,fstart,8000,datad);  % dummy
pw = langendijk( irw,ir2,bw,do,cp,s,bal,fstart,8000,dataw);  % warped
la=zeros(3,1);le=zeros(3,1);ci=zeros(3,2);
[la(1),le(1),ci(1,:)] = likelilangendijk( pb,pol,pol,targetb,responseb );
[la(2),le(2),ci(2,:)] = likelilangendijk( pd,pol,pol,targetd,responsed );
[la(3),le(3),ci(3,:)] = likelilangendijk( pd,pol,pol,targetw,responsew );

% plots
if plotflag==1
    plotlangendijk(pb,pol,pol,name,lat,delta,deltat,'baseline');
    hold on; h=plot( targetb, responseb, 'ko'); set(h,'MarkerFaceColor','w')
    plotlangendijk(pd,pol,pol,name,lat,delta,deltat,'dummy');
    hold on; h=plot( targetd, responsed, 'ko'); set(h,'MarkerFaceColor','w')
    plotlangendijk(pw,pol,pol,name,lat,delta,deltat,'warped');
    hold on; h=plot( targetw, responsew, 'ko'); set(h,'MarkerFaceColor','w')

    plotlikelilangendijk(la,le,ci,'spatstrat')
end