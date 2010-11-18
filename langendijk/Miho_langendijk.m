% localization model for Loca Photo Prestest according to Langendijk
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Robert Baumgartner, OEAW Acoustical Research Institute
% latest update: 2010-08-04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               SETTINGS                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global settings
plotflag=1;     % switch for plots
name='MMihocic';
lat=00;         % lateral angle of sagital plane
delta=0;        % lateral variability in degree
deltat=20;      % lateral variability of target angle in degree
% model settings
bw = 06;        % bandwidth of averagingfb as partial of an octave
do = 0;         % differential order
cp = 'std';   % comparison process; 'std' (default) or 'xcorr' 
s  = 2;         % standard deviation of transforming Gaussian function
bal= 1;         % balance of left to right channel; default: 1
fstart =2000;   % start frequency; minimum: 0,5kHz;   default: 2kHz
fend   =16000;   % end frequency;                      default: 16kHz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% good results for: cp='xcorr' & do=0 OR cp='std' & do=1

h = waitbar(0,'Please wait...');
b=load(['P:\Projects\HRTF_measurement\Loca Photo\Experiments\Miho Pretest\sources\hrtf_M_Mdtf.mat']);
a2=load(['P:\Projects\HRTF_measurement\Loca Photo\Experiments\Miho Pretest\sources\hrtf_M_Mdtf 8cc.mat']);
a3=load(['P:\Projects\HRTF_measurement\Loca Photo\Experiments\Miho Pretest\sources\hrtf_M_Cdtf MGD.mat']);
a4=load(['P:\Projects\HRTF_measurement\Loca Photo\Experiments\Miho Pretest\sources\hrtf_M_KEMAR normal pinna.mat']);
idx=find(b.posLIN(:,6)>=-(delta+0.01)/2+lat & b.posLIN(:,6)<=(delta+0.01)/2+lat); % median plane with lateral delta-variation; +0.01 because of rounding errors in posLIN
[pol,polidx]=unique(b.posLIN(idx,7));   % sorted polar angles
sagir1=double(b.hM(:,idx,:));           % unsorted impulse responses on sagital plane
ir1=sagir1(:,polidx,:);                 % sorted

sagir2=double(a2.hM(:,idx,:)); % unsorted
ir2=sagir2(:,polidx,:);        % sorted
sagir3=double(a3.hM(:,idx,:)); 
ir3=sagir3(:,polidx,:);        

idxK=find(a4.posLIN(:,4)>=-(delta+0.01)/2+lat & a4.posLIN(:,4)<=(delta+0.01)/2+lat); % median plane with lateral delta-variation; +0.01 because of rounding errors in posLIN
[polK,polidxK]=unique(a4.posLIN(idxK,5));   % sorted polar angles
sagirK=double(a4.hM(:,idxK,:));           % unsorted impulse responses on sagital plane                
ir4=sagirK(:,polidxK,:);

% Results of localization tests
test1=GetData(name, 'PreTest Mdtf', '..\..\..\Pilot Experiments\NH\MMihocic\Loca\Loca Photo\', 1);
test2=GetData(name, 'PreTest Mdtf 8cc', '..\..\..\Pilot Experiments\NH\MMihocic\Loca\Loca Photo\', 1);
test3=GetData(name, 'PreTest Cdtf MGD', '..\..\..\Pilot Experiments\NH\MMihocic\Loca\Loca Photo\', 1);
test4=GetData(name, 'PreTest KEMAR', '..\..\..\Pilot Experiments\NH\MMihocic\Loca\Loca Photo\', 1);

% Mdtf test
idxt1=find(test1(:,5)>=-(deltat+0.01)/2+lat & test1(:,5)<=(deltat+0.01)/2+lat);
[target1,polidxt1]=sort(test1(idxt1,6)); % sorted target angles of loc. test
response1=test1(idxt1,8); % unsorted response angles
response1=response1(polidxt1);
% Mdtf 8cc test
idxt2=find(test2(:,5)>=-(deltat+0.01)/2+lat & test2(:,5)<=(deltat+0.01)/2+lat);
[target2,polidxt2]=sort(test2(idxt2,6)); % sorted target angles of loc. test
response2=test2(idxt2,8); % unsorted response angles
response2=response2(polidxt2);
% Cdtf MGD test
idxt3=find(test3(:,5)>=-(deltat+0.01)/2+lat & test3(:,5)<=(deltat+0.01)/2+lat);
[target3,polidxt3]=sort(test3(idxt3,6)); % sorted target angles of loc. test
response3=test3(idxt3,8); % unsorted response angles
response3=response3(polidxt3);
% KEMAR test
idxt4=find(test4(:,5)>=-(deltat+0.01)/2+lat & test4(:,5)<=(deltat+0.01)/2+lat);
[target4,polidxt4]=sort(test4(idxt4,6)); % sorted target angles of loc. test
response4=test4(idxt4,8); % unsorted response angles
response4=response4(polidxt4);

waitbar(1/5)
p1 = langendijk( ir1,ir1,bw,do,cp,s,bal,fstart,fend);        % baseline
waitbar(2/5)
p2 = langendijk( ir2,ir1,bw,do,cp,s,bal,fstart,fend);  
waitbar(3/5)
p3 = langendijk( ir3,ir1,bw,do,cp,s,bal,fstart,fend); 
waitbar(4/5)
p4 = langendijk( ir4,ir1,bw,do,cp,s,bal,fstart,fend);
waitbar(5/5)
close(h)

% plots
if plotflag==1
    plotlangendijk(p1,pol,pol,name,lat,delta,deltat,'Mdtf (baseline)');
    hold on; h=plot( target1, response1, 'ko'); set(h,'MarkerFaceColor','w')
    plotlangendijk(p2,pol,pol,name,lat,delta,deltat,'Mdtf 8cc');
    hold on; h=plot( target2, response2, 'ko'); set(h,'MarkerFaceColor','w')
    plotlangendijk(p3,pol,pol,name,lat,delta,deltat,'Cdtf MGD');
    hold on; h=plot( target3, response3, 'ko'); set(h,'MarkerFaceColor','w')
    plotlangendijk(p4,polK,pol,name,lat,delta,deltat,'KEMAR');
    hold on; h=plot( target4, response4, 'ko'); set(h,'MarkerFaceColor','w')
end