function data = data_langendijk2002(flag)
%DATA_LANGENDIJK2002 Returns data points of response patterns in figures 7
%and 9 or of left-ear DTF plots in figure 11 of Langendijk's paper
%   Usage: data = data_langendijk2002(flag)
%
%   Output parameters:
%      data       - the data points from the given figure;
%                   in the case of response patterns (fig. 7&9) the first 
%                   row describes target position and the second one 
%                   belongs to the response position;
%                   in the case of DTF data (fig. 11) the first dimension
%                   of the data matrix describes frequency and the second one
%                   angle position, the FIRST COLUMN defines the actual angle
%                   positions
%
%   DATA_LANGENDIJK2002(flag) returns data points from the Langendijk
%   2002 paper. The flag may be one of:
%
%-    'P3-b'   - data from Fig.9; listener: P3, condition: 'baseline'.
%
%-    'P3-2o'  - data from Fig.9; listener: P3, condition: '2-oct'.
%
%-    'P3-1ol' - data from Fig.9; listener: P3, condition: '1-oct(low)'.
%
%-    'P3-1om' - data from Fig.9; listener: P3, condition: '1-oct(middle)'.
%
%-    'P3-1oh' - data from Fig.9; listener: P3, condition: '1-oct(high)'.
%
%-    'P6-b'   - data from Fig.9; listener: P6, condition: 'baseline'.
%
%-    'P6-2o'  - data from Fig.9; listener: P6, condition: '2-oct'.
%
%-    'P6-1ol' - data from Fig.9; listener: P6, condition: '1-oct(low)'.
%
%-    'P6-1om' - data from Fig.9; listener: P6, condition: '1-oct(middle)'.
%
%-    'P6-1oh' - data from Fig.9; listener: P6, condition: '1-oct(high)'.
%
%-    'P3-dtf' - DTF data from Fig.11; listener: P3. If you do not have the
%                bitmap of the JASA paper you will get a structure containing 
%                the dtf data in ARI format.
%
%-    'P6-dtf' - DTF data from Fig.11; listener: P6. If you do not have the
%                bitmap of the JASA paper you will get a structure containing 
%                the dtf data in ARI format.
%
%    If no flag is given, the function will print the list of valid flags.
%
%R  langendijk et al. 2002
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Robert Baumgartner, OEAW Acoustical Research Institute
% latest update: 2010-08-19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target=-55:29:235;
idb=round(0.5:0.1:11.4);
idc=round(0.5:0.2:11.4);
switch flag
    case 'P3-b'
        target=target(idb);
        response=zeros(1,110);
        response(1:10) =[-55,-55,-55,-55,-55,-55,-53,-54,-45,-45];
        response(11:20)=[-55,-40,-38,-34,-30,-29,-28,-22,-25,-15];
        response(21:30)=[-12,00,03,05,09,15,20,24,30,140];
        response(31:40)=[21,29,31,32,34,45,58,75,110,122];
        response(41:50)=[118,119,123,134,135,140,163,174,176,210];
        response(51:60)=[106,107,123,124,137,162,183,210,219,225];
        response(61:70)=[-55,95,113,119,120,121,126,161,181,214];
        response(71:80)=[119,130,142,143,170,184,190,199,225,235];
        response(81:90)=[146,153,154,155,157,161,180,181,181,210];
        response(91:100)=[93,94,170,180,182,183,184,191,230,231];
        response(101:110)=[230,231,235,235,235,235,235,235,235,234];
        data=[target;response];
    case 'P3-2o'
        target=target(idc);
        response=zeros(1,55);
        response(1:5)  =[-45,-12,165,181,182];
        response(6:10) =[-53,-30,-12,-11,192];
        response(11:15)=[-20,128,141,170,213];
        response(16:20)=[-55,-54,130,175,204];
        response(21:25)=[-30,72,172,180,202];
        response(26:30)=[164,181,219,235,234];
        response(31:35)=[-55,108,184,210,235];
        response(36:40)=[116,182,183,184,210];
        response(41:45)=[113,165,180,181,225];
        response(46:50)=[-54,-41,134,169,184];
        response(51:55)=[-55,234,234,235,235];
        data=[target;response];
    case 'P3-1ol'
        target=target(idc);
        response=zeros(1,55);
        response(1:5)  =[-54,-50,-49,-30,-25];
        response(6:10) =[-38,-37,-30,-27,-15];
        response(11:15)=[-28,-20,-5,3,13];
        response(16:20)=[-37,38,37,106,107];
        response(21:25)=[145,161,169,170,182];
        response(26:30)=[-50,-32,103,150,176];
        response(31:35)=[121,122,126,147,165];
        response(36:40)=[-55,110,165,168,175];
        response(41:45)=[144,152,154,180,187];
        response(46:50)=[170,175,176,212,223];
        response(51:55)=[-55,235,234,235,235];
        data=[target;response];
    case 'P3-1om'
        target=target(idc);
        response=zeros(1,55);
        response(1:5)  =[-55,-55,-54,-50,-48];
        response(6:10) =[-35,-28,-27,-25,-24];
        response(11:15)=[-51,-48,-25,-23,-18];
        response(16:20)=[-55,112,210,234,235];
        response(21:25)=[122,122,160,172,185];
        response(26:30)=[122,138,208,213,223];
        response(31:35)=[122,163,176,219,235];
        response(36:40)=[-55,-54,177,234,235];
        response(41:45)=[-55,128,163,180,224];
        response(46:50)=[-55,-54,150,151,183];
        response(51:55)=[-55,-54,228,234,235];
        data=[target;response];
    case 'P3-1oh'
        target=target(idc);
        response=zeros(1,55);
        response(1:5)  =[-55,-54,-50,-42,235];
        response(6:10) =[-55,-50,-40,-22,175];
        response(11:15)=[-29,3,60,145,150];
        response(16:20)=[14,130,166,168,180];
        response(21:25)=[40,122,123,148,180];
        response(26:30)=[122,123,128,208,215];
        response(31:35)=[119,133,138,145,175];
        response(36:40)=[122,175,180,181,182];
        response(41:45)=[110,175,175,176,190];
        response(46:50)=[-55,174,184,200,201];
        response(51:55)=[234,234,235,235,235];
        data=[target;response];
    case 'P6-b'
        target=target(idb);
        response=zeros(1,110);
        response(1:10) =[-55,-55,-55,-55,-55,-55,-53,-54,-45,-45];
        response(11:20)=[-50,-46,-41,-35,-36,-31,-32,-30,-25,-10];
        response(21:30)=[-25,-13,-8,-1,0,1,3,5,5,10];
        response(31:40)=[15,20,28,28,29,31,33,33,36,50];
        response(41:50)=[55,55,60,64,73,74,80,82,85,95];
        response(51:60)=[90,90,91,93,94,100,101,103,120,130];
        response(61:70)=[70,77,90,90,91,100,105,130,130,131];
        response(71:80)=[-50,-50,-40,0,95,100,101,105,162,167];
        response(81:90)=[145,146,150,160,165,169,170,171,180,181];
        response(91:100)=[185,186,195,196,200,209,210,211,215,220];
        response(101:110)=[210,215,220,225,230,231,235,235,234,234];
        data=[target;response];
    case 'P6-2o'
        target=target(idc);
        response=zeros(1,55);
        response(1:5)  =[-55,-54,-30,110,165];
        response(6:10) =[-55,-45,-40,-35,-20];
        response(11:15)=[-50,-35,-25,-20,125];
        response(16:20)=[-55,-50,-30,-31,180];
        response(21:25)=[-45,125,145,180,185];
        response(26:30)=[-55,-45,-30,-15,140];
        response(31:35)=[-55,-50,-45,-44,-20];
        response(36:40)=[-45,-40,-41,-35,170];
        response(41:45)=[-50,-35,100,180,181];
        response(46:50)=[-55,-54,-50,-45,140];
        response(51:55)=[-55,-55,-54,-54,-53];
        data=[target;response];
    case 'P6-1ol'
        target=target(idc);
        response=zeros(1,55);
        response(1:5)  =[-55,-54,-53,-49,-50];
        response(6:10) =[-35,-34,-30,-20,-19];
        response(11:15)=[-25,-10,-9,-11,5];
        response(16:20)=[25,30,31,35,40];
        response(21:25)=[65,70,95,96,100];
        response(26:30)=[60,85,90,91,110];
        response(31:35)=[75,85,90,100,105];
        response(36:40)=[-35,-25,-26,-20,80];
        response(41:45)=[0,185,184,186,190];
        response(46:50)=[175,195,205,210,215];
        response(51:55)=[225,230,233,234,235];
        data=[target;response];
    case 'P6-1om'
        target=target(idc);
        response=zeros(1,55);
        response(1:5)  =[-55,-54,-55,-54,-50];
        response(6:10) =[-55,-45,-44,-46,-35];
        response(11:15)=[-15,-10,5,6,10];
        response(16:20)=[5,15,20,30,31];
        response(21:25)=[45,60,61,65,66];
        response(26:30)=[-45,60,100,120,121];
        response(31:35)=[-55,80,105,175,176];
        response(36:40)=[-50,-40,-41,-10,10,];
        response(41:45)=[-55,-54,-50,-30,195];
        response(46:50)=[-50,100,105,106,135];
        response(51:55)=[-55,-55,-54,-54,235];
        data=[target;response];
    case 'P6-1oh'
        target=target(idc);
        response=zeros(1,55);
        response(1:5)  =[-55,-54,-55,-54,-35];
        response(6:10) =[-55,-54,-50,-51,-35];
        response(11:15)=[-15,-11,-10,-9,-5];
        response(16:20)=[10,25,35,55,105];
        response(21:25)=[40,45,50,51,75];
        response(26:30)=[-10,70,72,80,102];
        response(31:35)=[90,95,100,101,105];
        response(36:40)=[-55,-25,40,105,145];
        response(41:45)=[-55,-50,-45,-30,-25];
        response(46:50)=[-35,-34,-5,162,180];
        response(51:55)=[-55,-50,-54,235,234];
        data=[target;response];
    case 'P3-dtf'
        if exist('langendijk2002-dtfP3.bmp','file')==2
            [med,pol]=bmp2gr('langendijk2002-dtfP3');
            data=[pol;med];
        else
            data=load('hrtf_M_langendijk2002 P3');
        end
    case 'P6-dtf'
        if exist('langendijk2002-dtfP6.bmp','file')==2
            [med,pol]=bmp2gr('langendijk2002-dtfP6');
            data=[pol;med];
        else
            data=load('hrtf_M_langendijk2002 P6');
        end
    otherwise
        disp('You must specify one of the following flags: P3-b,P3-2o,P3-1ol,P3-1om,P3-1oh,P6-b,P6-2o,P6-1ol,P6-1om,P6-1oh,P3-dtf,P6-dtf')
end
end
            
            