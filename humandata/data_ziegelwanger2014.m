function data = data_ziegelwanger2014(varargin)
%DATA_ZIEGELWANGER2013  Data from Ziegelwanger and Majdak (2013)
%   Usage: data = data_ziegelwanger2014(flag)
%
%   `data_ziegelwanger2014(flag)` returns results for different HRTF
%   databases from Ziegelwanger and Majdak (2013).
%
%   The flag may be one of:
%  
%     'ARI'         ARI database. The output has the following
%                   fields: `data.results` and `data.subjects`.
%  
%     'CIPIC'       CIPIC database. The output has the following fields: 
%                   `data.results` and `data.subjects`.
%  
%     'LISTEN'      LISTEN database. The output has the following fields.
%                   `data.results` and `data.subjects`.
%  
%     'SPHERE_ROT'  HRTF sets for a rigid sphere placed in the center of
%                   the measurement setup and varying rotation. The
%                   output has the following fields: `data.results`,
%                   `data.subjects`, `data.phi`, `data.theta` and `data.radius`.
%  
%     'SPHERE_DIS'  HRTF sets for a rigid sphere with various positions in
%                   the measurement setup. The output has the following fields: 
%                   `data.results`, `data.subjects`, `data.xM`, `data.yM`,
%                   `data.zM` and `data.radius`.
%  
%     'NH89'        HRTF set of listener NH89 of the ARI database: The
%                   output has the following fields: `data.hM`,
%                   `data.meta` and `data.stimPar`.
%  
%     'reload'      Reload previously calculated results    
%  
%     'recalc'      Recalculate the results  
%
%   The fields are given by:
%
%     `data.results`     Results for all HRTF sets
%
%     `data.subjects`    IDs for HRTF sets
%
%     `data.phi`         Azimuth of ear position
%
%     `data.theta`       Elevation of ear position
%
%     `data.radius`      sphere radius
%
%     `data.xM`          x-coordinate of sphere center
%
%     `data.yM`          y-coordinate of sphere center
%
%     `data.zM`          z-coordinate of sphere center
%
%     `data`             SOFA object
% 
%   Examples:
%   ---------
% 
%   To get results from the ARI database, use::
%
%     data=data_ziegelwanger2014('ARI');
%
%   See also: ziegelwanger2014, ziegelwanger2014onaxis,
%   ziegelwanger2014offaxis, exp_ziegelwanger2014
%
%   References: ziegelwanger2014

% AUTHOR: Harald Ziegelwanger, Acoustics Research Institute, Vienna,
% Austria

%% ------ Check input options --------------------------------------------

% Define input flags
definput.flags.type = {'missingflag','ARI','CIPIC','LISTEN','SPHERE_DIS','SPHERE_ROT','NH89'};
definput.flags.results = {'reload','recalc'};

% Parse input options
[flags,keyvals]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
    flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
        sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
    error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
else
    hpath = which('hrtfinit');  % find local path of hrtf repository
    hpath = hpath(1:end-10);
  
%     if exist([hpath 'ziegelwanger2014'],'dir') ~= 7
%         fprintf([' Sorry! Before you can run this script, you have to download the HRTF Database from \n http://www.kfs.oeaw.ac.at/hrtf/database/amt/ziegelwanger2014.zip , \n unzip it, and move it into your HRTF repository \n ' hpath ' .\n' ' Then, press any key to quit pausing. \n'])
%         pause
%     end
    
    hpath = [hpath 'ziegelwanger2013' filesep];
end

%% ARI database
if flags.do_ARI
    
    if flags.do_recalc || ~exist([hpath 'ARI_results.mat'],'file')
        tmp=load([hpath 'info.mat']);
        data=tmp.info.ARI;
        for ii=1:length(data.subjects)
            disp(['Recalculate data for sujbect ' num2str(ii) '/' num2str(length(data.subjects)) ' of ARI database']);
            Obj=SOFAload([hpath 'ARI_' data.subjects{ii} '.sofa']);
            idx=find(mod(Obj.SourcePosition(:,2),10)==0);
            Obj.Data.IR=Obj.Data.IR(idx,:,:);
            Obj.SourcePosition=Obj.SourcePosition(idx,:);
            Obj.API.M=length(idx);
            Obj.MeasurementSourceAudioChannel=Obj.MeasurementSourceAudioChannel(idx,:);
            Obj.MeasurementAudioLatency=Obj.MeasurementAudioLatency(idx,:);
            Obj.DimSize.M=length(idx);
                
            [~,results(ii).MAX{1}]=ziegelwanger2014(Obj,1,0,1);
            [~,results(ii).CTD{1}]=ziegelwanger2014(Obj,2,0,1);
            [~,results(ii).AGD{1}]=ziegelwanger2014(Obj,3,0,1);
            [~,results(ii).MCM{1}]=ziegelwanger2014(Obj,4,0,1);

            [~,results(ii).MAX{2}]=ziegelwanger2014(Obj,1,2,1);
            [~,results(ii).CTD{2}]=ziegelwanger2014(Obj,2,2,1);
            [~,results(ii).AGD{2}]=ziegelwanger2014(Obj,3,2,1);
            [~,results(ii).MCM{2}]=ziegelwanger2014(Obj,4,2,1);

            [~,results(ii).MAX{3}]=ziegelwanger2014(Obj,1,3,1);
            [~,results(ii).CTD{3}]=ziegelwanger2014(Obj,2,3,1);
            [~,results(ii).AGD{3}]=ziegelwanger2014(Obj,3,3,1);
            [~,results(ii).MCM{3}]=ziegelwanger2014(Obj,4,3,1);

            [~,results(ii).MAX{4}]=ziegelwanger2014(Obj,1,1,1);
            [~,results(ii).CTD{4}]=ziegelwanger2014(Obj,2,1,1);
            [~,results(ii).AGD{4}]=ziegelwanger2014(Obj,3,1,1);
            [~,results(ii).MCM{4}]=ziegelwanger2014(Obj,4,1,1);
        end
        save([hpath 'ARI_results.mat'],'results');
        data.results=results;
    else
        tmp=load([hpath 'info.mat']);
        data=tmp.info.ARI;
        load([hpath 'ARI_results.mat']);
        data.results=results;
    end
    
end

%% CIPIC database
if flags.do_CIPIC
    
    if flags.do_recalc || ~exist([hpath 'CIPIC_results.mat'],'file')
        tmp=load([hpath 'info.mat']);
        data=tmp.info.CIPIC;
        for ii=1:length(data.subjects)
            disp(['Recalculate data for sujbect ' num2str(ii) '/' num2str(length(data.subjects)) ' of CIPIC database']);
            Obj=SOFAload([hpath 'CIPIC_' data.subjects{ii} '.sofa']);
                
            [~,results(ii).MAX{1}]=ziegelwanger2014(Obj,1,0,1);
            [~,results(ii).CTD{1}]=ziegelwanger2014(Obj,2,0,1);
            [~,results(ii).AGD{1}]=ziegelwanger2014(Obj,3,0,1);
            [~,results(ii).MCM{1}]=ziegelwanger2014(Obj,4,0,1);

            [~,results(ii).MAX{2}]=ziegelwanger2014(Obj,1,2,1);
            [~,results(ii).CTD{2}]=ziegelwanger2014(Obj,2,2,1);
            [~,results(ii).AGD{2}]=ziegelwanger2014(Obj,3,2,1);
            [~,results(ii).MCM{2}]=ziegelwanger2014(Obj,4,2,1);

            [~,results(ii).MAX{3}]=ziegelwanger2014(Obj,1,3,1);
            [~,results(ii).CTD{3}]=ziegelwanger2014(Obj,2,3,1);
            [~,results(ii).AGD{3}]=ziegelwanger2014(Obj,3,3,1);
            [~,results(ii).MCM{3}]=ziegelwanger2014(Obj,4,3,1);

            [~,results(ii).MAX{4}]=ziegelwanger2014(Obj,1,1,1);
            [~,results(ii).CTD{4}]=ziegelwanger2014(Obj,2,1,1);
            [~,results(ii).AGD{4}]=ziegelwanger2014(Obj,3,1,1);
            [~,results(ii).MCM{4}]=ziegelwanger2014(Obj,4,1,1);
        end
        save([hpath 'CIPIC_results.mat'],'results');
        data.results=results;
    else
        tmp=load([hpath 'info.mat']);
        data=tmp.info.CIPIC;
        load([hpath 'CIPIC_results.mat']);
        data.results=results;
    end
    
end

%% LISTEN database
if flags.do_LISTEN
    
    if flags.do_recalc || ~exist([hpath 'LISTEN_results.mat'],'file')
        tmp=load([hpath 'info.mat']);
        data=tmp.info.LISTEN;
        for ii=1:length(data.subjects)
            if ~strcmp(data.subjects{ii},'34')
                disp(['Recalculate data for subject ' num2str(ii) '/' num2str(length(data.subjects)) ' of LISTEN database']);
                Obj=SOFAload([hpath 'LISTEN_' data.subjects{ii} '.sofa']);
                Obj.Data.SamplingRate=48000;
                
                [~,results(ii).MAX{1}]=ziegelwanger2014(Obj,1,0,1);
                [~,results(ii).CTD{1}]=ziegelwanger2014(Obj,2,0,1);
                [~,results(ii).AGD{1}]=ziegelwanger2014(Obj,3,0,1);
                [~,results(ii).MCM{1}]=ziegelwanger2014(Obj,4,0,1);

                [~,results(ii).MAX{2}]=ziegelwanger2014(Obj,1,2,1);
                [~,results(ii).CTD{2}]=ziegelwanger2014(Obj,2,2,1);
                [~,results(ii).AGD{2}]=ziegelwanger2014(Obj,3,2,1);
                [~,results(ii).MCM{2}]=ziegelwanger2014(Obj,4,2,1);

                [~,results(ii).MAX{3}]=ziegelwanger2014(Obj,1,3,1);
                [~,results(ii).CTD{3}]=ziegelwanger2014(Obj,2,3,1);
                [~,results(ii).AGD{3}]=ziegelwanger2014(Obj,3,3,1);
                [~,results(ii).MCM{3}]=ziegelwanger2014(Obj,4,3,1);

                [~,results(ii).MAX{4}]=ziegelwanger2014(Obj,1,1,1);
                [~,results(ii).CTD{4}]=ziegelwanger2014(Obj,2,1,1);
                [~,results(ii).AGD{4}]=ziegelwanger2014(Obj,3,1,1);
                [~,results(ii).MCM{4}]=ziegelwanger2014(Obj,4,1,1);
            end
        end
        save([hpath 'LISTEN_results.mat'],'results');
        data.results=results;
    else
        tmp=load([hpath 'info.mat']);
        data=tmp.info.LISTEN;
        load([hpath 'LISTEN_results.mat']);
        data.results=results;
    end
    
end

%% SPHERE (Displacement) database
if flags.do_SPHERE_DIS
    
    if flags.do_recalc || ~exist([hpath 'Sphere_Displacement_results.mat'],'file')
        tmp=load([hpath 'info.mat']);
        data=tmp.info.Sphere_Displacement;
        results.p_onaxis=zeros(4,2,length(data.subjects));
        results.p_offaxis=zeros(7,2,length(data.subjects));
        for ii=1:length(data.subjects)
            disp(['Recalculate data for subject ' num2str(ii) '/' num2str(length(data.subjects)) ' of SPHERE_DIS database']);
            Obj=SOFAload([hpath 'Sphere_Displacement_' data.subjects{ii} '.sofa']);
                
            [~,results(ii).MAX{1}]=ziegelwanger2014(Obj,1,0,1);
            [~,results(ii).CTD{1}]=ziegelwanger2014(Obj,2,0,1);
            [~,results(ii).AGD{1}]=ziegelwanger2014(Obj,3,0,1);
            [~,results(ii).MCM{1}]=ziegelwanger2014(Obj,4,0,1);
                
            [~,results(ii).MAX{2}]=ziegelwanger2014(Obj,1,1,1);
            [~,results(ii).CTD{2}]=ziegelwanger2014(Obj,2,1,1);
            [~,results(ii).AGD{2}]=ziegelwanger2014(Obj,3,1,1);
            [~,results(ii).MCM{2}]=ziegelwanger2014(Obj,4,1,1);
        end
        save([hpath 'Sphere_Displacement_results.mat'],'results');
        data.results=results;
    else
        tmp=load([hpath 'info.mat']);
        data=tmp.info.Sphere_Displacement;
        load([hpath 'Sphere_Displacement_results.mat']);
        data.results=results;
    end
    
end

%% SPHERE (Rotation) database
if flags.do_SPHERE_ROT
    
    if flags.do_recalc || ~exist([hpath 'Sphere_Rotation_results.mat'],'file')
        tmp=load([hpath 'info.mat']);
        data=tmp.info.Sphere_Rotation;
        results.p=zeros(4,2,length(data.phi));
        for ii=1:length(data.subjects)
            disp(['Recalculate data for subject ' num2str(ii) '/' num2str(length(data.subjects)) ' of SPHERE_ROT database']);
            Obj=SOFAload([hpath 'Sphere_Rotation_' data.subjects{ii} '.sofa']);
                
            [~,results(ii).MAX{1}]=ziegelwanger2014(Obj,1,0,1);
            [~,results(ii).CTD{1}]=ziegelwanger2014(Obj,2,0,1);
            [~,results(ii).AGD{1}]=ziegelwanger2014(Obj,3,0,1);
            [~,results(ii).MCM{1}]=ziegelwanger2014(Obj,4,0,1);
                
            [~,results(ii).MAX{2}]=ziegelwanger2014(Obj,1,1,1);
            [~,results(ii).CTD{2}]=ziegelwanger2014(Obj,2,1,1);
            [~,results(ii).AGD{2}]=ziegelwanger2014(Obj,3,1,1);
            [~,results(ii).MCM{2}]=ziegelwanger2014(Obj,4,1,1);
        end
        save([hpath 'Sphere_Rotation_results.mat'],'results');
        data.results=results;
    else
        tmp=load([hpath 'info.mat']);
        data=tmp.info.Sphere_Rotation;
        load([hpath 'Sphere_Rotation_results.mat']);
        data.results=results;
    end
    
end

%% ARI database (NH89)
if flags.do_NH89

    data=SOFAload([hpath 'ARI_NH89.sofa']);
    
end