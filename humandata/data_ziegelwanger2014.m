function data = data_ziegelwanger2014(varargin)
%DATA_ZIEGELWANGER2014   Data from Ziegelwanger and Majdak (2014)
%   Usage: data = data_ziegelwanger2014(flag)
%
%   `data_ziegelwanger2014(flag)` returns results for different HRTF
%   databases from Ziegelwanger and Majdak (2014).
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
%     'Sphere'      HRTF set for a rigid sphere: The
%                   output has the following fields: `data.hM`,
%                   `data.meta` and `data.stimPar`.
%  
%     'SAT'         HRTF set for a rigid sphere combined with a neck and a
%                   torso: The output has the following fields: `data.hM`,
%                   `data.meta` and `data.stimPar`.
%  
%     'STP'         HRTF set for a rigid sphere combined with a neck, a
%                   torso and a pinna: The output has the following fields:
%                   `data.hM`,`data.meta` and `data.stimPar`.
%  
%     'NH89'        HRTF set of listener NH89 of the ARI database: The
%                   output has the following fields: `data.hM`,
%                   `data.meta` and `data.stimPar`.  
%  
%     'redo'        Recalculate the results  
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
%   Requirements: 
%   -------------
%
%   1) SOFA API from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
% 
%   2) Optimization Toolbox for Matlab
%
%   3) Data in hrtf/ziegelwanger2014
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
definput.flags.type = {'missingflag','ARI','CIPIC','LISTEN','SPHERE_DIS','SPHERE_ROT','NH89','Sphere','SAT','STP'};
definput.flags.redo = {'missingflag','redo'};

% Parse input options
[flags,~]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
    flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
        sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
    error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
else
    hpath = which('hrtfinit');  % find local path of hrtf repository
    hpath = hpath(1:end-10);
    hpath = [hpath 'ziegelwanger2014' filesep];
    
    if ~exist([hpath 'info.mat'],'file')
        urlwrite([SOFAdbURL filesep 'ziegelwanger2014/info.mat'],[hpath 'info.mat']);
    end
end

%% ARI database
if flags.do_ARI
    
    if flags.do_redo || ~exist([hpath '..' filesep '..' filesep 'experiments' filesep 'exp_ziegelwanger2014_ARI.mat'],'file')
        tmp=load([hpath 'info.mat']);
        data=tmp.info.ARI;
        for ii=1:length(data.subjects)
            disp(['Recalculate data for subject ' num2str(ii) '/' num2str(length(data.subjects)) ' (' data.subjects{ii} ') of ARI database']);
            Obj=SOFAload([hpath 'ARI_' data.subjects{ii} '.sofa']);
             
            if exist([hpath 'ARI_' data.subjects{ii} '.sofa.MCM.mat'],'file')
                load([hpath 'ARI_' data.subjects{ii} '.sofa.MCM.mat']);
            else
                [~,tmp]=ziegelwanger2014(Obj,4,0,0);
                toaEst=tmp.toa;
                save([hpath 'ARI_' data.subjects{ii} '.sofa.MCM.mat'],'toaEst');
            end
            [~,results(ii).MCM{1}]=ziegelwanger2014(Obj,toaEst,0,1e-8);
            [~,results(ii).MCM{2}]=ziegelwanger2014(Obj,toaEst,[0.05 0.01],1e-8);
        end
        save([hpath '..' filesep '..' filesep 'experiments' filesep 'exp_ziegelwanger2014_ARI.mat'],'results');
        data.results=results;
    else
        tmp=load([hpath 'info.mat']);
        data=tmp.info.ARI;
        load([hpath '..' filesep '..' filesep 'experiments' filesep 'exp_ziegelwanger2014_ARI.mat']);
        data.results=results;
    end
    
end

%% CIPIC database
if flags.do_CIPIC
    
    if flags.do_redo || ~exist([hpath '..' filesep '..' filesep 'experiments' filesep 'exp_ziegelwanger2014_CIPIC.mat'],'file')
        tmp=load([hpath 'info.mat']);
        data=tmp.info.CIPIC;
        for ii=1:length(data.subjects)
            disp(['Recalculate data for sujbect ' num2str(ii) filesep num2str(length(data.subjects)) ' of CIPIC database']);
            Obj=SOFAload([hpath 'CIPIC_' data.subjects{ii} '.sofa']);
             
            if exist([hpath 'CIPIC_' data.subjects{ii} '.sofa.MCM.mat'],'file')
                load([hpath 'CIPIC_' data.subjects{ii} '.sofa.MCM.mat']);
            else
                [~,tmp]=ziegelwanger2014(Obj,4,0,0);
                toaEst=tmp.toa;
                save([hpath 'CIPIC_' data.subjects{ii} '.sofa.MCM.mat'],'toaEst');
            end
            [~,results(ii).MCM{1}]=ziegelwanger2014(Obj,toaEst,0,1e-8);
            [~,results(ii).MCM{2}]=ziegelwanger2014(Obj,toaEst,[0.05 0.01],1e-8);
        end
        save([hpath '..' filesep '..' filesep 'experiments' filesep 'exp_ziegelwanger2014_CIPIC.mat'],'results');
        data.results=results;
    else
        tmp=load([hpath 'info.mat']);
        data=tmp.info.CIPIC;
        load([hpath '..' filesep '..' filesep 'experiments' filesep 'exp_ziegelwanger2014_CIPIC.mat']);
        data.results=results;
    end
    
end

%% LISTEN database
if flags.do_LISTEN
    
    if flags.do_redo || ~exist([hpath '..' filesep '..' filesep 'experiments' filesep 'exp_ziegelwanger2014_LISTEN.mat'],'file')
        tmp=load([hpath 'info.mat']);
        data=tmp.info.LISTEN;
        for ii=1:length(data.subjects)
            if ~strcmp(data.subjects{ii},'34')
                disp(['Recalculate data for subject ' num2str(ii) filesep num2str(length(data.subjects)) ' of LISTEN database']);
                Obj=SOFAload([hpath 'LISTEN_' data.subjects{ii} '.sofa']);
             
                if exist([hpath 'LISTEN_' data.subjects{ii} '.sofa.MCM.mat'],'file')
                    load([hpath 'LISTEN_' data.subjects{ii} '.sofa.MCM.mat']);
                else
                    [~,tmp]=ziegelwanger2014(Obj,4,0,0);
                    toaEst=tmp.toa;
                    save([hpath 'LISTEN_' data.subjects{ii} '.sofa.MCM.mat'],'toaEst');
                end
                [~,results(ii).MCM{1}]=ziegelwanger2014(Obj,toaEst,0,1e-8);
                [~,results(ii).MCM{2}]=ziegelwanger2014(Obj,toaEst,[0.05 0.01],1e-8);
            end
        end
        save([hpath '..' filesep '..' filesep 'experiments' filesep 'exp_ziegelwanger2014_LISTEN.mat'],'results');
        data.results=results;
    else
        tmp=load([hpath 'info.mat']);
        data=tmp.info.LISTEN;
        load([hpath '..' filesep '..' filesep 'experiments' filesep 'exp_ziegelwanger2014_LISTEN.mat']);
        data.results=results;
    end
    
end

%% SPHERE (Displacement) database
if flags.do_SPHERE_DIS
    
    if flags.do_redo || ~exist([hpath '..' filesep '..' filesep 'experiments' filesep 'exp_ziegelwanger2014_SPHERE_DIS.mat'],'file')
        tmp=load([hpath 'info.mat']);
        data=tmp.info.Displacement;
        results.p_onaxis=zeros(4,2,length(data.subjects));
        results.p_offaxis=zeros(7,2,length(data.subjects));
        for ii=1:length(data.subjects)
            disp(['Recalculate data for subject ' num2str(ii) filesep num2str(length(data.subjects)) ' of SPHERE_DIS database']);
            Obj=SOFAload([hpath 'Sphere_Displacement_' data.subjects{ii} '.sofa']);
             
            if exist([hpath 'Sphere_Displacement_' data.subjects{ii} '.sofa.MCM.mat'],'file')
                load([hpath 'Sphere_Displacement_' data.subjects{ii} '.sofa.MCM.mat']);
            else
                [~,tmp]=ziegelwanger2014(Obj,4,0,0);
                toaEst=tmp.toa;
                save([hpath 'Sphere_Displacement_' data.subjects{ii} '.sofa.MCM.mat'],'toaEst');
            end  
            [~,results(ii).MCM{1}]=ziegelwanger2014(Obj,toaEst,0,1e-8);
            [~,results(ii).MCM{2}]=ziegelwanger2014(Obj,toaEst,[0.05 0.01],1e-8);
        end
        save([hpath '..' filesep '..' filesep 'experiments' filesep 'exp_ziegelwanger2014_SPHERE_DIS.mat'],'results');
        data.results=results;
    else
        tmp=load([hpath 'info.mat']);
        data=tmp.info.Displacement;
        load([hpath '..' filesep '..' filesep 'experiments' filesep 'exp_ziegelwanger2014_SPHERE_DIS.mat']);
        data.results=results;
    end
    
end

%% SPHERE (Rotation) database
if flags.do_SPHERE_ROT
    
    if flags.do_redo || ~exist([hpath '..' filesep '..' filesep 'experiments' filesep 'exp_ziegelwanger2014_SPHERE_ROT.mat'],'file')
        tmp=load([hpath 'info.mat']);
        data=tmp.info.Rotation;
        results.p=zeros(4,2,length(data.phi));
        for ii=1:length(data.subjects)
            disp(['Recalculate data for subject ' num2str(ii) filesep num2str(length(data.subjects)) ' of SPHERE_ROT database']);
            Obj=SOFAload([hpath 'Sphere_Rotation_' data.subjects{ii} '.sofa']);
                
            if exist([hpath 'Sphere_Rotation_' data.subjects{ii} '.sofa.MAX.mat'],'file')
                load([hpath 'Sphere_Rotation_' data.subjects{ii} '.sofa.MAX.mat']);
            else
                [~,tmp]=ziegelwanger2014(Obj,1,0,0);
                toaEst=tmp.toa;
                save([hpath 'Sphere_Rotation_' data.subjects{ii} '.sofa.MAX.mat'],'toaEst');
            end
            [~,results(ii).MAX{1}]=ziegelwanger2014(Obj,toaEst,0,1e-8);
            
            if exist([hpath 'Sphere_Rotation_' data.subjects{ii} '.sofa.CTD.mat'],'file')
                load([hpath 'Sphere_Rotation_' data.subjects{ii} '.sofa.CTD.mat']);
            else
                [~,tmp]=ziegelwanger2014(Obj,2,0,0);
                toaEst=tmp.toa;
                save([hpath 'Sphere_Rotation_' data.subjects{ii} '.sofa.CTD.mat'],'toaEst');
            end
            [~,results(ii).CTD{1}]=ziegelwanger2014(Obj,toaEst,0,1e-8);
            
            if exist([hpath 'Sphere_Rotation_' data.subjects{ii} '.sofa.AGD.mat'],'file')
                load([hpath 'Sphere_Rotation_' data.subjects{ii} '.sofa.AGD.mat']);
            else
                [~,tmp]=ziegelwanger2014(Obj,3,0,0);
                toaEst=tmp.toa;
                save([hpath 'Sphere_Rotation_' data.subjects{ii} '.sofa.AGD.mat'],'toaEst');
            end
            [~,results(ii).AGD{1}]=ziegelwanger2014(Obj,toaEst,0,1e-8);
             
            if exist([hpath 'Sphere_Rotation_' data.subjects{ii} '.sofa.MCM.mat'],'file')
                load([hpath 'Sphere_Rotation_' data.subjects{ii} '.sofa.MCM.mat']);
            else
                [~,tmp]=ziegelwanger2014(Obj,4,0,0);
                toaEst=tmp.toa;
                save([hpath 'Sphere_Rotation_' data.subjects{ii} '.sofa.MCM.mat'],'toaEst');
            end
            [~,results(ii).MCM{1}]=ziegelwanger2014(Obj,toaEst,0,1e-8);
            [~,results(ii).MCM{2}]=ziegelwanger2014(Obj,toaEst,[0.05 0.01],1e-8);
        end
        save([hpath '..' filesep '..' filesep 'experiments' filesep 'exp_ziegelwanger2014_SPHERE_ROT.mat'],'results');
        data.results=results;
    else
        tmp=load([hpath 'info.mat']);
        data=tmp.info.Rotation;
        load([hpath '..' filesep '..' filesep 'experiments' filesep 'exp_ziegelwanger2014_SPHERE_ROT.mat']);
        data.results=results;
    end
    
end

%% ARI database (NH89)
if flags.do_NH89

    data=SOFAload([hpath 'ARI_NH89.sofa']);
    
    if exist([hpath 'ARI_NH89.sofa.MAX.mat'],'file')
        tmp=load([hpath 'ARI_NH89.sofa.MAX.mat']);
        data.Data.toaEst{1}=tmp.toaEst;
    else
        [~,tmp]=ziegelwanger2014(data,1,0,0);
        toaEst=tmp.toa;
        save([hpath 'ARI_NH89.sofa.MAX.mat'],'toaEst');
        data.Data.toaEst{1}=toaEst;
    end
    
    if exist([hpath 'ARI_NH89.sofa.CTD.mat'],'file')
        tmp=load([hpath 'ARI_NH89.sofa.CTD.mat']);
        data.Data.toaEst{2}=tmp.toaEst;
    else
        [~,tmp]=ziegelwanger2014(data,2,0,0);
        toaEst=tmp.toa;
        save([hpath 'ARI_NH89.sofa.CTD.mat'],'toaEst');
        data.Data.toaEst{2}=toaEst;
    end
    
    if exist([hpath 'ARI_NH89.sofa.AGD.mat'],'file')
        tmp=load([hpath 'ARI_NH89.sofa.AGD.mat']);
        data.Data.toaEst{3}=tmp.toaEst;
    else
        [~,tmp]=ziegelwanger2014(data,3,0,0);
        toaEst=tmp.toa;
        save([hpath 'ARI_NH89.sofa.AGD.mat'],'toaEst');
        data.Data.toaEst{3}=toaEst;
    end
    
    if exist([hpath 'ARI_NH89.sofa.MCM.mat'],'file')
        tmp=load([hpath 'ARI_NH89.sofa.MCM.mat']);
        data.Data.toaEst{4}=tmp.toaEst;
    else
        [~,tmp]=ziegelwanger2014(data,4,0,0);
        toaEst=tmp.toa;
        save([hpath 'ARI_NH89.sofa.MCM.mat'],'toaEst');
        data.Data.toaEst{4}=toaEst;
    end
    
end

%% Sphere
if flags.do_Sphere

    data=SOFAload([hpath 'Sphere.sofa']);
    
    if exist([hpath 'Sphere.sofa.MAX.mat'],'file')
        tmp=load([hpath 'Sphere.sofa.MAX.mat']);
        data.Data.toaEst{1}=tmp.toaEst;
    else
        [~,tmp]=ziegelwanger2014(data,1,0,0);
        toaEst=tmp.toa;
        save([hpath 'Sphere.sofa.MAX.mat'],'toaEst');
        data.Data.toaEst{1}=toaEst;
    end
    
    if exist([hpath 'Sphere.sofa.CTD.mat'],'file')
        tmp=load([hpath 'Sphere.sofa.CTD.mat']);
        data.Data.toaEst{2}=tmp.toaEst;
    else
        [~,tmp]=ziegelwanger2014(data,2,0,0);
        toaEst=tmp.toa;
        save([hpath 'Sphere.sofa.CTD.mat'],'toaEst');
        data.Data.toaEst{2}=toaEst;
    end
    
    if exist([hpath 'Sphere.sofa.AGD.mat'],'file')
        tmp=load([hpath 'Sphere.sofa.AGD.mat']);
        data.Data.toaEst{3}=tmp.toaEst;
    else
        [~,tmp]=ziegelwanger2014(data,3,0,0);
        toaEst=tmp.toa;
        save([hpath 'Sphere.sofa.AGD.mat'],'toaEst');
        data.Data.toaEst{3}=toaEst;
    end
    
    if exist([hpath 'Sphere.sofa.MCM.mat'],'file')
        tmp=load([hpath 'Sphere.sofa.MCM.mat']);
        data.Data.toaEst{4}=tmp.toaEst;
    else
        [~,tmp]=ziegelwanger2014(data,4,0,0);
        toaEst=tmp.toa;
        save([hpath 'Sphere.sofa.MCM.mat'],'toaEst');
        data.Data.toaEst{4}=toaEst;
    end
    
end

%% Sphere and Torso
if flags.do_SAT

    data=SOFAload([hpath 'SAT.sofa']);
    
    if exist([hpath 'SAT.sofa.MAX.mat'],'file')
        tmp=load([hpath 'SAT.sofa.MAX.mat']);
        data.Data.toaEst{1}=tmp.toaEst;
    else
        [~,tmp]=ziegelwanger2014(data,1,0,0);
        toaEst=tmp.toa;
        save([hpath 'SAT.sofa.MAX.mat'],'toaEst');
        data.Data.toaEst{1}=toaEst;
    end
    
    if exist([hpath 'SAT.sofa.CTD.mat'],'file')
        tmp=load([hpath 'SAT.sofa.CTD.mat']);
        data.Data.toaEst{2}=tmp.toaEst;
    else
        [~,tmp]=ziegelwanger2014(data,2,0,0);
        toaEst=tmp.toa;
        save([hpath 'SAT.sofa.CTD.mat'],'toaEst');
        data.Data.toaEst{2}=toaEst;
    end
    
    if exist([hpath 'SAT.sofa.AGD.mat'],'file')
        tmp=load([hpath 'SAT.sofa.AGD.mat']);
        data.Data.toaEst{3}=tmp.toaEst;
    else
        [~,tmp]=ziegelwanger2014(data,3,0,0);
        toaEst=tmp.toa;
        save([hpath 'SAT.sofa.AGD.mat'],'toaEst');
        data.Data.toaEst{3}=toaEst;
    end
    
    if exist([hpath 'SAT.sofa.MCM.mat'],'file')
        tmp=load([hpath 'SAT.sofa.MCM.mat']);
        data.Data.toaEst{4}=tmp.toaEst;
    else
        [~,tmp]=ziegelwanger2014(data,4,0,0);
        toaEst=tmp.toa;
        save([hpath 'SAT.sofa.MCM.mat'],'toaEst');
        data.Data.toaEst{4}=toaEst;
    end
    
end

%% Sphere, Torso and Pinna
if flags.do_STP

    data=SOFAload([hpath 'STP.sofa']);
    
    if exist([hpath 'STP.sofa.MAX.mat'],'file')
        tmp=load([hpath 'STP.sofa.MAX.mat']);
        data.Data.toaEst{1}=tmp.toaEst;
    else
        [~,tmp]=ziegelwanger2014(data,1,0,0);
        toaEst=tmp.toa;
        save([hpath 'STP.sofa.MAX.mat'],'toaEst');
        data.Data.toaEst{1}=toaEst;
    end
    
    if exist([hpath 'STP.sofa.CTD.mat'],'file')
        tmp=load([hpath 'STP.sofa.CTD.mat']);
        data.Data.toaEst{2}=tmp.toaEst;
    else
        [~,tmp]=ziegelwanger2014(data,2,0,0);
        toaEst=tmp.toa;
        save([hpath 'STP.sofa.CTD.mat'],'toaEst');
        data.Data.toaEst{2}=toaEst;
    end
    
    if exist([hpath 'STP.sofa.AGD.mat'],'file')
        tmp=load([hpath 'STP.sofa.AGD.mat']);
        data.Data.toaEst{3}=tmp.toaEst;
    else
        [~,tmp]=ziegelwanger2014(data,3,0,0);
        toaEst=tmp.toa;
        save([hpath 'STP.sofa.AGD.mat'],'toaEst');
        data.Data.toaEst{3}=toaEst;
    end
    
    if exist([hpath 'STP.sofa.MCM.mat'],'file')
        tmp=load([hpath 'STP.sofa.MCM.mat']);
        data.Data.toaEst{4}=tmp.toaEst;
    else
        [~,tmp]=ziegelwanger2014(data,4,0,0);
        toaEst=tmp.toa;
        save([hpath 'STP.sofa.MCM.mat'],'toaEst');
        data.Data.toaEst{4}=toaEst;
    end
    
end