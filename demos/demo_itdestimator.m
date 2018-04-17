function demo_itdestimator(varargin)
%   Usage: demo_itdestimator(flag) 
%
%   Input parameters:
% 
%       The following flags can be chosen:
%       'clean','noisy','stim' to create the data which was used
%                            &
%       'fig1,fig2,..,fig9 to create the plots shown in the paper 
%
% 
%   `demo_itdestimator` reproduces the results and figures of the technical report.
%
% 
%   Requirements: 
%   -------------
%
%   1) SOFA API from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
% 
%
%   Examples:
%   ---------
%
%   To display results for Fig.1 use :::
%
%     demo_itdestimator('fig1');
%
%   To display results for Fig.2 use :::
%
%     demo_itdestimator('fig2');
%
%   RB: CONTINUE ACCORDINGLY WITH ALL FIGURES...
% 
% 
% Copyright (C) 2009-2015 Piotr Majdak and the AMT team.
% This file is part of Auditory Modeling Toolbox (AMT) version 0.9.9
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% AUTHOR: Laurin Steidle


    
definput.import={'amt_cache'};
definput.flags.wat = {'clean','noisy','stim',...
    'fig1','fig2','fig3','fig4','fig5','fig6','fig7','fig8','fig9'};
definput.import={'amt_cache'}; % get the flags of amt_cache
[flags,~]=ltfatarghelper({},definput,varargin);

%% ---------------calculate-estimation-data-------------------------

% -------------------------------------------------------------------------
% calulates ITDs with all models for 'clean' synthesised IRs    
    if flags.do_clean

        modes = {'Threshold','Cen_e2','MaxIACCr', 'MaxIACCe',...
                'CenIACCr', 'CenIACCe', 'CenIACC2e', 'PhminXcor','IRGD'};
        n_modes = size(modes,2);

        
        tmp=amt_load('ziegelwanger2014','info.mat');
        % importing head parameter and unify units
        data=tmp.info.Rotation;  
        data.radius = data.radius*1e-3;
        data.phi = data.phi*pi/180;
        data.theta  = data.theta*pi/180;

        avg_headsize_idx = 15:15+13;
        n_heads = size(avg_headsize_idx,2);
        
        % initiate array
        msq = zeros(n_heads,n_modes); 
        
        for ii=avg_headsize_idx % goes thru all head geometries
            Obj=SOFAload(fullfile(SOFAdbPath, 'ziegelwanger2014',...
                [ 'Sphere_Rotation_' data.subjects{ii} '.sofa']));  
            
            %calculate theoretically expected ITDs
            p0_offaxis = [[data.radius(ii);0; 0; 0; 0;data.phi(ii)+pi/2;  data.theta(ii)*pi/180]...
                          [data.radius(ii);0; 0; 0; 0;data.phi(ii)-pi/2; -data.theta(ii)*pi/180]];
            disp(p0_offaxis)
            toa_sim_left    = ziegelwanger2014_offaxis(p0_offaxis(:,1),Obj.SourcePosition(:,1:2)*pi/180);
            toa_sim_right   = ziegelwanger2014_offaxis(p0_offaxis(:,2),Obj.SourcePosition(:,1:2)*pi/180);
            toadiff_sim     = toa_sim_left - toa_sim_right;  
            
            for kk = 1:n_modes % goes thru all estimation methods
                %calculate estimated ITDs
                toadiff = itdestimator(Obj,modes{kk},'bb'); 
                msq(ii-n_heads,kk) = sqrt(sum( (toadiff-toadiff_sim).^2 )/size(toadiff,1));         
            end
        end

        msq_summed = sum(msq,1)/size(msq,1);
        disp(modes);
        disp(msq_summed);
        save('clean_msq','msq')
    end



% -------------------------------------------------------------------------
    if flags.do_noisy
        
        modes = {'Threshold','Cen_e2','MaxIACCr', 'MaxIACCe',...
                'CenIACCr', 'CenIACCe', 'CenIACC2e', 'PhminXcor','IRGD'};
        n_modes = size(modes,2);
        
        % importing head parameter and unify units
        tmp=amt_load('ziegelwanger2014','info.mat');
        data=tmp.info.Rotation;  
        data.radius = data.radius*1e-3;
        data.phi = data.phi*pi/180;
        data.theta  = data.theta*pi/180;
        
        avg_headsize_idx = 15:15+13;
        n_heads = size(avg_headsize_idx,2);
        
        %define added noise levels of noise in dB of SNR
        noise_lvl = [50,40,30,20];
                
        % initiate array                
        msq = zeros(n_heads,size(noise_lvl,2),n_modes); 
        for ii=avg_headsize_idx % goes thru all head geometries
            for jj = 1:size(noise_lvl,2) % level of noise
                Obj=SOFAload(fullfile(SOFAdbPath, 'ziegelwanger2014',...
                            [ 'Sphere_Rotation_' data.subjects{ii} '.sofa']));

                % calculate theoretically expected ITDs
                p0_offaxis = [[data.radius(ii);0; 0; 0; 0;data.phi(ii)+pi/2;  data.theta(ii)*pi/180]...
                              [data.radius(ii);0; 0; 0; 0;data.phi(ii)-pi/2; -data.theta(ii)*pi/180]];
                disp(p0_offaxis)
                toa_sim_left    = ziegelwanger2014_offaxis(p0_offaxis(:,1),Obj.SourcePosition(:,1:2)*pi/180);
                toa_sim_right   = ziegelwanger2014_offaxis(p0_offaxis(:,2),Obj.SourcePosition(:,1:2)*pi/180);
                toadiff_sim     = toa_sim_left - toa_sim_right;  
                
                % adds noise
                for kk=1:Obj.API.M
                    for mm=1:Obj.API.R
                        ir_tmp = squeeze(Obj.Data.IR(kk,mm,:));
                        Obj.Data.IR(kk,mm,:) = ir_tmp + ...
                          setdbspl(noise(length(ir_tmp),1),dbspl(ir_tmp)-noise_lvl(jj));
%                         awgn(squeeze(Obj.Data.IR(kk,mm,:)),noise_lvl(jj));
                    end
                end

                for kk = 1:n_modes % goes thru all estimation methods
                    %calculate estimated ITDs
                    toadiff = itdestimator(Obj,modes{kk},'bb');
                    msq(ii-n_heads,jj,kk) = sqrt(sum( (toadiff-toadiff_sim).^2 )/size(toadiff,1));      
                end
            end
        end

        msq_summed = sum(msq,1)/size(msq,1);
        disp(modes);
        disp(msq_summed);
        save('noisy_msq','msq')
    end

    
% -------------------------------------------------------------------------
% creates data based on stimuli
% this is done to compare models found in itdestimator to the dietz model

    if flags.do_stim
        
        modes = {'MaxIACCr', 'MaxIACCe','IRGD','Dietz'};
        n_modes = size(modes,2);
        
        % importing head parameter and unify units
        tmp=amt_load('ziegelwanger2014','info.mat');
        data=tmp.info.Rotation;  
        data.radius = data.radius*1e-3;
        data.phi = data.phi*pi/180;
        data.theta  = data.theta*pi/180;

        avg_headsize_idx = 15:15+13;
        n_heads = size(avg_headsize_idx,2);
        
        % initiate array                    
        msq = zeros(n_heads,n_modes);

        for ii=avg_headsize_idx
            Obj=SOFAload(fullfile(SOFAdbPath, 'ziegelwanger2014',...
                        [ 'Sphere_Rotation_' data.subjects{ii} '.sofa']));
            % creates white noise and adjust object samlpe length
            fs  = Obj.Data.SamplingRate;
            white_noise = noise(fs/2);
            Obj.API.N = size(Obj.Data.IR,3) + size(white_noise,1) - 1;
       
            % adds white noise
            stim = zeros(Obj.API.M, Obj.API.R, Obj.API.N);
            for kk=1:Obj.API.M      %pos
                for mm=1:Obj.API.R  %ear
                    stim(kk,mm,:) = lconv(Obj.Data.IR(kk,mm,:),white_noise);
                end
            end
            Obj.Data.IR = stim;

            % calculate theoretically expected ITDs
            p0_offaxis = [[data.radius(ii);0; 0; 0; 0;data.phi(ii)+pi/2;  data.theta(ii)*pi/180]...
                          [data.radius(ii);0; 0; 0; 0;data.phi(ii)-pi/2; -data.theta(ii)*pi/180]];
            disp(p0_offaxis)
            toa_sim_left    = ziegelwanger2014_offaxis(p0_offaxis(:,1),Obj.SourcePosition(:,1:2)*pi/180);
            toa_sim_right   = ziegelwanger2014_offaxis(p0_offaxis(:,2),Obj.SourcePosition(:,1:2)*pi/180);
            toadiff_sim     = toa_sim_left - toa_sim_right;        

            % calculate estimated ITDs
            for kk = 1:n_modes    
               if strcmp(modes{kk},'Dietz')
               fprintf('Dietz \n')
                    for nn=1:Obj.API.M
                       insig = transpose(squeeze(stim(nn,:,:)));
                       dietz_out = dietz2011(insig,fs);
                       med = median(dietz_out.itd_lp);
                       toadiff(nn) = mean(med(1:8));
                       toadiff = transpose(toadiff);
                    end
                else
                   toadiff = itdestimator(Obj,modes{kk},'stim','lp');
                   
               end
                             
                msq(ii,kk) = sqrt(sum( (toadiff-toadiff_sim).^2 )/size(toadiff,1));
            end
            
            disp(modes);
            disp(msq_summed);
            save('stim_msq','msq')
        end
    end
    
%% ---------------figure-1------------------------------------------
% creates the threshold estmation model visualisation
    if flags.do_fig1
        nn = 115;
        col = jet(5);
        Obj = SOFAload(fullfile(SOFAdbPath,'baumgartner2017','hrtf b_nh14.sofa'));
        fs = Obj.Data.SamplingRate;
        
        fig = figure;

        [a1,~] = findpeaks(squeeze(Obj.Data.IR(nn,1,:)),'MinPeakProminence',0.005);
        th_value1 = max(a1)*0.2;
        [x1,~] = find(squeeze(Obj.Data.IR(nn,1,:))>th_value1,1);

        [a2,~] = findpeaks(squeeze(Obj.Data.IR(nn,2,:)),'MinPeakProminence',0.005);
        th_value2 = max(a2)*0.2;
        [x2,~] = find(squeeze(Obj.Data.IR(nn,2,:))>th_value2,1);

        ax1 = subplot(2,1,1);
        hold on
        plot((0:size(squeeze(Obj.Data.IR(nn,1,:)),1)-1)/fs*1e3,...
            squeeze(Obj.Data.IR(nn,1,:)),'k-')
        plot(ax1,[x1/fs*1e3; x1/fs*1e3],[-1; 1],'Color',col(5,:))
        plot(ax1,[x2/fs*1e3; x2/fs*1e3],[-1; 1],'Color',col(3,:))
        plot(ax1,[0; size(Obj.Data.IR(nn,1,:),3)/fs*1e3],[th_value1; th_value1],'Color',col(5,:))
        hold off
        title(ax1,'Left ear')
        xlim(ax1,[0 3])
        set(ax1,'Ytick',[-0.03,0,0.03])
        xlabel('Time (ms)')
        ylim(ax1,[-0.045 0.045])
        ylabel('Amplitude (a.u)')

        ax2 = subplot(2,1,2); 
        hold on
        plot((0:size(squeeze(Obj.Data.IR(nn,2,:)),1)-1)/fs*1e3,...
            squeeze(Obj.Data.IR(nn,2,:)),'k-')
        plot(ax2,[x2/fs*1e3; x2/fs*1e3],[-1; 1],'Color',col(5,:))
        plot(ax2,[0; size(Obj.Data.IR(nn,2,:),3)/fs*1e3],[th_value2; th_value2],'Color',col(5,:))
        hold off
        title(ax2,'Right ear')
        xlim(ax2,[0 3])
        set(ax2,'Ytick',[-0.03,0,0.03])
        xlabel('Time (ms)')
        ylim(ax2,[-0.045 0.045])
        ylabel('Amplitude (a.u.)')
        

    end

%% ---------------figure-2------------------------------------------
% creates cross correlation based estimatation visualisation
     if flags.do_fig2

        tmp=amt_load('ziegelwanger2014','info.mat');
        data=tmp.info.Rotation;  
        data.radius = data.radius*1e-3;
        data.phi = data.phi*pi/180;
        data.theta  = data.theta*pi/180;       
%         n_heads = size(data.subjects,2);
                
        fig = figure();
        hold on
        for ii=3:3 %n_heads
            jj = 115;
            Obj=SOFAload(fullfile(SOFAdbPath, 'ziegelwanger2014',...
                [ 'Sphere_Rotation_' data.subjects{ii} '.sofa']));

            plot(transpose(-239:239)/Obj.Data.SamplingRate*1e3,...
                envelope(xcorr(squeeze(Obj.Data.IR(jj,1,:)),squeeze(Obj.Data.IR(jj,2,:)))))
            plot(transpose(-239:239)/Obj.Data.SamplingRate*1e3,...
                xcorr(squeeze(Obj.Data.IR(jj,1,:)),squeeze(Obj.Data.IR(jj,2,:))))
        end
        %title( 'Examplary inter aural cross correlation' )
        xlabel('time delay (ms)')
        ylabel('correlation (a.u.)')
        legend('cross correlation','envelope of cc','Location','best')
        xlim([-1 1])
        hold off
        
     end
     
%% ---------------figure-3------------------------------------------
% creates group delay based estimatation visualisation
    if flags.do_fig3
        

        tmp=amt_load('ziegelwanger2014','info.mat');
        data=tmp.info.Rotation;  
%         n_heads = size(data.subjects,2);
                
        fig = figure();
        hold on
        for ii=5:5 %n_heads
            jj = 115;
            Obj=SOFAload(fullfile(SOFAdbPath, 'ziegelwanger2014',...
                [ 'Sphere_Rotation_' data.subjects{ii} '.sofa']));
            Ns  = Obj.API.N;
            fs  = Obj.Data.SamplingRate;
            [gd,w] = grpdelay(squeeze( Obj.Data.IR(jj,1,:) ),1,Ns,fs );
            plot(w,gd/fs*1e3)
            
        end
        title( 'Exemplary group delay' )
        xlabel('Frequency (Hz)')
        ylabel('Group delay (ms)')
        xlim([500,5000])
        ylim([1.35 1.8])
        hold off
        
    end

%% ---------------figure-4------------------------------------------
% -------------------------------------------------------------------------

    if flags.do_fig4

        Obj = SOFAload(fullfile(SOFAdbPath,'baumgartner2017','hrtf b_nh15.sofa'));
        white_noise = noise(258);

        % ---------------------- calculating toa ----------------------------------

        modes = {'Threshold','Cen_e2','MaxIACCr', 'MaxIACCe',...
                'CenIACCr', 'CenIACCe', 'CenIACC2e', 'PhminXcor','IRGD','Dietz'};
        n_modes = size(modes,2);
        cc = hsv(n_modes);

        plane_idx = find( Obj.SourcePosition(:,2) == 0 );
        plane_angle = Obj.SourcePosition(plane_idx,1);
        [n_angles,~] = size(plane_angle);

        for kk=1:n_modes
            if strcmp(modes{kk},'Dietz')
                fprintf('Dietz \n')
                for nn=1:n_angles
                     stimulus = SOFAspat(white_noise,Obj,plane_angle(nn),0);
                     dietz_out = dietz2011(stimulus,Obj.Data.SamplingRate);
                     med = median(dietz_out.itd);
                     toa_diff{kk}(nn) = mean(med(1:8));
                end
            else
            toa_diff{kk} = itdestimator(Obj,modes{kk});  
            end
        end
    
        fig = figure();
        hold on
        for kk=1:n_modes
            if strcmp(modes{kk},'Dietz') || strcmp(modes{kk},'Dietz lp')
                plot(plane_angle,toa_diff{kk},'color', cc(kk,:))
            else
                plot(plane_angle,toa_diff{kk}(plane_idx),'color', cc(kk,:))    
            end
        end
        
        xlabel('azimuthal angle (Deg)')
        ylabel('interaural time difference (s)')
        xlim([0,360])
        legend('Threshold','Cen_e2','MaxIACCr', 'MaxIACCe',...
                'CenIACCr', 'CenIACCe', 'CenIACC2e', 'PhminXcor',...
                'IRGD','Dietz','Location','eastoutside') 
        hold off
        
    end
%% ---------------figure-5------------------------------------------
% -------------------------------------------------------------------------

    
    if flags.do_fig5

    %load data
        data=data_ziegelwanger2014('SPHERE_ROT',flags.cachemode);
        for ii=1:length(data.results)
            p_onaxis{1}(:,:,ii)=data.results(ii).MAX{1}.p_onaxis;
            p_onaxis{2}(:,:,ii)=data.results(ii).CTD{1}.p_onaxis;
            p_onaxis{3}(:,:,ii)=data.results(ii).AGD{1}.p_onaxis;
            p_onaxis{4}(:,:,ii)=data.results(ii).MCM{1}.p_onaxis;
        end

    %plot figure
        fig = figure();

        lw2=2;%linewidth2
        ls='kx';%linestyle
        lc=[0 0 0];
        
        %phi
        subplot(211);
        var=[squeeze(p_onaxis{4}(2,1,:))/pi*180 ...
            squeeze(p_onaxis{4}(2,2,:))/pi*180 ...
            data.phi+ones(length(data.phi),1)*90 ...
            data.phi-ones(length(data.phi),1)*90];

        hold on
        for ch=1:size(p_onaxis{4},2)
            plot(15:23,var(15:23,2+ch),ls,'Linewidth',lw2,'color',lc)
            plot(24:28,var(24:28,2+ch),ls,'Linewidth',lw2,'color',lc)
        end

        clear var;
        set(gca,'YTick',[-90 90])
        set(gca,'xtick',1:42)
        xlim([15-1,15+14])
        set(gca,'Xticklabel',{''})
        ylabel('\phi (deg)')
        grid on
        set(gca,'GridLineStyle','--')
        
        
        %theta
        subplot(212);
        var=[squeeze(p_onaxis{4}(3,1,:))/pi*180 ...
            squeeze(p_onaxis{4}(3,2,:))/pi*180 ...
            data.theta ...
            -data.theta];

        hold on
        for ch=1:size(p_onaxis{4},2)
            plot(15:23,var(15:23,2+ch),ls,'Linewidth',lw2,'color',lc)
            plot(24:28,var(24:28,2+ch),ls,'Linewidth',lw2,'color',lc)           
        end

        clear var;
        ylabel('\theta (deg)')
        xlabel('Condition')
        ylim([-15,15])
        xlim([15-1,15+14])
        set(gca,'xtick',1:42)
        set(gca,'xticklabel',[])
        set(gca,'ytick',[-10 0 10])
        grid on
        set(gca,'GridLineStyle','--')
        
        %set(gcf,'color',[1 1 1]) 
    end


%% ---------------figure-6------------------------------------------
% -------------------------------------------------------------------------    
    if flags.do_fig6

        tmp=amt_load('ziegelwanger2014','info.mat');
        data=tmp.info.Rotation;  
        data.radius = data.radius*1e-3;
        data.phi = data.phi*pi/180;
        data.theta  = data.theta*pi/180;

        noise_lvl = [20,30,40];

        n_noise = size(noise_lvl,2);
        
        
        fig = figure();
        hold on
        for jj = 1:n_noise % level of noise
            Obj=SOFAload(fullfile(SOFAdbPath, 'ziegelwanger2014',...
                        [ 'Sphere_Rotation_' data.subjects{4} '.sofa']));

            IR = squeeze(Obj.Data.IR(777,1,:));
            IR = IR + setdbspl(noise(length(IR),1),dbspl(IR)-noise_lvl(jj));%%awgn(IR,noise_lvl(jj));

            plot(IR)

        end
        %title( 'Impulse response with added gaussian noise' )
        xlabel('Sample point')
        xlim([0 240])
        ylabel('Amplitude (a.u.)')
        legend('SNR 20dB','SNR 30dB','SNR 40dB','Location','best')
        hold off
        
    end
    
    
%% ---------------figure-7------------------------------------------
% tilted heads

    if flags.do_fig7
%         msq = amt_cache('get','clean_bb_msq',flags.cachemode);
%         if isempty(msq)
%           msq = demo_itdestimator('clean');
%           amt_cache('set','clean_bb_msq',msq);
%         end
        load('clean_bb_msq.mat')
        msq_bb_clean = msq(15:15+14-1,:);
        
        fig = figure();
        bar(msq_bb_clean(:,1,1)*1e6)
        xlabel('Heads as described in section 3.1')
        ylabel('ANR in \mus')
        colormap('lines')
        set(gca, 'YMinorGrid','on', 'YMinorGrid','on')
        
    end
    
    
%% ---------------figure-8------------------------------------------
% results1

    if flags.do_fig8
        load('clean_bb_msq.mat')
        msq_bb_clean = msq(15:15+14-1,:);
        load('clean_lp_msq.mat')
        msq_lp_clean = msq(15:15+14-1,:);
        
        load('noisy_bb_msq.mat')
        msq_bb_noisy = msq;
        load('noisy_lp_msq.mat')
        msq_lp_noisy = msq;
        
        msq_bb_clean = reshape(msq_bb_clean,14,1,9);
        msq_bb = cat(2,msq_bb_clean,msq_bb_noisy);
        msq_bb = msq_bb(1:end-5,:,:);
        msq_bb_std = squeeze( std(msq_bb) );
        msq_bb_sum = squeeze(sum(msq_bb,1)/size(msq_bb,1));
                
        msq_lp_clean = reshape(msq_lp_clean,14,1,9);
        msq_lp = cat(2,msq_lp_clean,msq_lp_noisy);
        msq_lp = msq_lp(1:end-5,:,:);
        msq_lp_std = squeeze( std(msq_lp )); 
        msq_lp_sum = squeeze(sum(msq_lp,1)/size(msq_lp,1));
        
        
        msq_int_sum = reshape(cat(2,cat(1,cat(1,msq_bb_sum,msq_lp_sum),zeros(5,9)),zeros(15,1)),5,30);
        msq_int_std = reshape(cat(2,cat(1,cat(1,msq_bb_std,msq_lp_std),zeros(5,9)),zeros(15,1)),5,30);
        M = size(msq_int_sum,1);
        N = size(msq_int_sum,2);
        K = numel(msq_int_sum);
        g = 0; k=0;

        bb_idx = 1:3:N-4;
        lp_idx = 2:3:N-4;
        zo_idx = [3:3:N-4,28,29,30];
%         len = size(bb_idx,2);
        
        fig = figure();
        hold on
        col = parula(10);
        for ii=1:M
            for jj = 1:N
                if any(bb_idx == jj)
                    handle(mod(g,9)+1) = barh(K,msq_int_sum(ii,jj)*1e6,1,'FaceColor',col(mod(g,9)+1,:));
%                     errorbar(K,msq_int_sum(ii,jj),msq_int_std(ii,jj),...
%                         'MarkerEdgeColor','k','MarkerFaceColor','w')
                    g = g+1;
                elseif any(lp_idx == jj)
                    barh(K,msq_int_sum(ii,jj)*1e6,1,'FaceColor',[1,1,1]*0.2)%col(mod(k,9)+1,:),'EdgeColor','r')
%                     errorbar(K,msq_int_sum(ii,jj),msq_int_std(ii,jj),...
%                         'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63])
                    k = k+1;
                end 
                K = K-1;
            end
        end
        xlim([5e0,2e3])
        xlabel('ANR in \mus')
        ylabel('SNR in dB')
        set(gca, 'XScale', 'log')
        set(gca, 'yTickLabel',{'20 dB','30 dB','40 dB','50 dB','clean'},'YTick',[18,48,78,108,138])
        set(gca, 'XMinorGrid','on', 'XMinorGrid','on')
               
        legend(handle, 'Threshold','Cen_e2','MaxIACCr', 'MaxIACCe',...
                    'CenIACCr', 'CenIACCe', 'CenIACC2e', 'PhminXcor','IRGD','Location','northeast')
        
       
    end
    
    
%% ---------------figure-9------------------------------------------
% results2

    if flags.do_fig9
        
        load('stim_msq_bb.mat')
        msq_stim_bb = msq(1:end-5,:);
        load('stim_msq_lp.mat')
        msq_stim_lp = msq(1:end-5,:);
        
        msq_stim = cat(2,msq_stim_bb,msq_stim_lp);
        msq_stim_sum = squeeze(sum(msq_stim,1)/size(msq_stim,1));
        msq_stim_std = std(msq_stim);
        
        
       
        fig = figure();
        hold on
        modes = {'MaxIACCr', 'MaxIACCe','IRGD','Dietz','MaxIACCr lp', 'MaxIACCe lp ','IRGD lp ','Dietz lp'};
        barh(1:size(msq_stim_sum,2),msq_stim_sum*1e6)
        errorbar(msq_stim_sum*1e6,1:size(msq_stim_sum,2),msq_stim_std*1e6,'horizontal','k.')
%                         'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63])
        xlim([0,400])
        xlabel('ANR in \mus')
        
        set(gca, 'YTickLabel',modes, 'YTick',1:size(msq_stim_sum,2))
        colormap('lines')
        set(gca, 'XMinorGrid','on', 'XGrid','on')
        hold off
        
                
    end 
end
