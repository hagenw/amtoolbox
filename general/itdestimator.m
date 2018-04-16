function [toa_diff,toa] = itdestimator(Obj,varargin)

%   Usage: itd = itd_estimator(data,[mode],[threshlvl],...
%                   [lowpass],[butterpoly],[upper_cutfreq]) 
%
%   Input parameters:
% 
%       data:       SOFA object
% 
%       mode:       (optional) Select one estimation methods
%                   ('Threshold','Cen_e2','MaxIACCr', 'MaxIACCe', 'CenIACCr',.. 
%                    'CenIACCe', 'CenIACC2e', 'PhminXcor','IRGD')
%             
%                   For details concerning estimation methods see:
%                   'http://asa.scitation.org/doi/10.1121/1.4996457'
%
%       threshlvl:  (optional)                  
%                   Set threshold level for 'Threshold' mode in dB
% 
%       lowpass:    (optional) Decide if lowpass shall be applied. 
%                   'lp' for lowpass, 'bb' for broadband
%
%       butterpoly: (optional) Select the order of the polynom
%                   applied in the butterworth filter. ( 2 =< i =< 6 )
% 
%       upper_cutfreq: (optional) Set frequency of lowpass cutoff in Hz      
% 
% 
%   Output parameters:
% 
%       itd:        interaural time difference in seconds
% 
% 
%   Purpose:
%   Estimates the ITD based on biaural impulse responses.
%   Several different estimaton methods can be chosen.
%   MaxIAACe is recommended.
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
%   Obj = SOFAload(fullfile(SOFAdbPath,'baumgartner2017','hrtf b_nh15.sofa'));
%   toa_diff = itdestimator(Obj,'MaxIACCe','lp','upper_cutfreq',3000)
%   
%   With these settings the estimator uses the MaxIAAce method and applies
%   a lowpass with a cut off frequency of 3kH.
%   
%   The output array is structured as the SOFA Data.IR
%   If you would like to select for example only data on the horrizontal
%   plane you could:
%
%   plane_idx = find( Obj.SourcePosition(:,2) == 0 );
%   plane_angle = Obj.SourcePosition(plane_idx,1);
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


    % ---------------------- ltfatarghelper -------------------------------

    definput.flags.mode = {'Threshold','Cen_e2','MaxIACCr', 'MaxIACCe',...
        'CenIACCr', 'CenIACCe', 'CenIACC2e', 'PhminXcor','IRGD'};
    definput.flags.lp = {'lp','bb'};
    definput.flags.peak = {'hp','fp'};
    definput.flags.toaguess = {'noguess','guesstoa'};
    
    definput.keyvals.threshlvl = -10;
    definput.keyvals.butterpoly = 6;
    definput.keyvals.upper_cutfreq = 3000;
    definput.keyvals.lower_cutfreq = 1000;
    definput.keyvals.avgtoa = 45;

    [flags,kv]=ltfatarghelper({},definput,varargin);

    
    % ---------------------- renaming input parameter ---------------------
    

    pos = Obj.API.M;
    ear = Obj.API.R;
    IR = Obj.Data.IR;
    fs  = Obj.Data.SamplingRate;
    Ns  = Obj.API.N;

    % ---------------------- initialising variables -----------------------

    toa = zeros(pos,ear);
    toa_diff = zeros(pos,1);


    % ---------------------- Applying low-pass ----------------------------

    if flags.do_lp
        fprintf('Applying Butterworth low pass \n')
        fprintf(strcat('Polynomial order of Butterworth filter: '...
            ,num2str(kv.butterpoly),'\n'))
        fprintf(strcat('Cut of frequency is: '...
            ,num2str(kv.upper_cutfreq),' Hz\n\n'))
        cut_off_freq_norm = norm_cut_off_freq(kv.upper_cutfreq,IR,fs);
        f_ir = zeros(pos,ear,Ns);
        for ii=1:pos
            for jj=1:ear  
                sir = squeeze( IR(ii,jj,:) );
                f_sir = butter_lp(sir, kv.butterpoly,...
                    cut_off_freq_norm );
                f_ir(ii,jj,:) = f_sir;
            end
        end
    
    else
        fprintf('No low pass filter is applied \n\n')
        f_ir = IR;
    end

    % ---------------------- estimating itd -------------------------------
    % ---------------------------------------------------------------------

    % ---------------------- Threshold ------------------------------------
    switch(flags.mode)
        case 'Threshold'
            fprintf('Threshold mode \n')
            fprintf(strcat('Threshold level is: ',...
                num2str(kv.threshlvl),'dB \n\n'))
            
            if flags.do_fp
                for ii=1:pos
                     for jj=1:ear
                        indB = 0.5*mag2db(squeeze(f_ir(ii,jj,:)).^2);
                        [~,B] = findpeaks(indB);
                        th_value = indB(B(1)) + kv.threshlvl;
                        toa(ii,jj) = find(indB>th_value,1);
                    end
                    toa_diff(ii) = toa(ii,1) - toa(ii,2);     
                end

            else
                for ii=1:pos
                    for jj=1:ear
                        indB = 0.5*mag2db(squeeze(f_ir(ii,jj,:)).^2);
                        th_value = indB(indB == max(indB)) + kv.threshlvl;
                        toa(ii,jj) = find(indB>th_value,1);
                    end
                    toa_diff(ii) = toa(ii,1) - toa(ii,2);     
                end
            end
                        

    % ---------------------- Cross-Correlation ----------------------------        
        case 'Cen_e2'
            fprintf('Cen-e2 mode \n')
            for ii=1:pos
                for jj = 1:ear
                    e_sir_sq = envelope(squeeze(f_ir(ii,jj,:))).^2;
                    toa(ii,jj) = centroid(transpose(1:Ns),e_sir_sq);
                end
                toa_diff(ii) = toa(ii,1) - toa(ii,2);     
            end        
            

        case 'MaxIACCr'
            fprintf('MaxIACCr mode \n')
            for ii=1:pos                
                cc = xcorr(squeeze(f_ir(ii,1,:)),squeeze(f_ir(ii,2,:)));
                [~,idx_lag] = max(abs(cc));
                toa_diff(ii) = idx_lag - Ns;                
            end
            if flags.do_guesstoa
                toa = guesstoa(toa_diff,toa, kv.avgtoa);
            end
            
        case 'MaxIACCe'
            fprintf('MaxIACCe mode \n')
            for ii=1:pos
                e_sir1 = envelope(squeeze(f_ir(ii,1,:)));
                e_sir2 = envelope(squeeze(f_ir(ii,2,:)));
                cc = xcorr(e_sir1,e_sir2);
                [~,idx_lag] = max(abs(cc));
                toa_diff(ii) = idx_lag - Ns;
            end
            if flags.do_guesstoa
                toa = guesstoa(toa_diff,toa, kv.avgtoa);
            end
            

        case 'CenIACCr'
            fprintf('CenIACCr mode \n')
            x = transpose(1:(Ns*2-1));
            for ii=1:pos                
                cc = xcorr(squeeze(f_ir(ii,1,:)),squeeze(f_ir(ii,2,:)));
                pos_cc = abs(cc);
                toa_diff(ii) = centroid(x,pos_cc)-Ns;
            end
            if flags.do_guesstoa
                toa = guesstoa(toa_diff,toa, kv.avgtoa);
            end
            

        case 'CenIACCe'
            fprintf('CenIACCe mode \n')
            x = transpose(1:(Ns*2-1));
            for ii=1:pos
                e_sir1 = envelope(squeeze(f_ir(ii,1,:)));
                e_sir2 = envelope(squeeze(f_ir(ii,2,:)));
                cc = xcorr(e_sir1,e_sir2);
                toa_diff(ii) = centroid(x,abs(cc))-Ns;
            end
            if flags.do_guesstoa
                toa = guesstoa(toa_diff,toa, kv.avgtoa);
            end
            

        case 'CenIACC2e'
            fprintf('CenIACC2e mode \n')
            x = transpose(1:(Ns*2-1));
            for ii=1:pos              
                e_sir1 = envelope(squeeze(f_ir(ii,1,:)));
                e_sir2 = envelope(squeeze(f_ir(ii,2,:)));           
                cc = xcorr(e_sir1,e_sir2).^2;
                toa_diff(ii) = centroid(x,abs(cc))-Ns;    
            end
            if flags.do_guesstoa
                toa = guesstoa(toa_diff,toa, kv.avgtoa);
            end
            

        case 'PhminXcor'
            fprintf('PhminXcor mode \n')
            ir_min=ARI_MinimalPhase(Obj);
            for ii=1:pos
                for jj=1:ear                    
                    cc = xcorr(squeeze(IR(ii,jj,:)),squeeze(ir_min(ii,jj,:)));
                    [~,toa(ii,jj)] = max(abs(cc));
                end
                toa_diff(ii) = toa(ii,1) - toa(ii,2);
            end
            
            
    % ---------------------- Groupdelay -----------------------------------
        case 'IRGD'
            fprintf('IRGD mode \n')
            for ii = 1:pos
                for jj = 1:ear
                    f_sir = squeeze( f_ir(ii,jj,:) );
                    [gd,w] = grpdelay(transpose(double(f_sir)),1,Ns,fs);
                    toa(ii,jj)=mean(gd(find(w>kv.lower_cutfreq): ...
                                        find(w>kv.upper_cutfreq)));
                end
                toa_diff(ii) = toa(ii,1) - toa(ii,2);
            end
    end
    toa_diff = toa_diff/fs;
    toa = toa/fs;
end



% -------------------------------------------------------------------------
% ---------------------- Functions ----------------------------------------
% -------------------------------------------------------------------------

% ---------------------- Centroid -----------------------------------------
function idx_cent = centroid(x,y)
    idx_cent = sum(x.*y)/sum(y);
end

% ---------------------- Butter low-pass ----------------------------------
function cut_off_freq_norm = norm_cut_off_freq(cut_off_freq,IR,fs)
    n = size(IR,1);
    n2=floor(n/2)+1;
    xr=(0:n2-1)*2/n;
    freq = xr*fs/2;
    cut_off_freq_norm = cut_off_freq/max(freq);
end

function x_lp = butter_lp(x, poly, cut_off_freq_norm )
    [a,b] = butter(poly,cut_off_freq_norm);
    x_lp = filter(a,b,x);
end

% ---------------------- guess toa ----------------------------------------
function toa = guesstoa(toa_diff,toa, avgtoa)
    toa(:,1) = toa(:,1) + avgtoa + toa_diff/2;
    toa(:,2) = toa(:,2) + avgtoa - toa_diff/2;
end
    
% ---------------------- Create minimal phase -----------------------------
% as used in ziegelwanger2014
function hMmin=ARI_MinimalPhase(Obj)
    hM=Obj.Data.IR;
    hMmin=hM;

    for jj=1:Obj.API.R
        for ii=1:Obj.API.M
            h=squeeze(hM(ii,jj,:));
            
            amp1=abs(fft(h));
            amp2=amp1;
            
            an2u=-imag(hilbert(log(amp1)));
            an2u=an2u(1:floor(length(h)/2)+1);
            
            an3u=[an2u; -flipud(an2u(2:end+mod(length(h),2)-1))];
            an3=an3u-round(an3u/2/pi)*2*pi;
            
            amp2=amp2(1:floor(length(h)/2)+1);
            amp3=[amp2; flipud(amp2(2:end+mod(length(h),2)-1))];
            
            h2=real(ifft(amp3.*exp(1i*an3)));
            hMmin(ii,jj,:)=h2(1:Obj.API.N);
        end
    end
end