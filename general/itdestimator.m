function [toa_diff,toa] = itdestimator(Obj,varargin)
%   Usage: itd = itd_estimator(data,[mode],[threshlvl],...
%                   [lowpass],[butterpoly],[upper_cutfreq]) 
%
%   Input parameters:
% 
%       data:       SOFA object or IR matrix with dimensions: 
%                     emitter x receiver x time
% 
%       fs:         sampling rate, used only if *data* provided as matrix
% 
%       mode:       (optional) Select one estimation methods
%                   ('Threshold' (default),'Cen_e2','MaxIACCr', 'MaxIACCe',
%                   'CenIACCr','CenIACCe', 'CenIACC2e', 'PhminXcor','IRGD')
%             
%                   For details concerning estimation methods see:
%                   'http://asa.scitation.org/doi/10.1121/1.4996457'
% 
%       lowpass:    (optional) Decide if lowpass shall be applied. 
%                   'lp' for lowpass (default), 'bb' for broadband
%
%       threshlvl:  (optional) Set threshold level for 'Threshold' mode in dB.        
%                   Default is -10 dB. 
%
%       butterpoly: (optional) Select the order of the polynom
%                   applied in the butterworth filter. ( 2 =< i =< 6 )
%                   Default is 6. 
% 
%       upper_cutfreq: (optional) Set frequency of lowpass cutoff in Hz.
%                      Default is 3000 Hz. 
% 
%       lower_cutfreq: (optional) Set frequency of highpass cutoff in Hz, 
%                      only used in IRGD mode. Default is 1000 Hz.    
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
%   a lowpass with a cut off frequency of 3 kHz.
%   
%   The output array is structured as the SOFA Data.IR
%   If you would like to select for example only data on the horizontal
%   plane you could:
%
%   plane_idx = find( Obj.SourcePosition(:,2) == 0 );
%   plane_angle = Obj.SourcePosition(plane_idx,1);

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

% AUTHOR: Laurin Steidle
% data matrix option: Robert Baumgartner


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
    definput.keyvals.fs = [];

    [flags,kv]=ltfatarghelper({},definput,varargin);

    
    % ---------------------- renaming input parameter ---------------------
    
    if isstruct(Obj)
      pos = Obj.API.M;
      ear = Obj.API.R;
      Ns  = Obj.API.N;
      IR = Obj.Data.IR;
      fs  = Obj.Data.SamplingRate;
    else
      pos = size(Obj,1);
      ear = size(Obj,2);
      Ns  = size(Obj,3);
      IR = Obj;
      if isempty(kv.fs)
        error('RB: No sampling rate (fs) provided.')
      end
      fs  = kv.fs;
    end

    % ---------------------- initialising variables -----------------------

    toa = zeros(pos,ear);
    toa_diff = zeros(pos,1);
    
    amt_disp('itdestimator \n','volatile')

    % ---------------------- Applying low-pass ----------------------------

    if flags.do_lp
        amt_disp('Applying Butterworth low pass \n','volatile')
        amt_disp(strcat('Polynomial order of Butterworth filter: '...
            ,num2str(kv.butterpoly),'\n'),'volatile')
        amt_disp(strcat('Cut-off frequency is: '...
            ,num2str(kv.upper_cutfreq),' Hz\n\n'),'volatile')
        cut_off_freq_norm = kv.upper_cutfreq/(fs/2);
        [lp_a,lp_b] = butter(kv.butterpoly,cut_off_freq_norm);
        f_ir = zeros(pos,ear,Ns);
        for ii=1:pos
            for jj=1:ear  
                sir = squeeze( IR(ii,jj,:) );
                f_sir = filter(lp_a,lp_b,sir);
                f_ir(ii,jj,:) = f_sir;
            end
        end
    
    else
        amt_disp('No low pass filter is applied \n\n','volatile')
        f_ir = IR;
    end

    % ---------------------- estimating itd -------------------------------
    % ---------------------------------------------------------------------

    % ---------------------- Threshold ------------------------------------
    switch(flags.mode)
        case 'Threshold'
            amt_disp('Threshold mode \n','volatile')
            amt_disp(strcat('Threshold level is: ',...
                num2str(kv.threshlvl),'dB \n\n'),'volatile')
            
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
            amt_disp('Cen-e2 mode \n','volatile')
            for ii=1:pos
                for jj = 1:ear
                    e_sir_sq = envelope(squeeze(f_ir(ii,jj,:))).^2;
                    toa(ii,jj) = centroid(transpose(1:Ns),e_sir_sq);
                end
                toa_diff(ii) = toa(ii,1) - toa(ii,2);     
            end        
            

        case 'MaxIACCr'
            amt_disp('MaxIACCr mode \n','volatile')
            for ii=1:pos                
                cc = xcorr(squeeze(f_ir(ii,1,:)),squeeze(f_ir(ii,2,:)));
                [~,idx_lag] = max(abs(cc));
                toa_diff(ii) = idx_lag - Ns;                
            end
            if flags.do_guesstoa
                toa = guesstoa(toa_diff,toa, kv.avgtoa);
            end
            
        case 'MaxIACCe'
            amt_disp('MaxIACCe mode \n','volatile')
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
            amt_disp('CenIACCr mode \n','volatile')
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
            amt_disp('CenIACCe mode \n','volatile')
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
            amt_disp('CenIACC2e mode \n','volatile')
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
            amt_disp('PhminXcor mode \n','volatile')
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
            amt_disp('IRGD mode \n','volatile','volatile')
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