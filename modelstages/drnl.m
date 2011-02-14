function outsig = drnl(insig,fc,fs,varargin)
%DRNL  Dual Resonance Nonlinear Filterbank
%   Usage: outsig = drnl(insig,fc,fs);
%
%   DRNL(insig,fc,fs) computes the Dual Resonance Non-Linear filterbank of
%   the input signal insig sampled at fs Hz with channels specified by the
%   center frequencies in fc. The DRNL is described in the paper
%   Lopez-Poveda and Meddis (2001). The DRNL models the basilar membrane
%   non-linearity.
%
%   The input to DRNL must be measured in stapes movement. The function
%   MIDDLEEARFILTER will do this, so this function must be called before
%   invoking DRNL.
%
%   See also: middlerearfilter, jepsen2008preproc
% 
%R   lopezpoveda2001hnc jepsen2008cmh

% AUTHOR: Morten Løve Jepsen
  
% Bugfix by Marton Marschall 9/2008
% Cleanup by Peter L. Soendergaard.
  
%DRNL for normal hearing, Morten 2007


definput.keyvals.lin_ngt = 2; % number of cascaded gammatone filters 
definput.keyvals.lin_nlp = 4; % no. of cascaded LP filters, orig 4
definput.keyvals.lin_fc = [-0.06762 1.01679]; %FOUND in 2001
definput.keyvals.lin_bw = [ .03728   .75     ]; %2001: [  .03728 .78563];
definput.keyvals.lin_gain = [4.20405 -.47909]; %FOUND in 2001
definput.keyvals.lin_lp_cutoff = [-0.06762 1.01 ]; %2001: [.06762 1.01679 ]

definput.keyvals.nlin_ngt_before = 2; % number of cascaded gammatone filters 
definput.keyvals.nlin_ngt_after = 2; % number of cascaded gammatone filters 
definput.keyvals.nlin_nlp = 1; % no. of cascaded LP filters in nlin path
definput.keyvals.nlin_fc_before = [-0.05252 1.01650]; %FOUND in 2001
definput.keyvals.nlin_fc_after  = [-0.05252 1.01650 ]; %FOUND in 2001
definput.keyvals.nlin_bw_before = [-0.03193 .7 ]; %2001: [-0.03193 .77426 ];
definput.keyvals.nlin_bw_after  = [-0.03193 .7 ]; %2001: [-0.03193 .77426 ];
definput.keyvals.nlin_lp_cutoff = [-0.05252 1.01 ];

definput.keyvals.nlin_a = [1.40298 .81916 ];  %FOUND in 2001
definput.keyvals.nlin_b = [1.61912 -.81867 ]; %FOUND in 2001
definput.keyvals.nlin_c = [-.60206 0]; % FOUND in 2001 c, compression coeff
definput.keyvals.nlin_d = 1;


[flags,kv]=ltfatarghelper({},definput,varargin);


do_testing=1;

% ---------------- main loop over center frequencies

% Code will fail for a row vector, FIXME
siglen = size(insig,1);
nsigs  = size(insig,2);
nfc    = length(fc);

outsig=zeros(siglen,nfc,nsigs,nsigs);

for ii=1:nfc

  % -------- Setup channel dependant definitions -----------------

  lin_fc        = polfun(kv.lin_fc,fc(ii));
  lin_bw        = polfun(kv.lin_bw,fc(ii));
  lin_lp_cutoff = polfun(kv.lin_lp_cutoff,fc(ii));
  lin_gain      = polfun(kv.lin_gain,fc(ii));
  
  nlin_fc_before = polfun(kv.nlin_fc_before,fc(ii));
  nlin_fc_after  = polfun(kv.nlin_fc_after,fc(ii));
  
  nlin_bw_before = polfun(kv.nlin_bw_before,fc(ii));
  nlin_bw_after  = polfun(kv.nlin_bw_after,fc(ii));
  
  nlin_lp_cutoff = polfun(kv.nlin_lp_cutoff,fc(ii));
  
  % a, the 1500 assumption is no good for compressionat low freq filters
  nlin_a = polfun(kv.nlin_a,min(fc(ii),1500));

  % b [(m/s)^(1-c)]
  nlin_b = polfun(kv.nlin_b,min(fc(ii),1500));
      
  nlin_c = polfun(kv.nlin_c,fc(ii));
  if do_testing
    [GTlin_b,GTlin_a] = coefGtDRNL(lin_fc,lin_bw,1,fs);
  else
    [GTlin_b,GTlin_a] = coefGtDRNL(lin_fc,lin_bw,kv.lin_ngt,fs);
  end;
  
  % Compute coefficients for the linear stage lowpass, use 2nd order
  % Butterworth.
  if do_testing
    [LPlin_b,LPlin_a] = coefLPDRNL(lin_lp_cutoff,fs);
  else
    [LPlin_b,LPlin_a] = butter(2,lin_lp_cutoff/(fs/2));
  end;

  if do_testing
    [GTnlin_b_before,GTnlin_a_before] = coefGtDRNL(nlin_fc_before,nlin_bw_before,...
                                                   1,fs);
    [GTnlin_b_after, GTnlin_a_after]  = coefGtDRNL(nlin_fc_after, nlin_bw_after,...
                                                   1,fs);
  else
    [GTnlin_b_before,GTnlin_a_before] = coefGtDRNL(nlin_fc_before,nlin_bw_before,...
                                                   kv.nlin_ngt_before,fs);
    [GTnlin_b_after, GTnlin_a_after]  = coefGtDRNL(nlin_fc_after, nlin_bw_after,...
                                                   kv.nlin_ngt_after,fs);    
  end;

  % Compute coefficients for the non-linear stage lowpass, use 2nd order
  % Butterworth.
  if do_testing
    [LPnlin_b,LPnlin_a] = coefLPDRNL(nlin_lp_cutoff,fs);
  else
    [LPnlin_b,LPnlin_a] = butter(2,nlin_lp_cutoff/(fs/2));
  end;

  % -------------- linear part --------------------------------

  % Apply linear gain
  y_lin = insig.*lin_gain; 
  
  % Now filtering.
  % Instead of actually perform multiply filtering, just convolve the
  % coefficients.      
  if do_testing
    for jj=1:kv.lin_ngt
      y_lin = filter(GTlin_b,GTlin_a,y_lin);
    end;
    
    for jj=1:kv.lin_nlp
      y_lin = filter(LPlin_b,LPlin_a,y_lin);
    end;

  else
    
    [blong,along]=convolveba(LPlin_b,LPlin_a,kv.lin_nlp);
    y_lin = filter(GTlin_b,GTlin_a,y_lin);    
    y_lin = filter(blong,along,y_lin);

    %y_lin = real(filter(conv(blong,GTlin_b),...
    %                    conv(along,GTlin_a),y_lin));    
  end;
  
  % -------------- Non-linear part ------------------------------
      
  % GT filtering before
  if do_testing
    y_nlin=insig;
    for jj=1:kv.nlin_ngt_before
      y_nlin = filter(GTnlin_b_before,GTnlin_a_before,y_nlin);
    end;
  else
    y_nlin = filter(GTnlin_b_before,GTnlin_a_before,insig);
  end;
  
  % Broken stick nonlinearity
  if kv.nlin_d~=1
    % Just to save some flops, make this optional.
    y_decide = [nlin_a*abs(y_nlin).^kv.nlin_d; ...
                nlin_b*(abs(y_nlin)).^nlin_c];
  else
    y_decide = [nlin_a*abs(y_nlin); ...
                nlin_b*(abs(y_nlin)).^nlin_c];    
  end;
  y_nlin = sign(y_nlin).* min(y_decide);
  
  % GT filtering after
  if do_testing
    for jj=1:kv.nlin_ngt_after
      y_nlin = filter(GTnlin_b_after,GTnlin_a_after,y_nlin);
    end;
  else
    y_nlin = filter(GTnlin_b_after,GTnlin_a_after,y_nlin);
  end;
  
  % then LP filtering
  if do_testing
    for jj=1:kv.nlin_nlp
      y_nlin = filter(LPnlin_b,LPnlin_a,y_nlin);
    end;
  else
    [blong,along]=convolveba(LPnlin_b,LPnlin_a,kv.nlin_nlp);
    y_nlin = filter(blong,along,y_nlin);
  end;
  
  outsig(:,ii,:) = reshape(y_lin + y_nlin,siglen,1,nsigs);    
  
end;
  
 
function outpar=polfun(par,fc)
  %outpar=10^(par(1)+par(2)*log10(fc));
  outpar=10^(par(1))*fc^par(2);