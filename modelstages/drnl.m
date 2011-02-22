function outsig = drnl(insig,fs,varargin)
%DRNL  Dual Resonance Nonlinear Filterbank
%   Usage: outsig = drnl(insig,fc,fs);
%
%   DRNL(insig,fc,fs) computes the Dual Resonance Non-Linear filterbank of
%   the input signal insig sampled at fs Hz with channels specified by the
%   center frequencies in fc. The DRNL is described in the paper
%   Lopez-Poveda and Meddis (2001). The DRNL models the basilar membrane
%   non-linearity.
%
%   The DRNL takes a lot of parameters which vary over frequency. Such a
%   parameter is described by a 1x2 vector [b a] and indicates that the
%   value of the parameter at frequency fc can be calculated by
%
%C     10^(b+alog10(fc));
%
%   The parameters are:
%
%-     'flow',flow - Set the lowest frequency in the filterbank to
%                    flow. Default value is 80 Hz.
%
%-     'fhigh',fhigh - Set the highest frequency in the filterbank to
%                    fhigh. Default value is 8000 Hz.
%
%-     'basef',basef - Ensure that the frequency basef is a center frequency
%                    in the filterbank. The default value of [] means
%                    no default.
%
%-     'middleear'  - Perform middleear filtering before the actual DRNL
%                     is applied using the middleear filter specified in
%                     Lopez-Poveda and Meddis (2001). This is the
%                     default.
%
%-     'nomiddleear' - No middle-ear filtering. Be carefull with this setting,
%                     as another scaling must then be perform to convert the
%                     input to stapes movement.
%
%-     'lin_ngt',n  - Number of cascaded gammatone filter in the linear
%                    part, default value is 2.
%
%-     'lin_nlp',n - Number of cascaded lowpass filters in the linear
%                    part, default value is 4
% 
%-     'lin_gain',g - Gain in the linear part, default value is [4.20405 ...
%                    .47909].
%
%-     'lin_fc',fc - Center frequencies of the gammatone filters in the
%                    linear part. Default value is [-0.06762 1.01679].
%
%-     'lin_bw',bw - Bandwidth of the gammatone filters in the linear
%                    part. Default value is [.03728  .78563]
%
%-     'lin_lp_cutoff',c - Cutoff frequency of the lowpass filters in the
%                    linear part. Default value is [-0.06762 1.01679 ]
%
%-     'nlin_ngt_before',n - Number of cascaded gammatone filters in the
%                    non-linear part before the broken stick
%                    non-linearity. Default value is 2.
%
%-     'nlin_ngt_after',n -  Number of cascaded gammatone filters in the
%                    non-linear part after the broken stick
%                    non-linearity. Default value is 2.
%
%-     'nlin_nlp',n - Number of cascaded lowpass filters in the
%                    non-linear part. Default value is 1.
%
%-     'nlin_fc_before',fc - Center frequencies of the gammatone filters in the
%                    non-linear part before the broken stick
%                    non-linearity. Default value is [-0.05252 1.01650].
%
%-     'nlin_fc_after',fc - Center frequencies of the gammatone filters in the
%                    non-linear part after the broken stick
%                    non-linearity. Default value is [-0.05252 1.01650].
%
%-     'nlin_bw_before',bw - Bandwidth of the gammatone filters in the
%                    non-linear part before the broken stick
%                    non-linearity. Default value is [-0.03193 .77426 ].
%
%-     'nlin_bw_after',w - Bandwidth of the gammatone filters in the
%                    non-linear part before the broken stick
%                    non-linearity. Default value is [-0.03193 .77426 ].
%
%-     'nlin_lp_cutoff',c - Cutoff frequency of the lowpass filters in the
%                    non-linear part. Default value is [-0.05252 1.01 ].
%
%-     'nlin_a',a - 'a' coefficient for the broken-stick non-linearity. Default
%                   value is [1.40298 .81916 ].
%
%-     'nlin_b',b = 'b' coefficient for the broken-stick non-linearity. Default
%                   value is [1.61912 -.81867
%
%-     'nlin_c',c = 'c' coefficient for the broken-stick non-linearity. Default
%                   value is [-.60206 0].
%
%-     'nlin_d',d = 'd' coefficient for the broken-stick non-linearity. Default
%                    value is 1.
%
%   See also: middlerearfilter, jepsen2008preproc
% 
%R   lopezpoveda2001hnc jepsen2008cmh

% AUTHOR: Morten Løve Jepsen
  
% Bugfix by Marton Marschall 9/2008
% Cleanup by Peter L. Soendergaard.
%
% In comparison to the oroginal code, the gammatone filter are computed
% by convolving the filter coefficients, instead of performing multiple
% runs of 'filter'. The lowpass filtering is still performed by mutliple
% runs through 'filter', as the linear-part lowpass filter turned out to
% be unstable when the coefficients was convolved.
  
% Import the parameters from the arg_drnl.m function.
definput.import={'drnl'};

[flags,kv]=ltfatarghelper({'flow','fhigh'},definput,varargin);

% find the center frequencies used in the filterbank, 1 ERB spacing
fc = erbspacebw(kv.flow, kv.fhigh, 1, kv.basef);

%% Apply the middle-ear filter
if flags.do_middleear
  
  me_fir = middleearfilter(fs);
  insig = filter(me_fir,1,insig);

end;

%% ---------------- main loop over center frequencies

% Code will fail for a row vector, FIXME
siglen = size(insig,1);
nsigs  = size(insig,2);
nfc    = length(fc);

outsig=zeros(siglen,nfc,nsigs,nsigs);

% Handle the compression limiting in the broken-stick non-linearity
if ~isempty(kv.compresslimit)
  fclimit=min(fc,kv.compresslimit);
else
  fclimit=fc;
end;

% Sanity checking, some center frequencies may go above the Nyquest
% frequency
% Happens for lin_lp_cutoff

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
  
  % Expand a and b using the (possibly) limited fc
  nlin_a = polfun(kv.nlin_a,fclimit(ii));

  % b [(m/s)^(1-c)]
  nlin_b = polfun(kv.nlin_b,fclimit(ii));
      
  nlin_c = polfun(kv.nlin_c,fc(ii));
  
  [GTlin_b,GTlin_a] = coefGtDRNL(lin_fc,lin_bw,kv.lin_ngt,fs);
  
  % Compute coefficients for the linear stage lowpass, use 2nd order
  % Butterworth.
  [LPlin_b,LPlin_a] = butter(2,lin_lp_cutoff/(fs/2));

  [GTnlin_b_before,GTnlin_a_before] = coefGtDRNL(nlin_fc_before,nlin_bw_before,...
                                                 kv.nlin_ngt_before,fs);
  [GTnlin_b_after, GTnlin_a_after]  = coefGtDRNL(nlin_fc_after, nlin_bw_after,...
                                                   kv.nlin_ngt_after,fs);    

  % Compute coefficients for the non-linear stage lowpass, use 2nd order
  % Butterworth.
  [LPnlin_b,LPnlin_a] = butter(2,nlin_lp_cutoff/(fs/2));

  % -------------- linear part --------------------------------

  % Apply linear gain
  y_lin = insig.*lin_gain; 
  
  % Gammatone filtering
  y_lin = filter(GTlin_b,GTlin_a,y_lin);    

  % Multiple LP filtering
  for jj=1:kv.lin_nlp
    y_lin = filter(LPlin_b,LPlin_a,y_lin);
  end;

  % -------------- Non-linear part ------------------------------
      
  % GT filtering before
  y_nlin = filter(GTnlin_b_before,GTnlin_a_before,insig);
  
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
  y_nlin = filter(GTnlin_b_after,GTnlin_a_after,y_nlin);

  % then LP filtering
  for jj=1:kv.nlin_nlp
    y_nlin = filter(LPnlin_b,LPnlin_a,y_nlin);
  end;
  
  outsig(:,ii,:) = reshape(y_lin + y_nlin,siglen,1,nsigs);    
  
end;
  
 
function outpar=polfun(par,fc)
  %outpar=10^(par(1)+par(2)*log10(fc));
  outpar=10^(par(1))*fc^par(2);