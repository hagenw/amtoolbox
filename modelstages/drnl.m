function outsig = drnl(insig,fc,fs)
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

n_lin_gt = 2; % number of cascaded gammatone filters 
n_lin_lp = 4; % no. of cascaded LP filters, orig 4

n_nlin_gt_before = 2; % number of cascaded gammatone filters 
n_nlin_gt_after = 2; % number of cascaded gammatone filters 
n_nlin_lp = 1; % no. of cascaded LP filters in nlin path

nlin_c = 10^(-.60206); % c, compression coeff

nlin_d = 1;

% ---------------- main loop over center frequencies

% Code will fail for a row vector, FIXME
siglen = size(insig,1);
nsigs  = size(insig,2);
nfc    = length(fc);

outsig=zeros(siglen,nfc,nsigs,nsigs);

for ii=1:nfc

  % -------- Setup channel dependant definitions -----------------
  fc_lin = 10^(-0.06762+1.01679*log10(fc(ii))); % Hz.
  bw_lin = 10^(.03728+.75*log10(fc(ii))); % Hz
  lp_lin_cutoff = 10^(-0.06762+1.01*log10(fc(ii)));
  
  lin_gain = 10^(4.20405 -.47909*log10(fc(ii)));
    
  fc_nlin_before = 10^(-0.05252+1.01650*log10(fc(ii))); % Hz
  fc_nlin_after  = 10^(-0.05252+1.01650*log10(fc(ii))); % Hz  
    
  bw_nlin_before = 10^(-0.03193+.7*log10(fc(ii))); % Hz  
  bw_nlin_after  = 10^(-0.03193+.7*log10(fc(ii))); % Hz  
  
  lp_nlin_cutoff = 10^(-0.05252+1.01*log10(fc(ii))); 
  
  if fc(ii)<=1000
    % a, the 1500 assumption is no good for compressionat low freq filters
    nlin_a = 10^(1.40298+.81916*log10(fc(ii))); 
    
    % b [(m/s)^(1-c)]
    nlin_b = 10^(1.61912-.81867*log10(fc(ii))); 
  else
    % a, the 1500 assumption is no good for compressionat low freq filters
    nlin_a = 10^(1.40298+.81916*log10(1500));
    
    % b [(m/s)^(1-c)]
    nlin_b = 10^(1.61912-.81867*log10(1500)); 
  end
  
  
  [GTlin_b,GTlin_a] = coefGtDRNL(fc_lin,bw_lin,n_lin_gt,fs);
  [LPlin_b,LPlin_a] = coefLPDRNL(lp_lin_cutoff,fs);

  [GTnlin_b_before,GTnlin_a_before] = coefGtDRNL(fc_nlin_before,bw_nlin_before,n_nlin_gt_before,fs);
  [GTnlin_b_after, GTnlin_a_after]  = coefGtDRNL(fc_nlin_after, bw_nlin_after,n_nlin_gt_after,fs);
  [LPnlin_b,LPnlin_a] = coefLPDRNL(lp_nlin_cutoff,fs);

  % -------------- linear part --------------------------------

  % Apply linear gain
  y_lin = insig.*lin_gain; 
  
  % Now filtering.
  % Instead of actually perform multiply filtering, just convolve the
  % coefficients.      
  %[blong,along]=convolveba(GTlin_b,GTlin_a,n_lin_gt);
  y_lin = real(filter(GTlin_b,GTlin_a,y_lin));

  % Same story: convolve the coefficients.
  [blong,along]=convolveba(LPlin_b,LPlin_a,n_lin_lp);
  y_lin = filter(blong,along,y_lin);
  
  % -------------- Non-linear part ------------------------------
      
  % GT filtering before
  y_nlin = real(filter(GTnlin_b_before,GTnlin_a_before,insig));
  
  % Broken stick nonlinearity
  if nlin_d~=1
    % Just to save some flops, make this optional.
    y_decide = [nlin_a*abs(y_nlin).^nlin_d; ...
                nlin_b*(abs(y_nlin)).^nlin_c];
  else
    y_decide = [nlin_a*abs(y_nlin); ...
                nlin_b*(abs(y_nlin)).^nlin_c];    
  end;
  y_nlin = sign(y_nlin).* min(y_decide);
  
  % GT filtering after
  y_nlin = real(filter(GTnlin_b_after,GTnlin_a_after,insig));
  
  % then LP filtering
  [blong,along]=convolveba(LPnlin_b,LPnlin_a,n_nlin_lp);
  y_nlin = real(filter(blong,along,y_nlin));
  
  outsig(:,ii,:) = reshape(y_lin + y_nlin,siglen,1,nsigs);    
  
end;
  
 
