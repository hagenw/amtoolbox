function [Pc_est SNRenv] = speechpercentcorrect(SNRenv_lin,material)
%SPEECHPERCENTCORRECT  Converts the overall SNRenv to percent correct
%   Usage: [Pc_est SNRenv ] = speechpercentcorrect(SNRenv_lin,material);
%
%   Input parameters:
%     SNRenv_lin :  matrix with the SNRenv values (not in dB),
%                   one for each SNR and processing combination.
% 
%     material   :  Specify the speech material used.
%                   Current: 'CLUE', 'DANTALE' or 'DANTALE II'
%
%   XXX Description is missing
%  
%   References: green1964effect

%  AUTHOR: Søren Jørgensen august 2010


if nargin < 2
    disp('you have to specify the spech material used: CLUE, DANTALE or DANTALEII')
end

% ---------- Determine the model parameters based on the speech material used
switch material
    case 'CLUE'
        m = 8000; 
        sigma_s = .6;
        
    case 'DANTALE'
        m = 8000 ; 
        sigma_s = .9;
     
    case 'DANTALEII'
        m = 50; 
        sigma_s = .9; 
    
end
% -------- The general conversion constants. They are the same for alle materials
    k = .82; 
    q = .5;

%--------  Only used for output:
   SNRenv = 10*log10(SNRenv_lin);

% ---------- Converting from SNRenv to d_prime  --------------
    d_prime = k*(SNRenv_lin).^q; %SNRm_lin; %

%----------- Converting from d_prime to Percent correct, Green and Birdsall (1964)----------
    Un = 1*norminv(1-(1/m)); 
    mn = Un + (.577 /Un);% F^(-1)[1/n] Basically gives the value that would be drawn from a normal destribution with probability p = 1/n.
    sig_n=  1.28255/Un; 
    Pc_est = normcdf(d_prime,mn,sqrt(sigma_s.^2+sig_n.^2))*100;  

        
