%DEMO_JOERGENSEN2013 Demo for the multi-resolution speech-based envelope spectrum model 
%
%   `demo_joergensen2013` computes the signal-to-noise envelope-power ratio
%   (SNRenv) for a speech sentence in noise at an SNR of -3 dB.
%   The SNRenv can be then used to predict speech intelligibility. 
%
%   Pcorrect indicates the probability of correctly understanding the
%   sentence based on results from Joergensen et al. (2013). 
%
%   See also: joergensen2013 joergensen2011
%
%   References: joergensen2013 joergensen2011predicting

% AUTHOR : Piotr Majdak

  % Load speech (from the CLUE corpus)
x=amtload('joergensen2011','Danish_CLUE_10sentence_samples_22kHz.mat');

speech  = x.sentenceArray{1};
sentenceFileLevel = -26.00; % The RMS level of all CLUE sentence files corresponds to ...
SPL = 65; % ... this sound pressure level
SNR = -3;
fs = 22050;
speech = speech*10^((SPL-sentenceFileLevel)/20);
N = length(speech);

  % Load speech-shaped noise from the CLUE corpus
noise = amtload('joergensen2011','SSN_CLUE_22kHz.wav');

% pick a random segment from the noise file
Nsegments = floor(length(noise)/N);
startIdx = randi(Nsegments-2 ,1)*N;
noise = noise(startIdx:startIdx+N -1)';
noise = noise./rms(noise)*10^((SPL-SNR)/20);
if size(noise) ~= size(speech), noise = noise'; end

  % Create a mixed signal
test = noise + speech;

  % Run the model without any priors
tmp = joergensen2013(test,noise,fs);
SNRenvs_noIOparameters = tmp.SNRenv

  % Run the model with parameters for the CLUE material from Joergensen et al., (2013).
IOparameters = [0.61 0.5 8000 0.6]; 
tmp = joergensen2013(test,noise,fs, IOparameters);
SNRenvs_withIOparameters = tmp.SNRenv
Pcorrect = tmp.P_correct
