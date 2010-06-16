%EXAMP_AUDSPECGRAM  Some auditory spectrograms.
%
%  Display some auditory spectrograms of the word 'greasy' sampled at 16
%  kHz.
%
%  FIGURE 1 Classic version
%
%    This figure shows a classic spectrogram using auditory filters of the
%    word 'greasy'
%
%  FIGURE 2 Full version
%
%    This figure shows the auditory spectrogram of the word 'greasy'. 
%
%  See also: audspecgram

disp('Type "help demo_audspecgram" to see a description of how this demo works.');

figure(1);
audspecgram(greasy,16000,'classic');

figure(2);
audspecgram(greasy,16000);





