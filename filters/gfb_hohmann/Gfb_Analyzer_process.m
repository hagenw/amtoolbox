function [output, analyzer] = Gfb_Analyzer_process(analyzer, input)
% [output, analyzer] = Gfb_Analyzer_process(analyzer, input)
%
% The analyzer processes the input data.
%
% PARAMETERS
% analyzer A Gfb_Analyzer struct created from Gfb_Analyzer_new. The
%          analyzer will be returned (with updated filter states) as
%          the second return parameter
% input   Either a row vector containing the input signal to process, or
%         a matrix containing different input signals for the different
%         bands.  Different rows correspond to different filter bands,
%         while different colums correspond to different times. 
% output  A matrix containing the analyzer's output signals.  The
%         rows correspond to the filter bands
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002, Sep 2003, Nov 2006, Jan 2007

% filename : Gfb_Analyzer_process.m


if (analyzer.fast)
  % use matlab extension for fast computation.
  [output, analyzer] = Gfb_Analyzer_fprocess(analyzer, input);
else
  number_of_bands = length(analyzer.center_frequencies_hz);
  output = zeros(number_of_bands, length(input));
  for band = [1:number_of_bands]
    [output(band,:), analyzer.filters(band)] = ...
        Gfb_Filter_process(analyzer.filters(band), ...
                           input );
  end
end

