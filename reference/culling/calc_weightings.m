function [weightings] = calc_weightings(fc)

weightings = zeros(1,length(fc));
bands = [0, 100, 200, 300, 400, 4400, 5300, 6400, 7700, 9500];
weights = [0, 0.0103, 0.0261, 0.0419, 0.0577, 0.0460, 0.0343, 0.0226, 0.0110, 0];
for n = 1:length(fc)
    if fc(n) >= 9500
        weightings(n) = 0;
    else
        i = find(bands > fc(n));
        weightings(n) = weights((i(1) - 1));
    end
end
weightings = weightings .* (1/sum(weightings));
