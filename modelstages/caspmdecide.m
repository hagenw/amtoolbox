function [detect,prob] = caspmdecide(mu,in_var,rule,numint)

% Y = 1, signal is detected

if in_var <= 0
    error('CASP:mdecide', 'in_var must be > 0');
end

if rule(1) > 1
    error('CASP:mdecide','Only x-down, 1-up procedure have been implemented')
end


switch numint
    case 2
        prob = 1 - (erfc((((mu/in_var)*0.707) - 0)     * sqrt(2)/2) / 2);
    case 3
        prob = 1 - (erfc((((mu/in_var)*0.765) - 0.423) * sqrt(2)/2) / 2);
    case 4
        prob = 1 - (erfc((((mu/in_var)*0.810) - 0.668) * sqrt(2)/2) / 2);
    otherwise
        error('CASP:mdecide', 'Only 2-, 3- and 4-AFC procedures are implemented');
end;

if rule(1) == 1 && max(prob) > (1 / (2 .^ (1/rule(2))))
    detect = 1; % signal heard
else
    detect = 0; % no signal heard
end;
