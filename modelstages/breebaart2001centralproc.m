function decision = breebaart2001centralproc(EI_map,monol,monor)

persistent template;
persistent templatesq;
persistent signaltemplate;
persistent maskernumber;
persistent signalnumber;
persistent template_ml;
persistent templatesq_ml;
persistent signaltemplate_ml;
persistent template_mr;
persistent templatesq_mr;
persistent signaltemplate_mr;

monofactor = 0.0003;
monol = cellfun(@(x) x*monofactor,monol,'un',0);
monor = cellfun(@(x) x*monofactor,monor,'un',0);

intnum = length(EI_map);
tempsize = size(EI_map{1});

intnoisevar = 1;
U_b = zeros(1,intnum);
U_ml = zeros(1,intnum);
U_mr = zeros(1,intnum);


    
if isempty(maskernumber)
    maskernumber = intnum-1;
    template = zeros(tempsize);
    templatesq = zeros(tempsize);
    template_ml = zeros(tempsize);
    templatesq_ml = zeros(tempsize);
    template_mr= zeros(tempsize);
    templatesq_mr = zeros(tempsize);
else
    maskernumber = maskernumber+intnum-1;
end

if isempty(signalnumber)
    signalnumber = 1;
    signaltemplate = zeros(tempsize);
    signaltemplate_ml = zeros(tempsize);
    signaltemplate_mr = zeros(tempsize);
else
    signalnumber = signalnumber + 1;
end

noisevar = templatesq - template.^2;
meandiff = signaltemplate - template;
weight = meandiff./(noisevar + intnoisevar);
Nuvar = intnoisevar*sum(sum(weight.^2));

noisevar_ml = templatesq_ml - template_ml.^2;
meandiff_ml = signaltemplate_ml - template_ml;
weight_ml = meandiff_ml./(noisevar_ml + intnoisevar);
Nuvar_ml = intnoisevar*sum(sum(weight_ml.^2));

noisevar_mr = templatesq_mr - template_mr.^2;
meandiff_mr = signaltemplate_mr - template_mr;
weight_mr = meandiff_mr./(noisevar_mr + intnoisevar);
Nuvar_mr = intnoisevar*sum(sum(weight_mr.^2));

for intcount = 1:intnum
    noise = randn;
    U_b(intcount) = sum(sum(weight.*(EI_map{intcount}-template))) + noise*sqrt(Nuvar);
    U_ml(intcount) = sum(sum(weight_ml.*(monol{intcount}-template_ml))) + noise*sqrt(Nuvar_ml);
    U_mr(intcount) = sum(sum(weight_mr.*(monor{intcount}-template_mr))) + noise*sqrt(Nuvar_mr);
end

if signalnumber > 1 && any(U_b) == 0 % no binaural decision and not first call
    U = mean([U_ml;U_mr]);
else
    U = mean([U_b;U_ml;U_mr]);
end

% If centralproc is called the first time, all the templates are empty
% In this case the results of U will be zeros. Therefore the response
% of max(U) will always be 1 = first occurence of 0.
[~,response] = max(U);

% binaural test
[~,response_b] = max(U_b);
if response_b ~= response && any(U_b) == 1
    fprintf('binaural decision is not the same: signalnumber %i, response = %i, binaural response = %i \n',signalnumber,response,response_b);
end

% update of templates
signaltemplate = ((signaltemplate*(signalnumber-1))+(EI_map{1}))./signalnumber;
signaltemplate_ml = ((signaltemplate_ml*(signalnumber-1))+(monol{1}))./signalnumber;
signaltemplate_mr = ((signaltemplate_mr*(signalnumber-1))+(monor{1}))./signalnumber;

adtemplate = zeros(tempsize);
adtemplatesq = zeros(tempsize);
adtemplate_ml = zeros(tempsize);
adtemplatesq_ml = zeros(tempsize);
adtemplate_mr = zeros(tempsize);
adtemplatesq_mr = zeros(tempsize);

for updatecounter = 2:intnum
    adtemplate = adtemplate + EI_map{updatecounter};
    adtemplatesq = adtemplatesq + EI_map{updatecounter}.^2;
    adtemplate_ml = adtemplate_ml + monol{updatecounter};
    adtemplatesq_ml = adtemplatesq_ml + monol{updatecounter}.^2;    
    adtemplate_mr = adtemplate_mr + monor{updatecounter};
    adtemplatesq_mr = adtemplatesq_mr + monor{updatecounter}.^2;
end

template = ((template*(maskernumber-(intnum-1)))+(adtemplate))./maskernumber;
templatesq = ((templatesq*(maskernumber-(intnum-1)))+(adtemplatesq).^2)./maskernumber; 
template_ml = ((template_ml*(maskernumber-(intnum-1)))+(adtemplate_ml))./maskernumber;
templatesq_ml = ((templatesq_ml*(maskernumber-(intnum-1)))+(adtemplatesq_ml).^2)./maskernumber; 
template_mr = ((template_mr*(maskernumber-(intnum-1)))+(adtemplate_mr))./maskernumber;
templatesq_mr = ((templatesq_mr*(maskernumber-(intnum-1)))+(adtemplatesq_mr).^2)./maskernumber; 


decision = response;