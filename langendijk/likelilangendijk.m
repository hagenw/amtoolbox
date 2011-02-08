function [ la,le,ci ] = likelilangendijk( p,rang,tang,target,response )
% LIKELILANGENDIJK Likelihood estimation for evaluating model performance 
% according to Langendijk et al. (2002)
% Usage:           [la,le,ci] = likelilangendijk(p,rang,tang,target,response)
% Input arguments:
%     p:           pdf matrix
%     rang:        polar angles of possible response angles
%     tang:        polar angles of possible target angles
%     target:      target polar angles of localization test
%     response:    response polar angles of localization test
% Output argument:
%     la:          actual likelihood
%     le:          expected likelihood
%     ci:          99% confidence interval for expected likelihood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR : Robert Baumgartner, OEAW Acoustical Research Institute
% latest update: 2010-08-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nt=length(target);

% pa represents pdf values of actual responses
pa=interp2([-90;rang(:);270],[-90;tang(:);270], [zeros(1,size(p,2)+2);[zeros(size(p,1),1),p,zeros(size(p,1),1)];zeros(1,size(p,2)+2)] ,target,response);
la=-2*sum(log(pa))*55/nt; % actual likelihood

% random generator
lex=zeros(100,1);
for ind=1:100
    pe=zeros(size(target));
    for ind1=1:nt 
        post=find(tang>=target(ind1),1); % target position
%         post=randi(size(p,2),1);
        posr = discreteinvrnd(p(:,post),1,1);
        pe(ind1)=p(posr,post);
    end
    lex(ind)=-2*sum(log(pe))*55/nt; 
end

le=mean(lex); % expected likelihood
err=2.58*std(lex);
ci=[le-err le+err]; % confidence interval

end