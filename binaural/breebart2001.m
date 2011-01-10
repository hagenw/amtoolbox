function decision=breebart2001(intervals,fs,tau,alpha)
%BREEBART2001  The Breebart2001 model
%
%  Input arguments:
%    Intervals  - Maxtrix of intervals. Dimensions: time, channel,
%                 interval no.
%
%
%  
  
% Compute the preprocessing of all the intervals
ir_all = breebart2001preproc(intervals,fs);

[siglen,nfreqchannels,naudiochannels,nifc] = size(ir_all);

ei_map = zeros(nifc, nfreqchannels, siglen);
for k=1:nifc
  for g=1:nfreqchannels
    ei_map(k,g,:) = eicell(squeeze(ir_all(:,g,:,k)),fs,tau,alpha);
  end
end




%%  ------- unprocess code ----------------

if 0
  %% Derive data for update or initialisation of templates
  % First iteration of the template of the target interval
  % g: number of channel
  
  TargetBin = zeros(length(s.ModelParameters.ChannelToLookAt), round(s.maskerlength*s.fs/1000));
  for g=1:length(s.ModelParameters.ChannelToLookAt),
    TargetBin(g,:) = squeeze(Bin(1,g,:))';
  end

  %% First iteration of the template of the reference interval
  % The average of the internal representation of all masker-alone
  % intervals
  
  ReferenceBin = zeros(length(s.ModelParameters.ChannelToLookAt), round(s.maskerlength*s.fs/1000));
  for g=1:length(s.ModelParameters.ChannelToLookAt),
    for j=2:s.nifc,
      ReferenceBin(g,:) = ReferenceBin(g,:) + squeeze(Bin(j,g,:))';
    end
    ReferenceBin(g,:) = ReferenceBin(g,:)./(s.nifc-1);
  end
  
  %% Test whether templates already exist
  % If no templates, create them
  % This is the first run store initial templates
  
  if ~isfield(s,'Ebar'),
    for g=1:length(s.ModelParameters.ChannelToLookAt),
      s.Ebar(g,:) = ReferenceBin(g,:);             % Masker-alone template
      s.Ebarsq(g,:) = ReferenceBin(g,:).^2 ;       % Same but squared
      s.EbarTarget(g,:) = TargetBin(g,:);      % Masker+signal template
    end
    s.gamma = 0;                           % Number of updates of template
  end
  

%% sigma2 is the variance in the internal reference representation
sigma2 = zeros(length(s.ModelParameters.ChannelToLookAt),round(s.maskerlength*s.fs/1000));
for g=1:length(s.ModelParameters.ChannelToLookAt),
    sigma2(g,:) = 1+s.Ebarsq(g,:)-(s.Ebar(g,:)).^2;
end



%% mu is difference of mean internal representation of reference and target
%% interval
mu = zeros(length(s.ModelParameters.ChannelToLookAt),round(s.maskerlength*s.fs/1000));
for g=1:length(s.ModelParameters.ChannelToLookAt),
    mu(g,:) = s.EbarTarget(g,:) - s.Ebar(g,:);
end

%% Nu noise eq to all internal noise
sigma2nuu = zeros(length(s.ModelParameters.ChannelToLookAt),1);
for g=1:length(s.ModelParameters.ChannelToLookAt),
    sigma2nuu(g) = sum((mu(g,:)./sigma2(g,:)).^2,2);
end

% sum accros channel.
sigma2nu=sum(sigma2nuu);

%% Distance
U = zeros(s.nifc,length(s.ModelParameters.ChannelToLookAt));
for k=1:s.nifc,
    for g=1:length(s.ModelParameters.ChannelToLookAt),
        U(k,g) = sum((mu(g,:)./sigma2(g,:)).*(squeeze(Bin(k,g,:))'-s.Ebar(g,:)));
    end
end
% sum across channel
if length(s.ModelParameters.ChannelToLookAt)>1
    V = sum(U,2);
else
    V=U;
end

% remove the variance of the noise on the target interval.
V(1) = V(1)-sqrt(sigma2nu);


%% Make decision from the model
[dummy, answer]=max(V);
% Translate x-index to interval index

if answer==1,       % the model pick up the good interval
    answer=corr;    % the place the listener would have picked up if he was correct sent
else                % if the model did not pick up the good interval
    %     if corr>1,      % one must picked up one of the other interval with a probabiloty of 0.5
    %         answer=1;
    %     else
    %         answer=2;
    %     end
        answer=corr;
        while (answer==corr)
           answer=ceil((s.nifc)*rand);
        end
%     if s.nifc==3
%         if corr>1,      % one must picked up one of the other interval with a probabiloty of 0.5
%             answer=1;
%         elseif corr<3
%             answer=3;
%         else
%             answer=2;
%         end
%     end

end




%% Update templates
if s.gamma>0,
    c1 = s.gamma/(s.gamma+1);
    c2 = 1/(s.gamma+1);
    for g=1:length(s.ModelParameters.ChannelToLookAt),
        s.Ebar(g,:) = c1*s.Ebar(g,:) + c2*ReferenceBin(g,:);
        s.EbarTarget(g,:) = c1*s.EbarTarget(g,:) + c2*TargetBin(g,:);
        s.Ebarsq(g,:) = c1*s.Ebarsq(g,:) + c2*(ReferenceBin(g,:).^2);
    end
    s.gamma = s.gamma + 1;
else
    s.gamma = 1;
end

end;

