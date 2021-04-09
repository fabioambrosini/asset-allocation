clear all

%load('Returns.mat') %AMZN - MSFT - ZM - ENPH - MRNA - CWBFX - BAFWX - MC.PA - 6758.T - 9988.HK - 0P0000ZZBQ
%equities=[Returns(:,1:5) Returns(:,7:end)]; %NO CWBFX

load('ReturnsNoAlibaba.mat')
equities=[ReturnsNoAlibaba(:,1:5) ReturnsNoAlibaba(:,7:end)]; %NO CWBFX

% remove outlier
boxplot(equities,'Labels',{'AMZN','MSFT','ZM','ENPH','MRNA','BAFWX','MC.PA','6758.T','0P0000ZZBQ'})
title('Historical data - boxplot')
%ZM: outlier: +57% non è un errore ma effettivamente registrato (21-8-2020 / 1-9-2020)
%ENPH: outlier: +43% non è un errore ma effettivamente registrato (18-2-2020 / 19-2-2020)
       %outlier: -25% (16-6-2020 / 17-6-2020) 
       %oulier: - 20% (seconda settimana marzo 2020)

m=mean(equities);
S=cov(equities);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. The IS optimal frontier is biased 
% 1. compute an estimation of the market parameters 
% 2. sample several times from the market and compute the EF
% 3. plot all the EF, their mean and the "true" frontier

T=size(equities,1);  %initial sample size
C=chol(S);
NSim=60; % # of samples from the market
NPort=25; % # of ptf along the frontier
NStock=size(equities,2); % # of different stocks

% risk, return and wgt of the portfolios
frisk=zeros(NPort, NSim);
frend=zeros(NPort, NSim);
fwgt=zeros(NPort,NStock, NSim);

for i=1:NSim %simulo r per ogni stock per ogni data NSim volte
    r=repmat(m,T,1)+ randn(T,length(m))*C; %normal returns
    SSim=cov(r); 
    mSim=mean(r); %m
    p = Portfolio('assetmean', mSim, 'assetcovar', SSim, 'lowerbudget', 1, 'upperbudget', 1, 'lowerbound', 0);
    fwgt(:,:,i) = estimateFrontier(p,NPort)'; % pesi portafogli sulla frontiera efficiente
    [frisk(:,i),frend(:,i)] = estimatePortMoments(p,fwgt(:,:,i)');
    plot(frisk(:,i),frend(:,i),'--','color',[.8,.8,.8])
    hold on
end

p = Portfolio('assetmean', m, 'assetcovar', S, 'lowerbudget', 1, 'upperbudget', 1, 'lowerbound', 0);
fwgtTrue = estimateFrontier(p,NPort);
[friskTrue,frendTrue] = estimatePortMoments(p,fwgtTrue);
plot(friskTrue,frendTrue,'r') 
friskMean=mean(frisk,2);
frendMean=mean(frend,2);
plot(friskMean,frendMean,'b') 





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. stimatori shrinkage
% same as before, but compare with the EF obtained from shrinkage
% estimators

T=size(equities,1); %SMALLER  
C=chol(S);
NSim=60; 
NPort=25; 
NStock=size(equities,2); 

% risk, return and wgt of the portfolios
frisk=zeros(NPort, NSim);
frend=zeros(NPort, NSim);
fwgt=zeros(NPort,NStock, NSim);

friskShr=zeros(NPort, NSim);
frendShr=zeros(NPort, NSim);
fwgtShr=zeros(NPort,NStock, NSim);

for i=1:NSim
    r=repmat(m,T,1)+ randn(T,length(m))*C;
    SSim=cov(r);
    mSim=mean(r);
    p = Portfolio('assetmean', mSim , 'assetcovar', SSim, 'lowerbudget', 1, 'upperbudget', 1, 'lowerbound', 0);
    fwgt(:,:,i) = estimateFrontier(p,NPort)';
    [frisk(:,i),frend(:,i)] = estimatePortMoments(p,fwgt(:,:,i)');

    SShr=lwcov(r);
    mtgt=mean(mean(r))*ones(1,NStock);
    mShr=jsmean(r,mtgt);
    q = Portfolio('assetmean', mShr , 'assetcovar', SShr, 'lowerbudget', 1, 'upperbudget', 1, 'lowerbound', 0);
    fwgtShr(:,:,i) = estimateFrontier(q,NPort)';
    [friskShr(:,i),frendShr(:,i)] = estimatePortMoments(q,fwgtShr(:,:,i)');
    plot(frisk(:,i),frend(:,i),'--','color',[.6,.6,.8]) 
    hold on
    plot(friskShr(:,i),frendShr(:,i),'--','color',[.8,.6,.6]) 
end

p = Portfolio('assetmean', m , 'assetcovar', S, 'lowerbudget', 1, 'upperbudget', 1, 'lowerbound', 0);
fwgtTrue = estimateFrontier(p,NPort)';
[friskTrue(:,i),frendTrue(:,i)] = estimatePortMoments(p,fwgtTrue');

figure()
plot(friskTrue,frendTrue,'r') 
friskMean=mean(frisk,2);
frendMean=mean(frend,2);
hold on
plot(friskMean,frendMean,'b') 
friskShrMean=mean(friskShr,2);
frendShrMean=mean(frendShr,2);
plot(friskShrMean,frendShrMean,'g') 

figure()
bar(fwgtShr(:,:,15),'stacked')
figure()
bar(fwgtTrue(:,:),'stacked')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2.5 Resampled and classical frontier
T=size(equities,1);
C=chol(S);
NRes=100;
NSim=60;
NPort=25;
NStock=size(equities,2);

frisk=zeros(NPort, NRes);
frend=zeros(NPort, NRes);
fwgt=zeros(NPort,NStock, NRes);

for i=1:NRes
    r=repmat(m,T,1)+ randn(T,length(m))*C;
    SSim=cov(r);
    mSim=mean(r);
    p = Portfolio('assetmean', mSim , 'assetcovar', SSim, 'lowerbudget', 1, 'upperbudget', 1, 'lowerbound', 0);
    fwgt(:,:,i) = estimateFrontier(p,NPort)';
    [frisk(:,i),frend(:,i)] = estimatePortMoments(p,fwgt(:,:,i)');
    %[frisk(:,i),frend(:,i),fwgt(:,:,i)]=frontcon(mSim*.1,SSim,NPort);
end
fwgtRE=mean(fwgt,3);
frendRE=fwgtRE*m';
friskRE=sqrt(diag(fwgtRE*S*fwgtRE'));
figure()
plot(friskRE,frendRE,'bx-')  
title('Efficient Frontier')
xlabel('Volatility')
ylabel('Return')
hold on

p = Portfolio('assetmean', m, 'assetcovar', S, 'lowerbudget', 1, 'upperbudget', 1, 'lowerbound', 0);
fwgtTrue = estimateFrontier(p,NPort)';
[friskTrue,frendTrue] = estimatePortMoments(p,fwgtTrue');
%[friskTrue,frendTrue,fwgtTrue]=frontcon(m,S,NPort);
plot(friskTrue,frendTrue,'r') 
legend('Resampled Frontier','True Frontier')
figure()                                                                 
bar(fwgtTrue,'stacked')
figure()
bar(fwgtRE,'stacked')
title('Portfolio Composition')
xlabel('Portfolios')
ylabel('Weights')
