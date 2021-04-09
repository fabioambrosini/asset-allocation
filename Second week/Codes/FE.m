clear all
clc

%Portfolio
load('PtfPrices.mat') %AMZN - MSFT - ZM - ENPH - MRNA - CWBFX - BAFWX - MC.PA - 6758.T  - 0P0000ZZBQ

%We do not consider the bond fund since its weight is fixed 
equities=[PtfPrices(:,1:5) PtfPrices(:,7:end)]; %NO CWBFX
returns_p = equities(1:end-1,:)./equities(2:end,:) - 1;
m=mean(returns_p);
S=cov(returns_p);

returns_CWBFX = (PtfPrices(1:end-1,6)./PtfPrices(2:end,6) - 1)'; %row vector

%Remove outliers
figure(1)
boxplot(returns_p,'Labels',{'AMZN','MSFT','ZM','ENPH','MRNA','BAFWX','MC.PA','6758.T','0P0000ZZBQ'})
title('Historical data - boxplot')


%Benchmark
load('BenchPrices.mat') %MIWO00000PUS - AGGH
returns_b = BenchPrices(1:end-1,:)./BenchPrices(2:end,:) - 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. The IS optimal frontier is biased 
% 1. compute an estimation of the market parameters 
% 2. sample several times from the market and compute the EF
% 3. plot all the EF, their mean and the "true" frontier

T=size(returns_p,1);  %initial sample size
C=chol(S);
NSim=60; % # of samples from the market
NPort=25; % # of ptf along the frontier
NStock=size(returns_p,2); % # of different stocks

% risk, return and wgt of the portfolios
frisk=zeros(NPort, NSim);
frend=zeros(NPort, NSim);
fwgt=zeros(NPort,NStock, NSim);

figure(2)
for i=1:NSim %simulo r per ogni stock per ogni data NSim volte
    r=repmat(m,T,1)+ randn(T,length(m))*C; %normal returns
    SSim=cov(r); 
    mSim=mean(r);
    p = Portfolio('assetmean', mSim, 'assetcovar', SSim, 'lowerbudget', 1, 'upperbudget', 1, 'lowerbound', 0);
    fwgt(:,:,i) = estimateFrontier(p,NPort)'; 
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


figure(3)
h=bar(fwgt(:,:,15),'stacked');
set(h(8),'facecolor',[0 1 1]);
set(h(9),'facecolor',[1 1 0]);
title('Portfolio composition  - IS')
xlabel('Portfolios')
ylabel('Weights')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. stimatori shrinkage
% same as before, but compare with the EF obtained from shrinkage
% estimators

T=size(returns_p,1); 
C=chol(S);
NSim=60; 
NPort=25; 
NStock=size(returns_p,2); 

% risk, return and wgt of the portfolios
frisk=zeros(NPort, NSim);
frend=zeros(NPort, NSim);
fwgt=zeros(NPort,NStock, NSim);

friskShr=zeros(NPort, NSim);
frendShr=zeros(NPort, NSim);
fwgtShr=zeros(NPort,NStock, NSim);

figure(4)
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

figure(5)
plot(friskTrue,frendTrue,'r') 
friskMean=mean(frisk,2);
frendMean=mean(frend,2);
hold on
plot(friskMean,frendMean,'b') 
friskShrMean=mean(friskShr,2);
frendShrMean=mean(frendShr,2);
plot(friskShrMean,frendShrMean,'g') 

figure(6)
h=bar(fwgtShr(:,:,15),'stacked');
set(h(8),'facecolor',[0 1 1]);
set(h(9),'facecolor',[1 1 0]);
title('Portfolio composition - Shrinkage estimator')
xlabel('Portfolios')
ylabel('Weights')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2.5 Resampled and classical frontier
T=size(returns_p,1);
C=chol(S);
NRes=100;
NSim=60;
NPort=25;
NStock=size(returns_p,2);

frisk=zeros(NPort, NRes);
frend=zeros(NPort, NRes);
fwgt=zeros(NPort,NStock, NRes);

figure(7)
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
plot(friskRE,frendRE,'c')  
title('Efficient Frontier')
xlabel('Standard Deviation')
ylabel('Return')
hold on

p = Portfolio('assetmean', m, 'assetcovar', S, 'lowerbudget', 1, 'upperbudget', 1, 'lowerbound', 0);
fwgtTrue = estimateFrontier(p,NPort)';
[friskTrue,frendTrue] = estimatePortMoments(p,fwgtTrue');
%[friskTrue,frendTrue,fwgtTrue]=frontcon(m,S,NPort);
plot(friskTrue,frendTrue,'r') 
legend('Resampled Frontier','True Frontier')
                                                                
[maxIR, index]  = IRcomputation(NPort,fwgtRE,returns_p,returns_CWBFX,returns_b);

figure(8)
h=bar(1,fwgtRE(index,:),0.25,'stacked');
set(h(8),'facecolor',[0 1 1]);
set(h(9),'facecolor',[1 1 0]);
title('Portfolio composition - Resampled frontier')
xlabel('Portfolios')
ylabel('Weights')
%% Final plots 
figure(9)
plot(friskTrue,frendTrue,'r')
hold on 
plot(friskMean,frendMean,'b')
plot(friskShrMean,frendShrMean,'g') 
plot(friskRE,frendRE,'k')
title('Efficient Frontiers')
legend('True','IS','Shrinkage estimator','Resampled frontier','location','SouthEast')
xlabel('Standard Deviation')
ylabel('Return')


ExcelWeights = (1/0.65)/100*[0, 9.269818, 1.771226, 5.251363, 3.848178, 19.193716, 10.515907, 8.408145, 6.741647]; %1 x 9
p = Portfolio('assetmean', m, 'assetcovar', S, 'lowerbudget', 1, 'upperbudget', 1, 'lowerbound', 0);
[risk,rend] = estimatePortMoments(p,ExcelWeights');

figure(10)
h=bar(1,ExcelWeights,0.25,'stacked');
set(h(8),'facecolor',[0 1 1]);
set(h(9),'facecolor',[1 1 0]);
title('Portfolio composition')
xlabel('Maximum IR portfolio')
ylabel('Weights')

figure(11)
plot(risk,rend,'bx','LineWidth',2,'MarkerSize',10)
hold on
plot(friskRE,frendRE,'k')
plot(friskRE(index),frendRE(index),'ro','LineWidth',2,'MarkerSize',10)
title('Efficient Frontier and maximization of IR')
xlabel('Standard Deviation')
ylabel('Return')
legend('maximum IR portfolio','Efficient Frontier','EF portfolio with maximum IR','location','SouthEast')
plot(friskRE,frendRE,'kx')  