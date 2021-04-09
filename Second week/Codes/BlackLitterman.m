clear all
clc

load('PtfPrices.mat') %AMZN - MSFT - ZM - ENPH - MRNA - CWBFX - BAFWX - MC.PA - 6758.T  - 0P0000ZZBQ

%We do not consider the bond fund since its weight is fixed 
equities=[PtfPrices(:,1:5) PtfPrices(:,7:end)]; %NO CWBFX
Returns= equities(1:end-1,:)./equities(2:end,:) - 1; %dalla piu recente alla piu vecchia
t=size(Returns,1);
% AMZN - MSFT - ZM - ENPH - MRNA - BAFWX - MC.PA - 6728.T - 0P0000ZZBQ
cap = [1.598e+12, 1.617e+12, 126.203e+09, 17.57e+09, 43.204e+09, 3.87e+09, 249.799e+09, 12.033e+12, 4.58e+09]; 
change = [0.8403, 0.8403, 0.8403, 0.8403, 0.8403, 0.8403, 1,0.0080, 1.1236];
capEUR = (cap.*change)'; %column vector 
%boxplot(Returns)

m=mean(Returns)';
S=cov(Returns);

returns_CWBFX = (PtfPrices(1:end-1,6)./PtfPrices(2:end,6) - 1)'; %row vector
load('BenchPrices.mat') %MIWO00000PUS - AGGH
returns_b = BenchPrices(1:end-1,:)./BenchPrices(2:end,:) - 1;
%% Black and Litterman
%Qualitative view - AMZN
%implied views
lambda=1.2;
w=capEUR/sum(capEUR); %column vector
mBL=2*lambda*S*w; %muBAR COLUMN VECTOR
%bar((1+[mBL m]).^250-1) %rendimento annuale

%qualitative view AMZN
v = zeros(1,size(Returns,2));
v(1,1) = 1; %AMZN
c = 1;
q = 2; %view su AMZN per riga sopra
mview = v*mBL + q.*sqrt(v*S*(v')); %mu_v
Sview = diag(c)^(1/2) *v*S*(v') *diag(c)^1/2;
mBLP=mBL+S*v'*((v*S*v'/t+Sview)\(mview-v*mBL))/t;
SBLP=(1+1/t)*S-S*v'*((v*S*v'/t+Sview)\v*S)/(t^2);
X = categorical({'AMZN', 'MSFT', 'ZM', 'ENPH', 'MRNA', 'BAFWX', 'MC.PA', '6758.T', '0P0000ZZBQ'});
X = reordercats(X,{'AMZN', 'MSFT', 'ZM', 'ENPH', 'MRNA', 'BAFWX', 'MC.PA', '6758.T', '0P0000ZZBQ'});
Y = (1+[mBLP mBL]).^250-1;
bar(X,Y);
%% Efficient Frontier
T = t; 
C=chol(S);
Cprior=(1+1/T)*C;
NSim=60; % # of samples from the market
NPort=25; % # of ptf along the frontier
NStock=size(Returns,2); % # of different stocks
% risk, return and wgt of the portfolios
frisk_prior=zeros(NPort, NSim);
frend_prior=zeros(NPort, NSim);
fwgt_prior=zeros(NPort,NStock, NSim);

figure(1)
for i=1:NSim %simulo r per ogni stock per ogni data NSim volte
    %prior distribution    
    r=repmat(mBL',T,1)+ randn(T,length(mBL'))*Cprior; %normal returns
    SSim=cov(r); 
    mSim=mean(r); 
    p = Portfolio('assetmean', mSim, 'assetcovar', SSim, 'lowerbudget', 1, 'upperbudget', 1, 'lowerbound', 0);
    fwgt_prior(:,:,i) = estimateFrontier(p,NPort)'; % pesi portafogli sulla frontiera efficiente
    [frisk_prior(:,i),frend_prior(:,i)] = estimatePortMoments(p,fwgt_prior(:,:,i)');
    plot(frisk_prior(:,i),frend_prior(:,i),'--','color',[.8,.8,.8])
    hold on
end
friskMean_prior=mean(frisk_prior,2);
frendMean_prior=mean(frend_prior,2);
plot(friskMean_prior,frendMean_prior,'b')
title('EF - Prior distribution')

figure(2)
h=bar(mean(fwgt_prior,3),'stacked');
set(h(8),'facecolor',[0 1 1]);
set(h(9),'facecolor',[1 1 0]);
title('Prior distribution')

%posterior distribution 
Cposterior=chol(SBLP);
NSim=60; % # of samples from the market
NPort=25; % # of ptf along the frontier
NStock=size(Returns,2); % # of different stocks
% risk, return and wgt of the portfolios
frisk_posterior=zeros(NPort, NSim);
frend_posterior=zeros(NPort, NSim);
fwgt_posterior=zeros(NPort,NStock, NSim);

figure(3)
for i=1:NSim %simulo r per ogni stock per ogni data NSim volte
    %posterior distribution    
    r=repmat(mBLP',T,1)+ randn(T,length(mBLP'))*Cposterior; %normal returns
    SSim=cov(r); 
    mSim=mean(r); 
    p = Portfolio('assetmean', mSim, 'assetcovar', SSim, 'lowerbudget', 1, 'upperbudget', 1, 'lowerbound', 0);
    fwgt_posterior(:,:,i) = estimateFrontier(p,NPort)'; % pesi portafogli sulla frontiera efficiente
    [frisk_posterior(:,i),frend_posterior(:,i)] = estimatePortMoments(p,fwgt_posterior(:,:,i)');
    plot(frisk_posterior(:,i),frend_posterior(:,i),'--','color',[.8,.8,.8])
    hold on
end
friskMean_posterior=mean(frisk_posterior,2);
frendMean_posterior=mean(frend_posterior,2);
plot(friskMean_posterior,frendMean_posterior,'r') 
title('Posterior distribution')
figure(4)
k=bar(mean(fwgt_posterior,3),'stacked');
set(k(8),'facecolor',[0 1 1]);
set(k(9),'facecolor',[1 1 0]);
title('Posterior distribution')

figure(5)
plot(friskMean_prior,frendMean_prior,'b')  
hold on
plot(friskMean_posterior,frendMean_posterior,'r') 
title('EF - Prior and Posterior distributions')
legend('Prior distribution','Posterior distribution','location','SouthEast')

[maxIR, index]  = IRcomputation(NPort,mean(fwgt_posterior,3),Returns,returns_CWBFX,returns_b);