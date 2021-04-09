function m=jsmean(sample,  mtgt)
% James Stein estimator for the mean
smean=mean(sample);
scov=cov(sample);
rho=max(eig(scov));
gamma=1/size(sample,1)*(sum(diag(scov)) -2*rho)/((smean-mtgt)*(smean-mtgt)');
gamma=max(0,min(1,gamma));
m=(1-gamma)*smean+gamma*mtgt;