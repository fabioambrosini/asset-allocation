function s = lwcov( sample )
%Ledoit-Wolf covariance estimator
scov=cov(sample);
covtgt=mean(diag(scov))*eye(size(scov,1));
num=0;
for i=1:size(sample,1)
num=num + trace((sample(i,:)'*sample(i,:)-scov)^2);
end
num=num/size(sample,1);
gamma=1/size(sample,1)*num/(sum(diag((scov-covtgt)^2)));
gamma=max(0,min(1,gamma));
s=(1-gamma)*scov+gamma*covtgt;
end

