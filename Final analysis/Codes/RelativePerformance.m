clear all
clc

initialDate = 0; %27/10 one day before start
secondDate = 10; %11/11 one day before first change 
thirdDate = 21; %26/11 one day before second change
endDate = 25; %02/12 last day


load('RendPtf'); 
load('RendBench');

date1 = initialDate+1;
date2 = secondDate;
covariance = cov(RendPtf(date1:date2),RendBench(date1:date2));
corr = covariance(1,2) / ( std(RendPtf(date1:date2))*std(RendBench(date1:date2)) ); 
R_squared = corr^2;

x = [1:1:25];
y = [RendPtf, RendBench];
bar(x,y)
title('Portfolio returns & Benchmark returns')
xlabel('Dates')
ylabel('Returns')
legend('portfolio','benchmark')

