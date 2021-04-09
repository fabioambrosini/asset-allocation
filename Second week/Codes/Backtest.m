clear all; 
clc;
load('DatesPrices.mat');

for ii = 1:size(Prices,2) %index of the first NOT NaN price for each ticker
    index = find(~isnan(Prices(:,ii)),1);
    
    id = Ticker(ii);
    quotes = Prices(index:end,ii); 
    dates = Dates(index+1:end,ii);
    Returns = quotes(2:end)./quotes(1:end-1) -1; %from the oldest date to the newest 

    SampleSize = size(Returns,1);
    TestWindowStart = 251; %lookback period of 1 year   
    TestWindow = TestWindowStart : SampleSize;
    EstimationWindowSize = 250; %yearly

    alpha = 0.99;
    VaR = zeros(length(TestWindow),1);

    for t = TestWindow
        jj = t - TestWindowStart + 1;
        EstimationWindow = t-EstimationWindowSize:t-1; %1 year backward
        X = Returns(EstimationWindow); 
        VaR(jj) = -quantile(X,1-alpha); 
    end

    figure;
    plot(dates(TestWindowStart:end),-VaR,'r')
    ylabel('VaR')
    xlabel('Date')
    hold on 
    plot(dates(TestWindowStart:end),Returns(TestWindow),'b') 
    legend({'99% Confidence Level'},'Location','NorthWest')
    title('VaR - Historical Simulation - ' + id)
    
    ReturnsTest = Returns(TestWindow);
    vbt = varbacktest(ReturnsTest,VaR,'PortfolioID',id,'VaRLevel', 0.99);
    disp('all period')
    summary(vbt)
    runtests(vbt)
    if(id~="ZM")
       StartCovidPeriod = find(dates >= '10-feb-2020',1); %Covid Index
       ReturnsTest = Returns(StartCovidPeriod:end);
       figure;
       plot(dates(StartCovidPeriod:end),-VaR(StartCovidPeriod-EstimationWindowSize:end),'r')
       ylabel('VaR')
       xlabel('Date')
       hold on 
       plot(dates(StartCovidPeriod:end),ReturnsTest,'b') 
       legend({'99% Confidence Level'},'Location','NorthWest')
       title('VaR - Historical Simulation - Covid period - ' +id)


       vbt = varbacktest(ReturnsTest,VaR(StartCovidPeriod-EstimationWindowSize:end),'PortfolioID',id,'VaRLevel', 0.99);
       disp('Covid period')
       summary(vbt)
       runtests(vbt) 
    end 
end 