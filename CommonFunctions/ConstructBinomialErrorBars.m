% binomial confidence interval 
function aa = ConstructBinomialErrorBars(data)

p=(data(:,3)+2)./(data(:,2)+4); % proportion testing positive (adjusted wald method)
delta=1.96*sqrt(((p.*(1-p))./(data(:,2)+4))); % +/- term given alpha=0.05

lower=max(p-delta,0); % lower bound of 95% CI
upper=min(p+delta,1); % upper bound of 95% CI

MfPrev=data(:,3)./data(:,2);
aa=[data(:,1) MfPrev MfPrev-lower upper-MfPrev]; % vector holding info for use with errorbar function
aa(:,[2 3 4]) = aa(:,[2 3 4])*100;

end