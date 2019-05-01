function mfSum = sumMfPrev(mfPrev,da,demog)

iAge = (1:da:length(mfPrev))';
sumPeople = sum(pi_PeopleFun(iAge/12,demog));
mfSum = sum(pi_PeopleFun(iAge/12.0,demog).*mfPrev)/sumPeople;

end