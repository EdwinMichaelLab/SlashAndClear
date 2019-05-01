%summing up the mf levels over ages for a population average
function mSum=sumMf(mVec,ageMthMax,da,demoX)
mSum0=sum(pi_PeopleFun(((1:da:ageMthMax)/12.0)',demoX));
mSum=sum(pi_PeopleFun(((1:da:ageMthMax)/12.0)',demoX).*mVec)/mSum0;
end
