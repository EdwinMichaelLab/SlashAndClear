% Function: calcL3Fun
%  
% Description: 
%       Function to calculate the value of L3 given the mf intensity
%       (mVec) and fly species. See Gambhir and Michael 2008.
%
% Inputs:
%       mVec,k0,kLin,r1,k1,beta,b1,sigma1,psi1,ageMthMax,da,k2,gam2,bCulex,demoX
%
% Outputs:
%       l3
% ________________________________________

function l3 = calcL3Fun_Oncho_exmort(mVec,beta,k0,kLin,k1,r1,sigma1,psi1,...
    b1,k2,gam2,tau,demoX,ageMthMax,da,ageLev,bCulex,larvaeDeath,excessMortality)

a2=(1:da:ageMthMax)';
sumPeople = sum(pi_PeopleFun(a2/12.0,demoX));
pI = pi_PeopleFun(a2/12.0,demoX)/sumPeople;

if bCulex ==1
    fI=1-densityDependenceFunCulex(mVec,...
        negBinshapeFun(mVec,k0,kLin),(r1/k1));
    pI = pi_PeopleFun(a2/12.0,demoX)/sumPeople;
    int1 =sum(pI.*fI);
    l3 = beta*k1*b1*int1/...
        (sigma1+beta*psi1+larvaeDeath+excessMortality*...
        Extra_mf_InducedMortality(sumPeople,demoX,da,ageMthMax,ageLev,mVec));
else
    fI=1-densityDependenceFunAnophWithoutTau(mVec,...
        negBinshapeFun(mVec,k0,kLin),k2,gam2,tau);
    pI = pi_PeopleFun(a2/12.0,demoX)/sumPeople;
    int1 = sum(pI.*fI);
    l3 = beta*k2*b1*int1/...
        (sigma1+beta*psi1+larvaeDeath+excessMortality*...
        Extra_mf_InducedMortality(sumPeople,demoX,da,ageMthMax,ageLev,mVec));
end

end

function f=densityDependenceFunCulex(m,k,gam)

f = (1+((m./k)*(1-exp(-gam)))).^(-k);
end

function f=densityDependenceFunAnophWithoutTau(m,k,k2,gam2,tau)

val2=gam2/k2;

fTerm1=(2*exp(val2*tau))./((1+ (m./k)*(1-exp(-val2))).^k);
fTerm2=(exp(2*val2*tau))./((1+ (m./k)*(1-exp(-2*val2))).^k);

f = (fTerm1-fTerm2);
end

function fVal0=Extra_mf_InducedMortality(sumPeople,demoX,da,ageMthMax,ageLev,mVec)
fVal0=0;
iAgeCa = 1;
for a2=1:da:ageMthMax
    pI = pi_PeopleFun(a2/12.0,demoX)/sumPeople;
    fVal0 = fVal0 + pI*min(a2/(12.0*ageLev),1)*mVec(iAgeCa);
    iAgeCa = iAgeCa + 1;
end
end
