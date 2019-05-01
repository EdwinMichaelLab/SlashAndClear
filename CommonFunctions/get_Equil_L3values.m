function [l3,flag,pVec,wVec,mVec,iVec] = get_Equil_L3values(bit,beta,alpha,k0,kLin,k1,r1,...
    sigma1,psi1,psi2s2,mu,gamma,b1,c,ageLev,k2,gam2,tau,immC,slopeC,...
    immCMin,larvaeDeath,PP,del,excessMortality,k0Spp,ageMthMax,da,bCulex,...
    demoX,tol1)

VoverH = bit/beta;
l3 = ones(size(VoverH)); % Initial guess value
l3VoverH = sum(l3.*VoverH);

% note that typically da = 1
gVec  = zeros((ageMthMax/da),1);   % used in the equation for wVec (gVec is never used but in equilibiumValuesOsStateVars_CFA, might be defined within that function)
pVec  = zeros((ageMthMax/da),1); % prepatent worm burden
wVec  = zeros((ageMthMax/da),1); % patent worm burden
mVec  = zeros((ageMthMax/da),1); % mf intensity
iVec  = zeros((ageMthMax/da),1); % measure of immunity

[pVec,wVec,mVec,iVec] = OnchoStateVarsEquilValues(beta,...
    k0,kLin,psi1,mu,alpha,gamma,c,ageLev,l3VoverH,psi2s2,immC,slopeC,...
    del,ageMthMax,da,immCMin,wVec,mVec,iVec,pVec,gVec,PP);

for iSpp = 1:length(VoverH)
    l3(iSpp) = calcL3Fun_Oncho_exmort(mVec,beta,k0Spp(iSpp),kLin,k1,r1,sigma1,psi1,...
        b1,k2,gam2,tau,demoX,ageMthMax,da,ageLev,bCulex(iSpp),larvaeDeath,excessMortality);
end
l3VoverH = sum(VoverH.*l3);
flag = 0;
while (flag < length(VoverH))
    [pVec,wVec,mVec,iVec] = OnchoStateVarsEquilValues(beta,...
        k0,kLin,psi1,mu,alpha,gamma,c,ageLev,l3VoverH,psi2s2,immC,slopeC,...
        del,ageMthMax,da,immCMin,wVec,mVec,iVec,pVec,gVec,PP);
    
    temp = l3;
    for iSpp = 1:length(VoverH)
        l3(iSpp) = calcL3Fun_Oncho_exmort(mVec,beta,k0Spp(iSpp),kLin,k1,r1,sigma1,psi1,...
            b1,k2,gam2,tau,demoX,ageMthMax,da,ageLev,bCulex(iSpp),larvaeDeath,excessMortality);
    end
    flag = 0;
    for iSpp = 1:length(l3)
        if abs(l3(iSpp)-temp(iSpp)) <=  tol1
            flag = flag+1;
        end
    end
    l3VoverH = sum(VoverH.*l3);
end
end
