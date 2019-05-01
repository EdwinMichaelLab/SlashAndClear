function [L3Values,mfPrevArray] = ...
    Calculate_EndemicEquil_L3_Mf(NParamVecs,...
    ageMthMax,bCulex,demoX,ParameterVectors,VoverH,da,toleranceX)

% initialize arrays to store equilibrium values
L3Values = zeros(NParamVecs,length(bCulex));
mfPrevArray = zeros(ageMthMax,NParamVecs,'single');

% loop through each parameter vector i
parfor i = 1:NParamVecs
    
    % access the parameters of index i
    [beta,~,~,alpha,k0,kLin,k1,r1,sigma1,psi1,psi2s2,mu,gamma,b1,c,...
        ageLev,k2,gam2,tau,immC,slopeC,k0_An,k0_Cx,immCMin,larvaeDeath,...
        PP,del,excessMortality] = get_theParameters(ParameterVectors,i);
    
    k0Spp=assign_k0Spp(bCulex,k0_Cx,k0_An);
    % given the sampled parameter vector, calculate the equilibrium values
    % of L3 and mf
    [l3Eq,~,~,~,mVec0,~] = get_Equil_L3values(VoverH(i,:)*beta,beta,alpha,k0,kLin,k1,r1,...
        sigma1,psi1,psi2s2,mu,gamma,b1,c,ageLev,k2,gam2,tau,immC,slopeC,...
        immCMin,larvaeDeath,PP,del,excessMortality,k0Spp,ageMthMax,da,bCulex,...
        demoX,toleranceX);
    
    % store equilibrium values given parameter vector i
    L3Values(i,:) = l3Eq; 
    mfPrevArray(:,i) = mfAgeprevFun(mVec0,negBinshapeFun(mVec0,k0,kLin));
    
end

end