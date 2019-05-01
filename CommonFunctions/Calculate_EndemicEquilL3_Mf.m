function [L3Values,mfPrevArray] = Calculate_EndemicEquilL3_Mf(NParamVecs,...
    ageMthMax,bCulex,demoX,ParameterVectors,da,toleranceX)

mfPrevArray = zeros(ageMthMax,NParamVecs);
L3Values=zeros(NParamVecs,length(bCulex));

gVec  = zeros((ageMthMax/da),1);
pVec  = zeros((ageMthMax/da),1);
wVec  = zeros((ageMthMax/da),1);
mVec  = zeros((ageMthMax/da),1);
iVec  = zeros((ageMthMax/da),1);
OvVec = zeros((ageMthMax/da),1);
IgVec = zeros((ageMthMax/da),1);

% loop through each parameter vector i
parfor i = 1:NParamVecs
    
    % access the parameters of index i
    [beta,HBIndex,GTCP,alpha,k1,r1,sigma1,psi1,psi2s2,mu,gamma,b1,c,...
        ageLev,k2,gam2,tau,immC,slopeC,k0_An,k0_Cx,immCMin,larvaeDeath,...
        PP,del,alpha2,gamma2,gamma3,k0,kLin,k0_Ag,kLin_Ag,k0_Ab,...
        kLin_Ab,alpha3max,h1,h2,VoverH] = get_theParameters(ParameterVectors,i);
    
    % given the sampled parameter vector, calculate the equilibrium values
    % of L3 and mf
    k0Spp=assign_k0Spp(bCulex,k0_Cx,k0_An);
    [l3Eq,~,~,~,mVec0,~,~,~] = get_Equil_L3values_30Mar2016(VoverH,beta,alpha,k0,kLin,...
        k1,r1,sigma1,psi1,psi2s2,mu,gamma,b1,c,ageLev,k2,gam2,tau,immC,slopeC,...
        immCMin,larvaeDeath,PP,del,alpha2,gamma2,alpha3max,h1,h2,gamma3,...
        k0Spp,ageMthMax,da,bCulex,demoX,toleranceX,gVec,pVec,...
        wVec,mVec,iVec,OvVec,IgVec);
    
    % store equilibrium values given parameter vector i
    L3Values(i,:)=l3Eq;
    mfPrevArray(:,i)=mfAgeprevFun(mVec0,negBinshapeFun(mVec0,k0,kLin));
    
    
end
end