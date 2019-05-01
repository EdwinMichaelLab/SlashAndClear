function [mfPrevIntv,MBRIntv,L3Intv,pIntv,wIntv,MfInt] = RunOnchoInterventionScenarios_slash(k,...
    ParameterVectors,L3Values,demoX,ageMthMax,da,bCulex,...
    AgeLimits,drugCoverage,VillageDrugEfficacy,...
    MDAInterval,NumYears,MultiVecMBR,VC0,slash_t,TBR)

mfPrevIntv = zeros(1+NumYears*12,length(k));
MBRIntv = zeros(1+NumYears*12,length(k));
L3Intv = zeros(1+NumYears*12,length(k));
pIntv = zeros(1+NumYears*12,length(k));
wIntv = zeros(1+NumYears*12,length(k));
MfInt = zeros(1+NumYears*12,length(k));

pVec = zeros((ageMthMax/da),1);
wVec = zeros((ageMthMax/da),1);
mVec = zeros((ageMthMax/da),1);
iVec = zeros((ageMthMax/da),1);
gVec = zeros((ageMthMax/da),1);

% slash parameters
% MBR_max, R_L, R_U, k1, k2, n, A, b1, b2, b3, b4
load('MBRFunParams.mat');
R = [R(end,:);R(1:end-1,:)]; %start sims in April at the beginning of rainy season

parfor i = 1:length(k)
    
    l3 = L3Values(i,:);
    
    [beta,~,~,alpha,k0,kLin,k1,r1,sigma1,psi1,psi2s2,mu,gamma,b1,c,...
        ageLev,k2,gam2,tau,immC,slopeC,k0_An,k0_Cx,immCMin,larvaeDeath,...
        PP,del,excessMortality] = get_theParameters(ParameterVectors,i);
    
    % adjust VoverH over time according to the application of larvicide
    % treatments
    VH = Temporal_Change_in_VoverH_due_to_VC_slash(VC0,...
        NumYears,beta,P(:,i),R(:,i),slash_t,MultiVecMBR(i));
    VoverH = VH(1,:);

    l3VoverH = sum(VoverH.*l3);
    
    [pVec0,wVec0,mVec0,iVec0] = ...
        OnchoStateVarsEquilValues(beta,...
        k0,kLin,psi1,mu,alpha,gamma,c,ageLev,l3VoverH,psi2s2,immC,slopeC,...
        del,ageMthMax,da,immCMin,wVec,mVec,iVec,pVec,gVec,PP);
    
    k0Spp = assign_k0Spp(bCulex,k0_Cx,k0_An);
    
    % model the given interventions
    [mfPrevArray,L3Values0,pArray,wArray,MfInt0] = ModelOnchoIntv(beta,...
        alpha,k0,kLin,k1,r1,sigma1,psi1,psi2s2,mu,gamma,b1,c,ageLev,k2,...
        gam2,tau,immC,immCMin,slopeC,PP,VH,del,ageMthMax,da,...
        bCulex,l3,larvaeDeath,excessMortality,gVec,pVec0,wVec0,mVec0,iVec0,demoX,...
        AgeLimits,drugCoverage,VillageDrugEfficacy(mod(i,100)+1,:),...
        MDAInterval,NumYears,k0Spp,TBR(i));
    
    mfPrevIntv(:,i) = mfPrevArray;
    MBRIntv(:,i) = VH.*beta;
    L3Intv(:,i) = L3Values0;
    pIntv(:,i) = pArray;
    wIntv(:,i) = wArray;
    MfInt(:,i) = MfInt0;
    
end
end

