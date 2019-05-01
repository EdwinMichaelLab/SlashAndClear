% Function: RunSelectParameterVectors
%
% Description: Runs through parameter vectors for baseline fitting and likelihood
%              calculation, adjusts for mf data quality.
% Inputs:
%   Site,NParamVecs,MfData,ABR,SIR_samples,demoA,demoB,ageMax,bCulex,
%       minPrev,maxPrev
%   Note: Not all inputs will be used depending on quality of the MfData.
%
% Outputs:
%   Writes to file, '../IO/OUT/ParamVectors%s.mat', Site.
%   Writes: bCulex, demog, ABR, ParameterVectors, L3Values, ageMthMax,
%           mfPrevArray.
% ________________________________________

function RunSelectParameterVectors(Site,NParamVecs,MfData,ABR,...
    SIR_samples,demoA,demoB,ageMax,bCulex,minPrev,maxPrev,ABRmax,...
    estABRflag)

ageMthMax = ageMax*12; % convert max age in years to months
demog = [1, 1, demoA, demoB, 1]; % set up demographic parameters

da = 1; % integration step-size of 1 month
toleranceX = 0.000001; % Used in endemic equilibrium of state variables

% Set random number stream
RandStream.setGlobalStream...
    (RandStream('mt19937ar','seed',rand*sum(10000*clock))); 

% Initialize Output Arrays
ParameterVectors = [];
L3Values1 = [];
mfPrevArray1 = [];
ABR1 = [];

% Loop to generate an initial sample of parameter vectors, run the model
% to equilibrium, and resample fits until the desired number of best 
% fits defined by SIR_samples is achieved. 

% Obtain the boolean variables needed to process MfData.
% Each of the following variables represent the quality of the MfData
[AgeMf_asInput, OverallMf_asInput, OverallMfRange_asInput] = get_dataQualityVars(MfData);

NumParam = 0;
numSamples = SIR_samples;

if OverallMfRange_asInput % No Baseline Mf, different loop condition in this case.
    numSamples = NParamVecs; 
end

warning('off','all'); % suppress warning messages related to parfor
options=optimset('MaxIter',25,'TolFun',1e-2,'TolX',1e-2); % fzero options for estimating ABR case

while ( NumParam < numSamples ) % numSamples = either SIR_samples or NParamVecs, depending on quality of MfData
      
    % calculate VoverH (this is a vector for one species, matrix for more than one species)
    % if ABR = NaN, no data is available so we will estimate from a prior
    % range defined in inputs
    % if ABR = value and estABRflag = 1 then we will estimate the ABR given an initial guess
    if isnan(ABR) || (~isnan(ABR) && (estABRflag)) 
        
        OverallMfPrev = sum(MfData(:,3))/sum(MfData(:,2));
        if isnan(ABR) % to estimate ABR from a prior range
            max_bit = ABRmax/12;
        elseif (estABRflag) % to estimate ABR given an initial guess to inform the prior range
            max_bit = 2*ABR/12;
        end
        
        [rangeParamVals,~,~] = Range_of_Parameters; % rangeParamVals used for determining how many parameters there are for allocating space in next line
        ParamVectors = zeros(length(rangeParamVals),NParamVecs);
        parfor i = 1:NParamVecs
            diff_mfprev=0;
            while abs(diff_mfprev - OverallMfPrev) <= OverallMfPrev
                ParamVectors0 = ParameterVectors_LHSDesign_based(1);
                [beta,~,~,alpha,k0,kLin,k1,r1,sigma1,psi1,psi2s2,mu,gamma,b1,c,...
                    ageLev,k2,gam2,tau,immC,slopeC,k0_An,k0_Cx,immCMin,larvaeDeath,...
                    PP,del,excessMortality] = get_theParameters(ParamVectors0,1);
                
                diff_mfprev = finding_VoverH_w_fzero(max_bit,OverallMfPrev,beta,...
                    alpha,k0,kLin,k1,r1,sigma1,psi1,psi2s2,mu,gamma,b1,c,ageLev,...
                    k2,gam2,tau,immC,slopeC,del,ageMthMax,da,demog,...
                    bCulex,toleranceX,k0_An,k0_Cx,immCMin,larvaeDeath,PP,excessMortality);
            end
            
            bitx = fzero(@(bitx)finding_VoverH_w_fzero(bitx,OverallMfPrev,beta,...
                alpha,k0,kLin,k1,r1,sigma1,psi1,psi2s2,mu,gamma,b1,c,ageLev,...
                k2,gam2,tau,immC,slopeC,del,ageMthMax,da,demog,...
                bCulex,toleranceX,k0_An,k0_Cx,immCMin,larvaeDeath,PP,excessMortality),...
                [0 max_bit],options);
            
            VoverH(i,:) = bitx/ParamVectors0(1); % VoverH = MBR/beta;
            ParamVectors(:,i) = ParamVectors0;
        end
    else
        % generate an initial sample of parameter vectors using Latin Hypercube
        % Sampling
        ParamVectors = ParameterVectors_LHSDesign_based(NParamVecs);
        [~,n] = size(ABR); % n is the number of species
        VoverH = (ABR/12)./repmat(ParamVectors(1,:)',[1 n]); % VoverH = MBR/beta
    end
    
    % Running to equilibrium
    [L3Values,mfPrevArray] = ...
        Calculate_EndemicEquil_L3_Mf(NParamVecs,...
        ageMthMax,bCulex,demog,ParamVectors,VoverH,da,toleranceX);
    
    % Resampling 
    % Resample the parameter vectors using 1) SIR algorithm and 2) pass/fail criteria
    if AgeMf_asInput % Mf age profile
        
        [kId1, ~] = get_modelFits_Mf_AgeData(MfData,...
            mfPrevArray,demog,da,ageMthMax,SIR_samples);
        kId = kId1; 
        
    elseif OverallMf_asInput % Overall Mf
        
        [kId1, ~] = get_modelFits_Mf_OverallData(MfData,...
            mfPrevArray,demog,da,ageMthMax,SIR_samples);
        kId = kId1;
        
    elseif OverallMfRange_asInput % No Baseline Mf
        
        OverallMf = zeros(length(mfPrevArray(1,:)),1);
        for i = 1:length(mfPrevArray(1,:))  %% loop through 200,000 parameters
            OverallMf(i) = getOverallPrev_fromAgeCurves(100*mfPrevArray(:,i),...
                demog,da,ageMthMax);
        end
        kId = find(OverallMf > minPrev & OverallMf < maxPrev);
        
    end
    
    % Store parameters, equilibrium L3, and an equilibrium mf values
    % corresponding to the resampled parameter vectors
    ParameterVectors(:,NumParam+1:NumParam+length(kId)) = ...
        ParamVectors(:,kId);
    
    L3Values1(NumParam+1:NumParam+length(kId),:) = ...
        L3Values(kId,:);
    
    mfPrevArray1(:,NumParam+1:NumParam+length(kId)) = ...
        mfPrevArray(:,kId);
    
    ABR1(NumParam+1:NumParam+length(kId),:) = ...
        12*VoverH(kId,:).*ParamVectors(1,kId)';
    
    % Update number of parameter vectors selected
    NumParam = NumParam + length(kId);
    fprintf(1,'Num selected kIds = %d, Total Num Params = %d\n',...
        length(kId),NumParam);
end

L3Values = L3Values1;
mfPrevArray = mfPrevArray1;
ABR = ABR1;

% save sampled best fitting outputs and basic descriptive parameters of
% the site
save(sprintf('../IO/OUT/ParamVectors%s.mat',Site,int2str(estABRflag)),'bCulex',...
    'demog','ABR','ParameterVectors','L3Values','ageMthMax',...
    'mfPrevArray');

end
