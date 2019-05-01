% Function: get_modelFits_Mf_OVerallData
%  
% Description: convert overall mf prevalence into expected age profiles
% patterns and calculate likelihoods using this theoretical data
%
% Inputs:
%   MfData0,mfPrevArray,demog,da,ageMthMax,SIR_samples
%
% Outputs:
%   kId1, kId2
% ________________________________________

function [kId1, kId2] = get_modelFits_Mf_OverallData(MfData0,...
    mfPrevArray,demog,da,ageMthMax,SIR_samples)

% Initialize the index and likelihood arrays
kId1 = []; % for indices of SIR resampled fits
kId2 = []; % for indices of pass/fail criteria fits

% set midage groups for constructed age profile
MidAge = [5 14.5 24.5 34.5 44.5 54.5 64.5];
OverallMfPrev = MfData0(:,3)/MfData0(:,2);
TotalMfSamples = MfData0(:,2);

% SIR_samples divided by how many types of age profiles
SIR_samples = round(SIR_samples/2);

for icurve = 1:2 % loop over three types of age curves
    
    LikArray1 = zeros(length(mfPrevArray(1,:)),1); % for likelihoods of SIR resampled fits
    LikArray2 = zeros(length(mfPrevArray(1,:)),1); % for likelihoods of pass/fail criteria fits
    
    % calculate theoretical age profile data
    MfData = getMfAgeProfile_fromOnchoCurves(TotalMfSamples,icurve,OverallMfPrev,ageMthMax/12,demog,MidAge);
    v=genvarname(sprintf('MfData%s',int2str(icurve)));
    eval([v '= MfData;']);
    
    % calculate 95% CI upper and lower bounds of each mf data point
    MfBounds = get_the95LU_bounds_agedata(MfData);
    
    % calculate likelihood of each sampled parameter vector using SIR and
    % pass/fail methods
    parfor i=1:length(mfPrevArray(1,:))
        LikArray1(i) = calculateLikelihoods(MfData,100*mfPrevArray(:,i),demog,da,ageMthMax);
        LikArray2(i) = likelihood_using_ranks_by_passFail_Mf_Only(mfPrevArray(:,i),MfBounds);
    end
    
    % resample by SIR algorithm, aggregate kIds over all age curves
    Indx1 = ReSampleMfPrevalence(LikArray1,SIR_samples);
    Indx = find( Indx1 > 0);
    if ~isempty(Indx)
        kId1 = [kId1; Indx1(Indx)];
    end
    
    % resample by pass/fail criteria, aggregate kIds over all age curves
    Indx2 = ReSampleMfPrevalence_passFail_Mf_only(LikArray2,SIR_samples);
    Indx = find( Indx2 > 0);
    if ~isempty(Indx)
        kId2 = [kId2; Indx2(Indx)];
    end
end

end
% #########################################################################
