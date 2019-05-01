% #########################################################################
function [kId1, kId2] = get_modelFits_Mf_AgeData(MfData,...
    mfPrevArray,demog,da,ageMthMax,SIR_samples)

% Initialize the index and likelihood arrays
kId1 = []; % for indices of SIR resampled fits
kId2 = []; % for indices of pass/fail criteria fits
LikArray1 = ones(length(mfPrevArray(1,:)),1); % for likelihoods of SIR resampled fits
LikArray2 = zeros(length(mfPrevArray(1,:)),1); % for likelihoods of pass/fail criteria fits

% calculate 95% CI upper and lower bounds of each mf data point
MfBounds = get_the95LU_bounds_agedata(MfData);

% calculate likelihood of each sampled parameter vector using SIR and
% pass/fail methods
parfor i=1:length(mfPrevArray(1,:))
    LikArray1(i)= calculateLikelihoods(MfData,100*mfPrevArray(:,i),demog,da,ageMthMax);
    LikArray2(i)= likelihood_using_ranks_by_passFail_Mf_Only(mfPrevArray(:,i),MfBounds);
end

% resample by SIR algorithm
Indx1 = ReSampleMfPrevalence(LikArray1,SIR_samples);
Indx = find( Indx1 > 0);
if ~isempty(Indx)
    kId1 = [kId1; Indx1(Indx)];
end

% resample by pass/fail criteria
Indx2 = ReSampleMfPrevalence_passFail_Mf_only(LikArray2,SIR_samples);
Indx = find( Indx2 > 0);
if ~isempty(Indx)
    kId2 = [kId2; Indx2(Indx)];
end

end
% #########################################################################
