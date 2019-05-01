function [mfPrevIntv,MBRIntv,L3Intv,pIntv,wIntv,MonthsMDA,MonthlyMDACov,MfInt] = Oncho_Intervention_MDA_VC(...
            ParameterVectors,L3Values,demog,ageMthMax,bCulex,...
            AgeLimits,MonthlyMDACov,RegimenEfficacy,...
            MDAInterval,NumYears,MultiVecMBR,VC,slash_t,EP_Th,TBR)

da = 1;
k = 1:length(L3Values);

[mfPrevIntv,MBRIntv,L3Intv,pIntv,wIntv,MfInt] = RunOnchoInterventionScenarios_slash(k,...
    ParameterVectors,L3Values,demog,ageMthMax,da,bCulex,...
    AgeLimits,MonthlyMDACov,RegimenEfficacy,...
    MDAInterval,NumYears,MultiVecMBR,VC,slash_t,TBR);

[MonthsMDA,~] = Time_toCross_below_WHO_Threshold(mfPrevIntv,EP_Th);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [timeX,tThresh] = Time_toCross_below_WHO_Threshold(mfPrevIntv,EP_Th)
[m,n] = size(mfPrevIntv);
timeX = [];
for i = 1:n
    id = find(mfPrevIntv(1:12:m,i) < EP_Th);
    if isempty(id)
        timeX = [timeX; NaN];
    else
        timeX = [timeX; 12*(id(1)-1)];
    end
end

tThresh = 1 + prctile(timeX,95); % 95% passes through threshold

end
