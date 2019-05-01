function mfPrevGroup = getOverallPrev_fromAgeCurves(mPrevVal,...
    demographic,da,ageMthMax)

warning('off','all');
ageGroup = ageMthMax/12;

mfPrevGroup = zeros(length(ageGroup),1);
iGroup = 1;
weightedMfPrev = 0;
numPeopleInGroup = 0;

iCheck =1;

for iAge = 1:da:ageMthMax
    
    if iGroup <= length(ageGroup)
        if iAge/12.0 < ageGroup(iGroup)
            pia=pi_PeopleFun(iAge/12.0,demographic);
            weightedMfPrev = weightedMfPrev + pia*mPrevVal(iCheck);
            numPeopleInGroup = numPeopleInGroup + pia;
        else
            mfPrevGroup(iGroup) = weightedMfPrev/(numPeopleInGroup*100);
            iGroup = iGroup + 1;
            weightedMfPrev = 0;
            numPeopleInGroup = 0;
        end
    end
    
    iCheck = iCheck + 1;
    
end

if iGroup == length(ageGroup)
    mfPrevGroup(iGroup) = weightedMfPrev/(numPeopleInGroup*100);
end

tol = 0.000001;

for iGroup = 1:length(ageGroup)
    if abs(mfPrevGroup(iGroup)) < tol
        mfPrevGroup(iGroup) = 0.0000001;
    elseif abs(1-mfPrevGroup(iGroup)) < tol
        mfPrevGroup(iGroup) = 0.999999;
    end
end
end


