function [mfAgeProfile] = getMfAgeProfile_fromOnchoCurves(TotalMfSamples,...
    curve,OverallMfPrev,ageMax,demog,MidAge)
% get the age-Group proportion from the demographics based calculation
AgeGroupProp = DemographicsBased_ageGroupProportions(demog,ageMax,MidAge);

%The first step is to evaluate the selected age curve theoretically for the
%country and TotalMfSamples given by the data

%Curve equations obtained from nlme fits to several data sets of the
%same curve type
MfPrevTheor = zeros(ageMax*12,1);

if curve == 1 %plateau
    for age = 0:ageMax*12
        MfPrevTheor(age+1,:) = (0.29*(age/12))/(1+((age/12)/3.29));
    end
else  %convex
    for age = 0:ageMax*12
        MfPrevTheor(age+1,:) = 0.08*(age/12)*exp(-0.03*(age/12));
    end
end

NoSampledTheor = round(TotalMfSamples.*AgeGroupProp);

NoPosTheor = zeros(length(MidAge),1);

Ages = 0:ageMax*12;
j = 1;
minAge = MidAge(j) - (MidAge(j+1)-MidAge(j))/2;
maxAge = MidAge(j) + (MidAge(j+1)-MidAge(j))/2;
id = find(minAge*12 <= Ages & Ages < maxAge*12);
NoPosTheor(j,:) = median(MfPrevTheor(id))*NoSampledTheor(j,:);
diffMidAge = (MidAge(3)-MidAge(2))/2;
for j = 2:length(MidAge)
    minAge = MidAge(j) - diffMidAge;
    maxAge = MidAge(j) + diffMidAge;
    id = find(minAge*12 <= Ages & Ages < maxAge*12);
    NoPosTheor(j,:) = median(MfPrevTheor(id))*NoSampledTheor(j,:);
end

OverallMfTheor = sum(NoPosTheor)/sum(NoSampledTheor);

%The second step is to construct the village age profile in proportion to
%the theoretical age profile constructed above
mfAgeProfile=zeros(length(MidAge),4);
mfAgeProfile(:,1)=MidAge;
mfAgeProfile(:,2)=NoSampledTheor;
mfAgeProfile(:,3)=min(round((OverallMfPrev/OverallMfTheor).*NoPosTheor),round(TotalMfSamples*AgeGroupProp));
mfAgeProfile(:,4)=[9; 19; 29; 39; 49; 59; 69];

end

