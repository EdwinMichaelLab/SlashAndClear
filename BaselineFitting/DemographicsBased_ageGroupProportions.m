function AgeGroupProp = DemographicsBased_ageGroupProportions(demog,...
    ageMax,MidAge)

pia = zeros(1+ageMax*12,1);
    for age = 0:ageMax*12
        pia(1+age,1) = demog(3)*exp(-demog(4)*age/12);
    end

Ages = 0:ageMax*12;

AgeGroupProp = zeros(length(MidAge),1);
j = 1;
minAge = MidAge(j) - (MidAge(j+1)-MidAge(j))/2;
maxAge = MidAge(j) + (MidAge(j+1)-MidAge(j))/2;
id = find(minAge*12 <= Ages & Ages < maxAge*12);
AgeGroupProp(j,1) = sum(pia(id))/sum(pia);

diffMidAge = (MidAge(3)-MidAge(2))/2;
for j = 2:length(MidAge)
    minAge = MidAge(j) - diffMidAge;
    maxAge = MidAge(j) + diffMidAge;
    id = find(minAge*12 <= Ages & Ages < maxAge*12);
    AgeGroupProp(j,1) = sum(pia(id))/sum(pia);
end
end
