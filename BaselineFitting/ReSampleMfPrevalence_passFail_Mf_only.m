% Resample based on pass/fail criteria to compensate for situations where
% likelihood-based method might fail 

function likArrayId=ReSampleMfPrevalence_passFail_Mf_only(likelihoodArray, NumRequired)

[m,n] = size(likelihoodArray);
likArrayId = [];
Num0 = 0;
iNum = 0;
likVals = unique(likelihoodArray(:,1));
likVals = likVals(end:-1:1);
[m1,n1] = size(likVals);

while Num0 < NumRequired && iNum < m1
    if iNum == 0
        id1 = find(likelihoodArray(:,1) >= likVals(iNum + 1));
    else
        id1 = find(likelihoodArray(:,1) >= likVals(iNum + 1) ...
            & likelihoodArray(:,1) < likVals(iNum));
    end
    if ~isempty(id1)
        likArrayId = [likArrayId; id1];
        Num0 = Num0 + length(id1);
    end
    iNum = iNum + 1;
end
end
