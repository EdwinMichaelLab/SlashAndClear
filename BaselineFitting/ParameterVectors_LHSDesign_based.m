% ParameterVectors_LHSDesign_based randomly samples 
% a set of NParamVecs Parameter vectors from their ranges

function [ ParameterVectors ] = ParameterVectors_LHSDesign_based( NParamVecs )

[rangeParamVals,~,minParamVals] = Range_of_Parameters;

ParameterVectors = zeros(length(minParamVals),NParamVecs);

for i=1:length(minParamVals)
    ParameterVectors(i,:) = minParamVals(i) + ...
        rangeParamVals(i)*lhsdesign(1,NParamVecs,'criterion','correlation');
end

ParameterVectors(16,:) = floor(ParameterVectors(16,:)+0.5); % ageLev should be non-fraction
ParameterVectors(26,:) = floor(ParameterVectors(26,:)+0.5); % PP should be non-fraction

end

