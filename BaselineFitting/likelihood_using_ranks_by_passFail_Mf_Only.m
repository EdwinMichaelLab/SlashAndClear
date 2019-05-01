function likelihood = likelihood_using_ranks_by_passFail_Mf_Only(mfPrev,Bounds)

likelihood = 0;

for j = 1:length(Bounds(:,1)) % loop through data of all age groups 
    
    % check if model-generated mf falls in 95% CI bounds of data
    id = find(100*mfPrev(Bounds(j,1)*12,1) >= Bounds(j,2) & ...
        100*mfPrev(Bounds(j,1)*12,1) <= Bounds(j,3));
    
    % count how many age groups "pass"
    if ~isempty(id)
        likelihood = likelihood + 1;
    end
end

end