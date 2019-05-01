% calculate the upper and lower prevalence values of 95% confidence interval 
function Bounds = get_the95LU_bounds_agedata(data)

aa=ConstructBinomialErrorBars(data);
Bounds = [aa(:,1), aa(:,2)-aa(:,3), aa(:,2)+aa(:,4)];

end
