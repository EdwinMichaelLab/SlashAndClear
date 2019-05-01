function maxTimeMonth = Time_toCross_below_Threshold(mfPrevIntv,EP_Th)
[m,n] = size(mfPrevIntv);
timeX = [];
for i = 1:n
    id = find(mfPrevIntv(1:12:m,i) < EP_Th);
    if isempty(id)
        timeX = [timeX; length(mfPrevIntv(:,1))];
    else
        timeX = [timeX; 12*(id(1)-1)];
    end
end

maxTimeMonth = prctile(timeX,95)/12; % 95% passes through threshold
% maxTimeMonth = mean(timeX)/12;

end
