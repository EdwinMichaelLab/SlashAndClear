Reg_min = [0.05 0.95 9];
Reg_max = [0.3 0.99 11];

RegimenEfficacy = Reg_min + (Reg_max-Reg_min).*lhsdesign(100,length(Reg_min));
RegimenEfficacy(:,3) = floor(RegimenEfficacy(:,3));

VillageDrugEfficacy = RegimenEfficacy;

save('RegEfficacy.mat','VillageDrugEfficacy');