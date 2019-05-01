function [xdemoA, xdemoB] = get_demoA_and_demoB(Country,Country_demo,...
    demoA,demoB)
for j=1:length(Country_demo)
    if strcmp(Country,char(Country_demo(j)))
        xdemoA=demoA(j);
        xdemoB=demoB(j);
        break;
    end
end
% fprintf(1,'%s,%f,%f\n',char(Country_demo(j)),xdemoA,xdemoB);
% pause;
end