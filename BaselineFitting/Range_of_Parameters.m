% min and max values for uniform prior distributions of model parameters 
function [rangeParamVals,maxParamVals,minParamVals] = Range_of_Parameters

%'beta';'HBIndex';'GTCP';'alpha';'k0';'kLin';'k1';'r1';'sigma1';'psi1';'psi2s2';
% 'mu';'gamma';'b1';'c';'ageLev';'k2';'gam2';'tau';'immC';'slopeC';'k0_An';'k0_Cx';'immCMin';
%'larvaeDeath';'PP';'del';excessMortality

minParamVals=[05;0.30;0.067;0.25;0.00036;0.000001;1.26;0.01471;1.5;0.124215;0.00004;0.0069;0.08;0.259;0.0001; 1;1.161;0.03687;0;0.5;0.10;0.00036;0.00036;0.0025;0.3333;9;0.00001;0.0208];
maxParamVals=[15;0.99;0.130;1.50;0.00440;0.022880;1.97;0.02253;8.5;0.703797;0.00400;0.0104;0.12;0.481;0.0010;25;1.537;0.04271;0;5.0;0.75;0.00440;0.00440;1.0000;1.1583;12;0.00010;0.0529];

rangeParamVals=maxParamVals-minParamVals;

end
