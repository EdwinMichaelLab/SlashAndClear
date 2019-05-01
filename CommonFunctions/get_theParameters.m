function [beta,HBIndex,GTCP,alpha,k0,kLin,k1,r1,sigma1,psi1,psi2s2,mu,gamma,b1,c,...
    ageLev,k2,gam2,tau,immC,slopeC,k0_An,k0_Cx,immCMin,larvaeDeath,...
    PP,del,excessMortality] = get_theParameters(ParameterVectors,i)

beta    = ParameterVectors(1,i);
HBIndex = ParameterVectors(2,i);
GTCP    = ParameterVectors(3,i);
beta    = HBIndex/GTCP;
alpha   = ParameterVectors(4,i);
k0  = ParameterVectors(5,i);
kLin   = ParameterVectors(6,i);
k1      = ParameterVectors(7,i);
r1      = ParameterVectors(8,i);
sigma1  = ParameterVectors(9,i);
psi1    = ParameterVectors(10,i);
psi2s2  = ParameterVectors(11,i);
mu      = ParameterVectors(12,i);
gamma   = ParameterVectors(13,i);
b1      = ParameterVectors(14,i);
c       = ParameterVectors(15,i);
ageLev  = ParameterVectors(16,i);
k2      = ParameterVectors(17,i);
gam2    = ParameterVectors(18,i);
tau     = ParameterVectors(19,i);
immC    = ParameterVectors(20,i);
slopeC  = ParameterVectors(21,i);
k0_An   = ParameterVectors(22,i);
k0_Cx   = ParameterVectors(23,i);
immCMin = ParameterVectors(24,i);
larvaeDeath = ParameterVectors(25,i);
PP      = ParameterVectors(26,i);
del     = ParameterVectors(27,i);
excessMortality = ParameterVectors(28,i);

end