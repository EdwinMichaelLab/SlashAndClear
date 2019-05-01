% Function to calculate the equilibrium values of worm-load mf-intensity
% and acquired immunity over age GIVEN the Equilibrium value of l3
% Changed on 21 April 2015
function [pVec,wVec,mVec,iVec]=OnchoStateVarsEquilValues(beta,...
    k0,kLin,psi1,mu,alpha,gamma,c,ageLev,l3VoverH,psi2s2,immC,slopeC,...
    del,ageMthMax,da,immCMin,wVec,mVec,iVec,pVec,gVec,PP)

SexRatio = 0.5;

pVec=pVec*0;
wVec=wVec*0;
mVec=mVec*0;
iVec=iVec*0;
gVec=gVec*0;

p1=0;
w1=0;
m1=0;
i1=0;
for a2=2:da:PP
    foi=min(a2/(12.0*ageLev),1);
    immSupp = (immCMin+immC*slopeC*(w1+p1))/(1+slopeC*(w1+p1));
    immFun = (1/(1+c*i1));
    gVec(a2)=foi*beta*psi1*psi2s2*immFun*immSupp*l3VoverH;
    Pworm=0; % Patent worm zero
    pVec(a2)=p1+(gVec(a2)-mu*p1-Pworm)*da;
    wVec(a2)=w1+(Pworm-mu*w1)*da;
    mVec(a2)=m1+(SexRatio*alpha*wormMatingProb(w1,...
        negBinshapeFun(w1,k0,kLin))*w1-gamma*m1)*da;
    iVec(a2)=i1+(p1+w1-del*i1)*da;
    
    p1 = pVec(a2);
    w1 = wVec(a2);
    m1 = mVec(a2);
    i1 = iVec(a2);
end
for a2=PP+1:da:ageMthMax
    foi=min(a2/(12.0*ageLev),1);
    immSupp = (immCMin+immC*slopeC*(w1+p1))/(1+slopeC*(w1+p1));
    immFun = (1/(1+c*i1));
    gVec(a2)=foi*beta*psi1*psi2s2*immFun*immSupp*l3VoverH;
    Pworm=gVec(a2-PP)*exp(-mu*PP); % Patent worm
    pVec(a2)=p1+(gVec(a2)-mu*p1-Pworm)*da;
    wVec(a2)=w1+(Pworm-mu*w1)*da;
    mVec(a2)=m1+(SexRatio*alpha*wormMatingProb(w1,...
        negBinshapeFun(w1,k0,kLin))*w1-gamma*m1)*da;
    iVec(a2)=i1+(p1+w1-del*i1)*da;
    
    p1 = pVec(a2);
    w1 = wVec(a2);
    m1 = mVec(a2);
    i1 = iVec(a2);
end

if (isreal(wVec)==0 || isreal(mVec)==0 || isreal(iVec)==0)
    pVec=pVec*0;
    wVec=wVec*0;
    mVec=mVec*0;
    iVec=iVec*0;
    gVec=gVec*0;
end

end


