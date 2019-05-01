function [MfPrevStore,L3Store,pArray,wArray,MfInt ] = ModelOnchoIntv(beta,...
    alpha,k0,kLin,k1,r1,sigma1,psi1,psi2s2,mu,gamma,b1,c,ageLev,k2,...
    gam2,tau,immC,immCMin,slopeC,PP,VoverH,del,ageMthMax,da,...
    bCulex,l3,larvaeDeath,excessMortality,gVec,pVec,wVec,mVec,iVec,demoX,...
    AgeLimits,drugCoverage,drug,...
    MDAInterval,NumYears,k0Spp,TBR)

SexRatio = 0.5;

p2 = pVec;
w2 = wVec;
m2 = mVec;
i2 = iVec;

MfPrevStore = zeros(12*NumYears+1,1,'single');
MfInt = zeros(12*NumYears+1,1,'single');
L3Store = zeros(12*NumYears+1,1,'single');
pArray = zeros(12*NumYears+1,1,'single');
wArray = zeros(12*NumYears+1,1,'single');

sumMf = sumMfPrev(m2,da,demoX);
MfPrevStore(1) = 100*mfAgeprevFun(sumMf,negBinshapeFun(sumMf,k0,kLin));
MfInt(1) = sumMf;
L3Store(1) = l3;

sumP = sumMfPrev(p2,da,demoX);
pArray(1) = 100*mfAgeprevFun(sumP,negBinshapeFun(sumP,k0,kLin));
sumW = sumMfPrev(w2,da,demoX);
wArray(1) = 100*mfAgeprevFun(sumW,negBinshapeFun(sumW,k0,kLin)); 

sinceTreat = 0;

for iTime = 1:da:NumYears*12
    
    % Initialize arrays
    p1 = 0;
    w1 = 0;
    m1 = 0;
    i1 = 0;
    
    pVec(1) = p1;
    wVec(1) = w1;
    mVec(1) = m1;
    iVec(1) = i1;
    
    l3VoverH = sum(VoverH(iTime,:).*l3);
    
    % children younger the the pre-patency period of the worms cannot
    % have any patent worms
    for a2 = 2:da:PP
        foi = min(a2/(12.0*ageLev),1);
        immSupp = (immCMin+immC*slopeC*(w1+p1))/(1+slopeC*(w1+p1));
        immFun = (1/(1+c*i1));
        gVec(a2,iTime) = foi*beta*psi1*psi2s2*immFun*immSupp*l3VoverH;
        Pworm = 0; % Patent worm
        pVec(a2) = max(p1 + (gVec(a2,iTime)-mu*p1),0)*da;
        wVec(a2) = w1 + (Pworm-mu*w1)*da;
        mVec(a2) = m1 + (SexRatio*alpha*wormMatingProb(w1,negBinshapeFun(w1,k0,kLin))*w1-gamma*m1)*da;
        iVec(a2) = i1 + (p1+w1-del*i1)*da;
        
        p1 = p2(a2);
        w1 = w2(a2);
        m1 = m2(a2);
        i1 = i2(a2);
    end
    
    % The main loop calculating the equilibrium (time-independent)
    % values of the arrays wVec, mVec and iVec for individuals older than
    % the worm pre-patency period starts below
    
    % 'If' option 1: Check whether it is time to administer chemo
    if (mod(iTime,MDAInterval(iTime)) == 1)
        sinceTreat = 1;
        for a2 = PP+1:da:ageMthMax
            foi = min(a2/(12.0*ageLev),1);
            gVec(a2,iTime) = foi*beta*psi1*psi2s2*immFun*immSupp*l3VoverH;
            if (a2 >= AgeLimits(1)*12 && a2 < AgeLimits(2)*12)
                gVec(a2,iTime) = (1-drug(1)*drugCoverage(iTime))*gVec(a2,iTime);
                pVec(a2) = (1-drug(1)*drugCoverage(iTime))*p1;
                wVec(a2) = (1-drug(1)*drugCoverage(iTime))*w1;
                mVec(a2) = (1-drug(2)*drugCoverage(iTime))*m1;
            else
                immSupp  = (immCMin+immC*slopeC*(w1+p1))/(1+slopeC*(w1+p1));
                immFun   = (1/(1+c*i1));
                Pworm    = gVec(a2-PP,max(iTime-PP,1))*exp(-mu*PP);
                pVec(a2) = max(p1+(gVec(a2,iTime)-mu*p1-gVec(a2-PP,max(iTime-PP,1))*exp(-mu*PP)),0)*da;
                wVec(a2) = w1+(Pworm-mu*w1)*da;
                mVec(a2) = m1+(SexRatio*alpha*wormMatingProb(w1,...
                    negBinshapeFun(w1,k0,kLin))*w1-gamma*m1)*da;
            end
            iVec(a2) = i1 + ((w1+p1) - del*i1)*da;
            p1 = p2(a2);
            w1 = w2(a2);
            m1 = m2(a2);
            i1 = i2(a2);
        end
        
    % 'If' option 2: If the current time is after one of the treatments but not before the
    % the drug has worn off, then the proportion of the
    % population that was treated will not produce Mf.

    % This treated proportion will have aged by a time that is
    % defined here as 'sinceTreat' below.
        
    elseif (sinceTreat < drug(3))
        sinceTreat = sinceTreat+1;
        for a2 = PP+1:da:ageMthMax
            foi = min(a2/(12.0*ageLev),1);
            immSupp = (immCMin+immC*slopeC*(w1+p1))/(1+slopeC*(w1+p1));
            immFun = (1/(1+c*i1));
            gVec(a2,iTime) = foi*beta*psi1*psi2s2*immFun*immSupp*l3VoverH;
            Pworm = gVec(a2-PP,max(iTime-PP,1))*exp(-mu*PP);
            pVec(a2) = max(p1+(gVec(a2,iTime)-mu*p1-gVec(a2-PP,max(iTime-PP,1))*exp(-mu*PP)),0)*da;
            wVec(a2) = w1+(Pworm-mu*w1)*da;
            if(a2 >= AgeLimits(1)*12+sinceTreat && a2 < ...
                    AgeLimits(2)*12+sinceTreat)
                % Only those untreated will produce Mf, hence a fraction
                % w*(1-coverage)
                mVec(a2) = m1+...
                    (SexRatio*alpha*wormMatingProb(w1,negBinshapeFun(w1,k0,kLin))*...
                    (1-drug(2)*drugCoverage(iTime))*w1-gamma*m1)*da;
            else
                mVec(a2) = m1+...
                    (SexRatio*alpha*wormMatingProb(w1,negBinshapeFun(w1,k0,kLin))*...
                    w1-gamma*m1)*da;
            end
            
            iVec(a2) = i1 + ((w1+p1) - del*i1)*da;
            
            p1 = p2(a2);
            w1 = w2(a2);
            m1 = m2(a2);
            i1 = i2(a2);
        end
        
    % 'If' option 3: If the current time doesn't lie in a
    % treatment/post-treatment-pre-waning period, just update the
    % variables normally
    
    else
        sinceTreat = 0;
        for a2 = PP+1:da:ageMthMax
            foi = min(a2/(12.0*ageLev),1);
            immSupp = (immCMin+immC*slopeC*(w1+p1))/(1+slopeC*(w1+p1));
            immFun = (1/(1+c*i1));
            gVec(a2,iTime) = foi*beta*psi1*psi2s2*immFun*immSupp*l3VoverH;
            Pworm = gVec(a2-PP,max(iTime-PP,1))*exp(-mu*PP);
            pVec(a2) = max(p1+(gVec(a2,iTime)-mu*p1-gVec(a2-PP,max(iTime-PP,1))*exp(-mu*PP)),0)*da;
            wVec(a2) = w1+(Pworm-mu*w1)*da;
            mVec(a2) = m1+(SexRatio*alpha*wormMatingProb(w1,...
                negBinshapeFun(w1,k0,kLin))*w1-gamma*m1)*da;
            iVec(a2) = i1 + ((w1+p1) - del*i1)*da;
            
            p1 = p2(a2);
            w1 = w2(a2);
            m1 = m2(a2);
            i1 = i2(a2);
            
        end
    end % end of 'if' block within each time step
    
    % recalculate the new value of l3
    % set l3 to zero if below threshold biting rate, no more transmission
    if VoverH(iTime+1,:) < TBR/12/beta
        l3 = 0;
    else
        for iSpp = 1:length(k0Spp)
            l3(iSpp) = calcL3Fun_Oncho_exmort(mVec,beta,k0Spp(iSpp),kLin,k1,r1,sigma1,psi1,...
                b1,k2,gam2,tau,demoX,ageMthMax,da,ageLev,bCulex(iSpp),larvaeDeath,excessMortality);
        end
    end
    
    % store current Mf and L3 values
    sumMf = sumMfPrev(mVec,da,demoX);
    MfPrevStore(iTime+1) = 100*mfAgeprevFun(sumMf,negBinshapeFun(sumMf,k0,kLin));
    MfInt(iTime+1) = sumMf;
    
    L3Store(iTime+1) = l3;

    sumP = sumMfPrev(pVec,da,demoX);
    pArray(iTime+1) = 100*mfAgeprevFun(sumP,negBinshapeFun(sumP,k0,kLin));
    sumW = sumMfPrev(wVec,da,demoX);
    wArray(iTime+1) = 100*mfAgeprevFun(sumW,negBinshapeFun(sumW,k0,kLin));
    
    % Set the dummy variables to store the current values so
    % that they can be used when the w,m,i arrays are updated at
    % the next timestep.
    
    p2 = pVec;
    w2 = wVec;
    m2 = mVec;
    i2 = iVec;
end
end
