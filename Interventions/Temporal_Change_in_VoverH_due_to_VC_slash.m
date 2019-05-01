function VH = Temporal_Change_in_VoverH_due_to_VC_slash(VC0,...
    NumYears,beta,P,R0,slash_t,MultiVecMBR)

MBRpeak = MultiVecMBR;
RL = P(1);
RU = P(2);
k1 = P(3);
k2 = P(4);
n = P(5);
A = P(6);

% rain described by cosinor function, loaded from MBRFunParams.mat
R = zeros(NumYears*12+1,1);
R(1:NumYears*12,1) = repmat(R0,NumYears,1);
R(end,1) = R0(1);

% MBR for given rainfall values
MBR_0 = zeros(length(R),1); % background MBR expected
for t = 1:length(R)
    MBR_0(t) = MBRpeak*(1-exp(-(R(t)/RL).^k1)).*exp(-(R(t)/RU).^k2); % MBR function of rain
    if t == 12
        while sum(MBR_0(1:t)) < MultiVecMBR*12
            MBRpeak = MBRpeak*1.01;
            for tt = 1:12
                MBR_0(tt) = MBRpeak*(1-exp(-(R(tt)/RL).^k1)).*exp(-(R(tt)/RU).^k2);
            end
        end
    end
end

% slash is implemented
if VC0 == 1
    MBR = zeros(length(R),1); % MBR given slash
    i = 1;
    for t = 1:length(R) % time since slash
        s = slash_t(i);
        ts = t-s;
        if ts < 1
            MBR(t) = MBR_0(t);
        else
            MBR(t) = MBR_0(t)*n*(1-exp(-A*ts)); % slash efficacy function 
        end
        if i < length(slash_t)
            if t == slash_t(i+1)
                i = i+1;
            end
        end
    end
else
    MBR = MBR_0;
end

VH = MBR/beta;

end