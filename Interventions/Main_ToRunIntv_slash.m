clear all; clc;

%% add CommonFunctions folder containing basic model functions and Baseline folder containing baseline fits to path
addpath('../CommonFunctions/');
addpath('../BaselineFitting/');

%% Load user defined variables within IO directory
load('../IO/IN/Inter_IN.mat');

% integration time step
da = 1; % 1 month

%% Intervention specifications

% Mass drug administration parameters
AgeLimits = [5 100]; % treatment for > 5 yrs old

load('RegEfficacy.mat'); % IVM worm kill rate, mf kill rate, sterilization period (months)
RegimenEfficacy = VillageDrugEfficacy;

setupParallelPool(); % Initialize the parallel workers

%% Run fitting procedure for each site

for iSites = 1:length(Sites)
    %% set up inputs
    
    % define data
    try
        MfData = eval(sprintf('%sMf',Sites{iSites}));
        NumYears = eval(sprintf('%sNumYears',Sites{iSites}));
        VC = eval(sprintf('%sVC',Sites{iSites}));
        MDAFreq = eval(sprintf('%sMDAFreq',Sites{iSites}));
        MDACov = eval(sprintf('%sMDACov',Sites{iSites}));
        slash_t = eval(sprintf('%sSlash_t',Sites{iSites}));
    catch ME
        % If the file is missing, complain but continue.
        if (strcmp(ME.identifier,'MATLAB:UndefinedFunction'))
            fprintf('\nWarning: Cannot find variables for %s\n', Sites{iSites});
            fprintf('Check PostIntv_data.m. Continuing...\n');
            continue;
            % Other error, throw it and die
        else
            rethrow(ME)
        end
    end
    
    %% load baseline fits
    % ABR, L3Values, ParameterVectors, ageMthMax, bCulex, demog, mfPrevArray are loaded.
    loadFile = sprintf('../IO/OUT/ParamVectors%s.mat',Sites{iSites});
    try
        load(loadFile);
    catch ME
        % If file not found, try the next one
        if (strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile'))
            fprintf('\nWarning: %s not found. Trying next baseline file...\n', loadFile);
            continue;
        else
            rethrow(ME) % If its not file not found, throw the error.
        end
    end
    
    MultiVecMBR = ABR/12;
    kId = 1:length(L3Values);
       
    %% model interventions
    
    loadFile = sprintf('../IO/OUT/Bpts_%s_ABR1.mat',Sites{iSites});
    try
        load(loadFile,'bRate');
    catch ME
        % If file not found, try the next one
        if (strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile'))
            fprintf('\nWarning: %s not found. Trying next baseline file...\n', loadFile);
            continue;
        else
            rethrow(ME) % If its not file not found, throw the error.
        end
    end
    TBR = bRate(:,1);
    
    % MDA Interval
    MDAInterval = zeros(1,NumYears*12);
    MDAInterval(1,1:end) = MDAFreq(1);
    
    % MDA Coverage
    for i = 1:length(MDACov)
        MonthlyMDACov = ones(NumYears*12,1)*MDACov(i);
    end
    
    % slash and clear intervention schedules
    for v = 1:4
        if v == 1
            slash_t = [1:12:50*12]; %annually month before peak observed biting (i.e. Apr slash, May peak biting)
            VC = 1;
        elseif v == 2
            slash_t = sort([[1:12:50*12],[3:12:50*12],[5:12:50*12],[7:12:50*12]]); % every other month during rainy (Apr, Jun, Aug, Oct)
            VC = 1;
        elseif v == 3
            slash_t = [1:1:50*12]; % monthly
            VC = 1;
        else
            VC = 0; % no slash and clear
        end
        
        % target threshold
        EP_Th = 1;
        
        [mfPrevIntv,MBRIntv,L3Intv,pIntv,wIntv,MonthsMDA,MonthlyMDACov,MfInt] ...
            = Oncho_Intervention_MDA_VC(...
            ParameterVectors,L3Values,demog,ageMthMax,bCulex,...
            AgeLimits,MonthlyMDACov,RegimenEfficacy,...
            MDAInterval,NumYears,MultiVecMBR,VC,slash_t,EP_Th,TBR);   
        
        save(sprintf('../IO/OUT/Intv%s_v%s.mat',char(Sites{iSites}),int2str(v)),...
            'mfPrevIntv','MBRIntv','L3Intv','pIntv','wIntv','RegimenEfficacy','MonthlyMDACov',...
            'MfData','EP_Th','MonthsMDA','slash_t','TBR');
    end
    
    
end


