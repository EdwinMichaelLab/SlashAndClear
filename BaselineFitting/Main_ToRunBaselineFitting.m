clear all; clc;
%% add CommonFunctions folder containing basic model functions to path
addpath('../CommonFunctions/');

%% User-defined inputs

% Parameters related to sampling and resampling
NParamVecs = 200000; % default: 200000
SIR_samples = 500; % default: 500

% try to load user-editable variables from IO/setup_Vars.m
% if no file exists, spit error message and die
try
    load('../IO/IN/Baseline_IN.mat');
catch ME
    if (strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile'))
        fprintf(1, '\nError: Could not read Baseline Input file from IO/.\nDid you run setup_Vars.m?\n');
        exit(2);  % Die with value of 2
    else
        rethrow(ME)
        exit(3); % Die with value of 3, unknown error
    end
end

% Import country-level age-demographic parameters
Import_CountryDemogParams;

% start parallel pool
setupParallelPool();

%% Run fitting procedure for each site

for iSites = 1:length(Sites)
    
    % define data
    MfData  = eval(sprintf('%sMf',Sites{iSites}));
    bCulex  = eval(sprintf('%sbCulex',Sites{iSites}));
    ABR     = eval(sprintf('%sABR',Sites{iSites}));
    
    % define demographic parameters demoA and demoB according to
    % country
    [demoA1, demoB1] = get_demoA_and_demoB(Countries{iSites},...
        Country_demo,demoA,demoB);
    
    % find maximum age in population
    if isnan(MfData)
        ageMax = 69;
    else
        ageMax = max(MfData(:,4));
    end
    
    % run fitting procedure and save results in IO/OUT
        RunSelectParameterVectors(Sites{iSites},NParamVecs,MfData,ABR,...
            SIR_samples,demoA1,demoB1,ageMax,bCulex,minPrev,maxPrev,...
            ABRmax,estABRflag)
    
end