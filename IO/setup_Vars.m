% Function: setup_Vars
%  
% Description: 
%   Instatiates all necessary variables for BaselineFitting Breakpoint
%   Calculations, and Interventions. This is the file which is edited by
%   the user for the location's data. This file will be ran which simply
%   instantiates the variables and saves them to a .mat file to be loaded.
% 
% Inputs:
%   None: Data is entered within this file.
%
% Outputs:
%   Variables are saved to files for: Baseline_IN.mat, Break_IN.mat, Inter_IN.mat
%_____________________________________________________________________________

addpath('../BaselineFitting/');
addpath('../Interventions/');
%{
 * BaselineFitting Variables
 * Orginally from Main_ToRunBaselineFitting.m
%}

% Site details: Relevant endemic regions, Country (pertains only Sub-Saharan Africa cases), Village

% country names should follow format of CountryDemogParams.csv
Countries = {'Uganda','Uganda','Uganda','Uganda'}; 
Sites = {'PalaurePacunaci', 'Masaloa', 'Nyimanji', 'OlimbuniAroga'};

% Site data: Mf prevalence, ABR, vector genera
% arranged in region-specific .m data files created by the user
[PalaurePacunaciMf, MasaloaMf, OlimbuniArogaMf, NyimanjiMf, ...
    PalaurePacunaciABR, MasaloaABR, OlimbuniArogaABR, NyimanjiABR, ...
    PalaurePacunacibCulex, MasaloabCulex, OlimbuniArogabCulex, NyimanjibCulex] = baseline_data;

% ABR range for use in sites without baseline ABR data
estABRflag = 1;
ABRmin = 100; 
ABRmax = 65000; 

% Mf range for use in sites without baseline Mf data 
minPrev = 0.05;
maxPrev = 0.99; 

save IN/Baseline_IN.mat Countries Sites ...
    PalaurePacunaciMf MasaloaMf OlimbuniArogaMf NyimanjiMf ...
    PalaurePacunaciABR MasaloaABR OlimbuniArogaABR NyimanjiABR ...
    PalaurePacunacibCulex MasaloabCulex OlimbuniArogabCulex NyimanjibCulex ...
    ABRmin ABRmax estABRflag minPrev maxPrev;


%{
 * Breakpoint Calculations Variables
 * Orginally from Main_ToCalcBreakpoints.m
%}

% Uses Regions, Countries, and sites  from above.
save('IN/Break_IN.mat', 'Countries', 'Sites');


%{
 * Interventions Variables
%}

% Uses Regions, Countries, and sites from above

% Site data: Mf prevalence, blood sample volume, MDA regimen, MDA frequency,...
% MDA coverage, number of years of treatment, vector control, switch year
% arranged in region-specific .m data files created by the user
[PalaurePacunaciMf, MasaloaMf, OlimbuniArogaMf, NyimanjiMf, ...
    PalaurePacunaciMDAFreq, MasaloaMDAFreq, OlimbuniArogaMDAFreq, NyimanjiMDAFreq, ...
    PalaurePacunaciSwitchYr, MasaloaSwitchYr, OlimbuniArogaSwitchYr, NyimanjiSwitchYr, ...
    PalaurePacunaciMDACov, MasaloaMDACov, OlimbuniArogaMDACov, NyimanjiMDACov, ...
    PalaurePacunaciNumYears, MasaloaNumYears, OlimbuniArogaNumYears, NyimanjiNumYears, ...
    PalaurePacunaciVC, MasaloaVC, OlimbuniArogaVC, NyimanjiVC,...
     PalaurePacunaciSlash_t, MasaloaSlash_t, OlimbuniArogaSlash_t, NyimanjiSlash_t] = Intv_data;

%save('IN/Inter_IN.mat', 'Countries', 'SSAsites', 'Regions');
save IN/Inter_IN.mat Countries Sites PalaurePacunaciMf MasaloaMf OlimbuniArogaMf NyimanjiMf ...
    PalaurePacunaciMDAFreq MasaloaMDAFreq OlimbuniArogaMDAFreq NyimanjiMDAFreq ...
    PalaurePacunaciSwitchYr MasaloaSwitchYr OlimbuniArogaSwitchYr NyimanjiSwitchYr ...
    PalaurePacunaciMDACov MasaloaMDACov OlimbuniArogaMDACov NyimanjiMDACov ...
    PalaurePacunaciNumYears MasaloaNumYears OlimbuniArogaNumYears NyimanjiNumYears ...
    PalaurePacunaciVC MasaloaVC OlimbuniArogaVC NyimanjiVC...
     PalaurePacunaciSlash_t MasaloaSlash_t OlimbuniArogaSlash_t NyimanjiSlash_t