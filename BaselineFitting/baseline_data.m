% baseline data for all sites

function [PalaurePacunaciMf, MasaloaMf, OlimbuniArogaMf, NyimanjiMf, ...
    PalaurePacunaciABR, MasaloaABR, OlimbuniArogaABR, NyimanjiABR, ...
    PalaurePacunacibCulex, MasaloabCulex, OlimbuniArogabCulex, NyimanjibCulex] = baseline_data

%% Baseline Mf Data
% Best data to provide would be age-stratified prevalence
% Variable names = {SiteName}Mf

% 1st column = mid-age of group; 2nd: Total number of samples; 3rd: Mf +ves;
% 4th: upper age of group

% If age-stratified data is not available, enter overall community
% prevalence in a single line following the column definitions above

PalaurePacunaciMf =[35 100 100 70];

MasaloaMf =[35 100 76 70];

OlimbuniArogaMf =[40 50 12 70];

NyimanjiMf =[35 160 93 70];


%% Baseline ABR data
% If not available for a particular site, enter value as NaN

PalaurePacunaciABR = NaN;
MasaloaABR = NaN;
OlimbuniArogaABR = NaN;
NyimanjiABR = NaN;

%% Black fly species flag based on dominant species 
% (0: no armature, S neavei, S damnosum, 1: armature)

PalaurePacunacibCulex = 0;
MasaloabCulex = 0;
OlimbuniArogabCulex = 0;
NyimanjibCulex = 0;

end
