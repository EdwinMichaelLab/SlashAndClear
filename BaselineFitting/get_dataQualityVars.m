% Function: get_dataQualityVars
%  
% Description: 
%   Tests the quality of the MF data and returns boolean variables corresponding to the
%   quality of the data for testing and conditional branching.
%  
% Inputs:
%   Mfdata
%
% Outputs:
%   Returns 3 boolean variables: AgeMf_asInput, OverallMf_asInput,
%   and OverallMfRange_asInput
%
% Date: Dec 19, 2017   CJK
% ________________________________________

function [AgeMf_asInput, OverallMf_asInput, OverallMfRange_asInput] = get_dataQualityVars(MfData)

% This is broken out into its own function to easily allow testing

% Setting variables to determine the quality of the data
% Essentially these if blocks are the same as the if statements within the previous 
% implementation of Main_ToRunBaselineFitting.m
if length(MfData(:,1)) > 1 % Mf age profile
    AgeMf_asInput = true;
else
    AgeMf_asInput = false; 
end
if length(MfData(:,1)) == 1 && ~isnan(MfData(1,1)) % Overall Mf
    OverallMf_asInput = true;
else
    OverallMf_asInput = false;
end
if isnan(MfData(1,1)) % No Baseline Mf
    OverallMfRange_asInput = true;
else
    OverallMfRange_asInput = false;
end

