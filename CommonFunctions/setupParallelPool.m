% Function: setupParallelPool.m
%
% Description: 
%   Initializes the parcluster and parpool used for all calculations.
%
% Inputs:
%   None.
%
% Outputs:
%   The parpool continues to exist after exiting function. No explicit output.
%
% Assumptions:
%   It is assumed the version of matlab used will be 2013 or more recent.
%   
% Date: Jan 17, 2018 CJK
% ________________________________________
function setupParallelPool()

NumWorkers = 4; % Change this variable depending on system capabilities.

p = gcp('nocreate');
if isempty(p)
    cc=parcluster();
    t=tempname();
    mkdir(t);
    cc.JobStorageLocation=t;
    parpool(cc, NumWorkers);
end
