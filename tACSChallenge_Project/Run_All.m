%% With this script all anylses or selected analyses can be run.
%% Script originally written by Anne-Fiona Griesfeller, CNRS Toulouse, in Jan 26
% H1 to H3 are always analyzed as a block. For H4-6 single hypotheses can be analyzed.  
% H1-3 must be run once before running H4-6, as these analyses require
% output from H1-2

clear
clc

%% Project paths
project_root   = 'C:\Users\griesfeller\Documents\tACSChallenge_Project';
data_root      = fullfile(project_root, 'data');
functions_root = fullfile(project_root, 'functions');
results_root   = fullfile(project_root, 'results');

addpath(genpath(functions_root));

%% Labs to include
lab_ids = [45]; % Change this to [18 19 20 ...] if further data is available. 

%% Switches
RunH1H2H3 = true;
RunH4H5H6 = false;

% Optional inside H4 H5 H6 
RunH4 = false;
RunH5 = false;
RunH6 = false;
MakePlotsH4toH6 = false;

%% H1 to H3 pooled output folder
outH1H2H3 = fullfile(results_root, 'H1H2H3', 'pooled');

%% H4 to H6 output folder
outH4H5H6 = fullfile(results_root, 'H4H5H6');

%% Run H1 to H3 pooled
if RunH1H2H3
    fprintf('\nRunning H1 H2 H3 pooled \n');
    OUT_H1H2H3 = tACSChallenge_runH1H2H3_Pooled( ...
        data_root, lab_ids, ...
        'OutDir', outH1H2H3, ...
        'perm', 100, ...
        'MakePlots', true);

    save(fullfile(outH1H2H3, sprintf('OUT_H1H2H3_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'))), 'OUT_H1H2H3', '-v7.3');
end

%% Run H4 to H6
if RunH4H5H6
    fprintf('\nRunning H4 H5 H6 \n');
    OUT_H4H5H6 = tACSChallenge_runH4H5H6( ...
        data_root, lab_ids, ...
        'OutDir', outH4H5H6, ...
        'RunH4', RunH4, ...
        'RunH5', RunH5, ...
        'RunH6', RunH6, ...
        'MakePlots', MakePlotsH4toH6);

    save(fullfile(outH4H5H6, sprintf('OUT_H4H5H6_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'))), 'OUT_H4H5H6', '-v7.3');
end

fprintf('\n Done \n');

