%% Run the tACS Challenge analyses, in full or in part.
%% Script originally written by Anne-Fiona Griesfeller, CNRS Toulouse
%
% H1 to H3 are always analysed as a block. Within H4 to H6, single hypotheses
% can be selected. H1 to H3 must be run once before H4 to H6, because those
% analyses read the PhaseMetrics.mat file that H1 to H3 writes into each lab
% folder.

clear
clc

%% Project paths
% The project root is derived from the location of this script, which sits in
% <project_root>/scripts. Every lab following the folder structure in the
% README can run this file unchanged. Otherwise change

project_root   = fileparts(fileparts(mfilename('fullpath')));
data_root      = fullfile(project_root, 'data');
functions_root = fullfile(project_root, 'functions');
results_root   = fullfile(project_root, 'results');

if ~isfolder(data_root)
    error(['Data folder not found at %s. Check that this script is in ' ...
        '<project_root>/scripts and that the folder structure matches the README, ' ...
        'or set project_root by hand above.'], data_root);
end

addpath(genpath(functions_root));

if isempty(which('circ_mtest'))
    error(['Circular Statistics Toolbox (CircStat) not found. Download it from ' ...
        'MATLAB Central File Exchange entry 10676 ' ...
        '(https://www.mathworks.com/matlabcentral/fileexchange/10676) ' ...
        'and unzip it into %s'], functions_root);
end

%% Labs to include
% Leave empty to analyse every lab folder in data directory.
% To restrict e.g., lab_ids = [3 18]
lab_ids = [];

if isempty(lab_ids)
    lab_dirs = dir(fullfile(data_root, 'L*'));
    lab_dirs = lab_dirs([lab_dirs.isdir]);
    lab_ids  = [];
    for iDir = 1:numel(lab_dirs)
        thisId = str2double(regexprep(lab_dirs(iDir).name, '^L', ''));
        if ~isnan(thisId)
            lab_ids(end+1) = thisId;
        end
    end
    lab_ids = sort(lab_ids);
    if isempty(lab_ids)
        error('No lab folders of the form L## were found in %s.', data_root);
    end
    fprintf('Labs found in %s: %s\n', data_root, mat2str(lab_ids));
end

%% Switches
RunH1H2H3 = true;
RunH4H5H6 = true;

% Optional, inside H4 to H6
RunH4 = true;
RunH5 = true;
RunH6 = true;
MakePlotsH4toH6 = true;

nPermutations = 1000;

% if true then PhaseMetrics.mat and AnalysisTable.mat are rebuild for all
forceRecompute = true;

%% Output folders
outH1H2H3 = fullfile(results_root, 'H1H2H3', 'pooled');
outH4H5H6 = fullfile(results_root, 'H4H5H6');

%% Run H1 to H3 pooled
if RunH1H2H3
    fprintf('\nRunning H1 H2 H3 pooled \n');
    OUT_H1H2H3 = tACSChallenge_runH1H2H3_Pooled( ...
        data_root, lab_ids, ...
        'OutDir', outH1H2H3, ...
        'perm', nPermutations, ...
        'MakePlots', true, ...
        'forceRecompute', forceRecompute);

    save(fullfile(outH1H2H3, sprintf('OUT_H1H2H3_%s.mat', ...
        datestr(now, 'yyyymmdd_HHMMSS'))), 'OUT_H1H2H3', '-v7.3');
end

%% Run H4 to H6
if RunH4H5H6

    if ~RunH1H2H3
        missing = {};
        for iLab = 1:numel(lab_ids)
            pm = fullfile(data_root, sprintf('L%02d', lab_ids(iLab)), 'PhaseMetrics.mat');
            if ~isfile(pm)
                missing{end+1} = sprintf('L%02d', lab_ids(iLab)); %#ok<SAGROW>
            end
        end
        if ~isempty(missing)
            warning(['PhaseMetrics.mat is missing for: %s. Set RunH1H2H3 = true and ' ...
                'run this script again'], strjoin(missing, ', '));
        end
    end

    fprintf('\nRunning H4 H5 H6 \n');
    OUT_H4H5H6 = tACSChallenge_runH4H5H6( ...
        data_root, lab_ids, ...
        'OutDir', outH4H5H6, ...
        'RunH4', RunH4, ...
        'RunH5', RunH5, ...
        'RunH6', RunH6, ...
        'MakePlots', MakePlotsH4toH6, ...
        'forceRecompute', forceRecompute);

    save(fullfile(outH4H5H6, sprintf('OUT_H4H5H6_%s.mat', ...
        datestr(now, 'yyyymmdd_HHMMSS'))), 'OUT_H4H5H6', '-v7.3');
end

fprintf('\n Done \n');
