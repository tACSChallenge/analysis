function OUT = tACSChallenge_runH4H5H6(data_root, lab_ids, varargin)
%% Wrapper for H4 H5 H6 global models.
%% script originally written by Anne-Fiona Griesfeller, CNRS Toulouse, in Jan 26
% Calls tACSChallenge_runGlobalModels and saves OUT.

p = inputParser;
p.addParameter('OutDir', '', @(x)ischar(x) || isstring(x));
p.addParameter('RunH4', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('RunH5', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('RunH6', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('MakePlots', true, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});

OutDir = string(p.Results.OutDir);
RunH4 = logical(p.Results.RunH4);
RunH5 = logical(p.Results.RunH5);
RunH6 = logical(p.Results.RunH6);
MakePlots = logical(p.Results.MakePlots);

if OutDir == ""
    OutDir = fullfile(string(data_root), 'results', 'H4H5H6');
end

if ~isfolder(OutDir)
    mkdir(OutDir);
end

% Safety check, PhaseMetrics must exist per lab for H4 H5 H6
for i = 1:numel(lab_ids)
    lab_path = fullfile(data_root, sprintf('L%d', lab_ids(i)));
    pm = fullfile(lab_path, 'PhaseMetrics.mat');
    if ~isfile(pm)
        warning('Missing PhaseMetrics.mat for lab %d at %s', lab_ids(i), pm);
    end
end

OUT = tACSChallenge_runGlobalModels( ...
    data_root, lab_ids, ...
    'RunH4', RunH4, ...
    'RunH5', RunH5, ...
    'RunH6', RunH6, ...
    'MakePlots', MakePlots, ...
    'OutDir', OutDir);

outFile = fullfile(OutDir, sprintf('GlobalModels_OUT_%s.mat', datestr(now, 'yyyymmdd_HHMMSS')));
save(outFile, 'OUT', '-v7.3');
fprintf('\nSaved H4 H5 H6 OUT struct to:\n  %s\n', outFile);

end
