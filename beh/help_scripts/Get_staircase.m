function results = Get_staircase(data_path, lab_num, sub_num, results_root)
%GET_STAIRCASE  Analyse the staircase logfile(s) of one tACS Challenge subject.
%
%   Get_staircase(data_path, lab_num, sub_num)
%   Get_staircase(data_path, lab_num, sub_num, results_root)
%   results = Get_staircase(...)
%
%   Combines the function tACSChallenge_AnalyzeStaircase.m (import / 
%   analysis / plot / export) into one function.
%
%   The function imports the logfiles of the staircase procedure and
%   calculates the threshold brightness (base brightness + target
%   brightness) which can then be used in the experiment. For every
%   staircase it writes an .xlsx table with all trials and a .png of the
%   staircase. The last value of the column "ThresholdBrightness" (= mean
%   final brightness of the last 10 trials) is the value to be used in the
%   experiment.
%
%   Trials in which the change in luminance results in no change (target
%   LED intensity equal to the base brightness) are ignored, as before.
%
% INPUT ------------------------------------------------------------------
%   data_path     Folder that is searched RECURSIVELY for the staircase
%                 logfiles, e.g.
%                 'O:\...\tACSChallenge_Project\data\'
%                 (a lab or subject folder also works, it is just faster)
%   lab_num       Lab number:      18, '18', 'L18'  are all accepted
%   sub_num       Subject number:  1, '01', 'S01'   are all accepted
%   results_root  (optional) folder in which the per-subject result folders
%                 are created. Default:
%                 <tACSChallenge_Project>\results\staircases
%                 (the project root is found by walking up from data_path)
%
% OUTPUT -----------------------------------------------------------------
%   results       1xN struct, one entry per staircase found, with fields
%                 label, file, subject, base_brightness,
%                 threshold_brightness, trials (table), data (raw signals),
%                 xlsx, png
%
% FILES WRITTEN ----------------------------------------------------------
%   into  <results_root>\<subject>\ , e.g. ...\results\staircases\sub-L18_S01\
%     one staircase:
%       sub-L18_S01_StaircaseResults_DATE_<dd_mm_yyyy_HH_MM_SS>.xlsx / .png
%     more than one staircase (chronological order, oldest = A):
%       sub-L18_S01_A_StaircaseResults_DATE_<...>.xlsx / .png
%       sub-L18_S01_B_StaircaseResults_DATE_<...>.xlsx / .png
%
% EXAMPLE ----------------------------------------------------------------
%   Get_staircase('O:\...\tACSChallenge_Project\data', 18, 1)
%   r = Get_staircase('O:\...\data', 'L18', 'S02');
%   r(1).threshold_brightness
%
% -------------------------------------------------------------------------
% Analysis adapted by Giuseppe Di Dona from tACSChallenge_ImportData
% (written by Florian Kasten), edited by Wahriman Andrade de Araújo.
% -------------------------------------------------------------------------

%% ------------------------------------------------------------------ input
if nargin < 3
    error('Get_staircase:notEnoughInput', ...
        'Syntax: Get_staircase(data_path, lab_num, sub_num[, results_root])');
end
if nargin < 4
    results_root = '';
end

data_path = char(data_path);
if ~exist(data_path, 'dir')
    error('Get_staircase:noDataPath', 'data_path does not exist:\n  %s', data_path);
end

labN = local_parseNumber(lab_num, 'lab_num');
subN = local_parseNumber(sub_num, 'sub_num');

% Patterns that identify this subject in a file name, tolerant to zero
% padding (L18 / L018) and to the separator used (sub-L18_S01, L18S1, ...).
% The (?![0-9]) prevents L1 from matching L18 and S1 from matching S12.
labPat = sprintf('L0*%d(?![0-9])', labN);
subPat = sprintf('S0*%d(?![0-9])', subN);

%% ---------------------------------------------------- where results go
if isempty(results_root)
    results_root = fullfile(local_findProjectRoot(data_path), 'results', 'staircases');
else
    results_root = char(results_root);
end

%% -------------------------------------------- find the staircase logfiles
files = local_findStaircaseFiles(data_path, labPat, subPat, results_root);

if isempty(files)
    error('Get_staircase:noFiles', ...
        ['No staircase logfile found for lab %d / subject %d in:\n  %s\n' ...
         'Expected a file such as sub-L%02d_S%02d_Staircase_<date>.tsv'], ...
        labN, subN, data_path, labN, subN);
end

% Subject code, taken from the file name exactly as before
% (everything in front of "_Staircase_"), e.g. sub-L18_S01
[~, firstName] = fileparts(files(1).name);
idx = regexpi(firstName, '_Staircase', 'once');
if isempty(idx) || idx == 1
    sj = sprintf('sub-L%02d_S%02d', labN, subN);   % fallback
else
    sj = firstName(1:idx-1);
end

outDir = fullfile(results_root, sj);
[ok, msg] = mkdir(outDir);
if ~ok
    error('Get_staircase:mkdir', 'Could not create output folder:\n  %s\n%s', outDir, msg);
end

% One timestamp for the whole run, so that the files of staircase A and B
% belong visibly together
datestring = local_timestamp();

fprintf('\n%s  |  lab %d, subject %d  |  %d staircase(s) found\n', sj, labN, subN, numel(files));
fprintf('Results -> %s\n', outDir);

%% ---------------------------------------------- analyse every staircase
nStair  = numel(files);
out     = struct([]);
labels  = arrayfun(@(k) local_letter(k), 1:nStair, 'UniformOutput', false);

for k = 1:nStair

    thisFile = fullfile(files(k).folder, files(k).name);

    % Letter only when the subject has more than one staircase
    if nStair > 1
        suffix = ['_' labels{k}];
    else
        suffix = '';
    end
    baseName = [sj, suffix, '_StaircaseResults_DATE_', datestring];

    fprintf('\n[%d/%d] %s\n', k, nStair, thisFile);

    try
        res = local_analyseStaircase(thisFile, sj, labels{k}, suffix, outDir, baseName);
    catch ME
        warning('Get_staircase:analysisFailed', ...
            'Staircase %s could not be analysed (%s):\n  %s\n  %s', ...
            labels{k}, ME.identifier, thisFile, ME.message);
        continue
    end

    if isempty(out)
        out = res;
    else
        out(end+1) = res; %#ok<AGROW>
    end
end

if isempty(out)
    error('Get_staircase:allFailed', ...
        'None of the %d staircase file(s) of %s could be analysed.', nStair, sj);
end

%% ------------------------------------------------------------- summary
fprintf('\n---------------------------------------------------------------\n');
for k = 1:numel(out)
    fprintf('%s  staircase %s : base = %.4f | threshold brightness = %.4f  (%d trials)\n', ...
        sj, out(k).label, out(k).base_brightness, out(k).threshold_brightness, out(k).n_trials);
end
fprintf('---------------------------------------------------------------\n\n');

if nargout > 0
    results = out;
end

end % ======================================================================


%% =========================== ANALYSIS ====================================
function res = local_analyseStaircase(filename, sj, label, suffix, outDir, baseName)
% Import one logfile, compute the threshold, write .xlsx and .png.

%% ------------------------------------------------------------ import
dataLines = [1, Inf];

opts = delimitedTextImportOptions("NumVariables", 14);
% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";
% Specify column names and types
opts.VariableNames = ["Var1", "Timeus", "TACSV", "BNCMode", "BNCInV", "BNCOutV",...
    "LeftButton", "RightButton", "LED0Bright", "LED1Bright", "LED2Bright", "LED3Bright", "LED4Bright", "LED5Bright"];
opts.SelectedVariableNames = ["Timeus", "TACSV", "BNCMode", "BNCInV", "BNCOutV",...
    "LeftButton", "RightButton", "LED0Bright", "LED1Bright", "LED2Bright", "LED3Bright", "LED4Bright", "LED5Bright"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double",...
    "double", "double", "double", "double", "double", "double", "double", "double"];
% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var1", "EmptyFieldRule", "auto");

% Import the data
tmp_matrix = readtable(filename, opts);

%% ------------------------------------------------- convert to output type
tmp_matrix = table2array(tmp_matrix);
tmp_matrix = tmp_matrix(2:end,:);

if size(tmp_matrix,1) < 2
    error('Get_staircase:emptyLog', 'The logfile contains no data rows.');
end

%% --------------------------------------------------------- extract data
% Delta T
data.raw_dt = diff(tmp_matrix(:,1));

% handle reset of processor clock in teensy during recording
% check if any dt is negative (step backwards in time)
if sum(data.raw_dt < 0) > 0

    % replace with the mean of all other dt
    data.raw_dt(data.raw_dt < 0) = mean(data.raw_dt(data.raw_dt > 0), 'omitnan');

end

data.dt = mean(data.raw_dt,'omitnan')/10^6;   % uS to Seconds
% Sampling Rate
data.Fs = 1/data.dt;
% tACS
data.tACS = tmp_matrix(:,2);
% BNC (V)
data.BNC_In  = tmp_matrix(:,4);
data.BNC_Out = tmp_matrix(:,5);
% Buttons
data.L_Button = tmp_matrix(:,6);
data.R_Button = tmp_matrix(:,7);
% LED
data.LEDs = tmp_matrix(:,8:13);

%% ------------------------------------------- onsets of button presses
% Fuse R_Button and L_Button
% takes the max resp over all L/R button so it keeps the value of the active
% button which is always superior to the other (one is 0 and the other is 1)
LR_Button = max([data.R_Button data.L_Button], [], 2);

% codes onsets (+1) and offsets (-1), a 0 is added in the beginning to
% correct the indexing
respOnsets = [0 ; diff(LR_Button)];

% removes offsets going from -1 to 0, only onsets remain (+1)
respOnsets(respOnsets < 0) = 0;

% this is when the subject pressed the button (in sample points)
RespLat = find(respOnsets > 0);

%% ----- leverage the inactive central LED to remove the intensity offset
% subtracts the value of the central led (1st column), transformed into a
% matrix, from all leds
data.LEDs = data.LEDs - repmat(data.LEDs(:,1), 1, size(data.LEDs,2));

%% ------------------------------------------- merge LED signals into one
% takes the max led over all leds so it keeps the value of the active LED,
% which is always superior to the others (always 0.2, now zero after the
% previous step)
LED = max(data.LEDs, [], 2);

%% ------------------------------------------------ onsets of target LED
LEDOnsets = [0; diff(LED)];   % same procedure used for the responses
LEDOnsets(LEDOnsets < 0) = 0;

% this is when the LEDs were on (in sample points)
LEDLat = find(LEDOnsets > 0);

if isempty(LEDLat)
    error('Get_staircase:noTrials', 'No LED onset found in this logfile.');
end

% this is the brightness of the LEDs in sample points -> TARGET BRIGHTNESS,
% which still needs the base brightness added to it
LEDBright = LED(LEDLat);

%% ------------------------------------------------- base brightness (a)
% Raw LED brightness columns: LED0Bright ... LED5Bright
data.LEDsRaw = tmp_matrix(:,8:13);
data.LEDs    = data.LEDsRaw;

% Dynamically estimate the base brightness from the first transition where
% the LEDs go from 0 to a non-zero value.
tol = 0.1;
baseShiftIdx = find( ...
    all(abs(data.LEDsRaw(1:end-1,:)) <= tol, 2) & ...
    any(abs(data.LEDsRaw(2:end,:))   >  tol, 2), ...
    1, 'first');

if isempty(baseShiftIdx)
    error('Get_staircase:noBaseBright', ...
        'Could not identify BASEBright: no 0-to-nonzero transition found in LED0Bright-LED5Bright.');
end

baseShiftIdx = baseShiftIdx + 1;
baseVals     = data.LEDsRaw(baseShiftIdx, abs(data.LEDsRaw(baseShiftIdx,:)) > tol);

BASEBright            = min(baseVals);
data.BASEBright       = BASEBright;
data.BASEBrightSample = baseShiftIdx;

%% ------------------------------------ matrix of behavioural responses
nTrials      = numel(LEDLat);
trials_stair = zeros(nTrials, 7);

for i = 1:nTrials

    curr_t = LEDLat(i);
    trials_stair(i,1) = curr_t;

    % the target is considered detected if a button press occurred within
    % 1.2 s (1200 samples)
    if any(RespLat > curr_t & RespLat < curr_t + 1200)

        trials_stair(i,2) = 1;                              % hit
        RT = RespLat(RespLat > curr_t & RespLat < curr_t + 1200);
        trials_stair(i,3) = min(RT) - curr_t;               % response time

    end

    trials_stair(i,4) = LEDBright(i);                       % target brightness (delta)
    trials_stair(i,5) = LEDBright(i) + BASEBright;          % final = target + base
    trials_stair(i,6) = i;                                  % trial number

end

%% ------------------- thresholded brightness from the last 10 trials
% The last value of column 7 is the thresholded brightness to be used.
nAvg = 10;
if nTrials < nAvg
    warning('Get_staircase:fewTrials', ...
        ['Only %d trials in this staircase - the threshold is averaged over ' ...
         'all available trials instead of the last %d.'], nTrials, nAvg);
    nAvg = nTrials;
end

for i = nAvg:nTrials
    trials_stair(i,7) = mean(trials_stair(i-nAvg+1 : i, 5));
end

threshold_bright = trials_stair(nTrials, 7);   % thresholded brightness

%% ------------------------------------------------ save the .xlsx table
trials_stair(:,8:9) = NaN;
trials_stair(1,8)   = BASEBright;
trials_stair(1,9)   = threshold_bright;

tab = array2table(trials_stair, 'VariableNames', ...
    {'Latency', 'Response', 'ISI', 'TargetBrightness_delta', 'FinalBrightness', ...
     'TrialNumber', 'ThresholdBrightness', 'Base_Brightness (a)', 'Target_Brightness (b)'});

xlsxFile = fullfile(outDir, [baseName '.xlsx']);
writetable(tab, xlsxFile);

%% ------------------------------------------------------ plot the result
markers    = {'v', '^'};
colors     = [1 0 0; 0 0.6 0];                 % Red  Green
markers_idx = trials_stair(:,2) + 1;           % 1 = miss (v, red), 2 = hit (^, green)

if isempty(suffix)
    figName = ['Staircase Subject ' sj];
else
    figName = ['Staircase Subject ' sj ' - staircase ' label];
end

fig = figure('Name', figName, 'Color', 'w');
hLine = plot(trials_stair(:,5), 'Color', [0.4 0.4 0.4]);
xlabel('Trials');
ylabel(sprintf('Final Brightness (Target + Base %g)', BASEBright));
title(figName, 'Interpreter', 'none');
hold on;

legH = hLine;
legL = {'Staircase'};
legendNames = {'Miss', 'Hit'};

for i = 1:numel(markers)
    idx = markers_idx == i;
    if ~any(idx), continue; end
    h = scatter(trials_stair(idx,6), trials_stair(idx,5), 30, 'filled', markers{i}, ...
        'MarkerEdgeColor', colors(i,:), 'MarkerFaceColor', colors(i,:));
    legH = [legH, h];             %#ok<AGROW>
    legL{end+1} = legendNames{i}; %#ok<AGROW>
end

% horizontal line over the trials the threshold was averaged from
hThr = plot([nTrials-nAvg+1 nTrials], [threshold_bright threshold_bright], ...
    'b-', 'LineWidth', 2);
legH = [legH, hThr];
legL{end+1} = sprintf('Threshold (last %d)', nAvg);

subStr = sprintf('Base Brightness (a) = %s and Target Brightness (b) = %s', ...
    num2str(BASEBright), num2str(threshold_bright));
try
    subtitle(subStr);                       % subtitle exists from R2020b on
catch
    title({figName, subStr}, 'Interpreter', 'none');
end

lgd = legend(legH, legL, 'Location', 'best');

% Make figure/axes light
ax = gca;
set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');

try
    % Title / labels / ticks text colour
    ax.Title.Color  = 'k';
    ax.XLabel.Color = 'k';
    ax.YLabel.Color = 'k';
    ax.ZLabel.Color = 'k';

    % Legend text + box
    lgd.TextColor = 'k';
    lgd.Color     = 'w';                    % legend background
    lgd.EdgeColor = 'k';                    % legend border
catch
    % older releases without these object properties - cosmetics only
end

% Save the result
pngFile = fullfile(outDir, [baseName '.png']);
saveas(fig, pngFile);

fprintf('   base = %.4f | threshold = %.4f | %d trials\n', BASEBright, threshold_bright, nTrials);
fprintf('   %s\n   %s\n', xlsxFile, pngFile);

%% ------------------------------------------------------------- output
res = struct( ...
    'label',                label, ...
    'file',                 filename, ...
    'subject',              sj, ...
    'n_trials',             nTrials, ...
    'base_brightness',      BASEBright, ...
    'threshold_brightness', threshold_bright, ...
    'trials',               tab, ...
    'data',                 data, ...
    'xlsx',                 xlsxFile, ...
    'png',                  pngFile);

end % local_analyseStaircase


%% ========================= FILE SELECTION ================================
function files = local_findStaircaseFiles(data_path, labPat, subPat, results_root)
% Recursively collect the staircase logfiles of one subject, oldest first.

% Only search the lab folder(s) when there are any - just for speed
searchRoots = {};
top = dir(data_path);
top = top([top.isdir]);
for k = 1:numel(top)
    if any(strcmp(top(k).name, {'.','..'})), continue; end
    if ~isempty(regexpi(top(k).name, labPat, 'once'))
        searchRoots{end+1} = fullfile(data_path, top(k).name); %#ok<AGROW>
    end
end
if isempty(searchRoots)
    searchRoots = {data_path};
end

all_files = struct('name', {}, 'folder', {}, 'datenum', {});
for k = 1:numel(searchRoots)
    all_files = [all_files, local_listFiles(searchRoots{k}, 0, 8)]; %#ok<AGROW>
end

if isempty(all_files), files = all_files; return; end

% Keep plausible logfiles only: *Staircase* , a text-like extension, not a
% result file written by this function, and not inside the results folder
keep = false(1, numel(all_files));
for k = 1:numel(all_files)
    [~, base, ext] = fileparts(all_files(k).name);
    if ~any(strcmpi(ext, {'.tsv', '.txt', '.csv', '.log', '.dat', ''})), continue; end
    if isempty(regexpi(base, 'Staircase', 'once')),           continue; end
    if ~isempty(regexpi(base, 'StaircaseResults', 'once')),   continue; end
    if ~isempty(results_root) && strncmpi(all_files(k).folder, results_root, numel(results_root))
        continue
    end
    keep(k) = true;
end
all_files = all_files(keep);
if isempty(all_files), files = all_files; return; end

% Match the subject on the file name; if nothing matches, try the full path
isSubject = false(1, numel(all_files));
for k = 1:numel(all_files)
    [~, base] = fileparts(all_files(k).name);
    isSubject(k) = ~isempty(regexpi(base, labPat, 'once')) && ...
                   ~isempty(regexpi(base, subPat, 'once'));
end
if ~any(isSubject)
    for k = 1:numel(all_files)
        p = fullfile(all_files(k).folder, all_files(k).name);
        isSubject(k) = ~isempty(regexpi(p, labPat, 'once')) && ...
                       ~isempty(regexpi(p, subPat, 'once'));
    end
end
files = all_files(isSubject);
if isempty(files), return; end

% Chronological order (oldest = staircase A). The date in the file name is
% used when present, the file date otherwise.
t = zeros(1, numel(files));
for k = 1:numel(files)
    [~, base] = fileparts(files(k).name);
    tok = regexp(base, '(\d{4})_(\d{2})_(\d{2})_(\d{2})_(\d{2})_(\d{2})', 'tokens', 'once');
    if numel(tok) == 6
        v = cellfun(@str2double, tok);
        t(k) = datenum(v(1), v(2), v(3), v(4), v(5), v(6));
    else
        t(k) = files(k).datenum;
    end
end
[~, order] = sort(t);
files = files(order);

end % local_findStaircaseFiles


function files = local_listFiles(root, depth, maxDepth)
% Recursive dir() that does not rely on the '**' wildcard (works in every
% MATLAB version). Hidden folders are skipped.

files = struct('name', {}, 'folder', {}, 'datenum', {});
if depth > maxDepth || ~exist(root, 'dir'), return; end

d = dir(root);
for k = 1:numel(d)

    nm = d(k).name;
    if strcmp(nm, '.') || strcmp(nm, '..') || (~isempty(nm) && nm(1) == '.')
        continue
    end

    if d(k).isdir
        files = [files, local_listFiles(fullfile(root, nm), depth+1, maxDepth)]; %#ok<AGROW>
    else
        files(end+1) = struct('name', nm, 'folder', root, 'datenum', d(k).datenum); %#ok<AGROW>
    end
end

end % local_listFiles


%% ============================ HELPERS ====================================
function root = local_findProjectRoot(data_path)
% Walk up from data_path to find the tACSChallenge_Project folder. If it is
% not there, the parent of a folder called "data" is used, otherwise
% data_path itself.

p = local_stripSep(local_absPath(data_path));

q = p;
for k = 1:10
    [parent, name] = fileparts(q);
    if strcmpi(name, 'tACSChallenge_Project')
        root = q;
        return
    end
    if isempty(parent) || strcmp(parent, q)
        break
    end
    q = parent;
end

[parent, name] = fileparts(p);
if strcmpi(name, 'data') && ~isempty(parent)
    root = parent;                 % ...\<project>\data  ->  ...\<project>
else
    root = p;
end

end % local_findProjectRoot


function p = local_absPath(p)
% Absolute path without changing the current folder.
p = char(p);
if exist(p, 'dir')
    d = dir(p);
    if ~isempty(d) && isfield(d, 'folder') && ~isempty(d(1).folder)
        p = d(1).folder;           % dir() returns the resolved absolute path
    end
end
end % local_absPath


function p = local_stripSep(p)
% Remove trailing file separators ('...\data\' -> '...\data')
while ~isempty(p) && (p(end) == '\' || p(end) == '/')
    p(end) = [];
end
end % local_stripSep


function n = local_parseNumber(v, name)
% Accept 18, '18', "18", 'L18', 'sub-L18' ... and return 18.
if isnumeric(v)
    if ~isscalar(v)
        error('Get_staircase:badInput', '%s must be a single number.', name);
    end
    n = round(v);
    return
end
s = char(v);
tok = regexp(s, '\d+', 'match', 'once');
if isempty(tok)
    error('Get_staircase:badInput', ...
        '%s must be a number or contain one (e.g. 18, ''18'', ''L18''), got ''%s''.', name, s);
end
n = str2double(tok);
end % local_parseNumber


function s = local_letter(k)
% 1 -> A, 2 -> B, ... 26 -> Z, 27 -> AA, ...
s = '';
while k > 0
    r = mod(k-1, 26);
    s = [char('A' + r) s]; %#ok<AGROW>
    k = floor((k-1)/26);
end
end % local_letter


function s = local_timestamp()
% dd_mm_yyyy_HH_MM_SS, like the datestr() call of the original script
try
    s = char(datetime('now', 'Format', 'dd_MM_yyyy_HH_mm_ss'));
catch
    s = datestr(now, 'dd_mm_yyyy_HH_MM_SS'); %#ok<DATST,TNOW1>
end
end % local_timestamp