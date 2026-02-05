function T = tACSChallenge_buildAnalysisTable(data_path, lab_id)
%% Build a unified analysis table for H4, H5, and H6 for one lab.
%% script originally written by Anne-Fiona Griesfeller, CNRS Toulouse, in Jan 26 

% INPUT
%   data_path : path to lab folder, e.g. 'C:\Users\griesfeller\Documents\Data\L18'
%   lab_id    : numeric lab id, e.g. 18
%
% OUTPUT
%   T         : MATLAB table (also saved as .mat and .csv in data_path)

if nargin < 2 || isempty(lab_id)
    % Try to infer lab id from folder name L##
    [~, labFolder] = fileparts(data_path);
    tok = regexp(labFolder, 'L(\d+)$', 'tokens', 'once');
    if ~isempty(tok)
        lab_id = str2double(tok{1});
    else
        lab_id = NaN;
    end
end

%% Find subject folders

dd = dir(fullfile(data_path, 'sub-*'));
dd = dd([dd.isdir]);
subjs = {dd.name};

if isempty(subjs)
    error('No subject folders found in %s', data_path);
end

%% Load PhaseMetrics (for H4 d_phase)

phaseFile = fullfile(data_path, 'PhaseMetrics.mat');
hasPhase  = isfile(phaseFile);

% Store d_phase as a numeric map: phaseBS.(subj).A/B/C
phaseBS = struct();

if hasPhase
    S = load(phaseFile);
    
    if isfield(S, 'subjs') && isfield(S, 'all_bs')
        subjsPM = S.subjs;
        all_bs  = S.all_bs;

        for si = 1:numel(subjsPM)
            sid = subjsPM{si};
            if ~ischar(sid) && ~isstring(sid)
                sid = char(string(sid));
            end

            % Map rows to condition labels
            % row 1 = A, row 2 = B, row 3 = C, row 4 = Sham
            key = matlab.lang.makeValidName(sid);  % converts sub-L18_S01 -> sub_L18_S01

phaseBS.(key) = struct();
phaseBS.(key).A = all_bs(1, si);
phaseBS.(key).B = all_bs(2, si);
phaseBS.(key).C = all_bs(3, si);

        end
    else
        warning('PhaseMetrics.mat found but does not contain expected variables subjs and all_bs. d_phase will be NaN.');
        hasPhase = false;
    end
end


%% Prepare output rows

rows = {}; 
row_i = 0;

conditions_patterns = {'*_A.*','*_B.*','*_C.*','*_Sh.*', ...
                       '*_A_*','*_B_*','*_C_*','*_Sh_*'};
condLabels = {'A','B','C','Sham'};

%% Loop subjects

for s = 1:numel(subjs)

    subj = subjs{s};
    subj_beh = fullfile(data_path, subj, 'beh');

    % Skip if no beh folder
    if ~isfolder(subj_beh)
        continue
    end

    % Determine Fs 
    anyFile = dir(fullfile(subj_beh, '*.tsv'));
    if isempty(anyFile)
        continue
    end
    d0 = tACSChallenge_ImportData(fullfile(anyFile(1).folder, anyFile(1).name));
    Fs = d0.Fs;
   
    %% Read metadata and compute covariates (subject level)
    
    subj_path = fullfile(data_path, subj);
    meta_file = fullfile(subj_path, 'metadata', [subj '_Meta_Data.tsv']);

    if ~isfile(meta_file)
        error('Meta data file not found for %s at %s', subj, meta_file);
    end

    meta = readtable(meta_file, 'FileType', 'text', 'Delimiter', '\t');
    m = meta(1,:);

    % Lab code
    if ismember('LabCode', m.Properties.VariableNames)
        lab_code_val = num2str(m.LabCode(1));
    else
        lab_code_val = num2str(lab_id);
    end

    % Gender
    if ismember('Gender', m.Properties.VariableNames)
        gender_val = string(m.Gender{1});
    else
        gender_val = "";
    end

    % Age (use Age_5_y_ if present, else Age if present)
    if ismember('Age_5_y_', m.Properties.VariableNames)
        age_val = m.Age_5_y_(1);
    elseif ismember('Age', m.Properties.VariableNames)
        age_val = m.Age(1);
    else
        age_val = NaN;
    end

    % EMLA yes or no to 0 or 1
    if ismember('EMLA_y_n_', m.Properties.VariableNames)
        emla_str = lower(string(m.EMLA_y_n_{1}));
        if emla_str == "yes"
            emla_val = 1;
        elseif emla_str == "no"
            emla_val = 0;
        else
            emla_val = NaN;
        end
    elseif ismember('EMLA', m.Properties.VariableNames)
        emla_val = m.EMLA(1);
    else
        emla_val = NaN;
    end

    % IAF and delta IAF
    if ismember('IAF', m.Properties.VariableNames)
        iaf_val = m.IAF(1);
    else
        error('IAF column not found in metadata for subject %s', subj);
    end
    delta_iaf_val = abs(iaf_val - 10);

    % Hemispheric impedance difference 
    left_pre   = mean([m.Imp1O1_Oz, m.Imp1CP5_O1, m.Imp1CP5_leftEye]);
    right_pre  = mean([m.Imp1O2_Oz, m.Imp1CP6_O2, m.Imp1CP6_rightEye]);
    left_post  = mean([m.Imp2O1_Oz, m.Imp2CP5_O1, m.Imp2CP5_leftEye]);
    right_post = mean([m.Imp2O2_Oz, m.Imp2CP6_O2, m.Imp2CP6_rightEye]);

    hemi_imp_diff_val = abs(mean([left_pre, left_post]) - mean([right_pre, right_post]));

    % Intensities for A B C
    intensities = [m.Intensity_A_, m.Intensity_B_, m.Intensity_C_];

%% Parse beh filenames once (needed for Sequence inference and offline sham)

ff = dir(fullfile(subj_beh, '*.tsv'));
if isempty(ff)
    continue
end

info = repmat(struct('blockNum',NaN,'cond',""), numel(ff), 1);
for i = 1:numel(ff)
    info(i) = tsc_parse_beh_filename(ff(i).name);
end

blockNum = [info.blockNum]';
cond     = string({info.cond})';

ok = ~isnan(blockNum) & cond ~= "";
ff = ff(ok);
blockNum = blockNum(ok);
cond = cond(ok);

[blockNum, ord] = sort(blockNum);
ff = ff(ord);
cond = cond(ord);


%% Infer fixed sequence (S1/S2/S3)

seqBlocks = [2 3 4 6 7 8 10 11 12];
seq1 = ["A","C","B","C","B","A","B","A","C"];
seq2 = ["B","A","C","A","C","B","C","B","A"];
seq3 = ["C","B","A","B","A","C","A","C","B"];

sequence_val = "Unknown";

obs = strings(1, numel(seqBlocks));
present = false(1, numel(seqBlocks));

for k = 1:numel(seqBlocks)
    b = seqBlocks(k);
    idx = find(blockNum == b, 1, 'first');
    if ~isempty(idx) && ismember(cond(idx), ["A","B","C"])
        obs(k) = cond(idx);
        present(k) = true;
    end
end

nPresent = sum(present);

if nPresent >= 3
    s1 = sum(obs(present) == seq1(present));
    s2 = sum(obs(present) == seq2(present));
    s3 = sum(obs(present) == seq3(present));

    scores = [s1 s2 s3];
    best = max(scores);

    if sum(scores == best) == 1
        if best == s1
            sequence_val = "S1";
        elseif best == s2
            sequence_val = "S2";
        else
            sequence_val = "S3";
        end
    else
        warning('Sequence tie for %s (scores S1=%d, S2=%d, S3=%d with %d blocks present).', ...
            subj, s1, s2, s3, nPresent);
    end
else
    warning('Too few active blocks to infer sequence for %s (only %d present).', subj, nPresent);
end



%% ONLINE performance for H5 (A,B,C) 

for c = 1:4
    if string(condLabels{c}) == "Sham"
        continue
    end

    thisLabel = string(condLabels{c});  

    condMatch = thisLabel;  % only A, B, C

    idxFiles = find(cond == condMatch);

    % Aggregate counters across blocks
    agg_nTargets = 0;
    agg_nHits = 0;
    agg_sumRT = 0;
    agg_nHitRT = 0;

    agg_nFA = 0;
    agg_nNonTargetWindows = 0;

    if ~isempty(idxFiles)
        for ii = 1:numel(idxFiles)
            fpath = fullfile(ff(idxFiles(ii)).folder, ff(idxFiles(ii)).name);
            data = tACSChallenge_ImportData(fpath);

            [LEDLat, RespLat, nSamplesBlock, Fs_block] = extractEventsFromBlock(data);

            % Use Fs from this block
            Fs = Fs_block;

            % If no targets detected in this block, skip it
            if isempty(LEDLat)
                continue
            end

            p = tsc_perf_from_sorted(LEDLat, RespLat, Fs, nSamplesBlock);

            agg_nTargets = agg_nTargets + p.nTargets;
            agg_nHits = agg_nHits + p.nHits;

            agg_sumRT = agg_sumRT + p.sumRT_secHits;
            agg_nHitRT = agg_nHitRT + p.nHitRT;

            agg_nFA = agg_nFA + p.nFA;
            agg_nNonTargetWindows = agg_nNonTargetWindows + p.nNonTargetWindows;
        end
    end

    % Compute aggregated rates and SDT from aggregated counts
    if agg_nTargets > 0
        hitRate = agg_nHits / agg_nTargets;
    else
        hitRate = NaN;
    end

    if agg_nHitRT > 0
        meanRT = agg_sumRT / agg_nHitRT;
    else
        meanRT = NaN;
    end

    if agg_nNonTargetWindows > 0
        faRate = agg_nFA / agg_nNonTargetWindows;
    else
        faRate = NaN;
    end

    % Edge corrections and SDT 
    hitRateCorr = edgeCorrect_local(hitRate, agg_nTargets);
    faRateCorr  = edgeCorrect_local(faRate, agg_nNonTargetWindows);

    if ~isnan(hitRateCorr) && ~isnan(faRateCorr)
        dprime = norminv(hitRateCorr) - norminv(faRateCorr);
        criterion = -0.5 * (norminv(hitRateCorr) + norminv(faRateCorr));
    else
        dprime = NaN;
        criterion = NaN;
    end

    % Write row
    row_i = row_i + 1;

    rows{row_i,1}  = subj;
    rows{row_i,2}  = lab_code_val;
    rows{row_i,3}  = "Online";
    rows{row_i,4}  = thisLabel;
    rows{row_i,5}  = "";

    rows{row_i,6}  = agg_nTargets;
    rows{row_i,7}  = agg_nHits;
    rows{row_i,8}  = hitRate;
    rows{row_i,9}  = meanRT;

    rows{row_i,10} = agg_nFA;
    rows{row_i,11} = agg_nNonTargetWindows;
    rows{row_i,12} = faRate;
    rows{row_i,13} = dprime;
    rows{row_i,14} = criterion;

    % H4 d_phase only for A,B,C online
    dph = NaN;
    if hasPhase && c <= 3
        key = matlab.lang.makeValidName(subj);
        cc  = condLabels{c}; 
        if isfield(phaseBS, key) && isfield(phaseBS.(key), cc)
            dph = phaseBS.(key).(cc);
        end
    end
    rows{row_i,15} = dph;

    % Covariates
    rows{row_i,16} = iaf_val;
    rows{row_i,17} = delta_iaf_val;
    rows{row_i,18} = intensities(c);
    rows{row_i,19} = hemi_imp_diff_val;
    rows{row_i,20} = gender_val;
    rows{row_i,21} = age_val;
    rows{row_i,22} = emla_val;
    rows{row_i,23} = sequence_val;
end


    
%% OFFLINE performance for H6 (sham only, labelled pre/postX)
   
    % Identify sham blocks
    isSh = cond == "Sh";
    shamBlocks = blockNum(isSh);
    shamFiles  = ff(isSh);

    % Expect sham blocks B1,B5,B9,B13 
    % Compute performance per sham block and label according to preceding active block
    for i = 1:numel(shamBlocks)

        bnum = shamBlocks(i);
        f = fullfile(shamFiles(i).folder, shamFiles(i).name);
        data = tACSChallenge_ImportData(f);

        [LEDLat, RespLat, nSamplesBlock, Fs_block] = extractEventsFromBlock(data);
        Fs = Fs_block;

if isempty(LEDLat)
    warning('No LED targets found in sham block %d for %s. Skipping.', bnum, subj);
    continue
end

p = tsc_perf_from_sorted(LEDLat, RespLat, Fs, nSamplesBlock);


        % Determine OfflineCondition label
        offLabel = "";
if bnum == 1
    offLabel = "pre";
elseif bnum == 5
    idx4 = find(blockNum == 4, 1, 'first');
    if ~isempty(idx4); offLabel = "post" + cond(idx4); else; offLabel = "postUnknown"; end
elseif bnum == 9
    idx8 = find(blockNum == 8, 1, 'first');
    if ~isempty(idx8); offLabel = "post" + cond(idx8); else; offLabel = "postUnknown"; end
elseif bnum == 13
    idx12 = find(blockNum == 12, 1, 'first');
    if ~isempty(idx12); offLabel = "post" + cond(idx12); else; offLabel = "postUnknown"; end
else
    offLabel = "postUnknown";
end

                row_i = row_i + 1;

        rows{row_i,1}  = subj;
        rows{row_i,2}  = lab_code_val;
        rows{row_i,3}  = "Offline";
        rows{row_i,4}  = "Sham";
        rows{row_i,5}  = offLabel;
        
        rows{row_i,6}  = p.nTargets;
        rows{row_i,7}  = p.nHits;
        rows{row_i,8}  = p.hitRate;
        rows{row_i,9}  = p.meanRT_sec;

        rows{row_i,10} = p.nFA;
        rows{row_i,11} = p.nNonTargetWindows;
        rows{row_i,12} = p.faRate;
        rows{row_i,13} = p.dprime;
        rows{row_i,14} = p.criterion;

        rows{row_i,15} = NaN;          % d_phase not defined offline

        rows{row_i,16} = iaf_val;
        rows{row_i,17} = delta_iaf_val;
        rows{row_i,18} = NaN;          % StimIntensity not defined for sham blocks
        rows{row_i,19} = hemi_imp_diff_val;
        rows{row_i,20} = gender_val;
        rows{row_i,21} = age_val;
        rows{row_i,22} = emla_val;
        rows{row_i,23} = sequence_val;
        
    end

end

%% 4) Convert to table and save

T = cell2table(rows, 'VariableNames', { ...
    'SubjectID','LabID','AnalysisWindow','Condition','OfflineCondition', ...
    'nTargets','nHits','hitRate','meanRT', ...
    'nFA','nNonTargetWindows','faRate','dprime','criterion', ...
    'd_phase', ...
    'IAF','deltaIAF','StimIntensity','HemiImpDiff','Gender','Age','EMLA','Sequence'});

% categories
T.SubjectID = categorical(string(T.SubjectID));
T.LabID = categorical(string(T.LabID));
T.AnalysisWindow = categorical(string(T.AnalysisWindow));
T.Condition = categorical(string(T.Condition));
T.OfflineCondition = categorical(string(T.OfflineCondition));
T.Gender = categorical(string(T.Gender));
T.Sequence = categorical(string(T.Sequence));


% Save
outMat = fullfile(data_path, 'AnalysisTable.mat');
outCsv = fullfile(data_path, 'AnalysisTable.csv');
save(outMat, 'T');
writetable(T, outCsv);

fprintf('Saved AnalysisTable with %d rows to:\n  %s\n  %s\n', height(T), outMat, outCsv);
end

%% Helper functions 

function [LEDLat, RespLat, nSamplesBlock, Fs] = extractEventsFromBlock(data)
% Extract target onsets and response onsets in samples from imported block.

Fs = data.Fs;

% block length in samples
nSamplesBlock = size(data.L_Button, 1);

% response onsets
respOnsetsL = diff(data.L_Button); respOnsetsL(respOnsetsL < 0) = 0;
RespLatL = find(respOnsetsL > 0);

respOnsetsR = diff(data.R_Button); respOnsetsR(respOnsetsR < 0) = 0;
RespLatR = find(respOnsetsR > 0);

RespLat = sort(unique([RespLatL; RespLatR]));

% target onsets
LEDs = data.LEDs - repmat(data.LEDs(:,1), 1, size(data.LEDs,2));
LED  = max(LEDs, [], 2);
LEDdiff = diff(LED); LEDdiff(LEDdiff < 0) = 0;
LEDLat = find(LEDdiff > 0);

LEDLat = LEDLat(:);
RespLat = RespLat(:);
end

function r = edgeCorrect_local(rIn, N)
    if isnan(rIn) || N <= 0
        r = NaN;
        return
    end

    if rIn >= 1
        r = (N - 1) / N;
    elseif rIn <= 0
        r = 1 / N;
    else
        r = rIn;
    end
end

