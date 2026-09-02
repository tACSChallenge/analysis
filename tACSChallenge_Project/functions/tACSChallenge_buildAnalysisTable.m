function T = tACSChallenge_buildAnalysisTable(data_path, lab_id)
% Analysis table is built per participant
%
% Blocks containing no LED targets are skipped.

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

fprintf('- Building analysis table for lab %d (%s) -\n', lab_id, data_path);

%% Find subject folders
dd = dir(fullfile(data_path, 'sub-*'));
dd = dd([dd.isdir]);
subjs = {dd.name};

if isempty(subjs)
    error('No subject folders found in %s', data_path);
end

%% Track why subjects get skipped
skipCounts = struct( ...
    'no_beh_folder',              0, ...
    'no_tsv_files',                0, ...
    'no_parsed_beh_files',        0);

metadataMissingCount = 0;

% Counters for incomplete blocks
nZeroRespActive   = 0;   % active blocks with targets but no responses
nZeroRespSham     = 0;   % sham blocks with targets but no responses
nNoTargetBlocks   = 0;   % blocks with no LED targets

%% Load PhaseMetrics (for H4 d_phase)
phaseFile = fullfile(data_path, 'PhaseMetrics.mat');
hasPhase  = isfile(phaseFile);

% Store d_phase as a numeric map
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
            
            % row 1 = A, row 2 = B, row 3 = C, row 4 = Sham
            key = matlab.lang.makeValidName(sid);
            phaseBS.(key)   = struct();
            phaseBS.(key).A = all_bs(1, si);
            phaseBS.(key).B = all_bs(2, si);
            phaseBS.(key).C = all_bs(3, si);
        end
    else
        warning(['PhaseMetrics.mat found but does not contain expected ' ...
            'variables subjs and all_bs. d_phase will be NaN.']);
        hasPhase = false;
    end
end

%% Prepare output rows
rows  = {};
row_i = 0;

condLabels = {'A','B','C','Sham'};

%% Loop subjects
for s = 1:numel(subjs)
    subj     = subjs{s};
    subj_beh = fullfile(data_path, subj, 'beh');

    % Skip if no beh folder
    if ~isfolder(subj_beh)
        skipCounts.no_beh_folder = skipCounts.no_beh_folder + 1;
        continue
    end

    % Determine Fs
    anyFile = dir(fullfile(subj_beh, '*.tsv'));
    if isempty(anyFile)
        skipCounts.no_tsv_files = skipCounts.no_tsv_files + 1;
        continue
    end

    d0 = tACSChallenge_ImportData(fullfile(anyFile(1).folder, anyFile(1).name));
    Fs = d0.Fs;

    subj_path = fullfile(data_path, subj);
    meta_file = fullfile(subj_path, 'metadata', [subj '_Meta_Data.tsv']);

    if ~exist(meta_file, 'file')
        metadataMissingCount = metadataMissingCount + 1;
        m = table();
    else
        meta = readtable(meta_file, 'FileType', 'text', 'Delimiter', '\t');
        m = meta(1,:);
    end

    if ismember('LabCode', m.Properties.VariableNames)
        lab_code_val = canonicalizeLabCode_local(m.LabCode(1), lab_id);
    else
        lab_code_val = canonicalizeLabCode_local(lab_id, lab_id);
    end

    % Gender
    if ismember('Gender', m.Properties.VariableNames)
        gender_val = toStringScalar_local(m.Gender);
    else
        gender_val = "";
    end

    % Age
    if ismember('Age_5_y_', m.Properties.VariableNames)
        age_val = toNumericScalar_local(m.Age_5_y_(1));
    elseif ismember('Age', m.Properties.VariableNames)
        age_val = toNumericScalar_local(m.Age(1));
    else
        age_val = NaN;
    end

    if ismember('Age_5_y_', m.Properties.VariableNames) || ismember('Age', m.Properties.VariableNames)
        if isnan(age_val)
            warning(['Age value could not be parsed as numeric for subject %s ' ...
                '(lab %d). Setting Age to NaN.'], subj, lab_id);
        end
    end

    % EMLA yes or no to 0 or 1
    if ismember('EMLA_y_n_', m.Properties.VariableNames)
        emla_str = lower(toStringScalar_local(m.EMLA_y_n_));
        if emla_str == "yes"
            emla_val = 1;
        elseif emla_str == "no"
            emla_val = 0;
        else
            emla_val = NaN;
        end
    elseif ismember('EMLA', m.Properties.VariableNames)
        emla_val = toNumericScalar_local(m.EMLA(1));
    else
        emla_val = NaN;
    end

    % IAF
    if ismember('IAF', m.Properties.VariableNames)
        iaf_val = toNumericScalar_local(m.IAF(1));
    else
        iaf_val = NaN;
    end

    delta_iaf_val = abs(iaf_val - 10);

    % Hemispheric impedance difference
    % Impedances, six electrode pairs, pre and post.
    imp1_O1_Oz        = toNumericScalar_local(readMetaField_local(m, 'Imp1O1_Oz',        lab_id, subj));
    imp1_CP5_O1       = toNumericScalar_local(readMetaField_local(m, 'Imp1CP5_O1',       lab_id, subj));
    imp1_CP5_leftEye  = toNumericScalar_local(readMetaField_local(m, 'Imp1CP5_leftEye',  lab_id, subj));
    imp1_O2_Oz        = toNumericScalar_local(readMetaField_local(m, 'Imp1O2_Oz',        lab_id, subj));
    imp1_CP6_O2       = toNumericScalar_local(readMetaField_local(m, 'Imp1CP6_O2',       lab_id, subj));
    imp1_CP6_rightEye = toNumericScalar_local(readMetaField_local(m, 'Imp1CP6_rightEye', lab_id, subj));

    imp2_O1_Oz        = toNumericScalar_local(readMetaField_local(m, 'Imp2O1_Oz',        lab_id, subj));
    imp2_CP5_O1       = toNumericScalar_local(readMetaField_local(m, 'Imp2CP5_O1',       lab_id, subj));
    imp2_CP5_leftEye  = toNumericScalar_local(readMetaField_local(m, 'Imp2CP5_leftEye',  lab_id, subj));
    imp2_O2_Oz        = toNumericScalar_local(readMetaField_local(m, 'Imp2O2_Oz',        lab_id, subj));
    imp2_CP6_O2       = toNumericScalar_local(readMetaField_local(m, 'Imp2CP6_O2',       lab_id, subj));
    imp2_CP6_rightEye = toNumericScalar_local(readMetaField_local(m, 'Imp2CP6_rightEye', lab_id, subj));

    left_pre   = mean([imp1_O1_Oz, imp1_CP5_O1, imp1_CP5_leftEye]);
    right_pre  = mean([imp1_O2_Oz, imp1_CP6_O2, imp1_CP6_rightEye]);
    left_post  = mean([imp2_O1_Oz, imp2_CP5_O1, imp2_CP5_leftEye]);
    right_post = mean([imp2_O2_Oz, imp2_CP6_O2, imp2_CP6_rightEye]);

    hemi_imp_diff_val = abs(mean([left_pre, left_post]) - mean([right_pre, right_post]));

    % Intensities for A B C
    intensities = [toNumericScalar_local(readMetaField_local(m, 'Intensity_A_', lab_id, subj)), ...
                   toNumericScalar_local(readMetaField_local(m, 'Intensity_B_', lab_id, subj)), ...
                   toNumericScalar_local(readMetaField_local(m, 'Intensity_C_', lab_id, subj))];

    %% Parse beh filenames once (needed for Sequence inference and offline sham)
    ff = dir(fullfile(subj_beh, '*.tsv'));

    if isempty(ff)
        skipCounts.no_parsed_beh_files = skipCounts.no_parsed_beh_files + 1;
        continue
    end

    info = repmat(struct('blockNum',NaN,'cond',""), numel(ff), 1);
    for i = 1:numel(ff)
        info(i) = tsc_parse_beh_filename(ff(i).name);
    end

    blockNum = [info.blockNum]';
    cond     = string({info.cond})';

    ok       = ~isnan(blockNum) & cond ~= "";
    ff       = ff(ok);
    blockNum = blockNum(ok);
    cond     = cond(ok);

    if isempty(ff)
        skipCounts.no_parsed_beh_files = skipCounts.no_parsed_beh_files + 1;
        continue
    end

    [blockNum, ord] = sort(blockNum);
    ff   = ff(ord);
    cond = cond(ord);

    %% Infer fixed sequence (S1/S2/S3)
    seqBlocks = [2 3 4 6 7 8 10 11 12];
    seq1 = ["A","C","B","C","B","A","B","A","C"];
    seq2 = ["B","A","C","A","C","B","C","B","A"];
    seq3 = ["C","B","A","B","A","C","A","C","B"];

    sequence_val = "Unknown";
    obs     = strings(1, numel(seqBlocks));
    present = false(1, numel(seqBlocks));

    for k = 1:numel(seqBlocks)
        b   = seqBlocks(k);
        idx = find(blockNum == b, 1, 'first');
        if ~isempty(idx) && ismember(cond(idx), ["A","B","C"])
            obs(k)     = cond(idx);
            present(k) = true;
        end
    end

    nPresent = sum(present);

    if nPresent >= 3
        s1 = sum(obs(present) == seq1(present));
        s2 = sum(obs(present) == seq2(present));
        s3 = sum(obs(present) == seq3(present));

        scores = [s1 s2 s3];
        best   = max(scores);

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
        idxFiles  = find(cond == condMatch);

        % Aggregate counters across blocks
        agg_nTargets          = 0;
        agg_nHits             = 0;
        agg_sumRT             = 0;
        agg_nHitRT            = 0;
        agg_nFA               = 0;
        agg_nNonTargetWindows = 0;
        nBlocksUsed = 0;

        if ~isempty(idxFiles)
            for ii = 1:numel(idxFiles)
                fpath = fullfile(ff(idxFiles(ii)).folder, ff(idxFiles(ii)).name);
                data  = tACSChallenge_ImportData(fpath);

                [LEDLat, RespLat, nSamplesBlock, Fs_block] = extractEventsFromBlock(data);

                % Use Fs from this block
                Fs = Fs_block;

                % No LED targets
                if isempty(LEDLat)
                    nNoTargetBlocks = nNoTargetBlocks + 1;
                    continue
                end

                % FULL block
                p = tsc_perf_from_sorted(LEDLat, RespLat, Fs, nSamplesBlock);

                if p.nHits == 0 && p.nFA == 0
                    nZeroRespActive = nZeroRespActive + 1;
                end

                nBlocksUsed = nBlocksUsed + 1;

                agg_nTargets          = agg_nTargets          + p.nTargets;
                agg_nHits             = agg_nHits             + p.nHits;
                agg_sumRT             = agg_sumRT             + p.sumRT_secHits;
                agg_nHitRT            = agg_nHitRT            + p.nHitRT;
                agg_nFA               = agg_nFA               + p.nFA;
                agg_nNonTargetWindows = agg_nNonTargetWindows + p.nNonTargetWindows;
            end
        end

        if ~isempty(idxFiles) && nBlocksUsed < numel(idxFiles)
            fprintf(['   [reduced condition] %s condition %s: %d of %d blocks used ' ...
                '(the rest registered no responses).\n'], ...
                subj, thisLabel, nBlocksUsed, numel(idxFiles));
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
            dprime    = norminv(hitRateCorr) - norminv(faRateCorr);
            criterion = -0.5 * (norminv(hitRateCorr) + norminv(faRateCorr));
        else
            dprime    = NaN;
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

        % H4 d_phase only for A,B,C
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
    isSh       = cond == "Sh";
    shamBlocks = blockNum(isSh);
    shamFiles  = ff(isSh);

    % Expects sham blocks B1,B5,B9,B13.
    % Compute performance per sham block and label according to preceding active block.
    for i = 1:numel(shamBlocks)
        bnum = shamBlocks(i);
        f    = fullfile(shamFiles(i).folder, shamFiles(i).name);
        data = tACSChallenge_ImportData(f);

        [LEDLat, RespLat, nSamplesBlock, Fs_block] = extractEventsFromBlock(data);
        Fs = Fs_block;

        if isempty(LEDLat)
            nNoTargetBlocks = nNoTargetBlocks + 1;
            continue
        end

        p = tsc_perf_from_sorted(LEDLat, RespLat, Fs, nSamplesBlock);

        if p.nHits == 0 && p.nFA == 0
            nZeroRespSham = nZeroRespSham + 1;
        end

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

%% Reports why subjects were skipped per lab
nSkippedTotal = skipCounts.no_beh_folder + skipCounts.no_tsv_files + ...
    skipCounts.no_parsed_beh_files;

if nSkippedTotal > 0 || row_i == 0
    fprintf('Lab %d (%s): %d of %d subject folders skipped entirely.\n', ...
        lab_id, data_path, nSkippedTotal, numel(subjs));
    fprintf('  no beh folder:            %d\n', skipCounts.no_beh_folder);
    fprintf('  no .tsv files in beh:     %d\n', skipCounts.no_tsv_files);
    fprintf('  no parsed beh files:      %d\n', skipCounts.no_parsed_beh_files);
end

if metadataMissingCount > 0
    fprintf(['Lab %d (%s): %d of %d subjects had no metadata file. ' ...
        'They were RETAINED (Condition/dprime/criterion/meanRT/LabID/Sequence ' ...
        'do not require metadata); IAF, Age, Gender, EMLA, HemiImpDiff, and ' ...
        'Intensity are NaN for them.\n'], ...
        lab_id, data_path, metadataMissingCount, numel(subjs));
end

%% Block summary for this lab
fprintf('Lab %d block summary:\n', lab_id);
fprintf('  blocks with no LED targets       : %d\n', nNoTargetBlocks);
fprintf('  active blocks with no responses  : %d\n', nZeroRespActive);
fprintf('  sham blocks with no responses    : %d\n', nZeroRespSham);

%% 4) Convert to table and save
varNames = {'SubjectID','LabID','AnalysisWindow','Condition','OfflineCondition', ...
    'nTargets','nHits','hitRate','meanRT', ...
    'nFA','nNonTargetWindows','faRate','dprime','criterion', ...
    'd_phase', ...
    'IAF','deltaIAF','StimIntensity','HemiImpDiff','Gender','Age','EMLA','Sequence'};

if row_i == 0
    % Typed empty table so concentration still works
    warning('tACSChallenge_buildAnalysisTable:NoRows', ...
        ['No usable rows were built for lab %d (%s). See the skip summary ' ...
         'printed above for the reason. Returning an empty table so other ' ...
         'labs can still be processed.'], lab_id, data_path);

    T = table( ...
        strings(0,1), strings(0,1), strings(0,1), strings(0,1), strings(0,1), ...
        zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
        zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
        zeros(0,1), ...
        zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), strings(0,1), zeros(0,1), zeros(0,1), strings(0,1), ...
        'VariableNames', varNames);
else
    T = cell2table(rows, 'VariableNames', varNames);
end

% categories
T.SubjectID        = categorical(string(T.SubjectID));
T.LabID            = categorical(string(T.LabID));
T.AnalysisWindow   = categorical(string(T.AnalysisWindow));
T.Condition        = categorical(string(T.Condition));
T.OfflineCondition = categorical(string(T.OfflineCondition));
T.Gender           = categorical(string(T.Gender));
T.Sequence         = categorical(string(T.Sequence));

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

% target onsets, detects LED increase and decrease
LEDs = data.LEDs - repmat(data.LEDs(:,1), 1, size(data.LEDs,2));
LED  = max(abs(LEDs), [], 2);

baseline = median(LED);
isTarget = LED ~= baseline;

LEDLat = find(diff(isTarget) == 1);

LEDLat  = LEDLat(:);
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

function v = toNumericScalar_local(x)
% Force a metadata field to a numeric scalar; NaN if unparseable.
    if iscell(x)
        x = x{1};
    end

    if isnumeric(x)
        v = double(x);
    elseif islogical(x)
        v = double(x);
    else
        v = str2double(string(x));
    end
end

function v = toStringScalar_local(col)
% Scalar string from a cell, string or categorical column.
    v = string(col(1));
end

function v = readMetaField_local(m, colName, lab_id, subj)
% Metadata field by name; NaN if the column is absent.
    if ismember(colName, m.Properties.VariableNames)
        v = m.(colName)(1);
    else
        warning('Metadata column "%s" not found for subject %s (lab %d). Setting %s to NaN.', ...
            colName, subj, lab_id, colName);
        v = NaN;
    end
end

function code = canonicalizeLabCode_local(rawVal, lab_id)
% Normalise a lab code to "L##" so one lab is not split into two levels.
    raw    = strtrim(string(rawVal));
    digits = regexp(raw, '\d+', 'match', 'once');

    if strlength(digits) == 0
        if ~isnan(lab_id)
            digits = string(lab_id);
        else
            code = raw;
            return
        end
    end

    code = "L" + digits;
end
