function OUT = tACSChallenge_runGlobalModels(root_path, lab_ids, varargin)
%% Run H4, H5, H6 mixed models across labs using unified AnalysisTable.
%% script originally written by Anne-Fiona Griesfeller, CerCo, Toulouse

p = inputParser;
p.addParameter('RunH4', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('RunH5', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('RunH6', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('MakePlots', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('OutDir', '', @(x) ischar(x) || isstring(x));
p.addParameter('forceRecompute', false, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});

RunH4 = logical(p.Results.RunH4);
RunH5 = logical(p.Results.RunH5);
RunH6 = logical(p.Results.RunH6);
MakePlots = logical(p.Results.MakePlots);
OutDir = string(p.Results.OutDir);
forceRecompute = logical(p.Results.forceRecompute);

if OutDir == ""
    OutDir = string(root_path);
end

% Create base output folder
if ~isfolder(OutDir)
    mkdir(OutDir);
end

% Create a subfolder per hypothesis for results.
outH4 = fullfile(OutDir, 'H4');
outH5 = fullfile(OutDir, 'H5');
outH6 = fullfile(OutDir, 'H6');

if ~isfolder(outH4); mkdir(outH4); end
if ~isfolder(outH5); mkdir(outH5); end
if ~isfolder(outH6); mkdir(outH6); end


OUT = struct();


%% Build or load per-lab tables, then combine them into one

Tall  = table();

for i = 1:numel(lab_ids)
    lab_id = lab_ids(i);
    lab_path = fullfile(root_path, sprintf('L%02d', lab_id));

    if ~isfolder(lab_path)
        warning('Lab folder not found: %s (skipping)', lab_path);
        continue
    end

    analysisFile = fullfile(lab_path, 'AnalysisTable.mat');

    % If lab already has a saved AnalysisTable.mat, and forcerecompute = false
    % it is loaded rather than rebuild again
    if isfile(analysisFile) && ~forceRecompute
        fprintf('Lab %d: found existing AnalysisTable.mat, loading (skipping recomputation).\n', lab_id);
        S = load(analysisFile, 'T');
        Tlab = S.T;
    else
        Tlab = tACSChallenge_buildAnalysisTable(lab_path, lab_id);
    end

    Tall = [Tall; Tlab];

    fprintf('Added lab %d (%d rows). Combined now %d rows.\n', lab_id, height(Tlab), height(Tall));
end

OUT.T_all     = Tall;

% Save combined table at root
outMat = fullfile(OutDir, 'AnalysisTable_allLabs.mat');
outCsv = fullfile(OutDir, 'AnalysisTable_allLabs.csv');
save(outMat, 'Tall');
writetable(Tall, outCsv);
fprintf('Saved combined table to:\n  %s\n  %s\n', outMat, outCsv);
%% Ensure variable types and references

Tall.SubjectID       = categorical(string(Tall.SubjectID));
Tall.LabID           = categorical(string(Tall.LabID));
Tall.AnalysisWindow  = categorical(string(Tall.AnalysisWindow));
Tall.Condition       = categorical(string(Tall.Condition));
Tall.OfflineCondition= categorical(string(Tall.OfflineCondition));
% Standardise Gender labels: "F"/"f"/"female" -> female, "M"/"m"/"male" -> male
g = lower(strtrim(string(Tall.Gender)));
g(g == "f") = "female";
g(g == "m") = "male";
Tall.Gender          = categorical(g);
Tall.Sequence        = categorical(string(Tall.Sequence));

% Reference levels
if any(categories(Tall.Condition) == "A")
    Tall.Condition = reordercats(Tall.Condition, ["A"; setdiff(categories(Tall.Condition), "A", 'stable')]);
end
if any(categories(Tall.OfflineCondition) == "pre")
    Tall.OfflineCondition = reordercats(Tall.OfflineCondition, ["pre"; setdiff(categories(Tall.OfflineCondition), "pre", 'stable')]);
end

    % Drop fixed-effect covariates that have no variability
    function [fixedTermsKept, dropped] = keepVaryingTerms(T, fixedTerms)
        fixedTermsKept = string(fixedTerms);
        dropped = strings(0,1);

        for k = 1:numel(fixedTermsKept)
            v = fixedTermsKept(k);

            if ~ismember(v, string(T.Properties.VariableNames))
                fixedTermsKept(k) = "";
                dropped(end+1,1) = v + " (missing column)";
                continue
            end

            x = T.(char(v));

            if iscategorical(x) || isstring(x) || ischar(x)
                x = categorical(string(x));
                x = x(~isundefined(x));
                nU = numel(categories(removecats(x)));
            else
                x = x(~isnan(x));
                nU = numel(unique(x));
            end

            if nU < 2
                fixedTermsKept(k) = "";
                dropped(end+1,1) = v + " (no variability)";
            end
        end

        fixedTermsKept = fixedTermsKept(fixedTermsKept ~= "");
    end


    % Choose random effects based on number of levels
    function re = randomPart(T, includeSeqAsRandom)
        reParts = strings(0,1);

        if numel(categories(T.SubjectID)) > 1
            reParts(end+1,1) = "(1|SubjectID)";
        end

        if numel(categories(T.LabID)) > 1
            reParts(end+1,1) = "(1|LabID)";
        end

        if includeSeqAsRandom && ismember('Sequence', T.Properties.VariableNames)
            if numel(categories(T.Sequence)) > 1
                reParts(end+1,1) = "(1|Sequence)";
            end
        end

        if isempty(reParts)
            re = "";
        else
            re = " + " + strjoin(cellstr(reParts), " + ");
        end
    end


%% H4 model (online conditions, outcome d_phase)

if RunH4
    fprintf('\nFitting H4 global mixed model (d_phase)...\n');
    reportH4 = fullfile(outH4, 'H4_model_report.txt');
    if isfile(reportH4)
        delete(reportH4);
    end

    T4 = Tall(Tall.AnalysisWindow == "Online" & ismember(string(Tall.Condition), ["A","B","C"]), :);
    T4 = rmmissing(T4, 'DataVariables', {'d_phase','deltaIAF'});

    % Drop phantom categories after subsetting
    vn = T4.Properties.VariableNames;
    for ii = 1:numel(vn)
        if iscategorical(T4.(vn{ii}))
            T4.(vn{ii}) = removecats(T4.(vn{ii}));
        end
    end

    fixedWanted = {'Condition','deltaIAF','StimIntensity','HemiImpDiff','Gender','Age','EMLA','IAF'};
    [fixedKept, dropped] = keepVaryingTerms(T4, fixedWanted);

    if ~isempty(dropped)
        fprintf('H4 dropping covariates: %s\n', strjoin(cellstr(dropped), ', '));
    end

    base = "d_phase ~ 1 + Condition*deltaIAF";

    addMain = setdiff(string(fixedKept), ["Condition","deltaIAF"], 'stable');
    if ~isempty(addMain)
        base = base + " + " + strjoin(cellstr(addMain), " + ");
    end

    re = randomPart(T4, true);
    formula = char(base + re);
    OUT.H4.formula_used = formula;

    try
        OUT.H4.lme = fitlme(T4, formula);
        disp(OUT.H4.lme);
        OUT.H4.anova = anova(OUT.H4.lme, 'DFMethod','satterthwaite');
        disp(OUT.H4.anova);
        appendModelToReport(reportH4, "H4 randomSeq", OUT.H4.lme, OUT.H4.anova, OUT.H4.formula_used);

        % Bayes factor for the deltaIAF by Condition interaction.
        baseRed4 = "d_phase ~ 1 + Condition + deltaIAF";
        if ~isempty(addMain)
            baseRed4 = baseRed4 + " + " + strjoin(cellstr(addMain), " + ");
        end
        OUT.H4.bayes = tACS_bayesFactorBIC(OUT.H4.lme, T4, char(baseRed4 + re));
        reportBayes(reportH4, "H4 deltaIAF x Condition interaction", OUT.H4.bayes);
    catch ME
        warning('H4 model failed: %s', ME.message);
        OUT.H4.lme = [];
        OUT.H4.anova = [];
        OUT.H4.bayes = [];
    end

    % H4 variant with Sequence as fixed effect
    try
        if ismember("Sequence", string(T4.Properties.VariableNames)) && numel(categories(T4.Sequence)) > 1

            fixedWanted2 = [cellstr(string(fixedKept)) {'Sequence'}];
            [fixedKept2, dropped2] = keepVaryingTerms(T4, fixedWanted2);
            if ~isempty(dropped2)
                fprintf('H4 fixedSeq dropping covariates: %s\n', strjoin(cellstr(dropped2), ', '));
            end

            base2 = "d_phase ~ 1 + Condition*deltaIAF";

            addMain2 = setdiff(string(fixedKept2), ["Condition","deltaIAF"], 'stable');
            if ~isempty(addMain2)
                base2 = base2 + " + " + strjoin(cellstr(addMain2), " + ");
            end

            re2 = randomPart(T4, false);
            formula2 = char(base2 + re2);

            OUT.H4.fixedSeq.formula_used = formula2;
            OUT.H4.fixedSeq.lme = fitlme(T4, formula2);
            disp(OUT.H4.fixedSeq.lme);
            OUT.H4.fixedSeq.anova = anova(OUT.H4.fixedSeq.lme, 'DFMethod','satterthwaite');
            disp(OUT.H4.fixedSeq.anova);
            appendModelToReport(reportH4, "H4 fixedSeq", OUT.H4.fixedSeq.lme, OUT.H4.fixedSeq.anova, OUT.H4.fixedSeq.formula_used);
        else
            OUT.H4.fixedSeq.lme = [];
            OUT.H4.fixedSeq.anova = [];
        end
    catch ME
        warning('H4 fixedSeq model failed: %s', ME.message);
        OUT.H4.fixedSeq.lme = [];
        OUT.H4.fixedSeq.anova = [];
    end

end


%% H5 models

if RunH5
    fprintf('\nFitting H5 global mixed models (dprime, criterion, meanRT)...\n');
    reportH5 = fullfile(outH5, 'H5_model_report.txt');
    if isfile(reportH5)
        delete(reportH5);
    end
    
    T5online = Tall(Tall.AnalysisWindow == "Online" & ismember(string(Tall.Condition), ["A","B","C"]), :);

    shamRows = Tall(Tall.AnalysisWindow == "Offline" & string(Tall.Condition) == "Sham", :);
    shamAgg  = table();
    subjList = unique(shamRows.SubjectID);
    for s = 1:numel(subjList)
        rows   = shamRows(shamRows.SubjectID == subjList(s), :);
        newRow = rows(1, :);
        newRow.dprime        = mean(rows.dprime,      'omitnan');
        newRow.criterion     = mean(rows.criterion,   'omitnan');
        newRow.meanRT        = mean(rows.meanRT,      'omitnan');
        newRow.HemiImpDiff   = mean(rows.HemiImpDiff, 'omitnan');
        newRow.StimIntensity = 0;
        shamAgg = [shamAgg; newRow];
    end

    T5 = [T5online; shamAgg];

    % Drop phantom categories
    vn = T5.Properties.VariableNames;
    for ii = 1:numel(vn)
        if iscategorical(T5.(vn{ii}))
            T5.(vn{ii}) = removecats(T5.(vn{ii}));
        end
    end

    % Make Sham the reference level
    if any(categories(T5.Condition) == "Sham")
        T5.Condition = reordercats(T5.Condition, ["Sham"; setdiff(categories(T5.Condition), "Sham", 'stable')]);
    end

    fixedWanted = {'Condition','StimIntensity','HemiImpDiff','Gender','Age','EMLA','IAF'};
    [fixedKept, dropped] = keepVaryingTerms(T5, fixedWanted);
    if ~isempty(dropped)
        fprintf('H5 dropping covariates: %s\n', strjoin(cellstr(dropped), ', '));
    end

    base = " ~ 1 + " + strjoin(cellstr(string(fixedKept)), " + ");
    re_rand = randomPart(T5, true);

    % Reduced fixed-effect set for the Bayes factor
    fixedRed5 = setdiff(string(fixedKept), "Condition", 'stable');
    if isempty(fixedRed5)
        baseRed5 = " ~ 1";
    else
        baseRed5 = " ~ 1 + " + strjoin(cellstr(fixedRed5), " + ");
    end

    % dprime
    try
        T5dp = rmmissing(T5, 'DataVariables', {'dprime'});
        formula = char("dprime" + base + re_rand);
        OUT.H5.dprime.randomSeq.formula_used = formula;
        OUT.H5.dprime.lme_randomSeq = fitlme(T5dp, formula);
        disp(OUT.H5.dprime.lme_randomSeq);
        OUT.H5.dprime.randomSeq.anova = anova(OUT.H5.dprime.lme_randomSeq, 'DFMethod','satterthwaite');
        disp(OUT.H5.dprime.randomSeq.anova);
        appendModelToReport(reportH5, "H5 dprime randomSeq", OUT.H5.dprime.lme_randomSeq, OUT.H5.dprime.randomSeq.anova, OUT.H5.dprime.randomSeq.formula_used);

        OUT.H5.dprime.bayes = tACS_bayesFactorBIC(OUT.H5.dprime.lme_randomSeq, T5dp, ...
            char("dprime" + baseRed5 + re_rand));
        reportBayes(reportH5, "H5 dprime: main effect of Condition", OUT.H5.dprime.bayes);
    catch ME
        warning('H5 dprime randomSeq model failed: %s', ME.message);
        OUT.H5.dprime.lme_randomSeq = [];
        OUT.H5.dprime.randomSeq.anova = [];
        OUT.H5.dprime.bayes = [];
    end

    % criterion
    try
        T5c = rmmissing(T5, 'DataVariables', {'criterion'});
        formula = char("criterion" + base + re_rand);
        OUT.H5.criterion.randomSeq.formula_used = formula;
        OUT.H5.criterion.lme_randomSeq = fitlme(T5c, formula);
        disp(OUT.H5.criterion.lme_randomSeq);
        OUT.H5.criterion.randomSeq.anova = anova(OUT.H5.criterion.lme_randomSeq, 'DFMethod','satterthwaite');
        disp(OUT.H5.criterion.randomSeq.anova);
        appendModelToReport(reportH5, "H5 criterion randomSeq", OUT.H5.criterion.lme_randomSeq, OUT.H5.criterion.randomSeq.anova, OUT.H5.criterion.randomSeq.formula_used);

        OUT.H5.criterion.bayes = tACS_bayesFactorBIC(OUT.H5.criterion.lme_randomSeq, T5c, ...
            char("criterion" + baseRed5 + re_rand));
        reportBayes(reportH5, "H5 criterion: main effect of Condition", OUT.H5.criterion.bayes);
    catch ME
        warning('H5 criterion randomSeq model failed: %s', ME.message);
        OUT.H5.criterion.lme_randomSeq = [];
        OUT.H5.criterion.randomSeq.anova = [];
        OUT.H5.criterion.bayes = [];
    end

    % meanRT
    try
        T5rt = rmmissing(T5, 'DataVariables', {'meanRT'});
        formula = char("meanRT" + base + re_rand);
        OUT.H5.meanRT.randomSeq.formula_used = formula;
        OUT.H5.meanRT.lme_randomSeq = fitlme(T5rt, formula);
        disp(OUT.H5.meanRT.lme_randomSeq);
        OUT.H5.meanRT.randomSeq.anova = anova(OUT.H5.meanRT.lme_randomSeq, 'DFMethod','satterthwaite');
        disp(OUT.H5.meanRT.randomSeq.anova);
        appendModelToReport(reportH5, "H5 meanRT randomSeq", OUT.H5.meanRT.lme_randomSeq, OUT.H5.meanRT.randomSeq.anova, OUT.H5.meanRT.randomSeq.formula_used);

        OUT.H5.meanRT.bayes = tACS_bayesFactorBIC(OUT.H5.meanRT.lme_randomSeq, T5rt, ...
            char("meanRT" + baseRed5 + re_rand));
        reportBayes(reportH5, "H5 meanRT: main effect of Condition", OUT.H5.meanRT.bayes);
    catch ME
        warning('H5 meanRT randomSeq model failed: %s', ME.message);
        OUT.H5.meanRT.lme_randomSeq = [];
        OUT.H5.meanRT.randomSeq.anova = [];
        OUT.H5.meanRT.bayes = [];
    end

    % Sequence as fixed effect instead of random
    if ismember("Sequence", string(T5.Properties.VariableNames)) && numel(categories(T5.Sequence)) > 1

        fixedWanted2 = [cellstr(string(fixedKept)) {'Sequence'}];
        [fixedKept2, dropped2] = keepVaryingTerms(T5, fixedWanted2);
        if ~isempty(dropped2)
            fprintf('H5 fixedSeq dropping covariates: %s\n', strjoin(cellstr(dropped2), ', '));
        end

        base2 = " ~ 1 + " + strjoin(cellstr(string(fixedKept2)), " + ");
        re_fixed = randomPart(T5, false);

        % dprime
        try
            T5dp = rmmissing(T5, 'DataVariables', {'dprime'});
            formula = char("dprime" + base2 + re_fixed);
            OUT.H5.dprime.fixedSeq.formula_used = formula;
            OUT.H5.dprime.lme_fixedSeq = fitlme(T5dp, formula);
            disp(OUT.H5.dprime.lme_fixedSeq);
            OUT.H5.dprime.fixedSeq.anova = anova(OUT.H5.dprime.lme_fixedSeq, 'DFMethod','satterthwaite');
            disp(OUT.H5.dprime.fixedSeq.anova);
            appendModelToReport(reportH5, "H5 dprime fixedSeq", OUT.H5.dprime.lme_fixedSeq, OUT.H5.dprime.fixedSeq.anova, OUT.H5.dprime.fixedSeq.formula_used);
        catch ME
            warning('H5 dprime fixedSeq model failed: %s', ME.message);
            OUT.H5.dprime.lme_fixedSeq = [];
            OUT.H5.dprime.fixedSeq.anova = [];
        end

        % criterion
        try
            T5c = rmmissing(T5, 'DataVariables', {'criterion'});
            formula = char("criterion" + base2 + re_fixed);
            OUT.H5.criterion.fixedSeq.formula_used = formula;
            OUT.H5.criterion.lme_fixedSeq = fitlme(T5c, formula);
            disp(OUT.H5.criterion.lme_fixedSeq);
            OUT.H5.criterion.fixedSeq.anova = anova(OUT.H5.criterion.lme_fixedSeq, 'DFMethod','satterthwaite');
            disp(OUT.H5.criterion.fixedSeq.anova);
            appendModelToReport(reportH5, "H5 criterion fixedSeq", OUT.H5.criterion.lme_fixedSeq, OUT.H5.criterion.fixedSeq.anova, OUT.H5.criterion.fixedSeq.formula_used);
        catch ME
            warning('H5 criterion fixedSeq model failed: %s', ME.message);
            OUT.H5.criterion.lme_fixedSeq = [];
            OUT.H5.criterion.fixedSeq.anova = [];
        end

        % meanRT
        try
            T5rt = rmmissing(T5, 'DataVariables', {'meanRT'});
            formula = char("meanRT" + base2 + re_fixed);
            OUT.H5.meanRT.fixedSeq.formula_used = formula;
            OUT.H5.meanRT.lme_fixedSeq = fitlme(T5rt, formula);
            disp(OUT.H5.meanRT.lme_fixedSeq);
            OUT.H5.meanRT.fixedSeq.anova = anova(OUT.H5.meanRT.lme_fixedSeq, 'DFMethod','satterthwaite');
            disp(OUT.H5.meanRT.fixedSeq.anova);
            appendModelToReport(reportH5, "H5 meanRT fixedSeq", OUT.H5.meanRT.lme_fixedSeq, OUT.H5.meanRT.fixedSeq.anova, OUT.H5.meanRT.fixedSeq.formula_used);
        catch ME
            warning('H5 meanRT fixedSeq model failed: %s', ME.message);
            OUT.H5.meanRT.lme_fixedSeq = [];
            OUT.H5.meanRT.fixedSeq.anova = [];
        end

    else
        OUT.H5.dprime.lme_fixedSeq = [];
        OUT.H5.criterion.lme_fixedSeq = [];
        OUT.H5.meanRT.lme_fixedSeq = [];
    end
end


%% H6 aftereffects model (offline sham blocks, outcomes dprime, criterion, meanRT)

if RunH6
    fprintf('\nFitting H6 global mixed models (offline sham)...\n');
    reportH6 = fullfile(outH6, 'H6_model_report.txt');
    if isfile(reportH6)
        delete(reportH6);
    end

    T6 = Tall(Tall.AnalysisWindow == "Offline" & string(Tall.Condition) == "Sham", :);
    T6 = T6(~isundefined(T6.OfflineCondition), :);

    % Drop phantom categories
    vn = T6.Properties.VariableNames;
    for ii = 1:numel(vn)
        if iscategorical(T6.(vn{ii}))
            T6.(vn{ii}) = removecats(T6.(vn{ii}));
        end
    end

    fixedWanted = {'OfflineCondition','StimIntensity','HemiImpDiff','Gender','Age','EMLA','IAF'};
    [fixedKept, dropped] = keepVaryingTerms(T6, fixedWanted);
    if ~isempty(dropped)
        fprintf('H6 dropping covariates: %s\n', strjoin(cellstr(dropped), ', '));
    end

    base = " ~ 1 + " + strjoin(cellstr(string(fixedKept)), " + ");
    re_rand = randomPart(T6, true);

    % Reduced fixed-effect set for the Bayes factor
    fixedRed6 = setdiff(string(fixedKept), "OfflineCondition", 'stable');
    if isempty(fixedRed6)
        baseRed6 = " ~ 1";
    else
        baseRed6 = " ~ 1 + " + strjoin(cellstr(fixedRed6), " + ");
    end

    % dprime
    try
        T6dp = rmmissing(T6, 'DataVariables', {'dprime'});
        formula = char("dprime" + base + re_rand);
        OUT.H6.dprime.randomSeq.formula_used = formula;
        OUT.H6.dprime.lme_randomSeq = fitlme(T6dp, formula);
        disp(OUT.H6.dprime.lme_randomSeq);
        OUT.H6.dprime.randomSeq.anova = anova(OUT.H6.dprime.lme_randomSeq, 'DFMethod','satterthwaite');
        disp(OUT.H6.dprime.randomSeq.anova);
        appendModelToReport(reportH6, "H6 dprime randomSeq", OUT.H6.dprime.lme_randomSeq, OUT.H6.dprime.randomSeq.anova, OUT.H6.dprime.randomSeq.formula_used);

        OUT.H6.dprime.bayes = tACS_bayesFactorBIC(OUT.H6.dprime.lme_randomSeq, T6dp, ...
            char("dprime" + baseRed6 + re_rand));
        reportBayes(reportH6, "H6 dprime: main effect of OfflineCondition", OUT.H6.dprime.bayes);
    catch ME
        warning('H6 dprime randomSeq model failed: %s', ME.message);
        OUT.H6.dprime.lme_randomSeq = [];
        OUT.H6.dprime.randomSeq.anova = [];
        OUT.H6.dprime.bayes = [];
    end

    % criterion
    try
        T6c = rmmissing(T6, 'DataVariables', {'criterion'});
        formula = char("criterion" + base + re_rand);
        OUT.H6.criterion.randomSeq.formula_used = formula;
        OUT.H6.criterion.lme_randomSeq = fitlme(T6c, formula);
        disp(OUT.H6.criterion.lme_randomSeq);
        OUT.H6.criterion.randomSeq.anova = anova(OUT.H6.criterion.lme_randomSeq, 'DFMethod','satterthwaite');
        disp(OUT.H6.criterion.randomSeq.anova);
        appendModelToReport(reportH6, "H6 criterion randomSeq", OUT.H6.criterion.lme_randomSeq, OUT.H6.criterion.randomSeq.anova, OUT.H6.criterion.randomSeq.formula_used);

        OUT.H6.criterion.bayes = tACS_bayesFactorBIC(OUT.H6.criterion.lme_randomSeq, T6c, ...
            char("criterion" + baseRed6 + re_rand));
        reportBayes(reportH6, "H6 criterion: main effect of OfflineCondition", OUT.H6.criterion.bayes);
    catch ME
        warning('H6 criterion randomSeq model failed: %s', ME.message);
        OUT.H6.criterion.lme_randomSeq = [];
        OUT.H6.criterion.randomSeq.anova = [];
        OUT.H6.criterion.bayes = [];
    end

    % meanRT
    try
        T6rt = rmmissing(T6, 'DataVariables', {'meanRT'});
        formula = char("meanRT" + base + re_rand);
        OUT.H6.meanRT.randomSeq.formula_used = formula;
        OUT.H6.meanRT.lme_randomSeq = fitlme(T6rt, formula);
        disp(OUT.H6.meanRT.lme_randomSeq);
        OUT.H6.meanRT.randomSeq.anova = anova(OUT.H6.meanRT.lme_randomSeq, 'DFMethod','satterthwaite');
        disp(OUT.H6.meanRT.randomSeq.anova);
        appendModelToReport(reportH6, "H6 meanRT randomSeq", OUT.H6.meanRT.lme_randomSeq, OUT.H6.meanRT.randomSeq.anova, OUT.H6.meanRT.randomSeq.formula_used);

        OUT.H6.meanRT.bayes = tACS_bayesFactorBIC(OUT.H6.meanRT.lme_randomSeq, T6rt, ...
            char("meanRT" + baseRed6 + re_rand));
        reportBayes(reportH6, "H6 meanRT: main effect of OfflineCondition", OUT.H6.meanRT.bayes);
    catch ME
        warning('H6 meanRT randomSeq model failed: %s', ME.message);
        OUT.H6.meanRT.lme_randomSeq = [];
        OUT.H6.meanRT.randomSeq.anova = [];
        OUT.H6.meanRT.bayes = [];
    end

    % Sequence as fixed effect instead of random
    if ismember("Sequence", string(T6.Properties.VariableNames)) && numel(categories(T6.Sequence)) > 1
        fixedWanted2 = [cellstr(string(fixedKept)) {'Sequence'}];
        [fixedKept2, dropped2] = keepVaryingTerms(T6, fixedWanted2);
        if ~isempty(dropped2)
            fprintf('H6 fixedSeq dropping covariates: %s\n', strjoin(cellstr(dropped2), ', '));
        end
        base2 = " ~ 1 + " + strjoin(cellstr(string(fixedKept2)), " + ");
        re_fixed = randomPart(T6, false);

        % dprime
        try
            T6dp = rmmissing(T6, 'DataVariables', {'dprime'});
            formula = char("dprime" + base2 + re_fixed);
            OUT.H6.dprime.fixedSeq.formula_used = formula;
            OUT.H6.dprime.lme_fixedSeq = fitlme(T6dp, formula);
            disp(OUT.H6.dprime.lme_fixedSeq);
            OUT.H6.dprime.fixedSeq.anova = anova(OUT.H6.dprime.lme_fixedSeq, 'DFMethod','satterthwaite');
            disp(OUT.H6.dprime.fixedSeq.anova);
            appendModelToReport(reportH6, "H6 dprime fixedSeq", OUT.H6.dprime.lme_fixedSeq, OUT.H6.dprime.fixedSeq.anova, OUT.H6.dprime.fixedSeq.formula_used);
        catch ME
            warning('H6 dprime fixedSeq model failed: %s', ME.message);
            OUT.H6.dprime.lme_fixedSeq = [];
            OUT.H6.dprime.fixedSeq.anova = [];
        end

        % criterion
        try
            T6c = rmmissing(T6, 'DataVariables', {'criterion'});
            formula = char("criterion" + base2 + re_fixed);
            OUT.H6.criterion.fixedSeq.formula_used = formula;
            OUT.H6.criterion.lme_fixedSeq = fitlme(T6c, formula);
            disp(OUT.H6.criterion.lme_fixedSeq);
            OUT.H6.criterion.fixedSeq.anova = anova(OUT.H6.criterion.lme_fixedSeq, 'DFMethod','satterthwaite');
            disp(OUT.H6.criterion.fixedSeq.anova);
            appendModelToReport(reportH6, "H6 criterion fixedSeq", OUT.H6.criterion.lme_fixedSeq, OUT.H6.criterion.fixedSeq.anova, OUT.H6.criterion.fixedSeq.formula_used);
        catch ME
            warning('H6 criterion fixedSeq model failed: %s', ME.message);
            OUT.H6.criterion.lme_fixedSeq = [];
            OUT.H6.criterion.fixedSeq.anova = [];
        end

        % meanRT
        try
            T6rt = rmmissing(T6, 'DataVariables', {'meanRT'});
            formula = char("meanRT" + base2 + re_fixed);
            OUT.H6.meanRT.fixedSeq.formula_used = formula;
            OUT.H6.meanRT.lme_fixedSeq = fitlme(T6rt, formula);
            disp(OUT.H6.meanRT.lme_fixedSeq);
            OUT.H6.meanRT.fixedSeq.anova = anova(OUT.H6.meanRT.lme_fixedSeq, 'DFMethod','satterthwaite');
            disp(OUT.H6.meanRT.fixedSeq.anova);
            appendModelToReport(reportH6, "H6 meanRT fixedSeq", OUT.H6.meanRT.lme_fixedSeq, OUT.H6.meanRT.fixedSeq.anova, OUT.H6.meanRT.fixedSeq.formula_used);
        catch ME
            warning('H6 meanRT fixedSeq model failed: %s', ME.message);
            OUT.H6.meanRT.lme_fixedSeq = [];
            OUT.H6.meanRT.fixedSeq.anova = [];
        end
    else
        OUT.H6.dprime.lme_fixedSeq = [];
        OUT.H6.criterion.lme_fixedSeq = [];
        OUT.H6.meanRT.lme_fixedSeq = [];
    end

end


%% Plots (H4, H5, H6)

if MakePlots
    fprintf('\nCreating plots for H4, H5, H6...\n');

    % H4: Scatter plots d_phase vs deltaIAF by condition (A/B/C)
   
    try
        T4p = Tall(Tall.AnalysisWindow == "Online" & ismember(string(Tall.Condition), ["A","B","C"]), :);
        T4p = T4p(~isnan(T4p.d_phase) & ~isnan(T4p.deltaIAF), :);

        if height(T4p) > 0
            fprintf('Generating global H4 scatter plots...\n');
            [condCols, condNames] = tACS_condColours(false);
            condKeys = {'A', 'B', 'C'};

            xAll = T4p.deltaIAF;
            yAll = T4p.d_phase;
            xPad = 0.04 * range(xAll);
            yPadLow  = 0.06 * range(yAll);
            yPadHigh = 0.16 * range(yAll);
            sharedX = [min(xAll) - xPad,     max(xAll) + xPad];
            sharedY = [min(yAll) - yPadLow,  max(yAll) + yPadHigh];

            letters    = {'A', 'B', 'C'};
            letterXoff = [-0.22, -0.10, -0.10];

            fig_h4 = figure('Color', 'w', 'Units', 'centimeters', ...
                'Position', [1 1 16.6 8.0]);
            set(fig_h4, 'Name', 'Global H4: d_phase vs deltaIAF by condition', ...
                'NumberTitle', 'off');
            tl4 = tiledlayout(fig_h4, 1, numel(condKeys), ...
                'TileSpacing', 'compact', 'Padding', 'loose');

            for i = 1:numel(condKeys)
                idx = string(T4p.Condition) == condKeys{i};

                ax4 = nexttile(tl4);
                tACS_scatterPanel(ax4, T4p.deltaIAF(idx), T4p.d_phase(idx), ...
                    condCols(i,:), condNames{i}, ...
                    '\DeltaIAF (Hz)', '{\itd}_{phase}', ...
                    'XLim', sharedX, 'YLim', sharedY, ...
                    'ShowYLabel', i == 1, ...
                    'PanelLetter', letters{i}, ...
                    'LetterXOffset', letterXoff(i));
            end

            saveas(fig_h4, fullfile(outH4, 'H4_scatter_dphase_deltaIAF.png'));
            savefig(fig_h4, fullfile(outH4, 'H4_scatter_dphase_deltaIAF.fig'));
            fprintf('Global H4 scatter plot created.\n');
        else
            fprintf('Skipping H4 plot (no valid rows).\n');
        end
    catch ME
        warning('Could not create H4 plot: %s', ME.message);
    end

    % H5: Online performance by condition (A/B/C) with Sham as reference
   
    try
     
        T5p_online = Tall(Tall.AnalysisWindow == "Online" & ismember(string(Tall.Condition), ["A","B","C"]), :);

        shamRows5 = Tall(Tall.AnalysisWindow == "Offline" & string(Tall.Condition) == "Sham", :);
        shamAgg5  = table();
        subjList5 = unique(shamRows5.SubjectID);
        for s = 1:numel(subjList5)
            rows5    = shamRows5(shamRows5.SubjectID == subjList5(s), :);
            newRow5  = rows5(1, :);
            newRow5.dprime    = mean(rows5.dprime,    'omitnan');
            newRow5.criterion = mean(rows5.criterion, 'omitnan');
            newRow5.meanRT    = mean(rows5.meanRT,    'omitnan');
            shamAgg5 = [shamAgg5; newRow5];
        end

        T5p = [T5p_online; shamAgg5];

        % Remove phantom categories
        vn5 = T5p.Properties.VariableNames;
        for ii = 1:numel(vn5)
            if iscategorical(T5p.(vn5{ii}))
                T5p.(vn5{ii}) = removecats(T5p.(vn5{ii}));
            end
        end

        if any(categories(T5p.Condition) == "Sham")
            presentCats = categories(T5p.Condition);
            newOrder = ["Sham"; setdiff(presentCats, "Sham", 'stable')];
            T5p.Condition = reordercats(T5p.Condition, newOrder);
        end

        if height(T5p) > 0
            fprintf('Generating H5 plots...\n');

            % Extract p-values for A, B, C vs Sham from fitted models
            pvals_H5 = struct('A_dp', NaN, 'B_dp', NaN, 'C_dp', NaN, ...
                              'A_cr', NaN, 'B_cr', NaN, 'C_cr', NaN, ...
                              'A_rt', NaN, 'B_rt', NaN, 'C_rt', NaN);

            if ~isempty(OUT.H5.dprime.lme_randomSeq)
                pvals_H5.A_dp = extractConditionP(OUT.H5.dprime.lme_randomSeq, "A");
                pvals_H5.B_dp = extractConditionP(OUT.H5.dprime.lme_randomSeq, "B");
                pvals_H5.C_dp = extractConditionP(OUT.H5.dprime.lme_randomSeq, "C");
            end
            if ~isempty(OUT.H5.criterion.lme_randomSeq)
                pvals_H5.A_cr = extractConditionP(OUT.H5.criterion.lme_randomSeq, "A");
                pvals_H5.B_cr = extractConditionP(OUT.H5.criterion.lme_randomSeq, "B");
                pvals_H5.C_cr = extractConditionP(OUT.H5.criterion.lme_randomSeq, "C");
            end
            if ~isempty(OUT.H5.meanRT.lme_randomSeq)
                pvals_H5.A_rt = extractConditionP(OUT.H5.meanRT.lme_randomSeq, "A");
                pvals_H5.B_rt = extractConditionP(OUT.H5.meanRT.lme_randomSeq, "B");
                pvals_H5.C_rt = extractConditionP(OUT.H5.meanRT.lme_randomSeq, "C");
            end

            fig = figure;
            set(fig, 'Name', 'Global H5: General performance by condition', 'NumberTitle', 'off');
            set(fig, 'Position', [100 100 1100 420]);

            tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

            ax1 = nexttile(tl, 1);
            plot_H5_panel(ax1, T5p, 'dprime', 'd''', ...
                [pvals_H5.A_dp, pvals_H5.B_dp, pvals_H5.C_dp]);

            ax2 = nexttile(tl, 2);
            plot_H5_panel(ax2, T5p, 'criterion', 'Criterion (c)', ...
                [pvals_H5.A_cr, pvals_H5.B_cr, pvals_H5.C_cr]);

            ax3 = nexttile(tl, 3);
            plot_H5_panel(ax3, T5p, 'meanRT', 'Mean RT (s)', ...
                [pvals_H5.A_rt, pvals_H5.B_rt, pvals_H5.C_rt]);

            title(tl, 'H5: General performance by stimulation condition', ...
                'FontSize', 11, 'FontWeight', 'bold');
            subtitle(tl, 'Error bars = \pm1 SEM  |  * p < .05 vs. Sham (LMM)', ...
                'FontSize', 8, 'Color', [0.4 0.4 0.4]);

            saveas(fig, fullfile(outH5, 'H5_overview.png'));
            savefig(fig, fullfile(outH5, 'H5_overview.fig'));
            fprintf('H5 plot saved.\n');
        else
            fprintf('Skipping H5 plots (no rows).\n');
        end
    catch ME
        warning('Could not create H5 plots: %s', ME.message);
    end

    % H6: After-effects, change from pre-stimulation sham (postA/B/C minus pre)
    try
        T6p = Tall(Tall.AnalysisWindow == "Offline" & string(Tall.Condition) == "Sham", :);
        T6p = T6p(~isundefined(T6p.OfflineCondition), :);

        % Remove phantom categories
        vn6 = T6p.Properties.VariableNames;
        for ii = 1:numel(vn6)
            if iscategorical(T6p.(vn6{ii}))
                T6p.(vn6{ii}) = removecats(T6p.(vn6{ii}));
            end
        end

        if height(T6p) > 0
            fprintf('Generating H6 plots...\n');

            % Compute per-participant change scores (post minus pre)
            outcomes  = {'dprime','criterion','meanRT'};
            postConds = {'postA','postB','postC'};

            subjList6 = unique(T6p.SubjectID);
            T6diff = table();

            for s = 1:numel(subjList6)
                subj     = subjList6(s);
                subjRows = T6p(T6p.SubjectID == subj, :);

                preRow = subjRows(string(subjRows.OfflineCondition) == "pre", :);
                if isempty(preRow); continue; end

                for pc = 1:numel(postConds)
                    postRow = subjRows(string(subjRows.OfflineCondition) == postConds{pc}, :);
                    if isempty(postRow); continue; end

                    newRow = postRow(1,:);
                    newRow.OfflineCondition = categorical({postConds{pc}});
                    for oc = 1:numel(outcomes)
                        ov = outcomes{oc};
                        newRow.(ov) = postRow.(ov)(1) - preRow.(ov)(1);
                    end
                    T6diff = [T6diff; newRow];
                end
            end

            % Remove phantom categories from diff table
            vn6d = T6diff.Properties.VariableNames;
            for ii = 1:numel(vn6d)
                if iscategorical(T6diff.(vn6d{ii}))
                    T6diff.(vn6d{ii}) = removecats(T6diff.(vn6d{ii}));
                end
            end

            % Extract p-values for postA, postB, postC vs pre from fitted models
            pvals_H6 = struct('A_dp', NaN, 'B_dp', NaN, 'C_dp', NaN, ...
                              'A_cr', NaN, 'B_cr', NaN, 'C_cr', NaN, ...
                              'A_rt', NaN, 'B_rt', NaN, 'C_rt', NaN);

            if ~isempty(OUT.H6.dprime.lme_randomSeq)
                pvals_H6.A_dp = extractOfflineP(OUT.H6.dprime.lme_randomSeq, "postA");
                pvals_H6.B_dp = extractOfflineP(OUT.H6.dprime.lme_randomSeq, "postB");
                pvals_H6.C_dp = extractOfflineP(OUT.H6.dprime.lme_randomSeq, "postC");
            end
            if ~isempty(OUT.H6.criterion.lme_randomSeq)
                pvals_H6.A_cr = extractOfflineP(OUT.H6.criterion.lme_randomSeq, "postA");
                pvals_H6.B_cr = extractOfflineP(OUT.H6.criterion.lme_randomSeq, "postB");
                pvals_H6.C_cr = extractOfflineP(OUT.H6.criterion.lme_randomSeq, "postC");
            end
            if ~isempty(OUT.H6.meanRT.lme_randomSeq)
                pvals_H6.A_rt = extractOfflineP(OUT.H6.meanRT.lme_randomSeq, "postA");
                pvals_H6.B_rt = extractOfflineP(OUT.H6.meanRT.lme_randomSeq, "postB");
                pvals_H6.C_rt = extractOfflineP(OUT.H6.meanRT.lme_randomSeq, "postC");
            end

            fig = figure;
            set(fig, 'Name', 'Global H6: After-effects relative to pre-stimulation', 'NumberTitle', 'off');
            set(fig, 'Position', [100 100 1100 420]);

            tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

            ax1 = nexttile(tl, 1);
            plot_H6_panel(ax1, T6p, 'dprime', 'd''', ...
                [pvals_H6.A_dp, pvals_H6.B_dp, pvals_H6.C_dp]);

            ax2 = nexttile(tl, 2);
            plot_H6_panel(ax2, T6p, 'criterion', 'Criterion (c)', ...
                [pvals_H6.A_cr, pvals_H6.B_cr, pvals_H6.C_cr]);

            ax3 = nexttile(tl, 3);
            plot_H6_panel(ax3, T6p, 'meanRT', 'Mean RT (s)', ...
                [pvals_H6.A_rt, pvals_H6.B_rt, pvals_H6.C_rt]);

            title(tl, 'H6: After-effects on general performance (sham blocks)', ...
                'FontSize', 11, 'FontWeight', 'bold');
            subtitle(tl, 'Error bars = \pm1 SEM  |  * p < .05 vs. Pre (LMM)  |  Dashed line = Pre mean', ...
                'FontSize', 8, 'Color', [0.4 0.4 0.4]);

            saveas(fig, fullfile(outH6, 'H6_overview.png'));
            savefig(fig, fullfile(outH6, 'H6_overview.fig'));
            fprintf('H6 plot saved.\n');
        else
            fprintf('Skipping H6 plots (no rows).\n');
        end
    catch ME
        warning('Could not create H6 plots: %s', ME.message);
    end

    fprintf('Finished creating plots.\n');
end



%% Append model output to a single report file

function appendModelToReport(reportFile, modelName, lmeObj, anovaTbl, formulaStr)
    if isempty(lmeObj)
        return
    end

    fid = fopen(reportFile, 'a');
    if fid == -1
        warning('Could not open report file for writing: %s', reportFile);
        return
    end

    fprintf(fid, '\n---\n');
    fprintf(fid, 'Model: %s\n', modelName);
    fprintf(fid, 'Formula:\n%s\n', formulaStr);
    fprintf(fid, 'N = %d observations, %d fixed-effect coefficients\n', ...
        lmeObj.NumObservations, lmeObj.NumCoefficients);

    if ~isempty(anovaTbl)
        try
            termNames = string(anovaTbl.Term);
            keyIdx = find(contains(termNames, "Condition") | ...
                          contains(termNames, "deltaIAF"));
            if ~isempty(keyIdx)
                fprintf(fid, '\nSummary of the terms of interest\n');
                fprintf(fid, '%-34s %10s %8s %10s %10s\n', ...
                    'Term', 'F', 'DF1', 'DF2', 'p');
                for kk = keyIdx(:)'
                    pv = anovaTbl.pValue(kk);
                    if pv < 0.0001
                        pTxt = '<.0001';
                    else
                        pTxt = sprintf('%.4f', pv);
                    end
                    fprintf(fid, '%-34s %10.3f %8d %10.1f %10s\n', ...
                        termNames(kk), anovaTbl.FStat(kk), anovaTbl.DF1(kk), ...
                        anovaTbl.DF2(kk), pTxt);
                end
            end
        catch
            
        end
    end

    fprintf(fid, '\nFull model output\n');
    lmeTxt = evalc('disp(lmeObj)');
    fprintf(fid, '%s\n', lmeTxt);

    if ~isempty(anovaTbl)
        fprintf(fid, '\nANOVA (Satterthwaite)\n');
        anovaTxt = evalc('disp(anovaTbl)');
        fprintf(fid, '%s\n', anovaTxt);
    end

    fclose(fid);
end

end



%%  Subfunctions


function plot_H5_panel(ax, Tp, yVar, yLabelText, pvals)

    condOrder = {'Sham','A','B','C'};
    [activeCols, activeNames] = tACS_condColours(false);
    condLabel = [{'Sham'}, activeNames];
    colours   = [0.50 0.50 0.50; activeCols];

    groups = cell(1, numel(condOrder));
    for j = 1:numel(condOrder)
        idx = string(Tp.Condition) == condOrder{j};
        groups{j} = Tp.(yVar)(idx);
    end

    tACS_raincloudPanel(ax, groups, condLabel, colours, ...
        'RefStars', pvals, 'Separator', true, ...
        'YLabel', yLabelText, 'XLabel', 'Condition');
end

function plot_H6_panel(ax, T6raw, yVar, yLabelText, pvals)

    condOrder = {'pre','postA','postB','postC'};
    [activeCols, activeNames] = tACS_condColours(false);
    condLabel = [{'Pre'}, cellfun(@(s) ['Post-' lower(s)], activeNames, 'UniformOutput', false)];
    colours   = [0.50 0.50 0.50; activeCols];

    groups = cell(1, numel(condOrder));
    for j = 1:numel(condOrder)
        idx = string(T6raw.OfflineCondition) == condOrder{j};
        groups{j} = T6raw.(yVar)(idx);
    end

    tACS_raincloudPanel(ax, groups, condLabel, colours, ...
        'RefStars', pvals, 'Separator', true, ...
        'YLabel', yLabelText, 'XLabel', 'Sham block');
end

function reportBayes(reportFile, effectName, bf)

    if isempty(bf) || ~isstruct(bf) || ~bf.ok
        return
    end

    lines = [ ...
        "";
        "Bayes factor, BIC approximation (Wagenmakers, 2007): " + effectName;
        "  BIC full / reduced : " + sprintf('%.2f / %.2f', bf.BIC_full, bf.BIC_reduced);
        "  Reduced model      : " + bf.formula_reduced;
        "  N                  : " + sprintf('%d', bf.nObs);
        "  BF01 / BF10        : " + sprintf('%.4g / %.4g', bf.BF01, bf.BF10);
        "  Evidence           : " + bf.evidence;
        ""];

    fprintf('\n%s\n', strjoin(cellstr(lines), newline));

    try
        fid = fopen(reportFile, 'a');
        if fid > 0
            fprintf(fid, '%s\n', strjoin(cellstr(lines), newline));
            fclose(fid);
        end
    catch
       
    end
end


function p = extractConditionP(lmeObj, condLabel)
% Extract p-value for a Condition contrast vs Sham from a fitted LMM (used by H5).
    p = NaN;
    if isempty(lmeObj); return; end
    try
        coeffTbl = lmeObj.Coefficients;
        rowNames = string(coeffTbl.Name);
        matchIdx = find(contains(rowNames, "Condition_" + condLabel) | ...
                       contains(rowNames, "Condition" + condLabel));
        if ~isempty(matchIdx)
            p = coeffTbl.pValue(matchIdx(1));
        end
    catch
    end
end


function p = extractOfflineP(lmeObj, condLabel)
% Extract p-value for an OfflineCondition contrast vs pre from a fitted LMM (used by H6).
    p = NaN;
    if isempty(lmeObj); return; end
    try
        coeffTbl = lmeObj.Coefficients;
        rowNames = string(coeffTbl.Name);
        matchIdx = find(contains(rowNames, "OfflineCondition_" + condLabel) | ...
                       contains(rowNames, "OfflineCondition" + condLabel));
        if ~isempty(matchIdx)
            p = coeffTbl.pValue(matchIdx(1));
        end
    catch
    end
end