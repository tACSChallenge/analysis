function OUT = tACSChallenge_runGlobalModels(root_path, lab_ids, varargin)
%% Run H4, H5, H6 mixed models across labs using unified AnalysisTable.
%% script originally written by Anne-Fiona Griesfeller, CNRS Toulouse, in Jan 26

% INPUT
%   root_path : folder that contains lab folders L##, e.g. 'C:\Users\...\Documents\Data'
%   lab_ids   : numeric vector of lab IDs to include, e.g. [18 28 31]
%
% If only one Hypothesis should be analyzed then an optional name-value can
% be added e.g. 'RunH4'(true/false). The default is always true.
%
% OUTPUT
%   OUT is a struct containing combined table and fitted models.

p = inputParser;
p.addParameter('RunH4', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('RunH5', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('RunH6', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('MakePlots', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('OutDir', '', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});

RunH4 = logical(p.Results.RunH4);
RunH5 = logical(p.Results.RunH5);
RunH6 = logical(p.Results.RunH6);
MakePlots = logical(p.Results.MakePlots);
OutDir = string(p.Results.OutDir);

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

Tall = table();

for i = 1:numel(lab_ids)
    lab_id = lab_ids(i);
    lab_path = fullfile(root_path, sprintf('L%d', lab_id));

    if ~isfolder(lab_path)
        warning('Lab folder not found: %s (skipping)', lab_path);
        continue
    end

    Tlab = tACSChallenge_buildAnalysisTable(lab_path, lab_id);
    Tall = [Tall; Tlab]; 
    fprintf('Added lab %d (%d rows). Combined now %d rows.\n', lab_id, height(Tlab), height(Tall));
end

OUT.T_all = Tall;

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
Tall.Gender          = categorical(string(Tall.Gender));
Tall.Sequence        = categorical(string(Tall.Sequence));

% Reference levels
if any(categories(Tall.Condition) == "A")
    Tall.Condition = reordercats(Tall.Condition, ["A"; setdiff(categories(Tall.Condition), "A", 'stable')]);
end
if any(categories(Tall.OfflineCondition) == "pre")
    Tall.OfflineCondition = reordercats(Tall.OfflineCondition, ["pre"; setdiff(categories(Tall.OfflineCondition), "pre", 'stable')]);
end

    % Helper: drop fixed-effect covariates that have no variability to(prevent "X must be full column rank" errors)

    function [fixedTermsKept, dropped] = keepVaryingTerms(T, fixedTerms)
    % Returns fixedTermsKept as a string array, with non varying or missing terms removed.
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


    % Helper: choose random effects based on number of levels
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
    T4 = rmmissing(T4, 'DataVariables', {'d_phase','deltaIAF'}); % must have both for interaction

    % Fixed effects for H4:
    fixedWanted = {'Condition','deltaIAF','StimIntensity','HemiImpDiff','Gender','Age','EMLA','IAF'};
    [fixedKept, dropped] = keepVaryingTerms(T4, fixedWanted);

    if ~isempty(dropped)
        fprintf('H4 dropping covariates: %s\n', strjoin(cellstr(dropped), ', '));
    end

    % Build formula:
    % d_phase ~ Condition*deltaIAF + StimIntensity + HemiImpDiff + Gender + Age + EMLA + random
    base = "d_phase ~ 1 + Condition*deltaIAF";

    % Add remaining main effects except Condition and deltaIAF
    addMain = setdiff(string(fixedKept), ["Condition","deltaIAF"], 'stable');
    if ~isempty(addMain)
        base = base + " + " + strjoin(cellstr(addMain), " + ");
    end

    re = randomPart(T4, true); % Sequence as random effect 
    formula = char(base + re);

    OUT.H4.formula_used = formula;

    try
        OUT.H4.lme = fitlme(T4, formula);
        disp(OUT.H4.lme);
        OUT.H4.anova = anova(OUT.H4.lme, 'DFMethod','satterthwaite');
        disp(OUT.H4.anova);
        appendModelToReport(reportH4, "H4 randomSeq", OUT.H4.lme, OUT.H4.anova, OUT.H4.formula_used);
    catch ME
        warning('H4 model failed: %s', ME.message);
        OUT.H4.lme = [];
        OUT.H4.anova = [];
    end
    

    % H4 variant with Sequence as fixed effect (instead of random)

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

        re2 = randomPart(T4, false); % no random Sequence
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


%% H5 models (online conditions, outcomes dprime, criterion, meanRT)

if RunH5
    fprintf('\nFitting H5 global mixed models (dprime, criterion, meanRT)...\n');
    reportH5 = fullfile(outH5, 'H5_model_report.txt');
if isfile(reportH5)
    delete(reportH5);
end

    T5 = Tall(Tall.AnalysisWindow == "Online" & ismember(string(Tall.Condition), ["A","B","C"]), :);

    fixedWanted = {'Condition','StimIntensity','HemiImpDiff','Gender','Age','EMLA','IAF'};
    [fixedKept, dropped] = keepVaryingTerms(T5, fixedWanted);
    if ~isempty(dropped)
        fprintf('H5 dropping covariates: %s\n', strjoin(cellstr(dropped), ', '));
    end

    base = " ~ 1 + " + strjoin(cellstr(string(fixedKept)), " + ");

    % First variant with Sequence as random effect
    re_rand = randomPart(T5, true);

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
    catch ME
        warning('H5 dprime randomSeq model failed: %s', ME.message);
        OUT.H5.dprime.lme_randomSeq = [];
        OUT.H5.dprime.randomSeq.anova = [];
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
    catch ME
        warning('H5 criterion randomSeq model failed: %s', ME.message);
        OUT.H5.criterion.lme_randomSeq = [];
        OUT.H5.criterion.randomSeq.anova = [];
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
    catch ME
        warning('H5 meanRT randomSeq model failed: %s', ME.message);
        OUT.H5.meanRT.lme_randomSeq = [];
        OUT.H5.meanRT.randomSeq.anova = [];
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


%% H6 aftereffects model (offline conditions (sham), outcomes dprime, criterion, meanRT)

if RunH6
    fprintf('\nFitting H6 global mixed models (offline sham)...\n');
reportH6 = fullfile(outH6, 'H6_model_report.txt');
if isfile(reportH6)
    delete(reportH6);
end

    T6 = Tall(Tall.AnalysisWindow == "Offline" & string(Tall.Condition) == "Sham", :);
    T6 = T6(~isundefined(T6.OfflineCondition), :);

    % H6 fixed effects:
    % OfflineCondition, StimIntensity (not defined here), HemiImpDiff, Gender, Age, EMLA, IAF
    
    fixedWanted = {'OfflineCondition','StimIntensity','HemiImpDiff','Gender','Age','EMLA','IAF'};
    [fixedKept, dropped] = keepVaryingTerms(T6, fixedWanted);
    if ~isempty(dropped)
        fprintf('H6 dropping covariates: %s\n', strjoin(cellstr(dropped), ', '));
    end

    base = " ~ 1 + " + strjoin(cellstr(string(fixedKept)), " + ");

    % Sequence as random effect
    re_rand = randomPart(T6, true);

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

    catch ME
        warning('H6 dprime randomSeq model failed: %s', ME.message);
        OUT.H6.dprime.lme_randomSeq = [];
        OUT.H6.dprime.randomSeq.anova = [];
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

    catch ME
        warning('H6 criterion randomSeq model failed: %s', ME.message);
        OUT.H6.criterion.lme_randomSeq = [];
        OUT.H6.criterion.randomSeq.anova = [];
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

    catch ME
        warning('H6 meanRT randomSeq model failed: %s', ME.message);
        OUT.H6.meanRT.lme_randomSeq = [];
        OUT.H6.meanRT.randomSeq.anova = [];
    end

    % Sequence as fixed effect instead of random effect
    if ismember("Sequence", string(T6.Properties.VariableNames)) && numel(categories(T6.Sequence)) > 1
        fixedWanted2 = [cellstr(string(fixedKept)) {'Sequence'}];
        % Make sure Sequence has variability, else drop
        [fixedKept2, dropped2] = keepVaryingTerms(T6, fixedWanted2);
        if ~isempty(dropped2)
            fprintf('H6 fixedSeq dropping covariates: %s\n', strjoin(cellstr(dropped2), ', '));
        end
        base2 = " ~ 1 + " + strjoin(cellstr(string(fixedKept2)), " + ");

        re_fixed = randomPart(T6, false); % no random Sequence

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

            conds = categories(T4p.Condition);

            figure;
            set(gcf, 'Name', 'Global H4: d_phase vs deltaIAF by condition', 'NumberTitle', 'off');

            for i = 1:numel(conds)
                thisCond = conds{i};
                idx      = T4p.Condition == thisCond;

                x = T4p.deltaIAF(idx);
                y = T4p.d_phase(idx);

                valid = ~isnan(x) & ~isnan(y);
                x = x(valid);
                y = y(valid);

                if numel(x) > 1
                    pp = polyfit(x, y, 1);
                    x_line = linspace(min(x), max(x), 100);
                    y_line = polyval(pp, x_line);
                else
                    pp = [NaN NaN];
                    x_line = [];
                    y_line = [];
                end

                subplot(1, numel(conds), i);
                plot(x, y, 'o', 'MarkerSize', 6);
                hold on;

                if ~isempty(x_line)
                    plot(x_line, y_line, '-', 'LineWidth', 1.5);
                end

                xlabel('\Delta IAF (Hz)');
                ylabel('d\_phase');
                title(sprintf('Condition %s', thisCond));
                grid on;

                if ~isnan(pp(1))
                    txt = sprintf('Slope = %.4f', pp(1));
                    xlim_current = xlim;
                    ylim_current = ylim;
                    text(xlim_current(1) + 0.05 * range(xlim_current), ...
                         ylim_current(2) - 0.1 * range(ylim_current), ...
                         txt, 'FontSize', 8, 'BackgroundColor', 'w');
                end
            end
            saveas(gcf, fullfile(outH4, 'H4_scatter_dphase_deltaIAF.png'));
            savefig(gcf, fullfile(outH4, 'H4_scatter_dphase_deltaIAF.fig'));

            fprintf('Global H4 scatter plot created.\n');
        else
            fprintf('Skipping H4 plot (no valid rows).\n');
        end
    catch ME
        warning('Could not create H4 plot: %s', ME.message);
    end

    
    % H5: Online performance by Condition
   
    try
        T5p = Tall(Tall.AnalysisWindow == "Online" & ismember(string(Tall.Condition), ["A","B","C"]), :);

        if height(T5p) > 0
            fprintf('Generating H5 plots...\n');
        fig = figure;
        set(fig, 'Name', 'Global H5 overview', 'NumberTitle', 'off');

        t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

        ax1 = nexttile(t, 1);
        plot_mean_sem_with_points_ax(ax1, T5p, 'Condition', 'dprime', 'Condition', 'dprime', 'dprime');

        ax2 = nexttile(t, 2);
        plot_mean_sem_with_points_ax(ax2, T5p, 'Condition', 'criterion', 'Condition', 'criterion', 'criterion');

        ax3 = nexttile(t, 3);
        T5p_rt = T5p(~isnan(T5p.meanRT), :);
        plot_mean_sem_with_points_ax(ax3, T5p_rt, 'Condition', 'meanRT', 'Condition', 'Mean RT (s)', 'Mean RT');

    saveas(fig, fullfile(outH5, 'H5_overview.png'));
    savefig(fig, fullfile(outH5, 'H5_overview.fig'));
   

        else
            fprintf('Skipping H5 plots (no rows).\n');
        end
    catch ME
        warning('Could not create H5 plots: %s', ME.message);
    end

    
    % H6: Offline aftereffects by OfflineCondition (pre/postA/postB/postC)
    
    try
        T6p = Tall(Tall.AnalysisWindow == "Offline" & string(Tall.Condition) == "Sham", :);
        T6p = T6p(~isundefined(T6p.OfflineCondition), :);

        if height(T6p) > 0
            fprintf('Generating H6 plots...\n');

           fig = figure;
            set(fig, 'Name', 'Global H6 overview', 'NumberTitle', 'off');

            t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

            ax1 = nexttile(t, 1);
        T6p_dp = T6p(~isnan(T6p.dprime), :);
        plot_mean_sem_with_points_ax(ax1, T6p_dp, 'OfflineCondition', 'dprime', 'Offline condition', 'dprime', 'dprime');

        ax2 = nexttile(t, 2);
        T6p_c = T6p(~isnan(T6p.criterion), :);
        plot_mean_sem_with_points_ax(ax2, T6p_c, 'OfflineCondition', 'criterion', 'Offline condition', 'criterion', 'criterion');

        ax3 = nexttile(t, 3);
        T6p_rt = T6p(~isnan(T6p.meanRT), :);
        plot_mean_sem_with_points_ax(ax3, T6p_rt, 'OfflineCondition', 'meanRT', 'Offline condition', 'Mean RT (s)', 'Mean RT');

    saveas(fig, fullfile(outH6, 'H6_overview.png'));
    savefig(fig, fullfile(outH6, 'H6_overview.fig'));
    
        else
            fprintf('Skipping H6 plots (no rows).\n');
        end
    catch ME
        warning('Could not create H6 plots: %s', ME.message);
    end

    fprintf('Finished creating plots.\n');
end

%% Helper: append model output to a single report file

function appendModelToReport(reportFile, modelName, lmeObj, anovaTbl, formulaStr)
    if isempty(lmeObj)
        return
    end

    fid = fopen(reportFile, 'a');
    if fid == -1
        warning('Could not open report file for writing: %s', reportFile);
        return
    end

    fprintf(fid, '\n====================\n');
    fprintf(fid, 'Model: %s\n', modelName);
    fprintf(fid, 'Formula:\n%s\n\n', formulaStr);

    % Capture printed output of disp into a string and write it
    lmeTxt = evalc('disp(lmeObj)');
    fprintf(fid, '%s\n', lmeTxt);

    if ~isempty(anovaTbl)
        fprintf(fid, '\nANOVA (Satterthwaite)\n');
        anovaTxt = evalc('disp(anovaTbl)');
        fprintf(fid, '%s\n', anovaTxt);
    end

    fclose(fid);
end

end  % end of tACSChallenge_runGlobalModels



%% Subfunction: mean + SEM by category with points

function plot_mean_sem_with_points_ax(ax, Tp, catVar, yVar, xLabelText, yLabelText, panelTitle)
    axes(ax);

    cats = categories(Tp.(catVar));
    nC = numel(cats);

    means = nan(nC,1);
    sems  = nan(nC,1);

    hold on;

    for j = 1:nC
        cj = cats{j};
        idx = Tp.(catVar) == cj;

        y = Tp.(yVar)(idx);
        y = y(~isnan(y));

        if isempty(y)
            continue
        end

        means(j) = mean(y);
        sems(j)  = std(y) / sqrt(numel(y));

        x0 = j * ones(sum(idx),1);
        xj = x0 + (rand(sum(idx),1) - 0.5) * 0.15;

        yy = Tp.(yVar)(idx);
        plot(xj, yy, 'o', 'MarkerSize', 5);
    end

    x = (1:nC)';
    errorbar(x, means, sems, 'LineStyle', 'none', 'LineWidth', 1.5);

    set(gca, 'XTick', 1:nC, 'XTickLabel', cats);
    xlabel(xLabelText);
    ylabel(yLabelText);
    title(panelTitle);
    grid on;
    box on;
end
