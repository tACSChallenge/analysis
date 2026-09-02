function tACS_exportBayesData(pooledFile, outDir)
% Helper for the non-registered / specified Bayesian analysis:
% reads the pooled metrics file written by tACSChallenge_runH1H2H3_Pooled and writes one long-format CSV that the
% accompanying R script consumes.
%
% The exported estimate is the bias-corrected phasic modulation per
% participant and condition:
%
%     delta = observed coefficient - mean of that participant's own
%             permutation distribution for that condition
%
% Subtracting each participant's own permuted mean removes the positive bias
% that the regression coefficient carries even under the null, so that a
% one-sample test of delta against zero is a meaningful test of phasic
% modulation.

    if nargin < 1 || isempty(pooledFile)
        here         = fileparts(mfilename('fullpath'));
        project_root = fileparts(here);   % bayes -> root
        pooledFile   = fullfile(project_root, 'results', 'H1H2H3', 'pooled', ...
                                'H1H2H3_pooled_metrics.mat');
    end
    if ~isfile(pooledFile)
        error(['Pooled metrics file not found: %s\nRun the H1 to H3 analysis ' ...
            'first, or pass the path explicitly.'], pooledFile);
    end
    if nargin < 2 || isempty(outDir)
        outDir = fileparts(pooledFile);
    end
    if ~isfolder(outDir); mkdir(outDir); end

    S = load(pooledFile, 'all_bs_allLabs', 'all_bs_perm_allLabs', ...
                         'subj_labels_allLabs');

    obs      = S.all_bs_allLabs;                          % nCond x nSubj
    permMean = mean(S.all_bs_perm_allLabs, 3, 'omitnan'); % nCond x nSubj
    labels   = S.subj_labels_allLabs;

    nSubj = size(obs, 2);
    if numel(labels) ~= nSubj
        warning(['Subject label count (%d) does not match the number of columns ' ...
            '(%d). Generic labels are used instead.'], numel(labels), nSubj);
        labels = strings(nSubj,1);
        for s = 1:nSubj
            labels(s) = sprintf('S%03d', s);
        end
    end

    condKeys = {'A', 'B', 'C'};
    nCond    = numel(condKeys);

    SubjectID = strings(nSubj*nCond, 1);
    Condition = strings(nSubj*nCond, 1);
    Amplitude = nan(nSubj*nCond, 1);
    PermMean  = nan(nSubj*nCond, 1);

    row = 0;
    for c = 1:nCond
        for s = 1:nSubj
            row = row + 1;
            SubjectID(row) = string(labels(s));
            Condition(row) = string(condKeys{c});
            Amplitude(row) = obs(c, s);
            PermMean(row)  = permMean(c, s);
        end
    end

    delta = Amplitude - PermMean;
    T = table(SubjectID, Condition, Amplitude, PermMean, delta);

    outFile = fullfile(outDir, 'bayes_H1H2_delta.csv');
    writetable(T, outFile);

    fprintf('Exported %d rows (%d participants x %d conditions) to:\n  %s\n', ...
        height(T), nSubj, nCond, outFile);
end
