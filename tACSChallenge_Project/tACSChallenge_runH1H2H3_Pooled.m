function OUT = tACSChallenge_runH1H2H3_Pooled(data_root, lab_ids, varargin)
%% Pooled H1 H2 H3 analysis across labs, while still writing PhaseMetrics.mat per lab.
%% script originally written by Anne-Fiona Griesfeller, CNRS Toulouse, in Jan 26
% Outputs only pooled results into OutDir.
%
% Required inputs
%   data_root: path to folder that contains L18, L19, etc
%   lab_ids: numeric vector of labs to include
%

p = inputParser;
p.addParameter('OutDir', '', @(x)ischar(x) || isstring(x));
p.addParameter('perm', 100, @(x)isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('conditions', {'*_A.*','*_B.*','*_C.*','*_Sh.*','*_A_*','*_B_*','*_C_*','*_Sh_*'}, @(x)iscell(x) && ~isempty(x));
p.addParameter('MakePlots', true, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});

OutDir = string(p.Results.OutDir);
perm = p.Results.perm;
conditions = p.Results.conditions;
MakePlots = logical(p.Results.MakePlots);

if OutDir == ""
    OutDir = fullfile(string(data_root), 'results', 'H1H2H3', 'pooled');
end

if ~isfolder(OutDir)
    mkdir(OutDir);
end

outH1 = fullfile(OutDir, 'H1');
outH2 = fullfile(OutDir, 'H2');
outH3 = fullfile(OutDir, 'H3');

if ~isfolder(outH1); mkdir(outH1); end
if ~isfolder(outH2); mkdir(outH2); end
if ~isfolder(outH3); mkdir(outH3); end

no_phases = 8;
n_cond = length(conditions) - 4;

% We pool across labs at the subject level by concatenating columns
all_ps_allLabs = [];
all_bs_allLabs = [];
all_pfs_allLabs = [];
all_bs_perm_allLabs = [];
all_hit_probs_allLabs = [];
subj_labels_allLabs = strings(0,1);

fprintf('\nRunning H1 H2 H3 pooled across %d labs.\n', numel(lab_ids));

for i = 1:numel(lab_ids)
    lab_id = lab_ids(i);
    lab_path = fullfile(data_root, sprintf('L%d', lab_id));

    if ~isfolder(lab_path)
        warning('Lab folder not found: %s (skipping)', lab_path);
        continue
    end

    subj_dirs = dir(fullfile(lab_path, 'sub-*'));
    subj_dirs = subj_dirs([subj_dirs.isdir]);
    subjs = {subj_dirs.name};

    fprintf('Lab %d: %d subjects.\n', lab_id, numel(subjs));

    if isempty(subjs)
        continue
    end

    % Per lab containers
    all_ps = zeros(n_cond, numel(subjs));
    all_bs = zeros(n_cond, numel(subjs));
    all_hit_probs = zeros(no_phases, n_cond, numel(subjs));
    all_bs_perm = zeros(n_cond, numel(subjs), perm);
    all_pfs = zeros(n_cond, numel(subjs));

    for s = 1:numel(subjs)
        curr_data = tACSChallenge_SortData(lab_path, subjs{s}, conditions);

        % Use a lab specific seed that stays stable across pooled runs
        seed = double(lab_id) * 1000 + s;

        [all_ps(:,s), all_bs(:,s), all_hit_probs(:,:,s), all_bs_perm(:,s,:), all_pfs(:,s)] = ...
            tACSChallenge_EvalData(curr_data, perm, seed);

        subj_labels_allLabs(end+1,1) = string(sprintf('L%d_%s', lab_id, subjs{s}));
    end

    % Save PhaseMetrics.mat inside the lab folder
    phase_results_file = fullfile(lab_path, 'PhaseMetrics.mat');
    save(phase_results_file, ...
        'all_bs', 'all_ps', 'all_pfs', 'all_hit_probs', 'all_bs_perm', ...
        'subjs', 'conditions', '-v7.3');
    fprintf('Saved PhaseMetrics to %s\n', phase_results_file);

    % Append into pooled containers
    all_ps_allLabs = [all_ps_allLabs, all_ps];
    all_bs_allLabs = [all_bs_allLabs, all_bs];
    all_pfs_allLabs = [all_pfs_allLabs, all_pfs];

    if isempty(all_bs_perm_allLabs)
        all_bs_perm_allLabs = all_bs_perm;
    else
        all_bs_perm_allLabs = cat(2, all_bs_perm_allLabs, all_bs_perm);
    end

    if isempty(all_hit_probs_allLabs)
        all_hit_probs_allLabs = all_hit_probs;
    else
        all_hit_probs_allLabs = cat(3, all_hit_probs_allLabs, all_hit_probs);
    end
end

% Duplicate first phase bin for plotting like in the original script
all_hit_probs_plot = all_hit_probs_allLabs;
all_hit_probs_plot(no_phases+1,:,:) = all_hit_probs_plot(1,:,:);

% Group level statistics
group_level_p    = zeros(n_cond,1);
group_level_comp = zeros(2,1);
pf_p             = zeros(2,1);

for cond = 1:n_cond
    currmean   = mean(all_bs_allLabs(cond,:), 2);
    currsurr   = squeeze(mean(all_bs_perm_allLabs(cond,:,:), 2));
    currsurrmu = mean(currsurr);
    currsurrsd = std(currsurr);
    [~, group_level_p(cond)] = ztest(currmean, currsurrmu, currsurrsd, 'tail', 'right');
end

[~, group_level_comp(1)] = ttest(all_bs_allLabs(1,:), all_bs_allLabs(2,:), 'tail', 'right');
[~, group_level_comp(2)] = ttest(all_bs_allLabs(1,:), all_bs_allLabs(3,:), 'tail', 'right');

pf_p(1) = circ_mtest(circ_dist(all_pfs_allLabs(1,:), all_pfs_allLabs(2,:)), 0);
pf_p(2) = circ_mtest(circ_dist(all_pfs_allLabs(1,:), all_pfs_allLabs(3,:)), 0);

% Write per hypothesis reports in their own folders
reportH1 = fullfile(outH1, 'H1_results.txt');
if isfile(reportH1); delete(reportH1); end
fid = fopen(reportH1, 'a');
fprintf(fid, 'H1: Phasic modulation per condition (z test vs permutation surrogate)\n');
for cond = 1:n_cond
    fprintf(fid, 'Condition %d (A=1, B=2, C=3, Sham=4): p = %.10f\n', cond, group_level_p(cond));
end
fclose(fid);

reportH2 = fullfile(outH2, 'H2_results.txt');
if isfile(reportH2); delete(reportH2); end
fid = fopen(reportH2, 'a');
fprintf(fid, 'H2: Contrasts (paired t tests, right tailed)\n');
fprintf(fid, 'A > B: p = %.10f\n', group_level_comp(1));
fprintf(fid, 'A > C: p = %.10f\n', group_level_comp(2));
fclose(fid);

reportH3 = fullfile(outH3, 'H3_results.txt');
if isfile(reportH3); delete(reportH3); end
fid = fopen(reportH3, 'a');
fprintf(fid, 'H3: Preferred phase differences (circular m tests)\n');
fprintf(fid, 'A vs B: p = %.10f\n', pf_p(1));
fprintf(fid, 'A vs C: p = %.10f\n', pf_p(2));
fclose(fid);

% Save pooled metrics
save(fullfile(OutDir, 'H1H2H3_pooled_metrics.mat'), ...
    'all_bs_allLabs', 'all_ps_allLabs', 'all_pfs_allLabs', 'all_hit_probs_allLabs', 'all_bs_perm_allLabs', ...
    'subj_labels_allLabs', 'conditions', 'lab_ids', ...
    'group_level_p', 'group_level_comp', 'pf_p', '-v7.3');

% Plots, pooled only
if MakePlots
    fig1 = figure;
    plot(-pi:pi/4:pi, squeeze(mean(all_hit_probs_plot,3)));
    xlabel('tACS phase'); ylabel('detection probability');
    legend(conditions(1:n_cond));
    saveas(fig1, fullfile(outH1, 'H1_hitprob_by_phase.png'));
    savefig(fig1, fullfile(outH1, 'H1_hitprob_by_phase.fig'));

    fig2 = figure;
    bar(mean(all_bs_allLabs,2));
    xlabel('conditions');
    xticklabels(conditions(1:n_cond));
    ylabel('modulation strength');
    saveas(fig2, fullfile(outH1, 'H1_modulation_strength.png'));
    savefig(fig2, fullfile(outH1, 'H1_modulation_strength.fig'));

    % This is the same figure but stored in H2 as well for convenience
    saveas(fig2, fullfile(outH2, 'H2_modulation_strength.png'));
    savefig(fig2, fullfile(outH2, 'H2_modulation_strength.fig'));
end

OUT = struct();
OUT.lab_ids = lab_ids;
OUT.conditions = conditions;
OUT.perm = perm;
OUT.subj_labels = subj_labels_allLabs;

OUT.all_ps = all_ps_allLabs;
OUT.all_bs = all_bs_allLabs;
OUT.all_pfs = all_pfs_allLabs;
OUT.all_bs_perm = all_bs_perm_allLabs;
OUT.all_hit_probs = all_hit_probs_allLabs;

OUT.group_level_p = group_level_p;
OUT.group_level_comp = group_level_comp;
OUT.pf_p = pf_p;
OUT.OutDir = OutDir;

end
