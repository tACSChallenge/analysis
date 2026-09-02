function OUT = tACSChallenge_runH1H2H3_Pooled(data_root, lab_ids, varargin)
%% Pooled H1 H2 H3 analysis across labs, plus writing PhaseMetrics.mat per lab.
%% script originally written by Anne-Fiona Griesfeller

p = inputParser;
p.addParameter('OutDir', '', @(x)ischar(x) || isstring(x));
p.addParameter('perm', 1000, @(x)isnumeric(x) && isscalar(x) && x > 0);

p.addParameter('conditions', {'*_A_*','*_B_*','*_C_*','*_Sh_*'}, ...
    @(x)iscell(x) && numel(x) == 4);
p.addParameter('MakePlots', true, @(x)islogical(x) || isnumeric(x));
p.addParameter('forceRecompute', false, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});

OutDir         = string(p.Results.OutDir);
perm           = p.Results.perm;
conditions     = p.Results.conditions;
MakePlots      = logical(p.Results.MakePlots);
forceRecompute = logical(p.Results.forceRecompute);

if OutDir == ""
    OutDir = fullfile(string(data_root), 'results', 'H1H2H3', 'pooled');
end

if ~isfolder(OutDir); mkdir(OutDir); end

outH1 = fullfile(OutDir, 'H1');
outH2 = fullfile(OutDir, 'H2');
outH3 = fullfile(OutDir, 'H3');

if ~isfolder(outH1); mkdir(outH1); end
if ~isfolder(outH2); mkdir(outH2); end
if ~isfolder(outH3); mkdir(outH3); end

no_phases = 8;
n_cond    = 4;

% Pooled containers
all_ps_allLabs        = [];
all_bs_allLabs        = [];
all_pfs_allLabs       = [];
all_bs_perm_allLabs   = [];
all_hit_probs_allLabs = [];
subj_labels_allLabs   = strings(0,1);

fprintf('\nRunning H1 H2 H3 pooled across %d labs.\n', numel(lab_ids));

for i = 1:numel(lab_ids)
    lab_id   = lab_ids(i);
    lab_path = fullfile(data_root, sprintf('L%02d', lab_id));

    if ~isfolder(lab_path)
        warning('Lab folder not found: %s (skipping)', lab_path);
        continue
    end

    phase_results_file = fullfile(lab_path, 'PhaseMetrics.mat');

    % if no forcerecompute, already saved phasemetrics is loaded
    recomputeThisLab = true;

    if isfile(phase_results_file) && ~forceRecompute

        fprintf('Lab %d: found existing PhaseMetrics.mat, checking it.\n', lab_id);

        S = load(phase_results_file);

        if isfield(S, 'perm')
            cachedPerm = S.perm;
        else
            cachedPerm = size(S.all_bs_perm, 3);
        end

        if cachedPerm ~= perm
            fprintf('Lab %d: cached PhaseMetrics.mat rejected (%d permutations, %d requested).\n', ...
                lab_id, cachedPerm, perm);
        else
            recomputeThisLab = false;
            fprintf('Lab %d: cache accepted (%d permutations), skipping recomputation.\n', ...
                lab_id, cachedPerm);
        end
    end

    if ~recomputeThisLab

        all_ps        = S.all_ps;
        all_bs        = S.all_bs;
        all_pfs       = S.all_pfs;
        all_hit_probs = S.all_hit_probs;
        all_bs_perm   = S.all_bs_perm;
        subjs         = S.subjs;

        for s = 1:numel(subjs)
            subj_labels_allLabs(end+1,1) = string(sprintf('L%d_%s', lab_id, subjs{s}));
        end

    else

        subj_dirs = dir(fullfile(lab_path, 'sub-*'));
        subj_dirs = subj_dirs([subj_dirs.isdir]);
        subjs     = {subj_dirs.name};

        fprintf('Lab %d: %d subjects.\n', lab_id, numel(subjs));

        if isempty(subjs); continue; end

        % Per-lab containers
        all_ps        = zeros(n_cond, numel(subjs));
        all_bs        = zeros(n_cond, numel(subjs));
        all_hit_probs = zeros(no_phases, n_cond, numel(subjs));
        all_bs_perm   = zeros(n_cond, numel(subjs), perm);
        all_pfs       = zeros(n_cond, numel(subjs));

        for s = 1:numel(subjs)
            curr_data = tACSChallenge_SortData(lab_path, subjs{s}, conditions);
            seed      = double(lab_id) * 1000 + s;

            [all_ps(:,s), all_bs(:,s), all_hit_probs(:,:,s), all_bs_perm(:,s,:), all_pfs(:,s)] = ...
                tACSChallenge_EvalData(curr_data, perm, seed);

            subj_labels_allLabs(end+1,1) = string(sprintf('L%d_%s', lab_id, subjs{s}));
        end

        % Save per-lab PhaseMetrics
        save(phase_results_file, ...
            'all_bs', 'all_ps', 'all_pfs', 'all_hit_probs', 'all_bs_perm', ...
            'subjs', 'conditions', 'perm', '-v7.3');
        fprintf('Saved PhaseMetrics to %s\n', phase_results_file);

    end

    % Append into pooled containers (runs for cached and recomputed)

    all_ps_allLabs  = [all_ps_allLabs,  all_ps];
    all_bs_allLabs  = [all_bs_allLabs,  all_bs];
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

% Duplicate first phase bin
all_hit_probs_plot = all_hit_probs_allLabs;
all_hit_probs_plot(no_phases+1,:,:) = all_hit_probs_plot(1,:,:);


%% Group level statistics

group_level_p    = zeros(n_cond,1);
group_level_comp = zeros(2,1);
pf_h             = false(2,1);

H1stats = struct('observed', {}, 'nullMean', {}, 'nullSD', {}, 'z', {}, ...
                 'critical', {}, 'p', {}, 'nSubjects', {});

for cond = 1:n_cond
    currmean   = mean(all_bs_allLabs(cond,:), 2);
    currsurr   = squeeze(mean(all_bs_perm_allLabs(cond,:,:), 2));
    currsurrmu = mean(currsurr);
    currsurrsd = std(currsurr);
    [~, group_level_p(cond)] = ztest(currmean, currsurrmu, currsurrsd, 'tail', 'right');

    H1stats(cond).observed  = currmean;
    H1stats(cond).nullMean  = currsurrmu;
    H1stats(cond).nullSD    = currsurrsd;
    H1stats(cond).z         = (currmean - currsurrmu) / currsurrsd;
    H1stats(cond).critical  = currsurrmu + 1.645 * currsurrsd;
    H1stats(cond).p         = group_level_p(cond);
    H1stats(cond).nSubjects = size(all_bs_allLabs, 2);
end

H2stats = struct('t', {}, 'df', {}, 'meanDiff', {}, 'p', {}, 'n', {});
ctrlIdx = [2 3];
for k = 1:2
    x = all_bs_allLabs(1,:);
    y = all_bs_allLabs(ctrlIdx(k),:);
    keep = ~isnan(x) & ~isnan(y);
    [~, pk, ~, st] = ttest(x(keep), y(keep), 'tail', 'right');
    group_level_comp(k)  = pk;
    H2stats(k).t         = st.tstat;
    H2stats(k).df        = st.df;
    H2stats(k).meanDiff  = mean(x(keep) - y(keep));
    H2stats(k).p         = pk;
    H2stats(k).n         = sum(keep);
end

H3stats = struct('h', {}, 'meanDiff', {}, 'ciLower', {}, 'ciUpper', {}, 'n', {});
for k = 1:2
    dphi = circ_dist(all_pfs_allLabs(1,:), all_pfs_allLabs(ctrlIdx(k),:));
    dphi = dphi(~isnan(dphi));
    [hk, muk, ulk, llk] = circ_mtest(dphi, 0);
    pf_h(k)            = logical(hk);
    H3stats(k).h        = logical(hk);
    H3stats(k).meanDiff = muk;
    H3stats(k).ciLower  = llk;
    H3stats(k).ciUpper  = ulk;
    H3stats(k).n        = numel(dphi);
end

pf_p = double(pf_h);


%% Reports

condNamesRep = {'Occipital', 'Retinal', 'Cutaneous', 'Sham'};

% H1
reportH1 = fullfile(outH1, 'H1_results.txt');
if isfile(reportH1); delete(reportH1); end
fid = fopen(reportH1, 'a');
fprintf(fid, 'H1: Phasic modulation, observed coefficient against the permuted null\n');
fprintf(fid, 'Labs: %s | Participants: %d | Permutations: %d\n', ...
    mat2str(lab_ids), size(all_bs_allLabs,2), perm);
fprintf(fid, 'One-tailed, critical z = 1.645\n\n');
fprintf(fid, '%-12s %10s %10s %10s %8s %10s %8s\n', ...
    'Condition', 'observed', 'null mean', 'null SD', 'z', 'crit', 'p');
for cond = 1:n_cond
    fprintf(fid, '%-12s %10.5f %10.5f %10.5f %8.3f %10.5f %8s\n', ...
        condNamesRep{cond}, H1stats(cond).observed, H1stats(cond).nullMean, ...
        H1stats(cond).nullSD, H1stats(cond).z, H1stats(cond).critical, ...
        pStr(H1stats(cond).p));
end
fclose(fid);

% H2
reportH2 = fullfile(outH2, 'H2_results.txt');
if isfile(reportH2); delete(reportH2); end
fid = fopen(reportH2, 'a');
fprintf(fid, 'H2: Occipital against each peripheral control\n');
fprintf(fid, 'Paired t tests, one-tailed. Bonferroni corrected alpha = 0.025 for two tests\n\n');
fprintf(fid, '%-24s %6s %10s %8s %8s %8s\n', ...
    'Contrast', 'n', 'mean diff', 't', 'df', 'p');
ctrlNamesRep = {'Retinal', 'Cutaneous'};
for k = 1:2
    fprintf(fid, '%-24s %6d %10.5f %8.3f %8d %8s\n', ...
        ['Occipital > ' ctrlNamesRep{k}], H2stats(k).n, H2stats(k).meanDiff, ...
        H2stats(k).t, H2stats(k).df, pStr(H2stats(k).p));
end
fclose(fid);

% H3
reportH3 = fullfile(outH3, 'H3_results.txt');
if isfile(reportH3); delete(reportH3); end
fid = fopen(reportH3, 'a');
fprintf(fid, 'H3: Difference in preferred phase, occipital against each control\n');
fprintf(fid, 'Circular one-sample test of the paired phase difference against 0\n\n');
fprintf(fid, 'circ_mtest returns h (0 = not rejected, 1 = rejected), not a p-value.\n\n');
fprintf(fid, '%-24s %6s %12s %12s %12s %6s\n', ...
    'Contrast', 'n', 'mean diff', 'CI lower', 'CI upper', 'h');
anyUndefined = false;
for k = 1:2
    if isnan(H3stats(k).ciLower) || isnan(H3stats(k).ciUpper)
        anyUndefined = true;
        fprintf(fid, '%-24s %6d %12.5f %12s %12s %6d\n', ...
            ['Occipital vs ' ctrlNamesRep{k}], H3stats(k).n, H3stats(k).meanDiff, ...
            'undefined', 'undefined', H3stats(k).h);
    else
        fprintf(fid, '%-24s %6d %12.5f %12.5f %12.5f %6d\n', ...
            ['Occipital vs ' ctrlNamesRep{k}], H3stats(k).n, H3stats(k).meanDiff, ...
            H3stats(k).ciLower, H3stats(k).ciUpper, H3stats(k).h);
    end
end
fprintf(fid, '\nAll values in radians.\n');

if anyUndefined
    fprintf(fid, ['\nCI undefined: resultant vector length too small for circ_confmean. ' ...
        'Where the CI is undefined, h = 0 means the test could not be carried out.\n']);
end
fclose(fid);


%% Save pooled metrics


save(fullfile(OutDir, 'H1H2H3_pooled_metrics.mat'), ...
    'all_bs_allLabs', 'all_ps_allLabs', 'all_pfs_allLabs', ...
    'all_hit_probs_allLabs', 'all_bs_perm_allLabs', ...
    'subj_labels_allLabs', 'conditions', 'lab_ids', ...
    'group_level_p', 'group_level_comp', 'pf_p', '-v7.3');

%% Plots

if MakePlots

    [cond_colours, cond_display, cond_xlabels] = tACS_condColours(true);

    n_h1_panels = 3;

    nulls_h1 = cell(1, n_h1_panels);
    dObs_h1  = nan(1, n_h1_panels);
    zVal_h1  = nan(1, n_h1_panels);

    for cond = 1:n_h1_panels
        surrDist       = squeeze(mean(all_bs_perm_allLabs(cond,:,:), 2));
        nulls_h1{cond} = surrDist(:);
        dObs_h1(cond)  = mean(all_bs_allLabs(cond,:), 2);
        zVal_h1(cond)  = (dObs_h1(cond) - mean(surrDist)) / std(surrDist);
    end

    allValues = [vertcat(nulls_h1{:}); dObs_h1(:)];
    xLo = min(allValues); xHi = max(allValues);
    pad = 0.05 * (xHi - xLo);
    edges = linspace(xLo - pad, xHi + pad, 31);

    maxProp = 0;
    for cond = 1:n_h1_panels
        maxProp = max(maxProp, max(histcounts(nulls_h1{cond}, edges, ...
            'Normalization', 'probability')));
    end
    yMaxH1 = ceil(maxProp * 1.15 * 100) / 100;

    fig_h1 = figure('Units', 'centimeters', 'Position', [2 2 16.6 6.5], 'Color', 'w');
    tiledlayout(1, n_h1_panels, 'TileSpacing', 'compact', 'Padding', 'loose');

    letters_h1  = {'A', 'B', 'C'};
    letterXoff  = [-0.20, -0.08, -0.08];

    for cond = 1:n_h1_panels
        ax = nexttile;
        hold(ax, 'on');

        histogram(ax, nulls_h1{cond}, 'BinEdges', edges, ...
            'Normalization', 'probability', ...
            'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 1);

        critVal = mean(nulls_h1{cond}) + 1.645 * std(nulls_h1{cond});
        xline(ax, critVal, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.2);
        xline(ax, dObs_h1(cond), '-', 'Color', cond_colours(cond,:), 'LineWidth', 2);

        xlim(ax, [edges(1) edges(end)]);
        ylim(ax, [0 yMaxH1]);

        xlabel(ax, cond_xlabels{cond});
        if cond == 1
            ylabel(ax, 'Probability');
        else
            ax.YTickLabel = [];
        end

        p_val = group_level_p(cond);
        if p_val < 0.001
            p_str = 'p < 0.001';
        else
            p_str = sprintf('p = %.3f', p_val);
        end
        title(ax, sprintf('%s\nz = %.2f, %s', cond_display{cond}, zVal_h1(cond), p_str), ...
            'FontName', 'Times New Roman', 'FontSize', 10, 'FontWeight', 'normal');

        set(ax, 'FontName', 'Times New Roman', 'FontSize', 10, 'Box', 'off', ...
            'TickDir', 'out', 'LineWidth', 0.75, 'Layer', 'top');

        text(ax, letterXoff(cond), 1.06, letters_h1{cond}, 'Units', 'normalized', ...
            'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Clipping', 'off');

        hold(ax, 'off');
    end

    saveas(fig_h1, fullfile(outH1, 'H1_observed_vs_null.png'));
    savefig(fig_h1, fullfile(outH1, 'H1_observed_vs_null.fig'));


    fig_phase = figure('Color', 'w');
    ax_ph = axes(fig_phase);
    hold(ax_ph, 'on');
    meanCurves = squeeze(mean(all_hit_probs_plot, 3));
    for cond = 1:3
        plot(ax_ph, -pi:pi/4:pi, meanCurves(:,cond), '-', ...
            'Color', cond_colours(cond,:), 'LineWidth', 1.5);
    end
    hold(ax_ph, 'off');
    xlabel(ax_ph, 'tACS phase (rad)');
    ylabel(ax_ph, 'Detection probability');
    legend(ax_ph, cond_display(1:3), 'Location', 'best');
    title(ax_ph, 'Detection probability by tACS phase', 'FontWeight', 'normal');
    set(ax_ph, 'FontName', 'Times New Roman', 'FontSize', 10, 'Box', 'off', ...
        'TickDir', 'out', 'LineWidth', 0.75);
    saveas(fig_phase, fullfile(outH1, 'H1_hitprob_by_phase.png'));
    savefig(fig_phase, fullfile(outH1, 'H1_hitprob_by_phase.fig'));

    %  H2 plot

    n_active      = 3;  % A, B, C only
    [cond_colors_h2, active_labels] = tACS_condColours(false);

    fig_h2 = figure('Units','normalized','Position',[0.15 0.2 0.45 0.55]);
    ax_h2  = axes(fig_h2);

    groups_h2 = cell(1, n_active);
    for k = 1:n_active
        groups_h2{k} = all_bs_allLabs(k,:).';
    end

    % Bonferroni corrected comparisons of the occipital condition against each
    % peripheral control
    brackets_h2(1).from  = 1;
    brackets_h2(1).to    = 2;
    brackets_h2(1).label = pval_label(group_level_comp(1));
    brackets_h2(2).from  = 1;
    brackets_h2(2).to    = 3;
    brackets_h2(2).label = pval_label(group_level_comp(2));

    rng(42);
    tACS_raincloudPanel(ax_h2, groups_h2, active_labels, cond_colors_h2, ...
        'Brackets', brackets_h2, ...
        'YLabel',   'Modulation strength (regression coefficient)', ...
        'XLabel',   'Condition');
    title(ax_h2, 'H2: Phasic modulation strength by condition', 'FontWeight', 'bold');

    saveas(fig_h2, fullfile(outH2, 'H2_modulation_strength.png'));
    savefig(fig_h2, fullfile(outH2, 'H2_modulation_strength.fig'));

end 


%% Output

OUT            = struct();
OUT.lab_ids    = lab_ids;
OUT.conditions = conditions;
OUT.perm       = perm;
OUT.subj_labels = subj_labels_allLabs;

OUT.all_ps       = all_ps_allLabs;
OUT.all_bs       = all_bs_allLabs;
OUT.all_pfs      = all_pfs_allLabs;
OUT.all_bs_perm  = all_bs_perm_allLabs;
OUT.all_hit_probs = all_hit_probs_allLabs;

OUT.group_level_p    = group_level_p;
OUT.group_level_comp = group_level_comp;
OUT.pf_p             = pf_p;   
OUT.H1stats          = H1stats;
OUT.H2stats          = H2stats;
OUT.H3stats          = H3stats;
OUT.OutDir           = OutDir;

end 

%% Helper functions


function s = pStr(p)
% Format a p-value at a precision that reflects what it can actually resolve.
    if isnan(p)
        s = 'NA';
    elseif p < 0.0001
        s = '<.0001';
    else
        s = sprintf('%.4f', p);
    end
end


function label = pval_label(p)
% Format a p-value as a significance label for plot brackets.
    if p < 0.001
        label = '***';
    elseif p < 0.01
        label = '**';
    elseif p < 0.025
        label = '*';
    else
        label = sprintf('p = %.3f', p);
    end
end
