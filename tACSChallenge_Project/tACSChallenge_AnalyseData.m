function [all_ps,all_bs,all_hit_probs,all_bs_perm,all_pfs,group_level_p,group_level_comp,pf_p] = tACSChallenge_AnalyseData(data_path, subjs, conditions, perm, varargin)
%% script originally written by Benedikt Zoefel, CNRS Toulouse, in October 2021
%% modified in April 22 and June 25 (minor fixes)
%% added permutations and preferred phases in April 24

%% requires circular statistics toolbox
%% Modified by Anne-Fiona Griesfeller in Jan 26

clc; close all;

p = inputParser;
p.addParameter('OutDir', '', @(x)ischar(x) || isstring(x));
p.addParameter('MakeOutputs', true, @(x)islogical(x) || isnumeric(x)); % plots + txt reports
p.parse(varargin{:});

OutBaseDir = string(p.Results.OutDir);
MakeOutputs = logical(p.Results.MakeOutputs);

if OutBaseDir == ""
    OutBaseDir = string(data_path);
end

outH1 = fullfile(OutBaseDir, 'H1');
outH2 = fullfile(OutBaseDir, 'H2');
outH3 = fullfile(OutBaseDir, 'H3');

if MakeOutputs
    if ~isfolder(outH1); mkdir(outH1); end
    if ~isfolder(outH2); mkdir(outH2); end
    if ~isfolder(outH3); mkdir(outH3); end
end


no_phases = 8;

% Number of true experimental conditions (A, B, C, Sham)
n_cond = length(conditions) - 4;

all_ps        = zeros(n_cond,length(subjs));
all_bs        = zeros(n_cond,length(subjs));
all_hit_probs = zeros(no_phases,n_cond,length(subjs));
all_bs_perm   = zeros(n_cond,length(subjs),perm);
all_pfs       = zeros(n_cond,length(subjs)); % preferred phases

for s = 1:length(subjs)
    curr_data = tACSChallenge_SortData(data_path, subjs{s}, conditions);
    [all_ps(:,s), all_bs(:,s), all_hit_probs(:,:,s), all_bs_perm(:,s,:), all_pfs(:,s)] = ...
        tACSChallenge_EvalData(curr_data,perm,s); 
end

all_hit_probs(no_phases+1,:,:) = all_hit_probs(1,:,:); % duplicate first phase bin

%% Plot detection probability vs phase (visualisation only)
if MakeOutputs
    figure
    plot(-pi:pi/4:pi, squeeze(mean(all_hit_probs,3)));
    xlabel('tACS phase'); ylabel('detection probability');
    legend(conditions(1:n_cond));  % only A,B,C,Sham

    saveas(gcf, fullfile(outH1, 'H1_hitprob_by_phase.png'));
    savefig(gcf, fullfile(outH1, 'H1_hitprob_by_phase.fig'));
end

%% Plot phasic modulation strength per condition

if MakeOutputs
    figure
    bar(mean(all_bs,2))
    xlabel('conditions');
    xticklabels(conditions(1:n_cond));
    ylabel('modulation strength');

    saveas(gcf, fullfile(outH1, 'H1_modulation_strength.png'));
    savefig(gcf, fullfile(outH1, 'H1_modulation_strength.fig'));

    saveas(gcf, fullfile(outH2, 'H2_modulation_strength.png'));
    savefig(gcf, fullfile(outH2, 'H2_modulation_strength.fig'));
end

% Group level statistics
group_level_p    = zeros(n_cond,1);
group_level_comp = zeros(2,1); % A vs B, A vs C
pf_p             = zeros(2,1);

%% H1: phasic effect per condition
for cond = 1:n_cond
    currmean   = mean(all_bs(cond,:),2);
    currsurr   = squeeze(mean(all_bs_perm(cond,:,:),2));
    currsurrmu = mean(currsurr);
    currsurrsd = std(currsurr);
    [~,group_level_p(cond)] = ztest(currmean,currsurrmu,currsurrsd,'tail','right');
end

%% H2: contrasts A > B and A > C (assumes order A=1, B=2, C=3, Sham=4)
[~,group_level_comp(1)] = ttest(all_bs(1,:),all_bs(2,:),'tail','right');
[~,group_level_comp(2)] = ttest(all_bs(1,:),all_bs(3,:),'tail','right');

%% H3: preferred phase differences (A vs B, A vs C)
pf_p(1) = circ_mtest(circ_dist(all_pfs(1,:),all_pfs(2,:)),0);
pf_p(2) = circ_mtest(circ_dist(all_pfs(1,:),all_pfs(3,:)),0);

%% Save hypothesis specific reports

if MakeOutputs
% H1 report
        reportH1 = fullfile(outH1, 'H1_results.txt');
    if isfile(reportH1); delete(reportH1); end
        fid = fopen(reportH1, 'a');
        fprintf(fid, 'H1: Phasic modulation per condition (z test vs permutation surrogate)\n');
    for cond = 1:n_cond
        fprintf(fid, 'Condition %d (A=1, B=2, C=3, Sham=4): p = %.10f\n', cond, group_level_p(cond));
    end
    fclose(fid);

% H2 report
        reportH2 = fullfile(outH2, 'H2_results.txt');
    if isfile(reportH2); delete(reportH2); end
        fid = fopen(reportH2, 'a');
        fprintf(fid, 'H2: Contrasts (paired t tests, right tailed)\n');
        fprintf(fid, 'A > B: p = %.10f\n', group_level_comp(1));
        fprintf(fid, 'A > C: p = %.10f\n', group_level_comp(2));
    fclose(fid);

% H3 report
        reportH3 = fullfile(outH3, 'H3_results.txt');
    if isfile(reportH3); delete(reportH3); end
        fid = fopen(reportH3, 'a');
        fprintf(fid, 'H3: Preferred phase differences (circular m tests)\n');
        fprintf(fid, 'A vs B: p = %.10f\n', pf_p(1));
        fprintf(fid, 'A vs C: p = %.10f\n', pf_p(2));
    fclose(fid);

end 

phase_results_file = fullfile(data_path, 'PhaseMetrics.mat');

save(phase_results_file, ...
     'all_bs', 'all_ps', 'all_pfs', 'all_hit_probs', 'all_bs_perm', ...
     'subjs', 'conditions');

fprintf('Saved phase metrics to %s\n', phase_results_file);
