function [ps,bs,hit_probs,bs_perm,pf] = tACSChallenge_EvalData(sorted_data,perm,subjno)
%% script originally written by Benedikt Zoefel, CNRS Toulouse, in October 2021
%% April 24: added permutations and preferred phase

%% sorted_data is the output from tACSChallenge_SortData

n_cond = length(sorted_data);% number of conditions to be analysed
ps = zeros(n_cond,1);
bs = zeros(n_cond,1);
bs_perm = zeros(n_cond,perm);
hit_probs = zeros(8,n_cond);
pf = zeros(n_cond,1); % preferred phase

rng(subjno);

for c = 1:n_cond
    
    curr_data = sorted_data{c};
    
    if size(curr_data,2) >=4 && any(curr_data(:,2))
        
    glm_input = zeros(size(curr_data,1),3);
    glm_input(:,1) = sin(curr_data(:,4));
    glm_input(:,2) = cos(curr_data(:,4));
    glm_input(:,3) = curr_data(:,2);
               
    g1 = fitglm(glm_input(:,1:2),glm_input(:,3),'distribution','binomial','link','logit');
    %p-value for phasic modulation of perception for single participant
    ps(c) = coefTest(g1);
    
    % extract the preferred phase fom the fit
    pf(c) = -1.*angle(g1.Coefficients.Estimate(3)+g1.Coefficients.Estimate(2)*1i);
    
    % regression coefficient
    curr_b = mnrfit(glm_input(:,1:2),glm_input(:,3)+1);
    % the higher the regression coefficient, the stronger the phasic
    % modulation of perception
    bs(c) = sqrt(curr_b(2)^2 + curr_b(3)^2);
    
    %% permutation procedure - shuffle phases
    for shuff = 1:perm
        shuffled_phases = curr_data(randperm(size(curr_data,1)),4);
        glm_input(:,1) = sin(shuffled_phases);
        glm_input(:,2) = cos(shuffled_phases);
        glm_input(:,3) = curr_data(:,2);
                  
        curr_b = mnrfit(glm_input(:,1:2),glm_input(:,3)+1);
        % the higher the regression coefficient, the stronger the phasic
        % modulation of perception
        bs_perm(c,shuff) = sqrt(curr_b(2)^2 + curr_b(3)^2);
    end
        
    %% divide detection probability into phase bins (mostly for visualisation) - number of phases fixed to 8 here
    for p = 2:8
        hit_probs(p,c) = sum(curr_data(curr_data(:,5)==p,2))/sum(curr_data(:,5)==p);
    end
    hit_probs(1,c) = (sum(curr_data(curr_data(:,5)==1,2))+sum(curr_data(curr_data(:,5)==9,2)))/ ...
        (sum(curr_data(:,5)==1)+sum(curr_data(:,5)==9));
    end
end



