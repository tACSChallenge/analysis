function bf = tACS_bayesFactorBIC(lmeFull, tbl, reducedFormula)
% Bayes factor for a nested pair of mixed models, BIC approximation
% (Wagenmakers, 2007): BF01 = exp((BIC_full - BIC_reduced)/2).
%
%   lmeFull         fitted LinearMixedModel
%   tbl             table the full model was fitted to
%   reducedFormula  formula for the nested reduced model

    bf = struct('BIC_full', NaN, 'BIC_reduced', NaN, 'BF01', NaN, 'BF10', NaN, ...
                'evidence', "not computed", 'formula_reduced', string(reducedFormula), ...
                'nObs', NaN, 'ok', false);

    if isempty(lmeFull)
        return
    end

    try
        % BIC is only comparable across different fixed effects under ML.
        if ~strcmpi(lmeFull.FitMethod, 'ML')
            warning('Model fitted by %s, not ML. No Bayes factor.', lmeFull.FitMethod);
            return
        end

        lmeReduced = fitlme(tbl, char(reducedFormula));

        if lmeReduced.NumObservations ~= lmeFull.NumObservations
            warning('Models used %d and %d observations. No Bayes factor.', ...
                lmeFull.NumObservations, lmeReduced.NumObservations);
            return
        end

        bf.BIC_full    = lmeFull.ModelCriterion.BIC;
        bf.BIC_reduced = lmeReduced.ModelCriterion.BIC;
        bf.nObs        = lmeFull.NumObservations;

        bf.BF01 = exp((bf.BIC_full - bf.BIC_reduced) / 2);
        bf.BF10 = 1 / bf.BF01;
        bf.evidence = evidenceLabel(bf.BF01, bf.BF10);
        bf.ok = true;

    catch ME
        warning('Bayes factor could not be computed: %s', ME.message);
    end
end


function label = evidenceLabel(BF01, BF10)
% Jeffreys scale, after Lee and Wagenmakers (2013).
    if BF01 >= 1
        b = BF01; direction = "for the null";
    else
        b = BF10; direction = "for the alternative";
    end

    if     b < 1;   strength = "no";
    elseif b < 3;   strength = "anecdotal";
    elseif b < 10;  strength = "moderate";
    elseif b < 30;  strength = "strong";
    elseif b < 100; strength = "very strong";
    else;           strength = "extreme";
    end

    label = strength + " evidence " + direction;
end
