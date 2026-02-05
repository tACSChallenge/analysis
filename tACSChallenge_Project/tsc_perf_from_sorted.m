function perf = tsc_perf_from_sorted(LEDLat, RespLat, Fs, nSamplesBlock)
%% Compute performance metrics from target and response onsets.
%% script originally written by Anne-Fiona Griesfeller, CNRS Toulouse, in Jan 26
% INPUTS (all in samples except Fs)
%   LEDLat        vector of target onset sample indices
%   RespLat       vector of response onset sample indices (both buttons combined)
%   Fs            sampling rate in Hz
%   nSamplesBlock number of samples in the block (used for last non target window)

% DEFINITIONS
%   Hit: first button press 150–1200 ms after target onset
%   Miss: no press in that window
%   False Alarm (FA): any button press not preceded by a target within previous 1200 ms
%   Non target response windows:
%     Pre window before first target
%     For each target i: (target_i + 1200 ms + 1) to next target
%     After last target: (last target + 1200 ms + 1) to end of block

% OUTPUT
%   nTargets, nHits, hitRate, meanRT_sec
%   nFA, nNonTargetWindows, faRate
%   dprime, criterion

    % Ensure column vectors and sorted
    LEDLat = sort(LEDLat(:));
    RespLat = sort(RespLat(:));

    hitStartSamp = round(0.150 * Fs);
    hitEndSamp   = round(1.200 * Fs);
    faLagSamp    = round(1.200 * Fs);

    % Initialize outputs
    perf = struct();
    perf.nTargets = numel(LEDLat);

    % If no targets, SDT is undefined under the preregistered hit rate definition
    if perf.nTargets == 0
        perf.nHits = 0;
        perf.hitRate = NaN;
        perf.meanRT_sec = NaN;

        perf.nFA = numel(RespLat);
        perf.nNonTargetWindows = 1; % whole block as a single non target window
        perf.faRate = perf.nFA / perf.nNonTargetWindows;

        perf.dprime = NaN;
        perf.criterion = NaN;
        return
    end

    
    %% Hits and RT (hits only)
   
    isRespUsedAsHit = false(size(RespLat));
    hitRT_samp = nan(perf.nTargets, 1);
    nHits = 0;

    for iT = 1:perf.nTargets
        t0 = LEDLat(iT);
        winStart = t0 + hitStartSamp;
        winEnd   = t0 + hitEndSamp;

        idx = find(RespLat >= winStart & RespLat <= winEnd & ~isRespUsedAsHit, 1, 'first');
        if ~isempty(idx)
            isRespUsedAsHit(idx) = true;
            nHits = nHits + 1;
            hitRT_samp(iT) = RespLat(idx) - t0;
        end
    end

    perf.nHits = nHits;
    perf.hitRate = nHits / perf.nTargets;

    if nHits > 0
    hitRT_sec = hitRT_samp(~isnan(hitRT_samp)) ./ Fs;
    perf.meanRT_sec = mean(hitRT_sec);

    % for aggregation across blocks
    perf.sumRT_secHits = sum(hitRT_sec);
    perf.nHitRT = numel(hitRT_sec);
else
    perf.meanRT_sec = NaN;

    % for aggregation across blocks
    perf.sumRT_secHits = 0;
    perf.nHitRT = 0;
end

 
    %% False alarms
    % Exclude presses used as hits, then apply preceding target rule
    
    candIdx = find(~isRespUsedAsHit);
    nFA = 0;

    for k = 1:numel(candIdx)
        pTime = RespLat(candIdx(k));

        lastT = find(LEDLat < pTime, 1, 'last');
        if isempty(lastT)
            nFA = nFA + 1;
        else
            if (pTime - LEDLat(lastT)) > faLagSamp
                nFA = nFA + 1;
            end
        end
    end

    perf.nFA = nFA;

    
    %% Non target response windows count
    
    nWindows = 0;

    % Pre window: start of block to just before first target
    if LEDLat(1) > 1
        nWindows = nWindows + 1;
    end

    % Between targets
    for iT = 1:(perf.nTargets - 1)
        wStart = LEDLat(iT) + faLagSamp + 1;
        wEnd   = LEDLat(iT + 1);
        if wStart <= wEnd
            nWindows = nWindows + 1;
        end
    end

    % After last target
    wStart = LEDLat(end) + faLagSamp + 1;
    wEnd   = nSamplesBlock;
    if wStart <= wEnd
        nWindows = nWindows + 1;
    end

    perf.nNonTargetWindows = nWindows;

    if nWindows > 0
        perf.faRate = nFA / nWindows;
    else
        perf.faRate = NaN;
    end

    
    %% Edge correction and SDT
    
    hitRateCorr = edgeCorrect(perf.hitRate, perf.nTargets);
    faRateCorr  = edgeCorrect(perf.faRate, perf.nNonTargetWindows);

    if ~isnan(hitRateCorr) && ~isnan(faRateCorr)
        zH = norminv(hitRateCorr);
        zF = norminv(faRateCorr);
        perf.dprime = zH - zF;
        perf.criterion = -0.5 * (zH + zF);
    else
        perf.dprime = NaN;
        perf.criterion = NaN;
    end
end

function r = edgeCorrect(rIn, N)
% Preregistered correction:
% If rate = 1, use (N−1)/N
% If rate = 0, use 1/N
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
