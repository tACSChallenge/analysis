function perf = tsc_perf_from_sorted(LEDLat, RespLat, Fs, nSamplesBlock, targetWindow)
%% Compute performance metrics from target and response onsets.
%% script originally written by Anne-Fiona Griesfeller

% DEFINITIONS
%   Hit: first button press 150-1200 ms after target onset
%   Miss: no press in that window
%   False alarm: any press not preceded by a target within 1200 ms
%   Non-target windows: before the first target; from target_i + 1200 ms to the
%   next target; from the last target + 1200 ms to the end of the block

    % Ensure column vectors and sorted
    LEDLat  = sort(LEDLat(:));
    RespLat = sort(RespLat(:));

    hitStartSamp = round(0.150 * Fs);
    hitEndSamp   = round(1.200 * Fs);
    faLagSamp    = round(1.200 * Fs);

    nT_all = numel(LEDLat);

    % Resolve target window
    if nargin < 5 || isempty(targetWindow)
        s = 1;
        e = nT_all;
    else
        s = max(1, targetWindow(1));
        e = min(nT_all, targetWindow(2));
    end

    perf = struct();

    % No targets: every press counts as a false alarm over one window.
    if nT_all == 0
        perf.nTargets = 0;
        perf.nHits = 0;
        perf.hitRate = NaN;
        perf.meanRT_sec = NaN;
        perf.sumRT_secHits = 0;
        perf.nHitRT = 0;
        perf.nFA = numel(RespLat);
        perf.nNonTargetWindows = 1;   % whole block as a single non-target window
        perf.faRate = perf.nFA / perf.nNonTargetWindows;
        perf.dprime = NaN;
        perf.criterion = NaN;
        return
    end

    % Requested window contains no targets.
    if s > e
        perf.nTargets = 0;
        perf.nHits = 0;
        perf.hitRate = NaN;
        perf.meanRT_sec = NaN;
        perf.sumRT_secHits = 0;
        perf.nHitRT = 0;
        perf.nFA = 0;
        perf.nNonTargetWindows = 0;
        perf.faRate = NaN;
        perf.dprime = NaN;
        perf.criterion = NaN;
        return
    end

    % Left bound: block start for the first window, else the onset of target s.
    if s == 1
        winStart = 1;
    else
        winStart = LEDLat(s);
    end
    % Right bound: just before target e+1, else the block end.
    if e < nT_all
        winEnd = LEDLat(e+1) - 1;
    else
        winEnd = nSamplesBlock;
    end

    % Restrict to the window and re-base sample indices to winStart.
    offset = winStart - 1;
    LEDLat  = LEDLat(s:e) - offset;
    RespLat = RespLat(RespLat >= winStart & RespLat <= winEnd) - offset;
    nSamplesBlock = winEnd - offset;

    perf.nTargets = numel(LEDLat);

    %% Hits and RT (hits only)
    isRespUsedAsHit = false(size(RespLat));
    hitRT_samp = nan(perf.nTargets, 1);
    nHits = 0;

    for iT = 1:perf.nTargets
        t0 = LEDLat(iT);
        winStartH = t0 + hitStartSamp;
        winEndH   = t0 + hitEndSamp;

        idx = find(RespLat >= winStartH & RespLat <= winEndH & ~isRespUsedAsHit, 1, 'first');
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
        perf.meanRT_sec    = mean(hitRT_sec);
        perf.sumRT_secHits = sum(hitRT_sec);
        perf.nHitRT        = numel(hitRT_sec);
    else
        perf.meanRT_sec    = NaN;
        perf.sumRT_secHits = 0;
        perf.nHitRT        = 0;
    end

    %% False alarms
    % Exclude presses used as hits, then apply the preceding-target rule.
    candIdx = find(~isRespUsedAsHit);
    nFA = 0;

    for k = 1:numel(candIdx)
        pTime = RespLat(candIdx(k));
        lastT = find(LEDLat < pTime, 1, 'last');
        if isempty(lastT)
            nFA = nFA + 1;
        elseif (pTime - LEDLat(lastT)) > faLagSamp
            nFA = nFA + 1;
        end
    end

    perf.nFA = nFA;

    %% Non target response windows count
    nWindows = 0;

    % Pre window: start of window to just before first target. After re-basing,
    % LEDLat(1) > 1 only for the first bin, so mid-block bins add no pre-window.
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

    % After last target in the window
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
    faRateCorr  = edgeCorrect(perf.faRate,  perf.nNonTargetWindows);

    if ~isnan(hitRateCorr) && ~isnan(faRateCorr)
        zH = norminv(hitRateCorr);
        zF = norminv(faRateCorr);
        perf.dprime    = zH - zF;
        perf.criterion = -0.5 * (zH + zF);
    else
        perf.dprime    = NaN;
        perf.criterion = NaN;
    end
end

function r = edgeCorrect(rIn, N)
% If rate = 1, use (N-1)/N
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
