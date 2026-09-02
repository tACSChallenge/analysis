function trials_sorted = tACSChallenge_SortData(data_path, subj, conditions)
%% script originally written by Florian Kasten, University of Oldenburg
%% modified by Benedikt Zoefel, CNRS Toulouse, in October 2021, April 2022 (clean tACS signal),
%% and June 2022 (correction for imperfect stimulation frequency)
%% modified by Florian Kasten in March 2022 (correction of target onsets)
%% modified by BZ in June 2025 (check for number of blocks)
%% modified by Anne-Fiona Griesfeller (left+right button presses, target detection)

if numel(conditions) ~= 4
    error(['conditions must contain exactly four filename patterns, one per ' ...
        'condition, in the order A, B, C, sham. Received %d.'], numel(conditions));
end

trials_sorted = cell(4,1);
phase_bins = 0:pi/4:2*pi;

for c = 1:4

    % subject folder plus beh subfolder
    subj_path = fullfile(data_path, subj, 'beh');

    curr_files = dir(fullfile(subj_path, conditions{c}));
    if ~isempty(curr_files)
        [~, ord] = sort(string({curr_files.name}));
        curr_files = curr_files(ord);
    end

    curr_cond = conditions{c};

    trialcounter = 0;
    trials_sorted{c} = [];
    
    isSham = contains(curr_cond, 'Sh');
    if ~isSham && numel(curr_files) ~= 3
        error('number of blocks for at least one condition is not three as expected');
    end
    if isSham && numel(curr_files) ~= 4
        error('number of sham blocks is not four as expected');
    end

    for b = 1:length(curr_files)

        %% import data from this block
        data = tACSChallenge_ImportData([curr_files(b).folder, filesep, curr_files(b).name]);

        Fs = data.Fs;

        %% get onsets of button presses from both buttons
       % left button
respOnsetsL = diff(data.L_Button);
respOnsetsL(respOnsetsL < 0) = 0;
RespLatL = find(respOnsetsL > 0);

% right button
respOnsetsR = diff(data.R_Button);
respOnsetsR(respOnsetsR < 0) = 0;
RespLatR = find(respOnsetsR > 0);

% combine both: any press on either button counts as a response
RespLat = sort(unique([RespLatL; RespLatR]));

        %% remove offset using central LED
        data.LEDs = data.LEDs - repmat(data.LEDs(:,1), 1, size(data.LEDs,2));

%% merge LED signals (handle increase AND decrease)
        LED = max(abs(data.LEDs), [], 2);

%% define baseline
        baseline = median(LED);

% detect deviation from baseline
        isTarget = LED ~= baseline;

% detect onsets only
        LEDLat = find(diff(isTarget) == 1);

% safety check
        if isempty(LEDLat)
        continue
    end

        %% tACS signal
        tACS = data.tACS;
        tACS = tACS - mean(tACS);
        
        % clean the signal
        start_idx = max(LEDLat(1) - round(1*Fs), 1);  % 1 s before first target, safety against negative indices
        tACS_to_fit = tACS(start_idx:end);

        % estimate its frequency (should be 9–11 Hz)
        [filt_b,filt_a] = butter(2,[9/(Fs/2) 11/(Fs/2)],'bandpass');
        tACS_filtered = filtfilt(filt_b,filt_a,tACS_to_fit);
        tACS_filtered_h = hilbert(tACS_filtered);
        ifq = instfreq(tACS_filtered, Fs);

        % reconstruct tACS sine
        t_fit = 0:1/Fs:(length(tACS_to_fit)-1)/Fs;
        % make a sine with the same average frequency and the same phase
        tACS_sine = cos(2*pi*t_fit*mean(ifq(4:end-3)) + angle(tACS_filtered_h(round(Fs)*10+1)));
        tACS_phase = angle(hilbert(tACS_sine));

        %% go through all trials in this block
        for i = 1:length(LEDLat)

            trialcounter = trialcounter + 1;
            curr_t = LEDLat(i);                         % target time in samples
            trials_sorted{c}(trialcounter,1) = curr_t;  % column 1: target sample index

            
            % Check for hits in preregistered hit window: 150–1200 ms after target
            winStart = curr_t + round(0.150 * Fs);
            winEnd   = curr_t + round(1.200 * Fs);
            idx_resp = RespLat >= winStart & RespLat <= winEnd;

            if any(idx_resp)
                % hit
                trials_sorted{c}(trialcounter,2) = 1;   % column 2: hit = 1
                thisRT = RespLat(find(idx_resp,1,'first')) - curr_t;
                trials_sorted{c}(trialcounter,3) = thisRT;  % column 3: RT in samples
            else
                % miss
                trials_sorted{c}(trialcounter,2) = 0;   % column 2: miss = 0
                trials_sorted{c}(trialcounter,3) = 0;   % column 3: RT = 0 for misses
            end

            % tACS phase at target time
            % index into tACS_phase: relative to start_idx and LEDLat(1)
            t_tACS = curr_t - start_idx + 1;
            if t_tACS < 1 || t_tACS > length(tACS_phase)
                % safety check: if out of bounds, set NaN
                trials_sorted{c}(trialcounter,4) = NaN;
                trials_sorted{c}(trialcounter,5) = NaN;
            else
                phi = tACS_phase(t_tACS);
                trials_sorted{c}(trialcounter,4) = phi;  % column 4: phase at target

                % nearest phase bin (for visualisation)
                [~, phase_bin] = min(abs(phi + pi - phase_bins));
                trials_sorted{c}(trialcounter,5) = phase_bin;  % column 5: phase bin index
            end

        end  % trials in block

    end  % blocks


end  % conditions
