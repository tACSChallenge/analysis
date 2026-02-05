%% write_eeg_meta.m
%
% This script generates .json and _channels.tsv files containing EEG metadata.
% 
% NOTE: This script is designed for BrainVision EEG files (.vhdr + .eeg).
% For other EEG formats, it may not work.
%
% Run this after you have finished running BIDSify.m (i.e., after files are
% already organised into sub-*/eeg/ folders).
%
% IMPORTANT: You MUST update the meta.* fields in Part 2 to match your EEG device.

%% -------- Part 1: Set paths / IDs --------
% Example:
data_path = "L03";
labnum = 3;
subnum = 1:20; 

% (Assume data_path, labnum, subnum already exist in your main script)

%% -------- Part 2: Update these meta fields --------
meta = struct();
meta.InstitutionName = "your institution";
meta.Manufacturer = "Brain Products";
meta.ManufacturersModelName = "BrainAmp MR plus";
meta.PowerLineFrequency = 50;
meta.RecordingType = "continuous";

meta.EEGReference = "Cz";
meta.EEGGround = "FPz";
meta.SoftwareFilters = "n/a";
meta.EEGChannelCount = 1;
meta.EOGChannelCount = 0;
meta.ECGChannelCount = 0;
meta.EMGChannelCount = 0;
meta.TriggerChannelCount = 0;
meta.MISCChannelCount = 0;

%% -------- Part 3: Loop over subjects and write JSON + channels.tsv --------
for s = subnum(:)'

    prefix = sprintf('sub-L%02d_S%02d', labnum, s);
    eeg_dir = fullfile(data_path, prefix, 'eeg');

    vhdr_list = dir(fullfile(eeg_dir, '*.vhdr'));

    for k = 1:numel(vhdr_list)

        vhdr_path = fullfile(vhdr_list(k).folder, vhdr_list(k).name);
        [~, baseName, ~] = fileparts(vhdr_path);

        % --- Read minimal header fields ---
        info = parse_vhdr_minimal(vhdr_path);

        % --- SamplingFrequency (Hz) from SamplingInterval (us) ---
        fs = 1e6 / info.SamplingInterval_us;

        % --- EEG file path (use DataFile field if possible) ---
        eeg_path = fullfile(fileparts(vhdr_path), char(info.DataFile));
        if ~exist(eeg_path, 'file')
            eeg_path = fullfile(fileparts(vhdr_path), [baseName '.eeg']);
        end

        % --- RecordingDuration (s) from file size ---
        bytesPerSample = bytes_per_sample_from_format(info.BinaryFormat);
        nBytes = double(dir(eeg_path).bytes);
        nSamples = floor(nBytes / (bytesPerSample * double(info.NumberOfChannels)));
        duration = nSamples / fs;

        % --- TaskName from filename: contains "pre" or "post" (case-insensitive) ---
        baseLower = lower(string(baseName));
        if contains(baseLower, "pre")
            prepost = "pre";
        else
            prepost = "post";
        end
        taskName = "EC resting-state recording " + prepost + " tACS challenge behaviour task";

        % --- Write JSON (same basename as .vhdr) ---
        J = meta;
        J.TaskName = taskName;
        J.SamplingFrequency = fs;
        J.RecordingDuration = duration;

        json_path = fullfile(eeg_dir, baseName + ".json");
        write_json(json_path, J);

        % --- Write channels.tsv (fixed content) ---
        ch_path = fullfile(eeg_dir, baseName + "_channels.tsv");
        write_channels_tsv(ch_path);
    end
end

disp('Done: JSON and channels.tsv generated.');

%% ================= Helper functions =================

function info = parse_vhdr_minimal(vhdr_path)
    txt = fileread(vhdr_path);
    lines = regexp(txt, '\r\n|\n|\r', 'split');

    SamplingInterval_us = NaN;
    NumberOfChannels    = NaN;
    DataFile            = "";
    BinaryFormat        = "";

    for i = 1:numel(lines)
        line = strtrim(lines{i});
        if line == "" || startsWith(line, ';') || startsWith(line, '[')
            continue;
        end
        parts = regexp(line, '=', 'split', 'once');
        if numel(parts) ~= 2, continue; end

        key = lower(strtrim(parts{1}));
        val = strtrim(parts{2});

        switch key
            case 'samplinginterval'
                SamplingInterval_us = str2double(val);
            case 'numberofchannels'
                NumberOfChannels = str2double(val);
            case 'datafile'
                DataFile = string(val);
            case 'binaryformat'
                BinaryFormat = string(val);
        end
    end

    info = struct();
    info.SamplingInterval_us = SamplingInterval_us;
    info.NumberOfChannels = NumberOfChannels;
    info.DataFile = DataFile;
    info.BinaryFormat = BinaryFormat;
end

function bps = bytes_per_sample_from_format(fmt)
    fmt = upper(strtrim(string(fmt)));
    switch fmt
        case {"INT_8","UINT_8"}
            bps = 1;
        case {"INT_16","UINT_16"}
            bps = 2;
        case {"INT_24","UINT_24"}
            bps = 3;
        case {"INT_32","UINT_32","IEEE_FLOAT_32"}
            bps = 4;
        case {"IEEE_FLOAT_64"}
            bps = 8;
        otherwise
            error('Unsupported BinaryFormat: %s', fmt);
    end
end

function write_json(json_path, J)
    try
        jsonText = jsonencode(J, 'PrettyPrint', true);
    catch
        jsonText = jsonencode(J);
    end
    fid = fopen(json_path, 'w', 'n', 'UTF-8');
    fprintf(fid, '%s\n', jsonText);
    fclose(fid);
end

function write_channels_tsv(ch_path)
    fid = fopen(ch_path, 'w', 'n', 'UTF-8');
    fprintf(fid, 'name\ttype\tunits\n');
    fprintf(fid, 'POz\tEEG\tµV\n');
    fclose(fid);
end
