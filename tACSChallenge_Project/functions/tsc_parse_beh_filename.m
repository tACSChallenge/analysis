function info = tsc_parse_beh_filename(fname)
%% Parse block number and condition label from a beh filename.
%% script originally written by Anne-Fiona Griesfeller


    info = struct('blockNum', NaN, 'cond', "");

    % Block number
    tokB = regexp(fname, '_B(\d+)_', 'tokens', 'once');
    if ~isempty(tokB)
        info.blockNum = str2double(tokB{1});
    end

    % Condition
    if contains(fname, '_Sh_') || contains(fname, '_Sh.')
        info.cond = "Sh";
    elseif contains(fname, '_A_')|| contains(fname, '_A.')
        info.cond = "A";
    elseif contains(fname, '_B_')|| contains(fname, '_B.')
        info.cond = "B";
    elseif contains(fname, '_C_')|| contains(fname, '_C.')
        info.cond = "C";
    else
        info.cond = "";
    end
end
