function info = tsc_parse_beh_filename(fname)
%% Parse block number and condition label from a beh filename.
%% script originally written by Anne-Fiona Griesfeller, CNRS Toulouse, in Jan 26
% OUTPUT:
%   info.blockNum   (double)
%   info.cond       ('A','B','C','Sh' or '' if unknown)

    info = struct('blockNum', NaN, 'cond', "");

    % Block number: look for _B<number>_
    tokB = regexp(fname, '_B(\d+)_', 'tokens', 'once');
    if ~isempty(tokB)
        info.blockNum = str2double(tokB{1});
    end

    % Condition: look for _A_ _B_ _C_ _Sh_ or _Sh.
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
