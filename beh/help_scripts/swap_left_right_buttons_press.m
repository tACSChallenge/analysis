function swap_left_right_buttons_press(data_path, lab_num, sub_num)
%SWAP_LEFT_RIGHT_BUTTONS_PRESS  Swap LeftButton/RightButton values in beh TSVs.
%
%   swap_left_right_buttons_press(DATA_PATH, LAB_NUM, SUB_NUM) swaps, row by
%   row, the values stored in the columns LeftButton and RightButton of every
%   behavioural TSV file (B1 ... B13) of the requested participants.
%
%   Nothing else is modified: the header, the column order, every other
%   column, the number formatting (0.0 stays 0.0), the line endings and the
%   file encoding are written back byte for byte. Files may have different
%   numbers of rows; each file is handled on its own.
%
%   DATA_PATH  folder of one lab, e.g. 'Z:\scripts\tACSChallenge_Project\data\L19'
%   LAB_NUM    lab number, e.g. 19   -> folders sub-L19_S01, sub-L19_S02, ...
%   SUB_NUM    vector of participants, e.g. 1:20
%
%   Example:
%       swap_left_right_buttons_press('Z:\scripts\tACSChallenge_Project\data\L19', 19, 1:20)
%
%   Printed to the Command Window:
%     - per file, how many values each column handed over to the other one.
%       The swap is unconditional - zeros are swapped exactly like any other
%       value - so this count always equals the number of data rows of the
%       file. A row-by-row check of the file as it was written to disk is
%       reported as OK / MISMATCH;
%     - a grand total over all processed files;
%     - a row-level sanity check on the LAST file that was actually modified
%       in the run, showing the LeftButton / RightButton values of individual
%       rows before and after the change.
%
%   Safety: the original file is copied to
%       <...>\beh\backup_before_button_swap\
%   before it is modified, and a file that already has a backup is SKIPPED,
%   so running this function twice cannot silently swap the data back. To
%   redo a file, delete its backup copy. Set MAKE_BACKUP = false to disable.
% ------------------------------------------------------------------------
% Script by Wahriman Andrade de Araújo.
% ------------------------------ settings --------------------------------
N_BLOCKS       = 13;                           % B1 ... B13
MAKE_BACKUP    = true;                         % keep a copy of the original
BACKUP_DIRNAME = 'backup_before_button_swap';
SANITY_NROWS   = 10;                           % rows printed in the sanity check
% ------------------------------------------------------------------------

TAB = char(9);
LF  = char(10);

nDone = 0; nSkip = 0; nBad = 0;
totRows = 0; totMovedL = 0; totMovedR = 0;
sanity = [];

fprintf('\n=== swap_left_right_buttons_press ===\n');
fprintf('data path : %s\n', data_path);
fprintf('lab       : L%02d\n', lab_num);
fprintf('subjects  : %s\n', mat2str(sub_num(:)'));
fprintf('"values moved" = rows whose value was handed over to the other column\n');
fprintf('every row is swapped, whatever the value, so it must equal "rows"\n');
fprintf('%s\n', repmat('-', 1, 78));

for s = sub_num(:)'

    subName = sprintf('sub-L%02d_S%02d', lab_num, s);
    behDir  = fullfile(data_path, subName, 'beh');

    fprintf('\n%s\n', subName);
    if exist(behDir, 'dir') ~= 7
        fprintf('   !! beh folder not found: %s  - skipped\n', behDir);
        nSkip = nSkip + N_BLOCKS;
        continue
    end

    for b = 1:N_BLOCKS

        % ---- locate the file of this block ------------------------------
        d = dir(fullfile(behDir, sprintf('*_B%d_*.tsv', b)));
        if isempty(d)
            d = dir(fullfile(behDir, sprintf('*_B%d.tsv', b)));
        end
        if ~isempty(d)
            d = d(~[d.isdir]);
        end
        if isempty(d)
            fprintf('   B%-2d : !! file not found - skipped\n', b);
            nSkip = nSkip + 1;  continue
        elseif numel(d) > 1
            fprintf('   B%-2d : !! %d files match - skipped (%s)\n', ...
                    b, numel(d), strjoin({d.name}, ', '));
            nSkip = nSkip + 1;  continue
        end
        fname = d(1).name;
        fpath = fullfile(behDir, fname);

        % ---- has this file already been swapped? -------------------------
        bakDir  = fullfile(behDir, BACKUP_DIRNAME);
        bakFile = fullfile(bakDir, fname);
        if MAKE_BACKUP && exist(bakFile, 'file') == 2
            fprintf('   B%-2d : !! backup already exists (file already swapped) - skipped\n', b);
            nSkip = nSkip + 1;  continue
        end

        % ---- read the file as raw bytes ---------------------------------
        raw = local_read(fpath);
        nl  = find(raw == LF, 1, 'first');
        if isempty(nl)
            fprintf('   B%-2d : !! no data rows - skipped\n', b);
            nSkip = nSkip + 1;  continue
        end
        head = raw(1:nl);          % header line, incl. its line break
        body = raw(nl+1:end);      % data rows, untouched bytes

        % ---- find the two columns in the header -------------------------
        hline = regexprep(raw(1:nl-1), '\r$', '');
        flds  = strsplit(hline, TAB, 'CollapseDelimiters', false);
        iL = find(strcmpi(flds, 'LeftButton'));
        iR = find(strcmpi(flds, 'RightButton'));
        if numel(iL) ~= 1 || numel(iR) ~= 1 || iR ~= iL + 1
            fprintf('   B%-2d : !! LeftButton/RightButton not found as adjacent columns - skipped\n', b);
            nSkip = nSkip + 1;  continue
        end

        % ---- swap column iL with column iL+1, row by row ----------------
        pat  = ['^(' sprintf('(?:[^\\t\\r\\n]*\\t){%d}', iL-1) ')([^\t\r\n]*)\t([^\t\r\n]*)'];
        rep  = ['$1$3' TAB '$2'];
        body2 = regexprep(body, pat, rep, 'lineanchors');

        % ---- values before / after --------------------------------------
        [vB, nRows] = local_cols(body,  iL);
        vA          = local_cols(body2, iL);

        if size(vB, 1) ~= nRows || size(vA, 1) ~= nRows
            fprintf('   B%-2d : !! could not read the button columns of every row - NOT written\n', b);
            nSkip = nSkip + 1;  continue
        end

        % ---- back up the original, then write it back --------------------
        if MAKE_BACKUP
            if exist(bakDir, 'dir') ~= 7, mkdir(bakDir); end
            [ok, msg] = copyfile(fpath, bakFile);
            if ~ok
                fprintf('   B%-2d : !! backup failed (%s) - NOT written\n', b, msg);
                nSkip = nSkip + 1;  continue
            end
        end
        local_write(fpath, [head body2]);
        nDone = nDone + 1;

        % ---- report -------------------------------------------------------
        % every row is swapped, whatever the value (0 included), so the number
        % of values each column handed over is always the number of data rows.
        movedL = sum(vA(:,2) == vB(:,1) | (isnan(vA(:,2)) & isnan(vB(:,1))));
        movedR = sum(vA(:,1) == vB(:,2) | (isnan(vA(:,1)) & isnan(vB(:,2))));
        totRows   = totRows   + nRows;
        totMovedL = totMovedL + movedL;
        totMovedR = totMovedR + movedR;

        if isequaln(vA, vB(:, [2 1]))
            flag = 'OK';
        else
            flag = 'MISMATCH !!';  nBad = nBad + 1;
        end

        fprintf('   B%-2d : Left -> %d values moved | Right -> %d values moved | %d rows | %s\n', ...
                b, movedL, movedR, nRows, flag);

        % ---- keep this file for the sanity check (last one wins) ----------
        sanity.sub   = subName;
        sanity.block = b;
        sanity.file  = fpath;
        sanity.vB    = vB;
        sanity.vA    = vA;
    end
end

% ------------------------------ totals ----------------------------------
fprintf('\n%s\n', repmat('-', 1, 78));
fprintf('files swapped: %d   files skipped: %d   mismatches: %d\n', nDone, nSkip, nBad);
if totRows > 0
    fprintf('all files    : Left -> %d values moved | Right -> %d values moved | %d rows\n', ...
            totMovedL, totMovedR, totRows);
end

% --------------------------- sanity check -------------------------------
fprintf('\n%s\n', repmat('=', 1, 78));
if isempty(sanity)
    fprintf('SANITY CHECK\n');
    fprintf('no file was modified in this run - nothing to show\n');
else
    fprintf('SANITY CHECK - last file modified: %s, block B%d\n', sanity.sub, sanity.block);
    fprintf('%s\n', sanity.file);
    vB = sanity.vB;  vA = sanity.vA;
    idx = find(vB(:,1) ~= 0 | vB(:,2) ~= 0);
    if isempty(idx)
        fprintf('no button press in this file - showing the first rows instead\n');
        idx = (1:size(vB,1))';
    else
        fprintf('first rows with a button press (%d press samples in total)\n', numel(idx));
    end
    idx = idx(1:min(SANITY_NROWS, numel(idx)));
    fprintf('%10s | %12s %12s | %12s %12s\n', ...
            'row', 'LeftButton', 'RightButton', 'LeftButton', 'RightButton');
    fprintf('%10s | %12s %12s | %12s %12s\n', '', 'BEFORE', 'BEFORE', 'AFTER', 'AFTER');
    fprintf('%s\n', repmat('-', 1, 68));
    for k = idx(:)'
        fprintf('%10d | %12g %12g | %12g %12g\n', k, vB(k,1), vB(k,2), vA(k,1), vA(k,2));
    end
    fprintf('%s\n', repmat('-', 1, 68));
    nS = size(vB, 1);
    fprintf('whole file : Left -> %d values moved | Right -> %d values moved | %d rows\n', ...
            sum(vA(:,2) == vB(:,1) | (isnan(vA(:,2)) & isnan(vB(:,1)))), ...
            sum(vA(:,1) == vB(:,2) | (isnan(vA(:,1)) & isnan(vB(:,2)))), nS);
    if isequaln(vA, vB(:, [2 1]))
        fprintf('every row: LeftButton(after) == RightButton(before) and vice versa -> OK\n');
    else
        fprintf('!! the columns were NOT exchanged correctly - check this file\n');
    end
end
fprintf('%s\n\n', repmat('=', 1, 78));

end % swap_left_right_buttons_press


% ======================= local helper functions =========================

function txt = local_read(f)
% Read the whole file as raw bytes (no encoding / line-ending translation).
fid = fopen(f, 'r');
if fid < 0, error('Cannot open %s', f); end
d = fread(fid, Inf, '*uint8');
fclose(fid);
txt = char(d(:)');
end


function local_write(f, txt)
% Write the raw bytes back.
fid = fopen(f, 'w');
if fid < 0, error('Cannot write %s', f); end
fwrite(fid, uint8(txt), 'uint8');
fclose(fid);
end


function [v, nRows] = local_cols(body, iL)
% Return the values of columns iL and iL+1 of every data row (nRows x 2).
LF = char(10);
nRows = numel(strfind(body, LF));
if ~isempty(body) && body(end) ~= LF, nRows = nRows + 1; end
pat = ['^' sprintf('(?:[^\\t\\r\\n]*\\t){%d}', iL-1) '([^\t\r\n]*)\t([^\t\r\n]*)[^\r\n]*'];
two = regexprep(body, pat, '$1 $2', 'lineanchors');
v   = sscanf(two, '%f %f', [2 Inf])';
end
