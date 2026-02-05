function data_array = tACSChallenge_ImportData(filename, dataLines)
%% script originally written by Florian Kasten, University of Oldenburg

%IMPORTFILE Import data from a text file
%  ERBLOCK120200724233639 = IMPORTFILE(FILENAME) reads data from text
%  file FILENAME for the default selection.  Returns the numeric data.
%
%  ERBLOCK120200724233639 = IMPORTFILE(FILE, DATALINES) reads data for
%  the specified row interval(s) of text file FILENAME. Specify
%  DATALINES as a positive scalar integer or a N-by-2 array of positive
%  scalar integers for dis-contiguous row intervals.

if nargin < 2
    dataLines = [1, Inf];
end

%%
opts = delimitedTextImportOptions("NumVariables", 14);
% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";
% Specify column names and types
opts.VariableNames = ["Var1", "Timeus", "TACSV", "BNCMode", "BNCInV", "BNCOutV",...
    "LeftButton", "RightButton", "LED0Bright", "LED1Bright", "LED2Bright", "LED3Bright", "LED4Bright", "LED5Bright"];
opts.SelectedVariableNames = ["Timeus", "TACSV", "BNCMode", "BNCInV", "BNCOutV",...
    "LeftButton", "RightButton", "LED0Bright", "LED1Bright", "LED2Bright", "LED3Bright", "LED4Bright", "LED5Bright"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double",...
    "double", "double", "double", "double", "double", "double", "double", "double"];
% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var1", "EmptyFieldRule", "auto");

% Import the data
tmp_matrix = readtable(filename, opts);

%% Convert to output type
tmp_matrix = table2array(tmp_matrix);
tmp_matrix = tmp_matrix(2:end,:);

%% Extract Data
% Delta T
data_array.raw_dt = diff(tmp_matrix(:,1));

% handle reset of processor clock in teensy during recording
% check if any dt is negative (step backwards in time)
if sum(data_array.raw_dt<0) > 0
    
    % replace with the mean of all other dt
    data_array.raw_dt(data_array.raw_dt < 0) =  mean(data_array.raw_dt(data_array.raw_dt > 0),'omitnan');
    
end

data_array.dt = mean(data_array.raw_dt,'omitnan')/10^6; %uS to Seconds
% Sampling Rate
data_array.Fs = 1/data_array.dt;
% tACS
data_array.tACS = tmp_matrix(:,2);
% BNC (V)
data_array.BNC_In = tmp_matrix(:,4);
data_array.BNC_Out = tmp_matrix(:,5);
% Buttons
data_array.L_Button = tmp_matrix(:,6);
data_array.R_Button = tmp_matrix(:,7);
% LED
data_array.LEDs = tmp_matrix(:,8:13);


end
