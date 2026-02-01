
                   
function [paf,psd,f] = IAF_estimate(data_path,labnum,subnum,sampling_rate)

%% Input: data_path: the path where we store the Data, e.g., 'Z\Data'
%         labnum: ID of the lab, must be numeric
%         subnum: ID of the subject, must be numeric
%         sampling_rate: EEG sampling rate in Hz, must be numeric

%% Output: paf: peak alpha frequency
%           where paf(1) is the PAF of pre eyes-closed EEG recording
%                 paf(2) is the PAF of pre eyes-open EEG recording
%                 paf(3) is the PAF of post eyes-closed EEG recording
%                 paf(4) is the PAF of post eyes-open EEG recording


%          psd: Power spectral density 
%            psd(1,:) is the PSD of pre eyes-closed EEG recording
%            psd(2,:) is the PSD of pre eyes-open EEG recording
%            ...
%            psd(4,:) is the PSD of of post eyes-open EEG recording             

%           f: frequencies corresponding to PSD 

% This function includes loading the EEG data part, preprocessing part and IAF estimation part. I primarily used 
% the function named 'restingIAF' from Corcoran (2018) in the IAF estimation part
% the script has been written by Yifan Shen (University of Glasogw, 3032397S@student.gla.ac.uk;
% and modified by Simon Hanslmayr (University of Glasgow, simon.hanslmayr@glasgow.ac.uk

%% Add necessary packages
addpath(genpath('restingIAF-master'));
addpath('eeglab2025.0.0');


%% Loading the EEG data

if labnum <10
labnum = strcat('0',num2str(labnum));
else
labnum = num2str(labnum);
end

if subnum<10
subnum = strcat('0',num2str(subnum));
else
subnum = num2str(subnum);
end

prefix = strcat('L',labnum,'_P',subnum);
prefix_S = strcat('L',labnum,'_S',subnum);

full_data_path = fullfile(data_path,prefix);

Pre_EC_file = strcat(prefix_S,'_Pre_EC');
Pre_EO_file = strcat(prefix_S,'_Pre_EO');
Post_EC_file =strcat(prefix_S,'_Post_EC');
Post_EO_file =strcat(prefix_S,'_Post_EO');

File_names = {Pre_EC_file,Pre_EO_file,Post_EC_file,Post_EO_file} ;

%% Preprocessing
eeglab nogui;

% Iteration for 4 conditions (pre_EC, pre_EO, post_EC, post_EO)
for i=1:4

% read the EEG data
EEG = pop_loadbv(full_data_path, strcat(File_names{i},'.vhdr'));

std_EEG = std(EEG.data);


% Filtering
EEG = pop_eegfiltnew(EEG, 1, 40);
EEG = pop_eegfiltnew(EEG, 49, 51, [], 1);  

% divide into 1s Epochs
epoch_length_sec = 1;
epoch_samples = EEG.srate * epoch_length_sec;

n_epochs = floor(EEG.pnts / epoch_samples);

EEG.event = []; 
for ii = 1:n_epochs
    EEG.event(ii).type = 'epoch_marker';
    EEG.event(ii).latency = (ii - 1) * epoch_samples + 1;
    EEG.event(ii).duration = 0;
end
EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = pop_epoch(EEG, {'epoch_marker'}, [0 1]);  
EEG = eeg_checkset(EEG);

% Select bad epochs

z_threshold = 3; 

bad_epochs = false(1, EEG.trials); 

for iiii = 1:n_epochs
    data = EEG.data(:, :, iiii);  
    

    data_vector = data(:);
    z_data = abs((data_vector) / std_EEG);  %
    if any(z_data > z_threshold) || any(z_data < -z_threshold)
        bad_epochs(iiii) = true; 
    end
end
% There was also a manual rejection step, but it may not be necessary since most signals appear to be clean.
% pop_eegplot(EEG, 1, 1, 0); 
% uiwait;
% 
%  bad_epoches_mannual = inputdlg({'bad_epoches'});
bad_epoches_mannual = {''};
%%
bad_epochs_str = strsplit(bad_epoches_mannual{1}, ',');  
bad_epochs_num = str2double(bad_epochs_str);

if ~isnan(bad_epochs_num)
bad_epochs(bad_epochs_num) = true; 
end


% record the bad_epochs
EEG = pop_rejepoch(EEG, bad_epochs, 0);

bad_epoch_indices{i} = find(bad_epochs);

% Save the bad epochs in preprocessing folder
folder_name = 'preprocessing_meta_data';
if ~exist(folder_name, 'dir')
    mkdir(folder_name);
end
save(fullfile(folder_name, strcat(File_names{i}, '_bad_epochs')), 'bad_epoch_indices');

% Convert epoched EEG to continuous

EEG = eeg_epoch2continuous(EEG);

% save the preprocessed data
output_folder = 'EEG_preprocessed';

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
save(fullfile(output_folder, File_names{i}), 'EEG');

%% IAF estimation
%[paf_sums,~,f] = restingIAF(double(EEG.data), 1, 1, [1, 40], sampling_rate,  [7 13], 11, 5);
[paf_sums,~,f] = restingIAF_sh(double(EEG.data), 1, 1, [1, 40], sampling_rate,  [7 13], 11, 5);

paf(i) = paf_sums.paf;
%if ~any(isnan(paf_sums.muSpec))
psd(i,:) = paf_sums.ps;
%else
%psd(i,:) = nan(1,161);    
%end

end



