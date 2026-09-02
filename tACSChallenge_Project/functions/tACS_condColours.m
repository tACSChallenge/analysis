function [colours, labels, xLabels] = tACS_condColours(includeSham)
% Condition colours and labels, shared by all figures.

    if nargin < 1 || isempty(includeSham); includeSham = false; end

    colours = [0.00 0.45 0.74;    % occipital
               0.85 0.33 0.10;    % retinal
               0.47 0.67 0.19;    % cutaneous
               0.50 0.50 0.50];   % sham

    labels  = {'Occipital', 'Retinal', 'Cutaneous', 'Sham'};

    xLabels = {'{\itd}_{occipital}', '{\itd}_{retinal}', ...
               '{\itd}_{cutaneous}', '{\itd}_{sham}'};

    if ~includeSham
        colours = colours(1:3, :);
        labels  = labels(1:3);
        xLabels = xLabels(1:3);
    end
end
