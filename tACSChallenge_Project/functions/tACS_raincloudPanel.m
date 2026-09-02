function tACS_raincloudPanel(ax, groups, labels, colours, varargin)
% Raincloud panel: half violin, jittered points, mean with SEM.

    q = inputParser;
    q.addParameter('RefStars',  [],    @(x) isnumeric(x));
    q.addParameter('Brackets',  [],    @(x) isstruct(x) || isempty(x));
    q.addParameter('Separator', false, @(x) islogical(x) || isnumeric(x));
    q.addParameter('YLabel',    '',    @(x) ischar(x) || isstring(x));
    q.addParameter('XLabel',    '',    @(x) ischar(x) || isstring(x));
    q.parse(varargin{:});

    refStars  = q.Results.RefStars;
    brackets  = q.Results.Brackets;
    separator = logical(q.Results.Separator);

    S  = styleDefaults();
    nC = numel(groups);
    hold(ax, 'on');

    for c = 1:nC
        drawOne(ax, c, groups{c}, colours(c,:), S);
    end

    allY = vertcat(groups{:});
    allY = allY(~isnan(allY));
    if isempty(allY); return; end
    yMin = min(allY); yMax = max(allY);
    yRange = yMax - yMin;
    if yRange == 0; yRange = 1; end

    topY = yMax;

    if ~isempty(refStars)
        refY = groups{1}; refY = refY(~isnan(refY));
        if ~isempty(refY)
            plot(ax, [0.4 nC+0.6], [mean(refY) mean(refY)], '--', ...
                'Color', S.refCol, 'LineWidth', 1.0);
        end
    end

    if separator
        plot(ax, [1.5 1.5], [yMin - yRange, yMax + yRange], ':', ...
            'Color', S.sepCol, 'LineWidth', 0.75);
    end

    if ~isempty(refStars)
        starY = yMax + 0.08*yRange;
        for c = 2:nC
            pj = refStars(c-1);
            if ~isnan(pj) && pj < 0.05
                text(ax, c, starY, '*', 'HorizontalAlignment', 'center', ...
                    'FontName', S.fontName, 'FontSize', S.letterSz+2, 'FontWeight', 'bold');
            end
        end
        topY = starY;
    end

    if ~isempty(brackets)
        for k = 1:numel(brackets)
            bY  = yMax + yRange * (0.10 + 0.14*k);
            tick = yRange * 0.035;
            x1 = brackets(k).from; x2 = brackets(k).to;
            plot(ax, [x1 x1 x2 x2], [bY-tick, bY, bY, bY-tick], ...
                'k-', 'LineWidth', 1.2);
            text(ax, mean([x1 x2]), bY + tick*0.6, brackets(k).label, ...
                'HorizontalAlignment', 'center', 'FontName', S.fontName, ...
                'FontSize', S.fontSize);
            topY = max(topY, bY + tick*2);
        end
    end

    set(ax, 'FontName', S.fontName, 'FontSize', S.fontSize, 'Box', 'off', ...
        'TickDir', 'out', 'LineWidth', 0.75, 'Layer', 'top');
    xlim(ax, [0.4 nC+0.6]);
    ylim(ax, [yMin - 0.1*yRange, topY + 0.12*yRange]);
    set(ax, 'XTick', 1:nC, 'XTickLabel', labels, 'XTickLabelRotation', 30);
    xlabel(ax, q.Results.XLabel, 'FontName', S.fontName, 'FontSize', S.fontSize);
    ylabel(ax, q.Results.YLabel, 'FontName', S.fontName, 'FontSize', S.fontSize);
end


function S = styleDefaults()
    S.fontName    = 'Times New Roman';
    S.fontSize    = 10;
    S.dotSize     = 9;      S.dotAlpha    = 0.30;
    S.meanSize    = 7;      S.barW        = 1.2;
    S.barCap      = 5;      S.black       = [0 0 0];
    S.violinW     = 0.36;   S.violinAlpha = 0.35;
    S.dotJit      = 0.10;   S.dotOff      = 0.14;
    S.sepCol      = [0.75 0.75 0.75];
    S.refCol      = [0.55 0.55 0.55];
    S.letterSz    = 12;
end


function drawOne(ax, xc, y, faceCol, S)
% Violin drawn only when n > 4.
    y = y(~isnan(y));
    if isempty(y); return; end

    if numel(y) > 4
        [fd, xi] = ksdensity(y);
        fd = fd / max(fd) * S.violinW;
        xv = [xc - fd, fliplr(repmat(xc, 1, numel(xi)))];
        yv = [xi, fliplr(xi)];
        patch(ax, xv, yv, faceCol, 'FaceAlpha', S.violinAlpha, ...
            'EdgeColor', faceCol*0.85, 'LineWidth', 0.5);
    end

    xj = xc + S.dotOff + (rand(numel(y),1)-0.5)*2*S.dotJit;
    scatter(ax, xj, y, S.dotSize, 'filled', 'MarkerFaceColor', faceCol, ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', S.dotAlpha);

    m = mean(y); sem = std(y)/sqrt(numel(y));
    errorbar(ax, xc, m, sem, 's', 'MarkerSize', S.meanSize, ...
        'MarkerFaceColor', S.black, 'MarkerEdgeColor', S.black, ...
        'Color', S.black, 'LineWidth', S.barW, 'CapSize', S.barCap);
end
