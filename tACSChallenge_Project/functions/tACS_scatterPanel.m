function tACS_scatterPanel(ax, x, y, faceCol, panelTitle, xLabelText, yLabelText, varargin)
% Scatter panel with least-squares fit, zero reference and faint grid.

    q = inputParser;
    q.addParameter('XLim',          [],    @isnumeric);
    q.addParameter('YLim',          [],    @isnumeric);
    q.addParameter('ShowYLabel',    true,  @(v) islogical(v) || isnumeric(v));
    q.addParameter('SlopeFormat',   'Slope = %.4f per Hz', @(v) ischar(v) || isstring(v) || isempty(v));
    q.addParameter('PanelLetter',   '',    @(v) ischar(v) || isstring(v));
    q.addParameter('LetterXOffset', -0.10, @isnumeric);
    q.parse(varargin{:});

    fontName = 'Times New Roman';
    fontSize = 10;

    valid = ~isnan(x) & ~isnan(y);
    x = x(valid); y = y(valid);

    hold(ax, 'on');

    xl = q.Results.XLim;
    if isempty(xl)
        if isempty(x), xl = [-1 1]; else, xl = [min(x) max(x)]; end
    end

    plot(ax, xl, [0 0], '-', 'Color', [0.80 0.80 0.80], 'LineWidth', 0.4, ...
        'HandleVisibility', 'off');

    plot(ax, x, y, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'none', ...
        'MarkerEdgeColor', faceCol, 'LineWidth', 0.8);

    slope = NaN;
    if numel(x) > 1
        xc = x - mean(x);
        yc = y - mean(y);
        slope = sum(xc .* yc) / sum(xc.^2);
        xFit  = linspace(xl(1), xl(2), 100);
        yFit  = mean(y) + slope * (xFit - mean(x));
        plot(ax, xFit, yFit, '-', 'Color', [0 0 0], 'LineWidth', 1.3);
    end

    hold(ax, 'off');

    set(ax, 'FontName', fontName, 'FontSize', fontSize, 'Box', 'off', ...
        'TickDir', 'out', 'LineWidth', 0.75, 'Layer', 'bottom');
    grid(ax, 'on');
    set(ax, 'GridAlpha', 0.15, 'GridColor', [0.4 0.4 0.4], ...
        'GridLineStyle', '-', 'MinorGridLineStyle', 'none');

    xlim(ax, xl);
    if ~isempty(q.Results.YLim); ylim(ax, q.Results.YLim); end

    title(ax, panelTitle, 'FontName', fontName, 'FontSize', fontSize, ...
        'FontWeight', 'normal');
    xlabel(ax, xLabelText, 'FontName', fontName, 'FontSize', fontSize);

    if logical(q.Results.ShowYLabel)
        ylabel(ax, yLabelText, 'FontName', fontName, 'FontSize', fontSize);
    else
        ylabel(ax, '');
        ax.YTickLabel = [];
    end

    fmt = q.Results.SlopeFormat;
    if ~isempty(fmt) && ~isnan(slope)
        text(ax, 0.035, 0.965, sprintf(char(fmt), slope), 'Units', 'normalized', ...
            'FontName', fontName, 'FontSize', 9, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    end

    if ~isempty(q.Results.PanelLetter)
        text(ax, q.Results.LetterXOffset, 1.06, char(q.Results.PanelLetter), ...
            'Units', 'normalized', 'FontName', fontName, 'FontSize', 12, ...
            'FontWeight', 'bold', 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'bottom', 'Clipping', 'off');
    end
end
