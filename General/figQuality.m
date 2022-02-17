function figQuality(figh, axh, dim, heat)
    %figh = figure handle
    %axh = axis handle
    %dim = [w h] in inches
    
    figh.Units = 'inches';
    figh.Position = [0 0 dim];
    figh.Renderer = 'Painters';
    figh.Color = [1 1 1];
    figh.PaperPositionMode = 'auto'; 
    if nargin>3
        if heat
            box on;
        else
            box off;
        end
    else
        box off;
    end
    %set(axh,'LooseInset',get(axh,'TightInset')) 
    
    allAxesInFigure = findall(figh,'type','axes');
    for i=1:size(allAxesInFigure,1)
        axh.XColor = 'k'; axh.YColor = 'k'; axh.ZColor = 'k';
        axh = allAxesInFigure(i);
        axh.FontSizeMode = 'manual';
        axh.FontSize = 10;
        plblsz = 10; %preferred label size
        axh.LabelFontSizeMultiplier = plblsz / axh.FontSize;
        axh.FontName = 'Arial';
        axh.LineWidth = .75;
        axh.TickDir = 'out';
    end
end
