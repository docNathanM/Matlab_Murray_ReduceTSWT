function [ fh,ah ] = nm_plotSpaceXY()
%nm_plotSpaceXY just sets a plot up the way I like.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set plot defaults:
    set(0,'DefaultLineLineWidth',1.5)
    set(0,'DefaultLineMarkerSize',12)
    set(0,'DefaultAxesFontSize',16)
    set(0,'DefaultFigureColor',[1,1,1])
    set(0,'DefaultTextFontSize',12)
    set(0,'DefaultTextInterpreter','latex')
    set(0,'DefaultTextFontName','Times-Roman')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     close all
    fh = figure;
    fh.Color = [1 1 1];
    fh.Position = [900 300 760 390];
    ah = axes;
    ah.XLabel.String = '$x$';
    ah.XLabel.Interpreter = 'latex';
    ah.XLabel.FontSize = 24;
    ah.YLabel.String = '$y$';
    ah.YLabel.Interpreter = 'latex';
    ah.YLabel.FontSize = 24;
    ah.FontSize = 20;
    ah.FontName = 'times';
    ah.Box = 'on';
    ah.NextPlot = 'add';
    
%     movegui(fh,'northeast');

end

