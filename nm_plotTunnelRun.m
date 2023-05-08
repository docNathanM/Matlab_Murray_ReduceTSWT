function [fh,ah,pH] = nm_plotTunnelRun(TRdata)
%% nm_plotTunnelRun() creates a plot of the TSWT run data in TRdata.
% Assumes that TRdata is loaded with nm_loadTR().
% Depends on nm_plotSpaceXY()

    % Create a switch to plot slightly differently for high Mach runs/
    if ( TRdata.Ptotal_psia_mn > 100 )
        plotCase = 'highMach';
    else
        plotCase = 'lowMach';
    end

    % The primary axis is going to be for P0 versus time ...
    [fh,ah(1)] = nm_plotSpaceXY();

    switch plotCase
        case 'highMach'
            fh.Position = [850 680 850 650];
        case 'lowMach'
            fh.Position = [850 680 850 500];
    end

    yPosA = ah(1).Position(2); % save this for later
    ah(1).Position(2) = ah(1).Position(2) + 0.20;
    ah(1).Position(4) = ah(1).Position(4) - 0.20;
    
    % This axis is for temperature and should overlay ah(1).
    ah(2) = axes;
    setMyFavoriteAxesProperties(ah(2));
    
    % Plot P_total and add entry to legend.
    pH.P0 = plot(ah(1),TRdata.time_sec,TRdata.Ptotal_psia_vals,'-k');
    pH.P0.LineWidth = 1.2;
    pH.P0.DisplayName = '$P_0$';
    legEntries = [pH.P0];

    % Plot P_static and add entry to legend.
    pH.Ps = plot(ah(1),TRdata.time_sec,TRdata.Pstatic_psia_vals,'-k');
    pH.Ps.LineWidth = 0.85;
    pH.Ps.Color = 0.5*[1 1 1];
    pH.Ps.DisplayName = '$P_s$';
    pH.Ps.Marker = "^";
    pH.Ps.MarkerSize = 4;
    pH.Ps.MarkerIndices = ...
        pH.Ps.MarkerIndices(TRdata.steadyIndexStartEnd(1):5:TRdata.steadyIndexStartEnd(2));
    
    legEntries = [legEntries pH.Ps];

    % If the 'steadyTime_sec' value is available, I want to plot vertical
    % lines to correspond to the start and end of this.
    if ( isfield(TRdata,'steadyTime_sec') )
        lh1 = line(ah(1),TRdata.steadyTime_sec(1)*[1 1],ah(1).YLim);
        lh1.LineStyle = '--';
        lh1.Color = [0 0 1];
        lh1.LineWidth = 0.5;
        lh1.DisplayName = '``\emph{Steady}"';
        
        legEntries = [legEntries lh1];

        lh2 = line(ah(1),TRdata.steadyTime_sec(2)*[1 1],ah(1).YLim);
        lh2.LineStyle = '--';
        lh2.Color = [0 0 1];
        lh2.LineWidth = 0.5;
    else
    end
    
    % For some testing we set a TTLdelay_sec value to start data
    % acquisition some time delta *after* the P0state changes to 7. This
    % state change corresponds to the tunnel controller turning on. If this
    % value is present in the data, I want to indicate that with a vertical
    % line.
    iX = find(TRdata.P0state == 7,1,'first');
    if ( isfield(TRdata,'TTLdelay_sec') && ~isempty(iX) )    
        TTLtime_sec = TRdata.time_sec(iX) + TRdata.TTLdelay_sec;
        lhTrig = line(ah(1),TTLtime_sec*[1 1],ah(1).YLim);
        lhTrig.LineStyle = '-';
        lhTrig.Color = [1 0 0];
        lhTrig.DisplayName = 'TTL On';
        
        legEntries = [legEntries lhTrig];
    else
    end
    
    % Compute and plot that tank and total temperatures ...
    Ttotal_degC = (TRdata.Ttotal_degF_vals - 32) .* (5/9);
    pH.TT210 = plot(ah(2),TRdata.time_sec,Ttotal_degC,'-.k');
    pH.TT210.DisplayName = '$T_0$';
    % Compute the tank temperature in degC from the TT100 value. The TT100
    % value is in degF in the raw data.    
    TankTemp_degC = (TRdata.WTCdataAll.TT100 - 32) .* (5/9);
    pH.TT100 = plot(ah(2),TRdata.time_sec,TankTemp_degC,':k');
    pH.TT100.LineWidth = 1.5;
    pH.TT100.Color = [0 0.65 0];
    pH.TT100.DisplayName = '$T_{tank}$';
    pH.TT100.Marker = '.';
    pH.TT100.MarkerSize = 6;
    pH.TT100.MarkerIndices = ...
        pH.TT100.MarkerIndices(TRdata.steadyIndexStartEnd(1):5:TRdata.steadyIndexStartEnd(2));
    
    % Add temperature handles to legend entries ...
    legEntries = [legEntries pH.TT210 pH.TT100];
   
    % Attempt to set smart limits for the plots ... for PRESSURE
    ah(1).YLim = [0 round(1.1*max(TRdata.Ptotal_psia_vals),-1)];
    if ( diff(nm_minmax( ah(1).YLim ))/10 < 4 )
        pDiff = 5;
    else
        switch plotCase
            case 'highMach'
                pDiff = 50;
            otherwise
                pDiff = 10;
        end
    end
    ah(1).YLim(1) = floor(ah(1).YLim(1)/pDiff) * pDiff;
    ah(1).YLim(2) = ceil(ah(1).YLim(2)/pDiff) * pDiff;
    ah(1).YTick = ah(1).YLim(1):pDiff:ah(1).YLim(2);
    nPdels = numel(ah(1).YTick)-1;
    % Mke sure the "steady" lines extend correctly.
    if ( isfield(TRdata,'steadyTime_sec') )
        lh1.YData = ah(1).YLim;
        lh2.YData = ah(1).YLim;
    end
    
    % Attempt to set smart limits for the plots ... for TEMPERATURE
    tMax = round(1.05*max(Ttotal_degC),0);
    tMin = round(0.95*min(Ttotal_degC),0);
    tRangeSet = round(diff([tMin tMax])/nPdels,0);
    if ( tRangeSet < 4 )
        nTdels = nPdels/2;
        tDiff = round(diff([tMin tMax])/nTdels,0);
    else
        nTdels = nPdels;
        tDiff = round(diff([tMin tMax])/nTdels,0);
    end
    ah(2).YLim(2) = tMax;
    ah(2).YLim(1) = tMax - nTdels*tDiff;    
    yTticks = ah(2).YLim(2):-tDiff:ah(2).YLim(1);
    ah(2).YTick = sort(yTticks);
    
    % Get TIME limits set and centered ...
    tOffSet = TRdata.time_sec(TRdata.steadyIndexStartEnd(1)) - ...
        TRdata.time_sec(1);
    ah(1).XLim = TRdata.steadyTime_sec + tOffSet*[-1 1];
    offLeft = TRdata.time_sec(1) - ah(1).XLim(1);
    offRight = ah(1).XLim(2) - TRdata.time_sec(end);
    diffOff = (offRight - offLeft)/2;
    ah(1).XLim(1) = ah(1).XLim(1) - diffOff;
    ah(1).XLim(2) = ah(1).XLim(2) - diffOff;

    % Apply TIME limits to ah(2) ...
    ah(2).XLim = ah(1).XLim;
    
    legH = legend(ah(1),legEntries);
    legH.Interpreter = 'latex';
    legH.Location = 'northoutside';
    legH.Orientation = 'horizontal';
    legH.FontSize = 18;
    legH.Title.String = sprintf('Tunnel Run %0.4d',TRdata.TR);

    ah(1).YLabel.String = 'Pressure (psia)';
    ah(2).YLabel.String = 'Temperature (deg.\ C)';
    ah(2).XLim = ah(1).XLim;
    ah(2).XTick = [];
    ah(2).YAxisLocation = 'right';
        
    switch plotCase
        case 'highMach'
    
            ah(3) = axes;
            setMyFavoriteAxesProperties(ah(3))
            ah(3).XLim = ah(1).XLim;
                        
            ah(1).XLabel.String = [];
            ah(3).XLabel.String = 'Time (sec)';
            ah(1).YLabel.String = '$P_0$ (psia)';
            ah(3).YLabel.String = '$P_{s}$ (psia)';
            
            ah(1).Position = [0.13 0.3677 0.7750 0.4923];
            ah(2).Position = [0.13 0.3677 0.7750 0.4923];
            ah(3).Position = [0.13 0.1100 0.7750 0.2000];
            
            % Copy Ps from ah(1) to ah(3) ... but don't delete it because I
            % need to keep the reference in the legend ... so just
            % eliminate the data from the line that exists in ah(1).
            pH.Ps2 = copyobj(pH.Ps,ah(3));
            pH.Ps.XData = [];
            pH.Ps.YData = [];

            if ( isfield(TRdata,'steadyTime_sec') )
                psYLim = ah(3).YLim;
                lh12 = copyobj(lh1,ah(3));
                lh12.YData = psYLim;
                lh22 = copyobj(lh2,ah(3));
                lh22.YData = psYLim;
            end
        
        otherwise
    end

    refresh(fh)
    
    movegui(fh,'north')

end

function setMyFavoriteAxesProperties(ahandle)
    ahandle.XLabel.Interpreter = 'latex';
    ahandle.XLabel.FontSize = 24;
    ahandle.YLabel.Interpreter = 'latex';
    ahandle.YLabel.FontSize = 24;
    ahandle.FontSize = 20;
    ahandle.FontName = 'times';
    ahandle.Box = 'on';
    ahandle.NextPlot = 'add';
    ahandle.Color = 'none';
end
