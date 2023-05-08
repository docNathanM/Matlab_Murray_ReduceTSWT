function [fh,ah,ph] = nm_plotTunnelRun(TRdata)

    [fh,ah] = nm_plotSpaceXY();

    ph = plot(TRdata.time_sec,TRdata.Ptotal_psia_vals,'-k');

    if ( isfield(TRdata,'steadyTime_sec') )
        lh1 = line(TRdata.steadyTime_sec(1)*[1 1],ah.YLim);
        lh1.LineStyle = '--';
        lh1.Color = [0 0 1];

        lh2 = line(TRdata.steadyTime_sec(2)*[1 1],ah.YLim);
        lh2.LineStyle = '--';
        lh2.Color = [0 0 1];
    else
    end

end

