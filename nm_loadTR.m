function [TRdata] = nm_loadTR(runNumIN,varargin)
%nm_loadTRcsv() => TRdata
% This script reads the tunnel run data and finds the steady portion of
% the run. The input to the function is the numeric run number of interest.

%% PARSE Input
p = inputParser;

defaultDoRebuild = false;
validInputNum = @(x) isnumeric(x);
defaultPathToCSV = 'DATA_WTC-0';

addRequired(p,'runNumIN',validInputNum);
addParameter(p,'reloadCSV',defaultDoRebuild,@(x) islogical(x));
addParameter(p,'pathToCSV',defaultPathToCSV,@(x) isfolder(x));

parse(p,runNumIN,varargin{:})

runNum = p.Results.runNumIN;
doRebuild = p.Results.reloadCSV;
pathToCSV = p.Results.pathToCSV;

%% Debug
clear variables
pathToCSV = 'DATA_WTC_runData';
doRebuild = true;
runNum = 2931;
%%

matFileName = sprintf('Run%d.mat',runNum);
pathToMATfiles = [pathToCSV '_MAT'];

if ( ~doRebuild && ...
        exist(fullfile(pathToMATfiles,matFileName),'file') == 2 )
    
    TRdata = load(fullfile(pathToMATfiles,matFileName));
    
else

csvFileName = sprintf('Run%d.csv',runNum);

%% Open the text file.
fileID = fopen(fullfile(pathToCSV,csvFileName),'r');

%% Determine Range of Data and the line where Configuration Data starts
% The CSV files have a break after the end of the log data with a listing
% that follows including the set points. First, I find the bounds of the
% log data. Then I use a switch-case structure to read the configuration
% information.

frewind(fileID);

tLine = fgetl(fileID);

if ( contains(tLine,'Run Time Limit') || ...
     contains(tLine,'Shutdown (PB010)') || ...
     contains(tLine,'Low Storage Pressure') || ...
     contains(tLine,'Po Valve Full Open') || ...
     contains(tLine,'MS Program Complete') )
    
TRdata.doTR = true;
    
frewind(fileID);

% Scan the file line by line to find the configuration data block ...
found = false;
cfgRowFound = false;
n = 1;
while not(found)
    tLine = fgetl(fileID);
    % NOTE: tLine = -1 if the line contains the eof marker
    if ( tLine == -1 )
        % EOF reached with no config block ... this is the case for some of
        % the earliest shake down runs.
%         endRow = n-1;
        found = true;
        TRdata.TR = runNum;
        TRdata.doTR = false;
        TRdata.RunMsg = 'No Config Block';
    elseif ( ~isempty(tLine) && strcmp(tLine(1:4),'NCPA') )
            found = true;
            cfgRowFound = true;
            endRow = n-2;
            cfgRow = n;   % NCPA W/T Test Configuration    
    else
        n = n+1;
    end
end

else
    
    TRdata.TR = runNum;
    TRdata.doTR = false;
    TRdata.RunMsg = tLine;
    
end

if ( TRdata.doTR )

startRow = 3; % Start row for the run data block is always 3.
delimiter = ',';

if ( cfgRowFound )
% Then, read the config information ...

frewind(fileID);
textscan(fileID,'%[^\n\r]', cfgRow(1)-1, 'ReturnOnError', false);
fgetl(fileID);
tLine = fgetl(fileID);
n=1;a={};
while ischar(tLine)
    a{n} = tLine;
    n=n+1;
    tLine = fgetl(fileID);
end

sa = strsplit(a{1},',');
TRdata.TR = str2double(sa{2});
TRdata.DateTime = datestr(datenum([sa{4} ' ' sa{6}],'mm/dd/yy HH:MM AM'));
TRdata.RunFile = sa{8};
sa = strsplit(a{2},',');
TRdata.RunMsg = sa{2};

for n=1:length(a)
    if ( contains(a{n},'Title') )
        a1 = strsplit(a{n},',');
        TRdata.tunnelConfig.Title = a1{2};
    end
    if ( contains(a{n},'M#') )
        a1 = strsplit(a{n},',');
        TRdata.tunnelConfig.MachSetPoint = a1{2};
    end
    if ( contains(a{n},'Po SP') )
        a1 = strsplit(a{n},',');
        TRdata.tunnelConfig.commandPzero_psi = a1{2};
    end
    switch a{n}
        case 'Choke Configuration'
            nn=n+1;
            a1 = strsplit(a{nn},',');
            TRdata.tunnelConfig.chokeFlap_percent = str2double(a1{2});
        case 'Plenum Exhaust Configuration'
            nn=n+1;
            a1 = strsplit(a{nn},',');
            TRdata.tunnelConfig.plenumExhaustValve_percent = ...
                str2double(a1{2});
        case 'Plenum Ejector Configuration'
            nn=n+1;
            a1 = strsplit(a{nn},',');
            TRdata.tunnelConfig.plenumEjectorSetPoint_psi = ...
                str2double(a1{2});
            nn=n+2;
            a1 = strsplit(a{nn},',');
            TRdata.tunnelConfig.plenumEjectorIntialPosition_percent = ...
                str2double(a1{2});
    end
end
 
else % Get minimal information from the first line in the file
    frewind(fileID);
    tLine = fgetl(fileID);
    
    sa = strsplit(tLine,' ');
    TRdata.TR = runNum;
    TRdata.DateTime = datestr(datenum([sa{2} ' ' sa{3} ' ' sa{4}],'mm/dd/yy HH:MM AM'));
    TRdata.RunFile = sa{1};
    TRdata.RunMsg = [sa{5:end}];
end
    
frewind(fileID);

%% Read the Columns of Run Data
% The run files don't all have the same number of columns, but the overall
% format is the same, and when a particular variable is available, its name
% is always the same. I can use the strsplit() function like I did above to
% find out how many columns are in the data file, and I can set the field
% names this way also.
%
% I already know where the startRow and endRow are located. The header
% values are on startRow-1 ...

for n = 1:startRow-1 % read each line and keep the last one read
    a = fgetl(fileID);
end

% Split into strings and trim whitespace
sa=strtrim(strsplit(a,','));
% Remove empty cells (sometimes the header row ends with a comma :(
empties = find(cellfun(@isempty,sa));
sa(empties) = [];

nCols = length(sa);

frewind(fileID);

% For more information, see the TEXTSCAN documentation.
formatSpec = [repmat('%s',1,nCols) '%[^\n\r]'];

% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, ...
    'Delimiter', delimiter, 'ReturnOnError', false);
% Remove the empties ...
dataArray(empties) = [];

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric strings to numbers.
% Replace non-numeric strings with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=1:nCols
    % Converts strings in the input cell array to numbers. Replaced
    % non-numeric strings with NaN.
    rawData = dataArray{col};
    for row = 1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric
        % prefixes and suffixes. This came from the Matlab import wizard,
        % so I have no idea what this regexpression does ...
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end

%% Replace non-numeric cells with NaN
% Find non-numeric cells
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); 
raw(R) = {NaN}; % Replace non-numeric cells

%% Create RAW Struct ...

rawC = cell(1,nCols);
for n = 1:nCols
    rawC{n} = cell2mat(raw(:,n));
end

dS = cell2struct(rawC,sa,2);
% Now, dS is a struct with the field names of the column labels, and
% each field is a row vector holding the run data for that column.

%% Construct the TRdata Struct as Desired

TRdata.time_sec = dS.Time;
TRdata.Mach_vals = dS.Mach;
TRdata.Ptotal_psia_vals = dS.Po;
TRdata.Pstatic_psia_vals = dS.Ps;
TRdata.Ttotal_degF_vals = dS.TT210;
TRdata.Re_10e6perFt_vals = dS.Re;

if ( isfield(dS,'Qpsi') )
    TRdata.Q_psi_vals = dS.Qpsi;
else
    TRdata.Q_psi_vals = 0.5 * 1.4 .* TRdata.Pstatic_psia_vals .* ...
        TRdata.Mach_vals.^2;
end

if ( isfield(dS,'Postate') )
    TRdata.P0state = dS.Postate;
end

if ( isfield(dS,'MSstate') )
    TRdata.MSstate = dS.MSstate;
end

if ( numel(unique(dS.MS_Pt)) > 1 )
    TRdata.AoA_SRV400 = dS.SRV400;
end

TRdata.WTCdataAll = dS;

%% Determine the Steady Portion of the Run
% I need to require TRdata.state.Postate = 7 if available. This was not
% available in the earlier version of the wind tunnel control system.
[rSt,rEnd,vC] = findSteady(TRdata);

if ( ~isempty(rSt) && numel(rSt) > 1 )
    TRdata.RunMsg = [TRdata.RunMsg ' Unsteady Run'];
    rSt = rSt(1);
    rEnd = rEnd(end);
end

if ( ~isempty(rSt) && ( rSt < rEnd ) ) % Was there a steady portion?
    % TRUE ... complete the TRdata calculations ...
    TRdata.steadyTime_sec = [rSt rEnd];
    TRdata.steadyIndexStartEnd = [...
        find(TRdata.time_sec > rSt,1), ...
        find(TRdata.time_sec < rEnd,1,'last')];
    TRdata.steadyCalc = vC;

    %% Add Mean Values
    % Based on steady time ...
    steadyRange = TRdata.steadyIndexStartEnd(1):TRdata.steadyIndexStartEnd(2);
    TRdata.Mach_mn = mean(TRdata.Mach_vals(steadyRange));
    TRdata.Ptotal_psia_mn = mean(TRdata.Ptotal_psia_vals(steadyRange));
    TRdata.Pstatic_psia_mn = mean(TRdata.Pstatic_psia_vals(steadyRange));
    TRdata.Ttotal_degF_mn = mean(TRdata.Ttotal_degF_vals(steadyRange));
    TRdata.Re_10e6perFt_mn = mean(TRdata.Re_10e6perFt_vals(steadyRange));
    TRdata.Q_psi_mn = mean(TRdata.Q_psi_vals(steadyRange));

    TRdata.TTR_eqn = '1 + ((1.4-1)/2)*Mach^2';
    TRdata.TTR = 1 + ((1.4-1)/2)*TRdata.Mach_mn^2;

    tempKofF = @(x) (x - 32)*(5/9) + 273.15;

    TRdata.Ttotal_K = tempKofF(TRdata.Ttotal_degF_mn);
    TRdata.Tstatic_K = TRdata.Ttotal_K / TRdata.TTR;
    TRdata.ca_mPerS_eqn = 'sqrt(1.4*287*TRdata.Tstatic_K)';
    TRdata.ca_mPerS = sqrt(1.4*287*TRdata.Tstatic_K);
    TRdata.U_mPerS = TRdata.ca_mPerS * TRdata.Mach_mn;

    alt_ft_fromP = @(x_mBar) ( 1 - (x_mBar/1013.25)^0.190284 ) * 145366.45;

    P_mBar = TRdata.Pstatic_psia_mn * 68.9475729;

    alt_ft = alt_ft_fromP(P_mBar);
    alt_km = alt_ft * 3.048e-4;

    TRdata.alt_km_eqn = 'alt_ft_fromP = @(x_mBar) ( 1 - (x_mBar/1013.25)^0.190284 ) * 145366.45';
    TRdata.alt_km = alt_km;
    
else % Change TRdata.doTR = false and stop ...
    
    TRdata.doTR = false;
    TRdata.RunMsg = 'NoSteadyPortionFound';
    
end

end

saveAsName = fullfile(pathToMATfiles,sprintf('Run%d.mat',TRdata.TR));
save(saveAsName,'-struct','TRdata')

end

end

function [runStart,runEnd,varCalc] = findSteady(tD)

% Looking for the longest steady portion of the run.
% First step ... get the wD-point variance in the Ptotal_psia

% 2023-04-27 NEM I was finding the edges of tT/max(tT), but that was not
% robust for large Pzero set points, so I am grabbing the set point from
% the tunnelConfig information and changing the edge finding to be based on
% tT./PsetPt ... see below ...
PsetPt = str2num(tD.tunnelConfig.commandPzero_psi);

wD = 5;
tT = zeros(1,length(tD.time_sec));
for n = ((wD-1)/2)+1:length(tD.time_sec)-((wD-1)/2)
    tT(n) = var(tD.Ptotal_psia_vals(n-((wD-1)/2):n+((wD-1)/2)));
end
% Step 2 find the longest portion
% b = [0 diff(tT/max(tT) < 0.0001)];
b = [0 diff(tT./PsetPt < 0.01)];
a = find(b);
if (b(a(1)) < 0)
    a = a(2:end);
    if (b(a(end)) > 0)
        a = a(1:end-1);
    end
end
c = a(2:2:end)-a(1:2:end); % segment lengths
d = find(c == max(c)); % longest segment
runStart = tD.time_sec(a(2*d-1)) + 0.5;
runEnd = tD.time_sec(a(2*d)) - 0.5;

if ( isfield(tD,'state') )
    runStart = max(runStart,min(tD.time_sec(tD.state.Postate==7)));
    runEnd = min(runEnd,max(tD.time_sec(tD.state.Postate==7)));
end

varCalc.movingVar = tT';
varCalc.diff = b';

end

