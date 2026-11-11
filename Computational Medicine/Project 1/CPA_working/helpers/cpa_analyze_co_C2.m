function results = cpa_analyze_co_C2(patientID,endHour,estimatorID)
%CPA_ANALYZE_CO_C2 Make a Figure-4-style analysis for any patient/estimator.
%
% This function performs the complete Question 3 workflow for a chosen
% patient, time limit, and allowed estimator:
%   1. Read the ABP waveform and detect pulse onsets with WABP.
%   2. Extract beat-to-beat features with ABPFEATURE.
%   3. Reject poor-quality beats with JSQI.
%   4. Estimate uncalibrated relative CO with the requested estimator.
%   5. Pair each intermittent thermodilution CO (TCO) measurement with the
%      good ABP beats in the immediately preceding one-minute window.
%   6. Use only the first paired TCO/relative-CO value for C2 calibration.
%   7. Plot calibrated CO, PP, MAP, and HR over the requested timeframe.
%
% Important distinction:
%   ESTIMATECO_V3 returns beat-to-beat relative CO and quality-filtered
%   features. It does not perform TCO pairing or C2 calibration. Those
%   operations are implemented explicitly below so they are easy to audit.
%
% Outputs:
%   RESULTS is a structure containing the calibration factor, continuous
%   one-minute trends, TCO measurements, and paired feature values.
%   The figure is saved in the workspace results folder as PNG and FIG.

if nargin < 1 || isempty(patientID), patientID = "s00020"; end
if nargin < 2 || isempty(endHour), endHour = 12; end
if nargin < 3 || isempty(estimatorID), estimatorID = 5; end
validateattributes(endHour,{'numeric'},{'scalar','positive','finite'}, ...
    mfilename,'endHour');
validateattributes(estimatorID,{'numeric'}, ...
    {'scalar','integer'},mfilename,'estimatorID');
assert(~ismember(estimatorID,[8 9 11]), ...
    'Estimators 8, 9, and 11 are not permitted for this assignment.');

%% 1. Locate files through portable relative paths
paths = cpa_setup_paths();
patientFiles = cpa_find_patient_files(patientID,paths);
patientID = patientFiles.patientID;
abpFile = patientFiles.abpFile;
trendFile = patientFiles.trendFile;
assert(~isempty(trendFile), ...
    'No trend/TCO file was found for patient %s.',patientID);
outputFolder = paths.resultsFolder;
algorithmName = estimatorName(estimatorID);

Fs = 125;                       % ABP samples per second
analysisEndSeconds = endHour*3600;
filterOrder = 1;                % no internal moving-average filter

%% 2. Read only the requested timeframe of ABP
% TIME is elapsed time in seconds; ABP is arterial pressure in mmHg.
% Reading one 12-hour block avoids loading the later part of the record.
fprintf('Reading patient %s through %.1f hours...\n',patientID,endHour);
[time, ABP] = readABPThroughTime(abpFile, analysisEndSeconds, Fs);

%% 3. Run the supplied ABP signal-processing functions
fprintf('Detecting pulse onsets and extracting ABP features...\n');
t_on = wabp(ABP);
feat = abpfeature(ABP, t_on);
[beatq, fractionGoodAllBeats] = jSQI(feat, t_on, ABP);

% FEAT columns used in this question:
%   column 2 = systolic pressure (Psys), mmHg
%   column 4 = diastolic pressure (Pdias), mmHg
%   column 5 = pulse pressure (PP = Psys-Pdias), mmHg
%   column 6 = mean arterial pressure (MAP), mmHg
%   column 7 = beat period, samples

%% 4. Obtain beat-to-beat uncalibrated relative CO
% ESTIMATECO_V3 removes beats for which BEATQ(:,1) is true. Its FEA output
% therefore remains aligned with CORELATIVE and COTIMEMINUTES.
[coRelative, coTimeMinutes, ~, goodFeatures] = ...
    estimateCO_v3(t_on, feat, beatq, estimatorID, filterOrder);

assert(~isempty(coRelative), 'estimateCO_v3 returned no CO estimates.');

% estimateCO_v3 expresses time relative to sample indices. Apply the small
% start-time correction so it aligns precisely with the elapsed-time column
% in the waveform and trend files.
coTimeMinutes = coTimeMinutes + (time(1) - 1/Fs)/60;

PP = goodFeatures(:,5);                  % pulse pressure, mmHg
MAP = goodFeatures(:,6);                 % mean arterial pressure, mmHg
beatPeriodSamples = goodFeatures(:,7);
HR = 60*Fs ./ beatPeriodSamples;         % heart rate, beats per minute

% HR meaning: if one beat lasts 100 samples at 125 Hz, its duration is
% 100/125 = 0.8 s and HR = 60/0.8 = 75 beats/min.

%% 5. Read episodic thermodilution CO measurements from the trend file
[trendTimeSeconds, trendCO] = readTrendCO(trendFile);
hasTCO = trendCO > 0 & trendTimeSeconds <= analysisEndSeconds;
tcoTimeMinutes = trendTimeSeconds(hasTCO)/60;
tcoMeasured = trendCO(hasTCO);            % measured TCO, L/min

assert(~isempty(tcoMeasured), ...
    'No nonzero TCO measurements were found in the requested timeframe.');

%% 6. Pair every TCO with the immediately preceding minute of good ABP
% Example: a TCO recorded at 12 min is paired with beats in (11,12] min.
% A median is used so one noisy but unflagged beat cannot dominate the pair.
nTCO = numel(tcoMeasured);
pairedRelativeCO = nan(nTCO,1);
pairedPP = nan(nTCO,1);
pairedMAP = nan(nTCO,1);
pairedHR = nan(nTCO,1);

for j = 1:nTCO
    precedingMinute = coTimeMinutes > tcoTimeMinutes(j)-1 & ...
                      coTimeMinutes <= tcoTimeMinutes(j);

    pairedRelativeCO(j) = median( ...
        coRelative(precedingMinute), 'omitnan');
    pairedPP(j) = median(PP(precedingMinute), 'omitnan');
    pairedMAP(j) = median(MAP(precedingMinute), 'omitnan');
    pairedHR(j) = median(HR(precedingMinute), 'omitnan');
end

assert(isfinite(pairedRelativeCO(1)) && pairedRelativeCO(1) > 0, ...
    'The first TCO has no valid preceding-minute relative CO estimate.');

%% 7. C2 calibration: use the first paired measurement only
% Do NOT use coRelative(1): that is the first beat in the record, not the
% beat data immediately preceding the first TCO measurement.
C2 = tcoMeasured(1) / pairedRelativeCO(1);
coCalibratedBeat = C2 .* coRelative;       % beat-to-beat calibrated CO
pairedCalibratedCO = C2 .* pairedRelativeCO;

%% 8. Build continuous one-minute trends for a readable Figure-4-style plot
% The paper associates each TCO with a one-minute ABP segment. We therefore
% summarize every consecutive minute using the same median operation. This
% creates one trend point per minute while retaining the
% beat-to-beat calibrated CO in RESULTS.coCalibratedBeat.
minuteEdges = 0:1:endHour*60;
minuteTimeMinutes = minuteEdges(2:end);  % plot at each window's end time
nMinutes = numel(minuteTimeMinutes);

coTrend = nan(nMinutes,1);
ppTrend = nan(nMinutes,1);
mapTrend = nan(nMinutes,1);
hrTrend = nan(nMinutes,1);

for m = 1:nMinutes
    inMinute = coTimeMinutes > minuteEdges(m) & ...
               coTimeMinutes <= minuteEdges(m+1);
    coTrend(m) = median(coCalibratedBeat(inMinute), 'omitnan');
    ppTrend(m) = median(PP(inMinute), 'omitnan');
    mapTrend(m) = median(MAP(inMinute), 'omitnan');
    hrTrend(m) = median(HR(inMinute), 'omitnan');
end

minuteTimeHours = minuteTimeMinutes/60;
tcoTimeHours = tcoTimeMinutes/60;

%% 9. Plot calibrated CO, PP, MAP, and HR
fig = figure('Color', 'w', 'Position', [80 50 1500 1000]);
layout = tiledlayout(fig, 4, 1, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
title(layout,sprintf(['Patient %s: %s continuous CO with ', ...
    'first-pair C2 calibration'],patientID,algorithmName), ...
    'FontWeight','bold','Color','k');

ax = gobjects(4,1);
ax(1) = nexttile(layout);
plot(ax(1), minuteTimeHours, coTrend, 'Color', [0.35 0.35 0.35], ...
    'LineWidth',1.0,'DisplayName',"Calibrated " + algorithmName + " CO");
hold(ax(1), 'on');
stem(ax(1), tcoTimeHours, tcoMeasured, 'k', 'filled', ...
    'LineWidth', 1.0, 'MarkerSize', 4, 'DisplayName', 'Measured TCO');
ylabel(ax(1), 'CO (L/min)');
legend(ax(1), 'Location', 'best');

ax(2) = nexttile(layout);
plot(ax(2), minuteTimeHours, ppTrend, 'Color', [0.35 0.35 0.35], ...
    'LineWidth', 1.0);
hold(ax(2), 'on');
stem(ax(2), tcoTimeHours, pairedPP, 'k', 'filled', ...
    'LineWidth', 1.0, 'MarkerSize', 4);
ylabel(ax(2), 'PP (mmHg)');

ax(3) = nexttile(layout);
plot(ax(3), minuteTimeHours, mapTrend, 'Color', [0.35 0.35 0.35], ...
    'LineWidth', 1.0);
hold(ax(3), 'on');
stem(ax(3), tcoTimeHours, pairedMAP, 'k', 'filled', ...
    'LineWidth', 1.0, 'MarkerSize', 4);
ylabel(ax(3), 'MAP (mmHg)');

ax(4) = nexttile(layout);
plot(ax(4), minuteTimeHours, hrTrend, 'Color', [0.35 0.35 0.35], ...
    'LineWidth', 1.0);
hold(ax(4), 'on');
stem(ax(4), tcoTimeHours, pairedHR, 'k', 'filled', ...
    'LineWidth', 1.0, 'MarkerSize', 4);
ylabel(ax(4), 'HR (bpm)');
xlabel(ax(4), 'Time (hours)');

for k = 1:numel(ax)
    ax(k).Color = 'w';
    ax(k).XColor = 'k';
    ax(k).YColor = 'k';
    ax(k).GridColor = [0.75 0.75 0.75];
    ax(k).GridAlpha = 0.4;
    ax(k).Toolbar.Visible = 'off';
    grid(ax(k), 'on');
    box(ax(k), 'on');
    xlim(ax(k),[0 endHour]);
end
linkaxes(ax, 'x');

baseName = sprintf('q3_%s_est%02d_C2_first%gh', ...
    patientID,estimatorID,endHour);
pngFile = fullfile(outputFolder,[baseName '.png']);
figFile = fullfile(outputFolder,[baseName '.fig']);

% Remove older outputs before writing. On Windows, EXPORTGRAPHICS can fail
% to overwrite an existing PNG even when the file is not read-only.
if isfile(pngFile)
    delete(pngFile);
end
if isfile(figFile)
    delete(figFile);
end
exportgraphics(fig, pngFile, 'Resolution', 300);
savefig(fig, figFile);

%% 10. Return useful intermediate and final results
results = struct();
results.patientID = patientID;
results.endHour = endHour;
results.estimatorID = estimatorID;
results.algorithmName = algorithmName;
results.C2 = C2;
results.fractionGoodAllBeats = fractionGoodAllBeats;
results.firstTCOminute = tcoTimeMinutes(1);
results.firstTCO = tcoMeasured(1);
results.firstPairedRelativeCO = pairedRelativeCO(1);
results.coTimeMinutes = coTimeMinutes;
results.coRelative = coRelative;
results.coCalibratedBeat = coCalibratedBeat;
results.goodFeatures = goodFeatures;
results.minuteTimeHours = minuteTimeHours;
results.coTrend = coTrend;
results.ppTrend = ppTrend;
results.mapTrend = mapTrend;
results.hrTrend = hrTrend;
results.tcoTimeHours = tcoTimeHours;
results.tcoMeasured = tcoMeasured;
results.pairedCalibratedCO = pairedCalibratedCO;
results.pairedPP = pairedPP;
results.pairedMAP = pairedMAP;
results.pairedHR = pairedHR;
results.pngFile = pngFile;
results.figFile = figFile;

fprintf('\nC2 calibration summary:\n');
fprintf('  First TCO time: %.2f min\n', results.firstTCOminute);
fprintf('  First measured TCO: %.3f L/min\n', results.firstTCO);
fprintf('  Preceding-minute median relative CO: %.6f\n', ...
    results.firstPairedRelativeCO);
fprintf('  C2 factor: %.6f\n', results.C2);
fprintf('  Overall good-beat fraction: %.1f%%\n', ...
    100*results.fractionGoodAllBeats);
fprintf('Saved figure to:\n  %s\n  %s\n', pngFile, figFile);
end

function [time, ABP] = readABPThroughTime(filename, endTimeSeconds, Fs)
%READABPTHROUGHTIME Read waveform rows from the start through a time limit.

fid = fopen(filename, 'rb');
assert(fid >= 0, 'Could not open ABP file: %s', filename);
cleanupObject = onCleanup(@() fclose(fid)); %#ok<NASGU>

firstLine = fgetl(fid);
firstValues = sscanf(firstLine, '%f%f');
assert(numel(firstValues) == 2, 'Could not parse the first ABP row.');
firstTime = firstValues(1);
frewind(fid);

numberOfRows = floor((endTimeSeconds-firstTime)*Fs) + 2;
dataCell = textscan(fid, '%f%f', numberOfRows, ...
    'CollectOutput', true, 'ReturnOnError', false);
data = dataCell{1};
data = data(data(:,1) <= endTimeSeconds + 0.5/Fs, :);

time = data(:,1);
ABP = data(:,2);
assert(time(end) >= endTimeSeconds-1/Fs, ...
    'The ABP file does not contain the full requested timeframe.');
assert(all(abs(diff(time)-1/Fs) < 1e-6), ...
    'ABP timestamps are not uniformly sampled at 125 Hz.');
end

function [elapsedSeconds, CO] = readTrendCO(filename)
%READTRENDCO Read elapsed time and thermodilution CO from the n.txt file.
%
% The trend file has 16 tab-separated columns and two header rows. Missing
% values are represented by "-"; STR2DOUBLE converts them to NaN.

fid = fopen(filename, 'rt');
assert(fid >= 0, 'Could not open trend file: %s', filename);
cleanupObject = onCleanup(@() fclose(fid)); %#ok<NASGU>

fgetl(fid);  % column-name header
fgetl(fid);  % unit header
formatSpec = repmat('%s', 1, 16);
raw = textscan(fid, formatSpec, 'Delimiter', '\t', ...
    'MultipleDelimsAsOne', false, 'ReturnOnError', false);

elapsedSeconds = str2double(raw{1});
CO = str2double(raw{16});
end

function name = estimatorName(estimatorID)
switch estimatorID
    case 1, name = "MAP baseline";
    case 2, name = "Pulse pressure (Windkessel)";
    case 3, name = "Systolic area";
    case 4, name = "Systolic area + Kouchoukos";
    case 5, name = "Liljestrand";
    case 6, name = "Herd";
    case 7, name = "Corrected impedance";
    case 10, name = "Diastolic decay";
    case 12, name = "Coalees";
    case 13, name = "Trivial ensemble";
    otherwise, name = "Estimator " + estimatorID;
end
end
