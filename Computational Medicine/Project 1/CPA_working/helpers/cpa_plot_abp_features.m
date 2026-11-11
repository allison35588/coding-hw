function results = cpa_plot_abp_features(patientID, targetHours)
%CPA_PLOT_ABP_FEATURES Plot ABP features for any downloaded patient.
%
% RESULTS = CPA_PLOT_ABP_FEATURES(PATIENTID) analyzes [10 11] hours.
% RESULTS = CPA_PLOT_ABP_FEATURES(PATIENTID, TARGETHOURS) analyzes any
% one or more requested times. PATIENTID can be 138, "138", or "s00138".
%
% This function reproduces traces similar to Figure 1 of Sun et al. (2009)
% for the first 20 complete ABP pulses beginning at each requested time.
% It uses the three supplied PhysioToolkit functions:
%   1. wabp       - detect the onset of each ABP pulse
%   2. abpfeature - calculate beat-level ABP features
%   3. jSQI       - identify pulses with abnormal signal quality
%
% Figure symbols follow the paper:
%   *  pulse onset
%   x  end systole estimated as 0.3*sqrt(beat period)
%   o  end systole estimated by the lowest-nonnegative-slope method
%
% The ABP text file is fixed-width and sampled at 125 Hz. Instead of
% loading the entire file, READABPWINDOW seeks directly to short windows
% surrounding the requested hours. A 30-second preroll lets WABP initialize its
% detector, and 60 seconds after the requested hour provides enough beats
% for JSQI (which requires at least 30 detected onsets).
%
% Output filenames include the normalized patient ID and requested times.
%
% Reference:
% Sun JX et al. Cardiac output estimation using arterial blood pressure
% waveforms. Crit Care Med. 2009;37(1):72-80.

%% File locations and analysis settings
if nargin < 1 || isempty(patientID)
    patientID = "s00020";
end

if nargin < 2 || isempty(targetHours)
    targetHours = [10, 11];
end
validateattributes(targetHours, {'numeric'}, ...
    {'vector', 'real', 'finite', 'nonnegative'}, mfilename, 'targetHours');
targetHours = targetHours(:).';  % Use a row vector consistently.

paths = cpa_setup_paths();
patientFiles = cpa_find_patient_files(patientID,paths);
patientID = patientFiles.patientID;
abpFile = patientFiles.abpFile;
outputFolder = paths.resultsFolder;

Fs = 125;                   % ABP sampling frequency (Hz)
nPulses = 20;
prerollSeconds = 30;
postrollSeconds = 60;

%% Detect ABP features in the two requested windows
results = repmat(struct(), numel(targetHours), 1);

for k = 1:numel(targetHours)
    targetSeconds = targetHours(k) * 3600;
    windowStart = targetSeconds - prerollSeconds;
    windowEnd = targetSeconds + postrollSeconds;

    [timeSeconds, abp] = readABPWindow( ...
        abpFile, windowStart, windowEnd, Fs);

    % Run the supplied PhysioToolkit ABP processing sequence.
    onset = wabp(abp);
    features = abpfeature(abp, onset);
    [beatQ, fractionGood] = jSQI(features, onset, abp);

    assert(~isempty(beatQ), ...
        'jSQI needs at least 30 detected onsets at hour %g.', targetHours(k));

    onsetTime = timeSeconds(onset);
    firstBeat = find(onsetTime >= targetSeconds, 1, 'first');
    assert(~isempty(firstBeat), ...
        'No ABP pulse was detected after hour %g.', targetHours(k));

    beatIndex = firstBeat:(firstBeat + nPulses - 1);
    assert(beatIndex(end) + 1 <= numel(onset), ...
        'The extracted window does not contain 20 complete pulses.');
    assert(beatIndex(end) <= size(features, 1), ...
        'Feature extraction did not return all 20 requested pulses.');

    % Include the next onset as the right boundary so pulse 20 is complete.
    plotSamples = onset(beatIndex(1)):onset(beatIndex(end) + 1);

    results(k).targetHour = targetHours(k);
    results(k).targetSeconds = targetSeconds;
    results(k).timeSeconds = timeSeconds;
    results(k).abp = abp;
    results(k).onset = onset;
    results(k).features = features;
    results(k).beatQ = beatQ;
    results(k).fractionGoodInAnalysisWindow = fractionGood;
    results(k).beatIndex = beatIndex;
    results(k).selectedFractionGood = ...
        mean(~beatQ(beatIndex, 1));
    results(k).plotSamples = plotSamples;
end

%% Plot the two sets of 20 pulses
figureHeight = max(500, 420 * numel(targetHours));
fig = figure('Color', 'w', 'Position', [100 100 1500 figureHeight]);
layout = tiledlayout(fig, numel(targetHours), 1, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
title(layout, sprintf('Patient %s: ABP waveform features', patientID), ...
    'FontWeight', 'bold', 'Color', 'k');

axesHandles = gobjects(numel(targetHours), 1);
for k = 1:numel(targetHours)
    ax = nexttile(layout);
    axesHandles(k) = ax;
    ax.Color = 'w';
    ax.XColor = 'k';
    ax.YColor = 'k';
    ax.GridColor = [0.75 0.75 0.75];
    ax.GridAlpha = 0.45;
    ax.Toolbar.Visible = 'off';
    hold(ax, 'on');

    targetSeconds = results(k).targetSeconds;
    timeSeconds = results(k).timeSeconds;
    abp = results(k).abp;
    onset = results(k).onset;
    features = results(k).features;
    beatQ = results(k).beatQ;
    beatIndex = results(k).beatIndex;
    plotSamples = results(k).plotSamples;

    xWave = timeSeconds(plotSamples) - targetSeconds;
    plot(ax, xWave, abp(plotSamples), 'k-', 'LineWidth', 1.1, ...
        'DisplayName', 'ABP waveform');

    onsetSamples = onset(beatIndex);
    plot(ax, timeSeconds(onsetSamples) - targetSeconds, ...
        abp(onsetSamples), 'k*', 'MarkerSize', 7, 'LineWidth', 1.0, ...
        'DisplayName', 'Pulse onset');

    % ABPFEATURE column 9: end systole by 0.3*sqrt(beat period).
    endSystoleHeuristic = round(features(beatIndex, 9));
    plot(ax, timeSeconds(endSystoleHeuristic) - targetSeconds, ...
        abp(endSystoleHeuristic), 'kx', 'MarkerSize', 7, ...
        'LineWidth', 1.4, ...
        'DisplayName', 'End systole: 0.3 sqrt(beat period)');

    % ABPFEATURE column 11: end systole by lowest nonnegative slope.
    endSystoleSlope = round(features(beatIndex, 11));
    plot(ax, timeSeconds(endSystoleSlope) - targetSeconds, ...
        abp(endSystoleSlope), 'ko', 'MarkerSize', 6, ...
        'LineWidth', 1.2, 'MarkerFaceColor', 'none', ...
        'DisplayName', 'End systole: lowest nonnegative slope');

    % JSQI marks abnormal pulses. Red diamonds are an added annotation and
    % are not one of the three symbols in the original Figure 1.
    badBeatIndex = beatIndex(beatQ(beatIndex, 1));
    if ~isempty(badBeatIndex)
        badOnsets = onset(badBeatIndex);
        plot(ax, timeSeconds(badOnsets) - targetSeconds, ...
            abp(badOnsets), 'rd', 'MarkerSize', 8, 'LineWidth', 1.2, ...
            'DisplayName', 'jSQI-rejected pulse');
    end

    nGood = sum(~beatQ(beatIndex, 1));
    title(ax, sprintf(['First %d pulses beginning at %g h ', ...
        '(%d/%d pass jSQI)'], nPulses, targetHours(k), nGood, nPulses), ...
        'Color', 'k');
    xlabel(ax, sprintf('Time after %g h (s)', targetHours(k)), ...
        'Color', 'k');
    ylabel(ax, 'Arterial blood pressure (mmHg)', 'Color', 'k');
    xlim(ax, [xWave(1), xWave(end)]);
    grid(ax, 'on');
    box(ax, 'on');
    ax.FontSize = 11;
end

% Add one legend key for the red jSQI marker even when rejections occur
% only in the second panel (as they do for patient s00020).
hasRejectedPulse = any(arrayfun(@(s) ...
    any(s.beatQ(s.beatIndex, 1)), results));
if hasRejectedPulse
    plot(axesHandles(1), NaN, NaN, 'rd', 'MarkerSize', 8, ...
        'LineWidth', 1.2, 'DisplayName', 'jSQI-rejected pulse');
end

legendHandle = legend(axesHandles(1), 'Location', 'northoutside', ...
    'Orientation', 'horizontal', 'NumColumns', 3);
legendHandle.Color = 'w';
legendHandle.TextColor = 'k';
legendHandle.EdgeColor = [0.5 0.5 0.5];

hourTags = arrayfun(@(h) sprintf('%gh', h), targetHours, ...
    'UniformOutput', false);
fileTag = strrep(strjoin(hourTags, '_'), '.', 'p');
baseFilename = [char(patientID) '_ABP_features_' fileTag];
pngFile = fullfile(outputFolder, [baseFilename '.png']);
figFile = fullfile(outputFolder, [baseFilename '.fig']);
if isfile(pngFile), delete(pngFile); end
if isfile(figFile), delete(figFile); end
exportgraphics(fig, pngFile, 'Resolution', 300);
savefig(fig, figFile);

fprintf('Saved figure to:\n  %s\n  %s\n', pngFile, figFile);
end

function [timeSeconds, abp] = readABPWindow(filename, startTime, endTime, Fs)
%READABPWINDOW Read a short interval from a large two-column ABP file.
%
% The first value on each line is elapsed time in seconds and the second is
% ABP in mmHg. A timestamp-guided binary seek jumps near the requested time
% without loading or scanning the preceding hours. This works with both
% fixed-width and occasionally variable-width rows in the patient files.

fid = fopen(filename, 'rb');
assert(fid >= 0, 'Could not open ABP file: %s', filename);
cleanupObject = onCleanup(@() fclose(fid)); %#ok<NASGU>

firstLine = fgetl(fid);
assert(ischar(firstLine), 'The ABP file is empty.');
firstValues = sscanf(firstLine, '%f%f');
assert(numel(firstValues) == 2, 'Could not parse the first ABP row.');
firstTime = firstValues(1);
assert(startTime >= firstTime, ...
    ['The requested analysis window begins at %.3f h, but this patient''s ', ...
     'ABP file begins at %.3f h. Choose a later target hour.'], ...
    startTime/3600, firstTime/3600);

fseek(fid, 0, 'eof');
fileBytes = ftell(fid);
recordStart = findRecordAtOrAfterTime(fid, startTime, fileBytes);
seekStatus = fseek(fid, recordStart, 'bof');
assert(seekStatus == 0, 'Could not seek to the requested ABP interval.');

% Read two extra samples to absorb floating-point rounding at the edges,
% then retain only timestamps inside the requested interval.
numberOfRows = ceil((endTime - startTime) * Fs) + 2;

dataCell = textscan(fid, '%f%f', numberOfRows, ...
    'CollectOutput', true, 'ReturnOnError', false);
data = dataCell{1};
data = data(data(:, 1) <= endTime + 0.5/Fs, :);
assert(size(data, 1) >= floor((endTime - startTime) * Fs), ...
    'ABP file ended before the requested interval was completely read.');

timeSeconds = data(:, 1);
abp = data(:, 2);

% Confirm that binary seeking found the first sample at or after STARTTIME.
assert(timeSeconds(1) >= startTime - 0.5/Fs && ...
       timeSeconds(1) <= startTime + 1.5/Fs, ...
    'Timestamp seek did not land at the requested ABP interval.');
assert(all(abs(diff(timeSeconds) - 1/Fs) < 1e-6), ...
    'ABP timestamps are not uniformly sampled at 125 Hz.');
end

function recordStart = findRecordAtOrAfterTime(fid, targetTime, fileBytes)
%FINDRECORDATORAFTERTIME Locate a timestamp using byte-position bisection.

lowByte = 0;
highByte = fileBytes - 1;

% Each probe discards a possible partial row, reads the next complete row,
% and uses its timestamp to choose the earlier or later half of the file.
while highByte - lowByte > 4096
    probeByte = floor((lowByte + highByte) / 2);
    fseek(fid, probeByte, 'bof');
    if probeByte > 0
        fgetl(fid);  % Discard the partial row at the probe position.
    end

    probeStart = ftell(fid);
    probeLine = fgetl(fid);
    if ~ischar(probeLine)
        highByte = probeByte;
        continue
    end

    probeValues = sscanf(probeLine, '%f%f');
    if numel(probeValues) < 2
        highByte = probeByte;
        continue
    end

    if probeValues(1) < targetTime
        lowByte = ftell(fid);
    else
        highByte = probeStart;
    end
end

% Start slightly before the lower bound, align to a complete row, and scan
% the small remaining interval to find the first timestamp >= TARGETTIME.
scanByte = max(0, lowByte - 4096);
fseek(fid, scanByte, 'bof');
if scanByte > 0
    fgetl(fid);
end

while true
    candidateStart = ftell(fid);
    candidateLine = fgetl(fid);
    assert(ischar(candidateLine), ...
        'Requested time occurs after the end of the ABP file.');
    candidateValues = sscanf(candidateLine, '%f%f');
    if numel(candidateValues) == 2 && candidateValues(1) >= targetTime
        recordStart = candidateStart;
        return
    end
end
end
