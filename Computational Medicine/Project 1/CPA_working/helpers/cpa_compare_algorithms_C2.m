function results = cpa_compare_algorithms_C2( ...
    patientID,endHour,estimatorIDs,outputTag)
%CPA_COMPARE_ALGORITHMS_C2 Compare calibrated CO-from-ABP algorithms.
%
% RESULTS = CPA_COMPARE_ALGORITHMS_C2(PATIENTID, ENDHOUR, ESTIMATORIDS)
% processes the patient waveform from its start through ENDHOUR, applies
% the same ABP preprocessing and jSQI filtering to every estimator, and
% independently calibrates each algorithm using the first TCO pairing (C2).
%
% Default comparison:
%   ID 5 - Liljestrand
%   ID 1 - mean arterial pressure (MAP) baseline
%   ID 7 - Wesseling corrected impedance
%   ID 14 - Parlikar beat-averaged Windkessel estimator
%
% Performance metrics exclude the first TCO point because that point was
% used for calibration and therefore is not an independent validation.

if nargin < 2 || isempty(endHour)
    endHour = 12;
end
if nargin < 3 || isempty(estimatorIDs)
    estimatorIDs = [5, 1, 7, 14];
end
if nargin < 4 || strlength(string(outputTag))==0
    outputTag = "q4";
end
outputTag = lower(string(outputTag));
assert(isscalar(outputTag) && all(isstrprop(char(outputTag),'alphanum')), ...
    'outputTag must contain only letters and numbers, such as q4 or q5.');
validateattributes(endHour, {'numeric'}, ...
    {'scalar', 'positive', 'finite'}, mfilename, 'endHour');
validateattributes(estimatorIDs, {'numeric'}, ...
    {'vector', 'integer'}, mfilename, 'estimatorIDs');
assert(~any(ismember(estimatorIDs, [8, 9, 11])), ...
    'Estimators 8, 9, and 11 are not permitted for this assignment.');

paths = cpa_setup_paths();
patientFiles = cpa_find_patient_files(patientID,paths);
patientID = patientFiles.patientID;
Fs = 125;
endTimeSeconds = endHour*3600;

%% Locate patient waveform and trend/TCO files
abpFile = patientFiles.abpFile;
trendFile = patientFiles.trendFile;
assert(~isempty(trendFile), ...
    'No trend/TCO file was found for patient %s.',patientID);

%% Run the shared ABP signal-processing pipeline once
fprintf('Reading %s through %.1f hours...\n', patientID, endHour);
[time, ABP] = readABPThroughTime(abpFile, endTimeSeconds, Fs);

fprintf('Running wabp, abpfeature, and jSQI for %s...\n', patientID);
t_on = wabp(ABP);
feat = abpfeature(ABP, t_on);
[beatq, fractionGood] = jSQI(feat, t_on, ABP);

%% Read intermittent TCO reference measurements
[trendTimeSeconds, trendCO] = readTrendCO(trendFile);
hasTCO = trendCO > 0 & ...
         trendTimeSeconds >= time(1) & ...
         trendTimeSeconds <= endTimeSeconds;
tcoTimeMinutes = trendTimeSeconds(hasTCO)/60;
tcoMeasured = trendCO(hasTCO);
assert(numel(tcoMeasured) >= 2, ...
    ['At least two overlapping TCO measurements are required to ', ...
     'calibrate and independently compare performance for %s.'], patientID);

%% Estimate and independently C2-calibrate each algorithm
nAlgorithms = numel(estimatorIDs);
algorithmNames = strings(nAlgorithms,1);
C2 = nan(nAlgorithms,1);
MAE = nan(nAlgorithms,1);
RMSE = nan(nAlgorithms,1);
Bias = nan(nAlgorithms,1);
NRMSEpct = nan(nAlgorithms,1);
NValidation = zeros(nAlgorithms,1);
CoveragePct = nan(nAlgorithms,1);
pairedPredictions = nan(numel(tcoMeasured), nAlgorithms);
minuteTrends = cell(nAlgorithms,1);
algorithmDiagnostics = cell(nAlgorithms,1);

% PP, MAP, and HR describe the shared ABP waveform. They do not depend on
% which CO estimator is selected, so calculate them once and reuse them in
% the Figure-4-style comparison rather than drawing duplicate copies.
pairedPP = nan(numel(tcoMeasured),1);
pairedMAP = nan(numel(tcoMeasured),1);
pairedHR = nan(numel(tcoMeasured),1);

minuteEdges = 0:1:endHour*60;
minuteTimeMinutes = minuteEdges(2:end);
minuteTimeHours = minuteTimeMinutes/60;

%% Build PP, MAP, and HR trends once from the shared good ABP beats
% These features describe the ABP waveform, not a particular CO estimator.
nFeatureBeats = size(feat,1);
goodFeatureBeat = ~logical(beatq(1:nFeatureBeats,1));
featureTimeMinutes = t_on(1:nFeatureBeats)/(Fs*60) + ...
    (time(1)-1/Fs)/60;
PP = feat(:,5);
MAP = feat(:,6);
HR = 60*Fs./feat(:,7);

ppTrend = nan(numel(minuteTimeMinutes),1);
mapTrend = nan(numel(minuteTimeMinutes),1);
hrTrend = nan(numel(minuteTimeMinutes),1);
for m = 1:numel(minuteTimeMinutes)
    inMinute = goodFeatureBeat & ...
        featureTimeMinutes>minuteEdges(m) & ...
        featureTimeMinutes<=minuteEdges(m+1);
    ppTrend(m) = median(PP(inMinute),'omitnan');
    mapTrend(m) = median(MAP(inMinute),'omitnan');
    hrTrend(m) = median(HR(inMinute),'omitnan');
end
for j = 1:numel(tcoMeasured)
    precedingMinute = goodFeatureBeat & ...
        featureTimeMinutes>tcoTimeMinutes(j)-1 & ...
        featureTimeMinutes<=tcoTimeMinutes(j);
    pairedPP(j) = median(PP(precedingMinute),'omitnan');
    pairedMAP(j) = median(MAP(precedingMinute),'omitnan');
    pairedHR(j) = median(HR(precedingMinute),'omitnan');
end

for a = 1:nAlgorithms
    estimatorID = estimatorIDs(a);
    algorithmNames(a) = estimatorName(estimatorID);

    if estimatorID==14
        % Parlikar requires the original onset pressures ABP(T_ON), which
        % are not accepted by the supplied ESTIMATECO_V3 interface.
        [relativeCO,coTimeMinutes,~,parlikarDiagnostics] = ...
            cpa_est14_parlikar(ABP,t_on,feat,beatq,Fs,361);
        parlikarDiagnostics.timeMinutesAllBeats = ...
            parlikarDiagnostics.timeMinutesAllBeats + ...
            (time(1)-1/Fs)/60;
        algorithmDiagnostics{a} = parlikarDiagnostics;
    else
        [relativeCO,coTimeMinutes,~,~] = ...
            estimateCO_v3(t_on,feat,beatq,estimatorID,1);
        algorithmDiagnostics{a} = struct();
    end
    coTimeMinutes = coTimeMinutes + (time(1)-1/Fs)/60;
    CoveragePct(a) = 100*numel(relativeCO)/nFeatureBeats;

    % Pair every TCO with the preceding minute of the current estimator.
    pairedRelative = nan(numel(tcoMeasured),1);
    for j = 1:numel(tcoMeasured)
        precedingMinute = ...
            coTimeMinutes > tcoTimeMinutes(j)-1 & ...
            coTimeMinutes <= tcoTimeMinutes(j);
        pairedRelative(j) = median( ...
            relativeCO(precedingMinute), 'omitnan');
    end

    assert(isfinite(pairedRelative(1)) && pairedRelative(1) > 0, ...
        'The first TCO has no valid preceding-minute estimate for %s.', ...
        algorithmNames(a));

    % Each estimator has a different arbitrary relative scale, so each gets
    % its own C2 factor from the same first measured TCO.
    C2(a) = tcoMeasured(1)/pairedRelative(1);
    calibratedBeatCO = C2(a).*relativeCO;
    pairedPredictions(:,a) = C2(a).*pairedRelative;
    if estimatorID==14
        parlikarDiagnostics.C2 = C2(a);
        parlikarDiagnostics.compliance_mL_per_mmHg = C2(a)/0.06;
        parlikarDiagnostics.coCalibratedLminAllBeats = ...
            C2(a).*parlikarDiagnostics.coRawAllBeats;
        algorithmDiagnostics{a} = parlikarDiagnostics;
    end

    % One median per consecutive minute gives a readable continuous trend.
    trend = nan(numel(minuteTimeMinutes),1);
    for m = 1:numel(minuteTimeMinutes)
        inMinute = coTimeMinutes > minuteEdges(m) & ...
                   coTimeMinutes <= minuteEdges(m+1);
        trend(m) = median(calibratedBeatCO(inMinute), 'omitnan');
    end
    minuteTrends{a} = trend;

end

%% Calculate performance on one common set of validation TCO points
% Point #1 is excluded because it defines C2. If any algorithm lacks a
% preceding-minute estimate at a later TCO, that TCO is excluded for every
% algorithm. This keeps the numerical comparison apples-to-apples.
validationRows = (2:numel(tcoMeasured)).';
commonValid = isfinite(tcoMeasured(validationRows)) & ...
    all(isfinite(pairedPredictions(validationRows,:)),2);
commonValidationRows = validationRows(commonValid);
assert(~isempty(commonValidationRows), ...
    'No common post-calibration TCO points are available for comparison.');

measuredCommon = tcoMeasured(commonValidationRows);
for a = 1:nAlgorithms
    predicted = pairedPredictions(commonValidationRows,a);
    errors = predicted-measuredCommon;
    NValidation(a) = numel(errors);
    MAE(a) = mean(abs(errors));
    RMSE(a) = sqrt(mean(errors.^2));
    Bias(a) = mean(errors);
    NRMSEpct(a) = 100*RMSE(a)/mean(measuredCommon);
end

metricsTable = table(estimatorIDs(:),algorithmNames,C2,MAE,RMSE, ...
    Bias,NRMSEpct,NValidation,CoveragePct, ...
    'VariableNames', {'EstimatorID','Algorithm','C2','MAE_Lmin', ...
    'RMSE_Lmin','Bias_Lmin','NRMSE_pct','NValidation','Coverage_pct'});

%% Plot continuous trends, paired predictions, and validation errors
colors = lines(nAlgorithms);
tcoTimeHours = tcoTimeMinutes/60;
fig = figure('Color','w','Position',[80 80 1500 900]);
layout = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
title(layout, sprintf(['Patient %s: C2-calibrated CO algorithm ', ...
    'comparison (first %.0f h)'], patientID, endHour), ...
    'FontWeight','bold','Color','k');

ax = gobjects(3,1);
ax(1) = nexttile(layout);
hold(ax(1),'on');
for a = 1:nAlgorithms
    plot(ax(1), minuteTimeHours, minuteTrends{a}, ...
        'Color',colors(a,:),'LineWidth',1.2, ...
        'DisplayName',algorithmNames(a));
end
stem(ax(1),tcoTimeHours,tcoMeasured,'k','filled', ...
    'LineWidth',1,'MarkerSize',4,'DisplayName','Measured TCO');
ylabel(ax(1),'CO (L/min)');
legend(ax(1),'Location','best','NumColumns',2);

ax(2) = nexttile(layout);
hold(ax(2),'on');
plot(ax(2),tcoTimeHours,tcoMeasured,'ko-','LineWidth',1.2, ...
    'MarkerFaceColor','k','DisplayName','Measured TCO');
for a = 1:nAlgorithms
    plot(ax(2),tcoTimeHours,pairedPredictions(:,a),'o-', ...
        'Color',colors(a,:),'LineWidth',1.1,'MarkerSize',5, ...
        'DisplayName',algorithmNames(a) + " paired estimate");
end
ylabel(ax(2),'Paired CO (L/min)');

ax(3) = nexttile(layout);
hold(ax(3),'on');
yline(ax(3),0,'k:','LineWidth',1);
for a = 1:nAlgorithms
    errors = nan(size(tcoMeasured));
    errors(commonValidationRows) = ...
        pairedPredictions(commonValidationRows,a)-measuredCommon;
    plot(ax(3),tcoTimeHours,errors,'o-', ...
        'Color',colors(a,:),'LineWidth',1.1,'MarkerSize',5, ...
        'DisplayName',algorithmNames(a));
end
ylabel(ax(3),'Error (L/min)');
xlabel(ax(3),'Time (hours)');

for k = 1:numel(ax)
    ax(k).Color='w';
    ax(k).XColor='k';
    ax(k).YColor='k';
    ax(k).GridColor=[0.75 0.75 0.75];
    ax(k).GridAlpha=0.4;
    ax(k).Toolbar.Visible='off';
    grid(ax(k),'on');
    box(ax(k),'on');
    xlim(ax(k),[0 endHour]);
end
linkaxes(ax,'x');

baseName = sprintf('%s_%s_algorithm_comparison_C2_first%gh', ...
    outputTag,patientID,endHour);
pngFile = fullfile(paths.resultsFolder,[baseName '.png']);
figFile = fullfile(paths.resultsFolder,[baseName '.fig']);
csvFile = fullfile(paths.resultsFolder,[baseName '_metrics.csv']);
if isfile(pngFile), delete(pngFile); end
if isfile(figFile), delete(figFile); end
exportgraphics(fig,pngFile,'Resolution',300);
savefig(fig,figFile);
writetable(metricsTable,csvFile);

%% Make a Figure-4-style plot with CO, PP, MAP, and HR
% The top panel compares all independently C2-calibrated CO
% algorithms are placed on the same axes. The three lower panels provide
% the physiologic ABP context requested in Q3. Because every estimator uses
% the same accepted beats, PP, MAP, and HR are common to all algorithms.
clinicalFig = figure('Color','w','Position',[100 60 1500 1000]);
clinicalLayout = tiledlayout(clinicalFig,4,1, ...
    'TileSpacing','compact','Padding','compact');
title(clinicalLayout,sprintf( ...
    'Patient %s: calibrated CO algorithms and ABP features (first %.0f h)', ...
    patientID,endHour),'FontWeight','bold','Color','k');

clinicalAx = gobjects(4,1);
clinicalAx(1) = nexttile(clinicalLayout);
hold(clinicalAx(1),'on');
for a = 1:nAlgorithms
    plot(clinicalAx(1),minuteTimeHours,minuteTrends{a}, ...
        'Color',colors(a,:),'LineWidth',1.1, ...
        'DisplayName',algorithmNames(a));
end
stem(clinicalAx(1),tcoTimeHours,tcoMeasured,'k','filled', ...
    'LineWidth',1,'MarkerSize',4,'DisplayName','Measured TCO');
ylabel(clinicalAx(1),'CO (L/min)');
legend(clinicalAx(1),'Location','best','NumColumns',2);

clinicalAx(2) = nexttile(clinicalLayout);
plot(clinicalAx(2),minuteTimeHours,ppTrend, ...
    'Color',[0.35 0.35 0.35],'LineWidth',1);
hold(clinicalAx(2),'on');
stem(clinicalAx(2),tcoTimeHours,pairedPP,'k','filled', ...
    'LineWidth',0.9,'MarkerSize',3);
ylabel(clinicalAx(2),'PP (mmHg)');

clinicalAx(3) = nexttile(clinicalLayout);
plot(clinicalAx(3),minuteTimeHours,mapTrend, ...
    'Color',[0.35 0.35 0.35],'LineWidth',1);
hold(clinicalAx(3),'on');
stem(clinicalAx(3),tcoTimeHours,pairedMAP,'k','filled', ...
    'LineWidth',0.9,'MarkerSize',3);
ylabel(clinicalAx(3),'MAP (mmHg)');

clinicalAx(4) = nexttile(clinicalLayout);
plot(clinicalAx(4),minuteTimeHours,hrTrend, ...
    'Color',[0.35 0.35 0.35],'LineWidth',1);
hold(clinicalAx(4),'on');
stem(clinicalAx(4),tcoTimeHours,pairedHR,'k','filled', ...
    'LineWidth',0.9,'MarkerSize',3);
ylabel(clinicalAx(4),'HR (bpm)');
xlabel(clinicalAx(4),'Time (hours)');

for k = 1:numel(clinicalAx)
    clinicalAx(k).Color='w';
    clinicalAx(k).XColor='k';
    clinicalAx(k).YColor='k';
    clinicalAx(k).GridColor=[0.75 0.75 0.75];
    clinicalAx(k).GridAlpha=0.4;
    clinicalAx(k).Toolbar.Visible='off';
    grid(clinicalAx(k),'on');
    box(clinicalAx(k),'on');
    xlim(clinicalAx(k),[0 endHour]);
end
linkaxes(clinicalAx,'x');

clinicalBaseName = sprintf('%s_%s_figure4_style_C2_first%gh', ...
    outputTag,patientID,endHour);
clinicalPngFile = fullfile(paths.resultsFolder, ...
    [clinicalBaseName '.png']);
clinicalFigFile = fullfile(paths.resultsFolder, ...
    [clinicalBaseName '.fig']);
if isfile(clinicalPngFile), delete(clinicalPngFile); end
if isfile(clinicalFigFile), delete(clinicalFigFile); end
exportgraphics(clinicalFig,clinicalPngFile,'Resolution',300);
savefig(clinicalFig,clinicalFigFile);

%% Return results
results = struct();
results.patientID = patientID;
results.endHour = endHour;
results.outputTag = outputTag;
results.estimatorIDs = estimatorIDs(:);
results.algorithmNames = algorithmNames;
results.C2 = C2;
results.fractionGood = fractionGood;
results.minuteTimeHours = minuteTimeHours;
results.minuteTrends = minuteTrends;
results.tcoTimeHours = tcoTimeHours;
results.tcoMeasured = tcoMeasured;
results.pairedPredictions = pairedPredictions;
results.commonValidationRows = commonValidationRows;
results.commonValidationTimeHours = ...
    tcoTimeHours(commonValidationRows);
results.metricsTable = metricsTable;
results.algorithmDiagnostics = algorithmDiagnostics;
parlikarIndex = find(estimatorIDs==14,1,'first');
if isempty(parlikarIndex)
    results.parlikar = struct();
else
    results.parlikar = algorithmDiagnostics{parlikarIndex};
    results.parlikar.estimatorIndex = parlikarIndex;
end
results.ppTrend = ppTrend;
results.mapTrend = mapTrend;
results.hrTrend = hrTrend;
results.pairedPP = pairedPP;
results.pairedMAP = pairedMAP;
results.pairedHR = pairedHR;
results.pngFile = pngFile;
results.figFile = figFile;
results.csvFile = csvFile;
results.figure4StylePngFile = clinicalPngFile;
results.figure4StyleFigFile = clinicalFigFile;

fprintf('\n%s algorithm comparison (validation excludes calibration point):\n', ...
    patientID);
disp(metricsTable);
fprintf(['Saved comparison outputs to:\n  %s\n  %s\n  %s\n', ...
    '  %s\n  %s\n'],pngFile,figFile,csvFile, ...
    clinicalPngFile,clinicalFigFile);
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
    case 14, name = "Parlikar";
    otherwise, name = "Estimator " + estimatorID;
end
end

function [time,ABP] = readABPThroughTime(filename,endTimeSeconds,Fs)
fid = fopen(filename,'rb');
assert(fid>=0,'Could not open ABP file: %s',filename);
cleanupObject = onCleanup(@() fclose(fid)); %#ok<NASGU>
firstLine = fgetl(fid);
firstValues = sscanf(firstLine,'%f%f');
assert(numel(firstValues)==2,'Could not parse first ABP row.');
firstTime = firstValues(1);
assert(firstTime < endTimeSeconds, ...
    'ABP waveform starts after the requested comparison timeframe.');
frewind(fid);
numberOfRows = floor((endTimeSeconds-firstTime)*Fs)+2;
dataCell = textscan(fid,'%f%f',numberOfRows, ...
    'CollectOutput',true,'ReturnOnError',false);
data = dataCell{1};
data = data(data(:,1)<=endTimeSeconds+0.5/Fs,:);
time = data(:,1);
ABP = data(:,2);
assert(time(end)>=endTimeSeconds-1/Fs, ...
    'ABP file does not extend through the requested timeframe.');
end

function [elapsedSeconds,CO] = readTrendCO(filename)
fid = fopen(filename,'rt');
assert(fid>=0,'Could not open trend file: %s',filename);
cleanupObject = onCleanup(@() fclose(fid)); %#ok<NASGU>
fgetl(fid);
fgetl(fid);
raw = textscan(fid,repmat('%s',1,16),'Delimiter','\t', ...
    'MultipleDelimsAsOne',false,'ReturnOnError',false);
elapsedSeconds = str2double(raw{1});
CO = str2double(raw{16});
end
