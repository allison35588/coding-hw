function outputs = run_Q5_parlikar_comparison( ...
    patientIDs,endHours,estimatorIDs)
%RUN_Q5_PARLIKAR_COMPARISON Run Q5 for any list of downloaded patients.
%
% Examples:
%   outputs = run_Q5_parlikar_comparison();
%   outputs = run_Q5_parlikar_comparison([20 214],12);
%   outputs = run_Q5_parlikar_comparison( ...
%       ["s00020" "s00214" "s00321"],[12 12 12],[5 1 7 14]);
%
% ENDHOURS may be one scalar shared by all patients or one value per
% patient. Every estimator within a patient uses the same timeframe,
% preprocessing, TCO pairing, and first-pair C2 calibration.

if nargin<1 || isempty(patientIDs)
    patientIDs = [20 214];
end
if nargin<2 || isempty(endHours)
    endHours = 12;
end
if nargin<3 || isempty(estimatorIDs)
    estimatorIDs = [5 1 7 14];
end

if isnumeric(patientIDs)
    patientIDs = string(patientIDs(:));
else
    patientIDs = string(patientIDs(:));
end
nPatients = numel(patientIDs);

validateattributes(endHours,{'numeric'}, ...
    {'vector','positive','finite'},mfilename,'endHours');
if isscalar(endHours)
    endHours = repmat(endHours,nPatients,1);
else
    endHours = endHours(:);
    assert(numel(endHours)==nPatients, ...
        'endHours must be scalar or have one value per patient.');
end
validateattributes(estimatorIDs,{'numeric'}, ...
    {'vector','integer'},mfilename,'estimatorIDs');

codeFolder = fileparts(mfilename('fullpath'));
addpath(fullfile(codeFolder,'helpers'));
paths = cpa_setup_paths();

%% Run the identical comparison workflow for every requested patient
patientResults = cell(nPatients,1);
patientLabels = strings(nPatients,1);
combinedMetrics = table();
for p = 1:nPatients
    fprintf('\n=== Q5: patient %s ===\n',patientIDs(p));
    patientResults{p} = cpa_compare_algorithms_C2( ...
        patientIDs(p),endHours(p),estimatorIDs,"q5");
    patientLabels(p) = patientResults{p}.patientID;

    patientTable = addvars(patientResults{p}.metricsTable, ...
        repmat(patientLabels(p),height(patientResults{p}.metricsTable),1), ...
        'Before',1,'NewVariableNames','Patient');
    if p==1
        combinedMetrics = patientTable;
    else
        combinedMetrics = [combinedMetrics;patientTable]; %#ok<AGROW>
    end
end

%% Macro averages give every patient equal weight
% This avoids letting a patient with many more TCO measurements dominate
% the overall conclusion.
nAlgorithms = numel(estimatorIDs);
algorithmNames = patientResults{1}.algorithmNames;
macroMAE = nan(nAlgorithms,1);
macroRMSE = nan(nAlgorithms,1);
macroAbsBias = nan(nAlgorithms,1);
macroNRMSE = nan(nAlgorithms,1);
for a = 1:nAlgorithms
    rows = combinedMetrics.EstimatorID==estimatorIDs(a);
    macroMAE(a) = mean(combinedMetrics.MAE_Lmin(rows),'omitnan');
    macroRMSE(a) = mean(combinedMetrics.RMSE_Lmin(rows),'omitnan');
    macroAbsBias(a) = mean(abs(combinedMetrics.Bias_Lmin(rows)),'omitnan');
    macroNRMSE(a) = mean(combinedMetrics.NRMSE_pct(rows),'omitnan');
end
macroMetrics = table(estimatorIDs(:),algorithmNames,macroMAE,macroRMSE, ...
    macroAbsBias,macroNRMSE,'VariableNames', ...
    {'EstimatorID','Algorithm','MacroMAE_Lmin','MacroRMSE_Lmin', ...
     'MacroAbsBias_Lmin','MacroNRMSE_pct'});

bestByRMSE = strings(nPatients,1);
for p = 1:nPatients
    [~,bestIndex] = min(patientResults{p}.metricsTable.RMSE_Lmin);
    bestByRMSE(p) = patientResults{p}.metricsTable.Algorithm(bestIndex);
end
[~,bestMacroIndex] = min(macroMetrics.MacroRMSE_Lmin);
bestTable = table(patientLabels,endHours,bestByRMSE, ...
    'VariableNames',{'Patient','EndHour','BestAlgorithmByRMSE'});

combinedCsv = fullfile(paths.resultsFolder, ...
    'q5_all_patients_metrics_C2.csv');
macroCsv = fullfile(paths.resultsFolder, ...
    'q5_macro_average_metrics_C2.csv');
bestCsv = fullfile(paths.resultsFolder, ...
    'q5_best_algorithm_by_patient_C2.csv');
writetable(combinedMetrics,combinedCsv);
writetable(macroMetrics,macroCsv);
writetable(bestTable,bestCsv);

%% Flexible summary chart: one tile per patient
figureHeight = max(650,420*ceil(nPatients/2));
fig = figure('Color','w','Position',[100 80 1500 figureHeight]);
layout = tiledlayout(fig,'flow','TileSpacing','compact','Padding','compact');
title(layout,['Q5 validation performance on common TCO points after ', ...
    'first-pair C2 calibration'],'FontWeight','bold');

for p = 1:nPatients
    ax = nexttile(layout);
    values = [patientResults{p}.metricsTable.MAE_Lmin, ...
              patientResults{p}.metricsTable.RMSE_Lmin];
    bar(ax,values,'grouped');
    ax.XTick = 1:nAlgorithms;
    ax.XTickLabel = cellstr(algorithmNames);
    ax.XTickLabelRotation = 20;
    ylabel(ax,'Error (L/min)');
    title(ax,sprintf('Patient %s (n = %d common validation TCOs)', ...
        patientLabels(p),patientResults{p}.metricsTable.NValidation(1)));
    legend(ax,{'MAE','RMSE'},'Location','best');
    grid(ax,'on');
    box(ax,'on');
end

summaryPng = fullfile(paths.resultsFolder, ...
    'q5_parlikar_performance_summary_C2.png');
summaryFig = fullfile(paths.resultsFolder, ...
    'q5_parlikar_performance_summary_C2.fig');
if isfile(summaryPng), delete(summaryPng); end
if isfile(summaryFig), delete(summaryFig); end
exportgraphics(fig,summaryPng,'Resolution',300);
savefig(fig,summaryFig);

%% Return both indexed and named access to patient results
outputs = struct();
outputs.patientIDs = patientLabels;
outputs.endHours = endHours;
outputs.estimatorIDs = estimatorIDs(:);
outputs.patientResults = patientResults;
for p = 1:nPatients
    fieldName = matlab.lang.makeValidName(char(patientLabels(p)));
    outputs.(fieldName) = patientResults{p};
end
outputs.combinedMetrics = combinedMetrics;
outputs.macroMetrics = macroMetrics;
outputs.bestByPatient = bestTable;
outputs.bestMacroByRMSE = macroMetrics.Algorithm(bestMacroIndex);
outputs.combinedCsv = combinedCsv;
outputs.macroCsv = macroCsv;
outputs.bestCsv = bestCsv;
outputs.summaryPng = summaryPng;
outputs.summaryFig = summaryFig;
outputs.q5ForQ6Mat = fullfile(paths.resultsFolder, ...
    'q5_outputs_for_q6.mat');
save(outputs.q5ForQ6Mat,'outputs','-v7.3');

fprintf('\n=== Q5 best algorithm by patient (common validation TCOs) ===\n');
disp(bestTable);
fprintf('\n=== Q5 macro-average performance (patients equally weighted) ===\n');
disp(macroMetrics);
fprintf('Best macro-average RMSE: %s\n',outputs.bestMacroByRMSE);
fprintf(['Saved Q5 summary files to:\n  %s\n  %s\n  %s\n  %s\n', ...
    'Saved Q5-to-Q6 data to:\n  %s\n'], ...
    combinedCsv,macroCsv,bestCsv,summaryPng,outputs.q5ForQ6Mat);
end
