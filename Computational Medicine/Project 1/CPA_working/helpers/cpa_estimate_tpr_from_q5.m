function results = cpa_estimate_tpr_from_q5(q5PatientResult)
%CPA_ESTIMATE_TPR_FROM_Q5 Estimate continuous TPR from saved Parlikar output.
%
% RESULTS = CPA_ESTIMATE_TPR_FROM_Q5(Q5PATIENTRESULT) consumes one patient
% result returned by CPA_COMPARE_ALGORITHMS_C2. It does not rerun WABP,
% feature extraction, or the Parlikar estimator.

assert(isstruct(q5PatientResult) && isfield(q5PatientResult,'parlikar'), ...
    'Input must be one patient result returned by the Q5 comparison.');
p = q5PatientResult.parlikar;
requiredFields = {'C2','tauSeconds','coCalibratedLminAllBeats', ...
    'MAPmmHg','deltaPmmHg','beatPeriodSeconds','timeMinutesAllBeats', ...
    'validOutputMask'};
for k = 1:numel(requiredFields)
    assert(isfield(p,requiredFields{k}), ...
        'Q5 Parlikar output is missing field %s. Rerun Q5.', ...
        requiredFields{k});
end

%% Recover arterial compliance from the Q5 calibration factor
% Q5 maps [mmHg/s] to [L/min]. Multiplication by 1000/60 converts that
% scale to [mL/s], leaving compliance in [mL/mmHg].
C = p.C2*1000/60;  % same as p.C2/0.06

tau = p.tauSeconds;
MAP = p.MAPmmHg;
deltaP = p.deltaPmmHg;
T = p.beatPeriodSeconds;
coLmin = p.coCalibratedLminAllBeats;
co_mLs = coLmin*1000/60;
beatTimeHours = p.timeMinutesAllBeats/60;

%% TPR from R=tau/C and independent verification with paper equation 11
TPR = tau./C;
denominator = co_mLs-C.*deltaP./T;
TPRcheck = MAP./denominator;

valid = p.validOutputMask & isfinite(TPR) & isfinite(TPRcheck) & ...
    isfinite(beatTimeHours) & tau>0 & C>0 & denominator>0;
TPR(~valid) = NaN;
TPRcheck(~valid) = NaN;

agreementError = TPRcheck-TPR;
maxAgreementError = max(abs(agreementError),[],'omitnan');

%% One-minute median trend for a readable continuous plot
endHour = q5PatientResult.endHour;
minuteEdges = 0:endHour*60;
minuteTimeHours = minuteEdges(2:end)/60;
tprTrend = nan(numel(minuteTimeHours),1);
for m = 1:numel(minuteTimeHours)
    inMinute = p.timeMinutesAllBeats>minuteEdges(m) & ...
        p.timeMinutesAllBeats<=minuteEdges(m+1);
    tprTrend(m) = median(TPR(inMinute),'omitnan');
end

%% Plot the one-minute median TPR used for presentation
% Beat-to-beat TPR remains available in RESULTS.TPR for auditing. The plot
% uses only non-overlapping one-minute medians so tens of thousands of beat
% values and isolated outliers do not obscure the physiologic trend.
fig = figure('Color','w','Position',[100 100 1500 650]);
ax = axes(fig);
plot(ax,minuteTimeHours,tprTrend,'b-','LineWidth',1.4, ...
    'DisplayName','One-minute median TPR');
xlabel(ax,'Time (hours)');
ylabel(ax,'TPR [mmHg/(mL/s)]');
title(ax,sprintf(['Patient %s: Parlikar TPR (C = %.3f mL/mmHg; ', ...
    'first-pair C2 calibration)'],q5PatientResult.patientID,C));
xlim(ax,[0 endHour]);
grid(ax,'on');
box(ax,'on');
legend(ax,'Location','best');

% Prevent a few isolated numerical spikes from flattening the useful trace.
finiteTPR = sort(TPR(isfinite(TPR)));
if numel(finiteTPR)>=20
    upperIndex = max(1,round(0.995*numel(finiteTPR)));
    upperLimit = finiteTPR(upperIndex);
    if isfinite(upperLimit) && upperLimit>0
        ylim(ax,[0 1.1*upperLimit]);
    end
end

paths = cpa_setup_paths();
baseName = sprintf('q6_%s_parlikar_TPR_first%gh', ...
    q5PatientResult.patientID,endHour);
pngFile = fullfile(paths.resultsFolder,[baseName '.png']);
figFile = fullfile(paths.resultsFolder,[baseName '.fig']);
if isfile(pngFile), delete(pngFile); end
if isfile(figFile), delete(figFile); end
exportgraphics(fig,pngFile,'Resolution',300);
savefig(fig,figFile);

results = struct();
results.patientID = q5PatientResult.patientID;
results.endHour = endHour;
results.C2 = p.C2;
results.compliance_mL_per_mmHg = C;
results.beatTimeHours = beatTimeHours;
results.tauSeconds = tau;
results.TPR = TPR;
results.TPRcheck = TPRcheck;
results.validMask = valid;
results.validFraction = mean(valid);
results.maxEquation11AgreementError = maxAgreementError;
results.minuteTimeHours = minuteTimeHours;
results.tprTrend = tprTrend;
results.pngFile = pngFile;
results.figFile = figFile;

fprintf(['%s TPR: C = %.4f mL/mmHg, valid beats = %.1f%%, ', ...
    'max Eq.11 check difference = %.3g\n'],q5PatientResult.patientID,C, ...
    100*results.validFraction,maxAgreementError);
fprintf('Saved Q6 TPR figure to:\n  %s\n  %s\n',pngFile,figFile);
end
