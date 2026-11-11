function [coRelative,timeMinutes,goodFeatures,diagnostics] = ...
    cpa_est14_parlikar(ABP,t_on,feat,beatq,Fs,windowBeats)
%CPA_EST14_PARLIKAR Parlikar beat-averaged Windkessel CO estimator.
%
% [CORELATIVE,TIMEMINUTES,GOODFEATURES,DIAGNOSTICS] =
% CPA_EST14_PARLIKAR(ABP,T_ON,FEAT,BEATQ,FS,WINDOWBEATS) implements
% equations 6, 8, 9, and 10 of Parlikar et al. (2007).
%
% The returned CORELATIVE is intentionally uncalibrated. Its scale includes
% the unknown arterial compliance. The caller must use the same TCO/C2
% calibration used for the other CO-from-ABP estimators.
%
% Inputs use the same beat definitions as ABPFEATURE:
%   FEAT has N rows, while T_ON has N+1 pulse onsets.
%   FEAT(:,4) = DAP, FEAT(:,6) = MAP, FEAT(:,7) = period [samples].
%
% A 361-beat window is the default. The paper describes an odd centered
% window but reports 360-beat windows in its results table; 361 preserves
% an exactly centered estimate with 180 beats on either side.

if nargin < 5 || isempty(Fs), Fs = 125; end
if nargin < 6 || isempty(windowBeats), windowBeats = 361; end
validateattributes(Fs,{'numeric'},{'scalar','positive','finite'}, ...
    mfilename,'Fs');
validateattributes(windowBeats,{'numeric'}, ...
    {'scalar','integer','positive'},mfilename,'windowBeats');
assert(mod(windowBeats,2)==1, ...
    'windowBeats must be odd so the window has one center beat.');

N = size(feat,1);
assert(numel(t_on)>=N+1, ...
    'T_ON must contain at least N+1 onsets for N feature rows.');
assert(size(beatq,1)>=N, ...
    'BEATQ must contain at least one quality row per feature beat.');

ABP = ABP(:);
t_on = t_on(:);

%% Quantities in equations 6, 8, and 9
% Beat n runs from onset t_n to t_(n+1), so every vector below has N rows.
onsetPressure = ABP(t_on(1:N+1));
deltaP = diff(onsetPressure);            % Eq. 6: P(t_n+1)-P(t_n)
MAP = feat(:,6);                         % cycle-averaged pressure
DAP = feat(:,4);                         % diastolic pressure
T = diff(t_on(1:N+1))/Fs;                % beat duration [s]

alpha = 2;
estimatedPP = alpha.*(MAP-DAP);          % Eq. 9

% Rearrange Eq. 8 as MAP*beta = (estimatedPP-deltaP)/T, beta=1/tau.
x = MAP;
y = (estimatedPP-deltaP)./T;
goodBeat = ~logical(beatq(1:N,1));
validEquation = goodBeat & isfinite(x) & isfinite(y) & ...
    isfinite(T) & T>0;

%% Centered-window least-squares estimate of beta = 1/tau
% For a scalar no-intercept regression x*beta=y, the least-squares result
% is sum(x*y)/sum(x^2). MOVSUM computes that exact solution for every
% centered window without running a separate backslash operation per beat.
halfWindow = floor(windowBeats/2);
xy = zeros(N,1);
xx = zeros(N,1);
xy(validEquation) = x(validEquation).*y(validEquation);
xx(validEquation) = x(validEquation).^2;

sumXY = movsum(xy,[halfWindow halfWindow]);
sumXX = movsum(xx,[halfWindow halfWindow]);
validCount = movsum(double(validEquation),[halfWindow halfWindow]);

beta = nan(N,1);
hasFullWindow = false(N,1);
if N>=windowBeats
    hasFullWindow(halfWindow+1:N-halfWindow) = true;
end
enoughValid = validCount>=ceil(0.80*windowBeats);
canEstimate = hasFullWindow & enoughValid & sumXX>0;
betaCandidate = sumXY(canEstimate)./sumXX(canEstimate);
keepBeta = isfinite(betaCandidate) & betaCandidate>0;
candidateIndex = find(canEstimate);
beta(candidateIndex(keepBeta)) = betaCandidate(keepBeta);

%% Equation 10 without the unknown arterial compliance C_n
% The quantity in parentheses has pressure/time units. C2 calibration later
% supplies the patient/algorithm-specific scale needed to express L/min.
coRaw = deltaP./T + MAP.*beta;

% Preserve the original chronological sequence during window estimation;
% only now remove bad center beats and beats without a valid time constant.
validOutput = goodBeat & isfinite(coRaw);
coRelative = coRaw(validOutput);
timeAllMinutes = t_on(1:N)/(Fs*60);
timeMinutes = timeAllMinutes(validOutput);
goodFeatures = feat(validOutput,:);

%% Diagnostics used to document reproducibility and data coverage
diagnostics = struct();
diagnostics.windowBeats = windowBeats;
diagnostics.alpha = alpha;
diagnostics.beta = beta;
diagnostics.tauSeconds = 1./beta;
diagnostics.coRawAllBeats = coRaw;
diagnostics.timeMinutesAllBeats = timeAllMinutes;
diagnostics.deltaPmmHg = deltaP;
diagnostics.beatPeriodSeconds = T;
diagnostics.MAPmmHg = MAP;
diagnostics.DAPmmHg = DAP;
diagnostics.estimatedPPmmHg = estimatedPP;
diagnostics.goodBeatMask = goodBeat;
diagnostics.validOutputMask = validOutput;
diagnostics.validEquationFraction = mean(validEquation);
diagnostics.outputFractionAllBeats = mean(validOutput);
diagnostics.nInputBeats = N;
diagnostics.nOutputBeats = nnz(validOutput);
diagnostics.nNonpositiveRawCO = nnz(validOutput & coRaw<=0);
end
