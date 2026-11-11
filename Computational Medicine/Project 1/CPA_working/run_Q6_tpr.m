function outputs = run_Q6_tpr(q5Input)
%RUN_Q6_TPR Generate Parlikar TPR figures from Q5 outputs.
%
% Recommended same-session use:
%   q5 = run_Q5_parlikar_comparison([20 214],12);
%   q6 = run_Q6_tpr(q5);
%
% A saved Q5 MAT file can also be used:
%   q6 = run_Q6_tpr("results/q5_outputs_for_q6.mat");

codeFolder = fileparts(mfilename('fullpath'));
addpath(fullfile(codeFolder,'helpers'));
paths = cpa_setup_paths();

if nargin<1 || isempty(q5Input)
    q5Input = fullfile(paths.resultsFolder,'q5_outputs_for_q6.mat');
end
if ischar(q5Input) || (isstring(q5Input) && isscalar(q5Input))
    saved = load(q5Input,'outputs');
    assert(isfield(saved,'outputs'), ...
        'The MAT file does not contain a Q5 variable named outputs.');
    q5Input = saved.outputs;
end

assert(isstruct(q5Input), ...
    'q5Input must be Q5 outputs or the path to a saved Q5 MAT file.');
if isfield(q5Input,'patientResults')
    q5PatientResults = q5Input.patientResults;
else
    q5PatientResults = {q5Input};
end

nPatients = numel(q5PatientResults);
patientResults = cell(nPatients,1);
patientIDs = strings(nPatients,1);
for p = 1:nPatients
    patientResults{p} = cpa_estimate_tpr_from_q5(q5PatientResults{p});
    patientIDs(p) = patientResults{p}.patientID;
end

outputs = struct();
outputs.patientIDs = patientIDs;
outputs.patientResults = patientResults;
for p = 1:nPatients
    fieldName = matlab.lang.makeValidName(char(patientIDs(p)));
    outputs.(fieldName) = patientResults{p};
end
end
