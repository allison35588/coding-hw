function results = run_patient20_Q4_presentation()
%RUN_PATIENT20_Q4_PRESENTATION Compare three CO algorithms for s00020.

codeFolder = fileparts(mfilename('fullpath'));
addpath(fullfile(codeFolder, 'helpers'));
cpa_setup_paths();

% 5 = Liljestrand, 1 = MAP baseline, 7 = corrected impedance.
results = cpa_compare_algorithms_C2("s00020", 12, [5, 1, 7]);
end
