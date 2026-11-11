function outputs = run_patient214_presentation()
%RUN_PATIENT214_PRESENTATION Generate Problems 2 and 4 for s00214.

codeFolder = fileparts(mfilename('fullpath'));
addpath(fullfile(codeFolder, 'helpers'));
cpa_setup_paths();

fprintf('\n=== Problem 2: patient s00214 trace at 10 hours ===\n');
outputs.problem2 = cpa_plot_abp_features("s00214",10);

fprintf('\n=== Problem 4: patient s00214 algorithm comparison ===\n');
outputs.problem4 = cpa_compare_algorithms_C2( ...
    "s00214", 12, [5, 1, 7]);
end
