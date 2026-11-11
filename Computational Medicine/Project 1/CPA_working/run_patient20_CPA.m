% RUN_PATIENT20_CPA  Single script for submitted Problems 1 and 3.
%
% TEAM INFORMATION (replace before submission)
% Team #: TODO
% Team members: TODO
%
% This is the ONE file instructors should press "Run" on for the Coding Progress
% Assignment. It generates only the required patient-s00020 results:
%   Problem 1: 20-pulse feature plots at 10 and 11 hours.
%   Problem 3: first-12-hour Liljestrand CO with first-pair C2 calibration.
%
% The helper functions locate the project folders automatically, so this
% script does not depend on MATLAB's current working folder.

codeFolder = fileparts(mfilename('fullpath'));
addpath(fullfile(codeFolder, 'helpers'));
cpa_setup_paths();

fprintf('\n=== Problem 1: patient s00020 ABP feature traces ===\n');
CPA_outputs.problem1 = cpa_plot_abp_features("s00020",[10 11]);

fprintf('\n=== Problem 3: patient s00020 Liljestrand C2 analysis ===\n');
CPA_outputs.problem3 = cpa_analyze_co_C2("s00020",12,5);

fprintf('\nPatient s00020 CPA results are complete.\n');
