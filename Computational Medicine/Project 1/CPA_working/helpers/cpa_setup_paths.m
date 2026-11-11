function paths = cpa_setup_paths()
%CPA_SETUP_PATHS Build portable project paths relative to this code folder.
%
% No user-specific absolute path appears here. The same CPA_WORKING folder
% can therefore be moved to another computer. Supported layouts are:
%
%   project/CPA_working/                 (this code)
%   project/PhysioToolkit.../            (supplied toolkit; current layout)
%   project/s00020/                      (patient data; current layout)
%
% or a self-contained submission layout:
%
%   CPA_working/PhysioToolkit.../
%   CPA_working/data/s00020/
%
% Explicit ADDPATH calls are intentional. GENPATH is avoided because the
% downloaded project contains alternate estimateCO files with the same name.

helpersFolder = fileparts(mfilename('fullpath'));
codeFolder = fileparts(helpersFolder);
parentFolder = fileparts(codeFolder);

toolkitCandidates = { ...
    fullfile(codeFolder,'PhysioToolkitCardiacOutput_MatlabCode'), ...
    fullfile(parentFolder,'PhysioToolkitCardiacOutput_MatlabCode')};
toolkitFolder = '';
for k = 1:numel(toolkitCandidates)
    if isfolder(toolkitCandidates{k})
        toolkitFolder = toolkitCandidates{k};
        break
    end
end
assert(~isempty(toolkitFolder), [ ...
    'Could not find PhysioToolkitCardiacOutput_MatlabCode relative to ', ...
    'the CPA_working folder. Keep the supplied toolkit either inside ', ...
    'CPA_working or beside it.']);

% Put our functions first. Add supplied folders at the end so a duplicate
% helper cannot silently override the documented CPA implementation.
addpath(codeFolder,'-begin');
addpath(helpersFolder,'-begin');
addpath(parentFolder,'-end'); % supports the current estimateCO_v3 location
addpath(fullfile(toolkitFolder,'2analyze'),'-end');
addpath(fullfile(toolkitFolder,'3estimate'),'-end');
modifiedEstimateFolder = fullfile(toolkitFolder,'Modified estimateCO');
if isfolder(modifiedEstimateFolder)
    addpath(modifiedEstimateFolder,'-end');
end

paths = struct();
paths.helpersFolder = helpersFolder;
paths.codeFolder = codeFolder;
paths.parentFolder = parentFolder;
paths.workspaceFolder = parentFolder; % backward-compatible field name
paths.toolkitFolder = toolkitFolder;

% Patient folders may be under CPA_working/data, directly under
% CPA_working, or beside CPA_working. CPA_FIND_PATIENT_FILES checks these
% roots in order without depending on the MATLAB current folder.
paths.dataRoots = {fullfile(codeFolder,'data'),codeFolder,parentFolder};

% Generated files stay with the portable code folder.
paths.resultsFolder = fullfile(codeFolder,'results');
if ~isfolder(paths.resultsFolder)
    mkdir(paths.resultsFolder);
end
end
