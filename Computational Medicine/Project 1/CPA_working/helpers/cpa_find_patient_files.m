function files = cpa_find_patient_files(patientID,paths)
%CPA_FIND_PATIENT_FILES Locate ABP and trend files under relative data roots.
%
% FILES = CPA_FIND_PATIENT_FILES(PATIENTID,PATHS) accepts values such as
% 20, "20", or "s00020". It never searches from the Windows Desktop or
% MATLAB's current folder; all candidate roots come from CPA_SETUP_PATHS.

patientID = normalizePatientID(patientID);
abpFiles = strings(0,1);
trendFiles = strings(0,1);

for k = 1:numel(paths.dataRoots)
    patientFolder = fullfile(paths.dataRoots{k},char(patientID));
    if ~isfolder(patientFolder)
        continue
    end

    abpMatch = dir(fullfile(patientFolder,'**', ...
        [char(patientID) '*_ABP.txt']));
    for j = 1:numel(abpMatch)
        abpFiles(end+1,1) = string(fullfile( ...
            abpMatch(j).folder,abpMatch(j).name)); %#ok<AGROW>
    end

    trendMatch = dir(fullfile(patientFolder,'**', ...
        [char(patientID) '*n.txt']));
    for j = 1:numel(trendMatch)
        trendFiles(end+1,1) = string(fullfile( ...
            trendMatch(j).folder,trendMatch(j).name)); %#ok<AGROW>
    end
end

abpFiles = unique(abpFiles,'stable');
trendFiles = unique(trendFiles,'stable');
assert(numel(abpFiles)==1, ...
    'Expected exactly one ABP file for %s; found %d.', ...
    patientID,numel(abpFiles));

files = struct();
files.patientID = patientID;
files.abpFile = char(abpFiles(1));
files.trendFile = '';
if numel(trendFiles)==1
    files.trendFile = char(trendFiles(1));
elseif numel(trendFiles)>1
    error('Expected at most one trend/TCO file for %s; found %d.', ...
        patientID,numel(trendFiles));
end
end

function patientID = normalizePatientID(patientID)
if isnumeric(patientID)
    patientNumber = patientID;
else
    text = erase(lower(strtrim(string(patientID))),"s");
    patientNumber = str2double(text);
end
assert(isscalar(patientNumber) && isfinite(patientNumber) && ...
    patientNumber>=0 && patientNumber==fix(patientNumber), ...
    'patientID must look like 20 or s00020.');
patientID = string(sprintf('s%05d',patientNumber));
end
