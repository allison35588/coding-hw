# MATLAB Code Usage

Run the main files from the `CPA_working` folder. Main files call the helper
functions in `helpers` automatically.

## Main code

### `run_patient20_CPA.m`

Runs the submitted patient-s00020 analyses:

- Problem 1: 20 ABP pulses beginning at 10 and 11 hours.
- Problem 3: first 12 hours of estimator 5 with first-pair C2 calibration.

```matlab
run("run_patient20_CPA.m")
```

Results are placed in the workspace variable:

```matlab
CPA_outputs.problem1
CPA_outputs.problem3
```

### `run_patient20_Q4_presentation.m`

Compares estimators 5, 1, and 7 over the first 12 hours of s00020.

```matlab
q4Patient20 = run_patient20_Q4_presentation();
```

### `run_patient214_presentation.m`

Generates the 10-hour ABP feature trace and the Q4 estimator comparison for
s00214.

```matlab
patient214 = run_patient214_presentation();
```

Returned fields:

```matlab
patient214.problem2
patient214.problem4
```

### `run_Q5_parlikar_comparison.m`

Runs Parlikar estimator 14 and the selected reference estimators for any list
of patients.

```matlab
q5 = run_Q5_parlikar_comparison(patientIDs,endHours,estimatorIDs);
```

Examples:

```matlab
q5 = run_Q5_parlikar_comparison([20 214],12,[5 1 7 14]);

q5 = run_Q5_parlikar_comparison( ...
    [20 214 321],[12 12 10],[5 1 7 14]);
```

`endHours` may be one value for all patients or one value per patient.

Important returned fields:

```matlab
q5.patientResults
q5.s00020
q5.s00214
q5.combinedMetrics
q5.macroMetrics
q5.bestByPatient
```

Q5 also saves its complete output for later Q6 use:

```matlab
results/q5_outputs_for_q6.mat
```

### `run_Q6_tpr.m`

Calculates and plots one-minute median Parlikar TPR using the tau and C2
values already produced by Q5.

Same MATLAB session:

```matlab
q5 = run_Q5_parlikar_comparison([20 214],12,[5 1 7 14]);
q6 = run_Q6_tpr(q5);
```

Later MATLAB session:

```matlab
q6 = run_Q6_tpr("results/q5_outputs_for_q6.mat");
```

If the default MAT file exists, this also works:

```matlab
q6 = run_Q6_tpr();
```

Important returned fields:

```matlab
q6.patientResults
q6.s00020.TPR
q6.s00020.TPRcheck
q6.s00020.tauSeconds
q6.s00020.tprTrend
q6.s00020.minuteTimeHours
```

## Helper functions

### `cpa_setup_paths`

Adds the required code folders and returns the paths used by the other
functions. Normally called automatically by the main code.

```matlab
paths = cpa_setup_paths();
```

Useful returned fields:

```matlab
paths.codeFolder
paths.helpersFolder
paths.toolkitFolder
paths.dataRoots
paths.resultsFolder
```

### `cpa_find_patient_files`

Finds one patient's ABP and trend/TCO files.

```matlab
paths = cpa_setup_paths();
files = cpa_find_patient_files("s00020",paths);
```

Returned fields:

```matlab
files.patientID
files.abpFile
files.trendFile
```

Patient IDs may be entered as:

```matlab
20
"20"
"s00020"
```

### `cpa_plot_abp_features`

Plots 20 ABP pulses and the onset/end-systole features for any patient and
one or more requested hours.

```matlab
results = cpa_plot_abp_features(patientID,targetHours);
```

Examples:

```matlab
results = cpa_plot_abp_features(20,[10 11]);
results = cpa_plot_abp_features("s00214",10);
```

### `cpa_analyze_co_C2`

Generates a Figure-4-style continuous CO, PP, MAP, and HR analysis for one
patient and one estimator using first-pair C2 calibration.

```matlab
results = cpa_analyze_co_C2(patientID,endHour,estimatorID);
```

Example:

```matlab
results = cpa_analyze_co_C2(20,12,5);
```

Important returned fields:

```matlab
results.C2
results.coTrend
results.ppTrend
results.mapTrend
results.hrTrend
results.tcoMeasured
results.pairedCalibratedCO
```

### `cpa_compare_algorithms_C2`

Compares multiple estimators for one patient using identical preprocessing,
first-pair C2 calibration, and common post-calibration TCO validation points.

```matlab
results = cpa_compare_algorithms_C2( ...
    patientID,endHour,estimatorIDs,outputTag);
```

Examples:

```matlab
q4 = cpa_compare_algorithms_C2(20,12,[5 1 7],"q4");
q5 = cpa_compare_algorithms_C2(20,12,[5 1 7 14],"q5");
```

Important returned fields:

```matlab
results.metricsTable
results.minuteTrends
results.pairedPredictions
results.commonValidationRows
results.parlikar
```

Estimator IDs used in this project:

```text
1   MAP baseline
5   Liljestrand
7   Corrected impedance
14  Parlikar
```

### `cpa_est14_parlikar`

Calculates uncalibrated beat-level Parlikar CO. This helper does not perform
TCO pairing or C2 calibration.

```matlab
[coRelative,timeMinutes,goodFeatures,diagnostics] = ...
    cpa_est14_parlikar(ABP,t_on,feat,beatq,Fs,windowBeats);
```

Typical call:

```matlab
[coRelative,timeMinutes,goodFeatures,diagnostics] = ...
    cpa_est14_parlikar(ABP,t_on,feat,beatq,125,361);
```

Useful diagnostic fields:

```matlab
diagnostics.beta
diagnostics.tauSeconds
diagnostics.coRawAllBeats
diagnostics.timeMinutesAllBeats
diagnostics.deltaPmmHg
diagnostics.beatPeriodSeconds
diagnostics.MAPmmHg
diagnostics.validOutputMask
```

### `cpa_estimate_tpr_from_q5`

Calculates TPR for one patient from that patient's Q5 result. It does not
rerun ABP processing or the Parlikar estimator.

```matlab
tpr = cpa_estimate_tpr_from_q5(q5PatientResult);
```

Example:

```matlab
q5 = run_Q5_parlikar_comparison(20,12,[5 1 7 14]);
tpr = cpa_estimate_tpr_from_q5(q5.s00020);
```

Important returned fields:

```matlab
tpr.compliance_mL_per_mmHg
tpr.tauSeconds
tpr.TPR
tpr.TPRcheck
tpr.tprTrend
tpr.minuteTimeHours
tpr.validFraction
```
