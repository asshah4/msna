% Peformance Analysis of MSNA compared to Manual Annotation
%
% 	Variables that can be modified
% 		filterSize = # of samples for medfilt, default 100 samples
% 		scalingSample = # of peaks to account for scaling to max size of ~100
% 		baselineFunction = mean versus median for baseline
% 		noiseLevelFunction = how to calculate noise (e.g. IQR, SD)
% 		thresholdFunction = calculate types of threshold	
% 		qualitySQI = value for RR signal quality, default 0.9
% 		lowerLatency = lower bounds of RR latency
% 		upperLatency = upper bounds of RR latency
% 		peakDist = value in seconds
% 		peakWidth = minimum of the width of the peaks	
% 		peakAmplitude = choose between static or dynamic thresholds
% 		peakSlope = choice of what slope to use
%
% Will require a series of functions to return data in a certain manner 
 	

%% Pseudocode

% Read all the data into a vector
% Identify the global variables, including sampling frequency
% Call peak detection algorithm, which labels all peaks
% Identify all ECG and BP signal to generate a approximate RR series
% Clean / filter the peak data

%% Global variables and settings

% Clear workspace
clear; clc; close all;

% Add to path, need to be in Signal Processing / MSNA folderA
addpath(genpath(pwd));

% Folder holding data
raw_folder = [pwd filesep 'raw_data'];

% Target folder for patient data (required by HRV tooblox
proc_folder = [pwd filesep 'proc_data'];

% Identify all MSNA data in folder
files = dir(fullfile(raw_folder, '*.txt'));
patients = regexprep({files.name}, '.txt', '');
numsub = length(patients);

% Sampling frequency
fs = 1000;


%% Write this into file for each patient

% move into raw_data folder
cd(raw_folder);

parfor i = 1:numsub
	names(i) = patients(i); % name of patient being analyzed
	[timestamp, rawMSNA, ecg, bp, t, N] = ...
		ExtractRawSignal(names{i}, fs); % Extract raw data
	[rr] = HeartRateAnalysis(names{i}, ecg, bp, t, fs); % RR ints
	
	% Temporary latency value
	latency = 1.3;
	[procMSNA, bandMSNA, pks, locs] = ...
		FindAllPeaks(rawMSNA, t, fs, rr, latency); % Peak analysis
	
	% Calculate summary data
	x = length(t)/fs/60;
	y = length(locs)/x;
	z = length(locs)/length(rr)*100;
	
	% Data output
	freq(i) = round(y);
	incidence(i) = round(z);
end

% Write the data into an excel files
patid = names';
freq = freq';
incidence = incidence';
T = table(patid, freq, incidence);
writetable(T, [pwd filesep 'matlab_burst_algorithm.csv']);

%% Troubleshooting

fileName = 'NVX_CKD_020_010918EOS';

% Extract raw data
[timestamp, rawMSNA, ecg, bp, t, N] = ...
	ExtractRawSignal(fileName, fs); 

% RR intervals
[rr] = HeartRateAnalysis(fileName, ecg, bp, t, fs);

% Peak analysis
% Define latency in seconds
latency = 1.3;
[procMSNA, bandMSNA, pks, locs] = ...
	FindAllPeaks(rawMSNA, t, fs, rr, latency); 


% Visualize the plot, amplify the ECG for visualization
plot(t, bandMSNA, ...
	locs, pks, 'o', ...
	t, ecg*-100)
ylim([-100 100]);


% Calculate summary data
x = length(t)/fs/60;
y = length(locs)/x;
z = length(locs)/length(rr)*100;

% Data output
freq = round(y);
incidence = round(z);
freq
incidence
