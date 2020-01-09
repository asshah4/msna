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

%% Limit peaks to match for heart rate

% Plot corrected peaks
plot(t, bandMSNA, locs, pks, 'o')
xlim([0 t(end)]), ylim([-50, 50])
legend('Signal', 'ECG-limited peaks')

% Name of patient
fprintf(['Name: %s\n'], fileName)

% Number of bursts over entire time series
x = length(t)/fs/60;
y = length(locs)/x;
fprintf(['Bandpassed burst rate (burst/minute) \n= %d bursts/minute, ' ...
	'\ntotal time is %d minutes\n'], round(y), round(x))

% Bursts incidence = bursts / 100 beats
x = length(locs)/length(rr)*100;
fprintf(['Bandpassed burst incidence (bursts/100beats) \n= %d bursts/100beats, ' ...
	'\ntotal beats is %d, average HR = %d\n'], ...
	round(x), length(rr), round(60/mean(diff(rr))))

%% Write this into file for each patient

for i = 1:numsub
	names(i) = patients(i); % name of patient being analyzed
	[timestamp, rawMSNA, ecg, bp, t, N] = ...
		ExtractRawSignal(names{i}, fs); % Extract raw data
	[rr] = HeartRateAnalysis(names{i}, ecg, bp, t, fs); % RR ints
	[procMSNA, bandMSNA, pks, locs] = ...
		FindAllPeaks(rawMSNA, t, fs, rr); % Peak analysis
	
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
writetable(T, [pwd filesep 'burst_algorithm.csv']);

%% Troubleshooting

fileName = 'NVX_CKD_018_073117';

% Extract raw data
[timestamp, rawMSNA, ecg, bp, t, N] = ...
	ExtractRawSignal(fileName, fs); 

% RR intervals
[rr] = HeartRateAnalysis(fileName, ecg, bp, t, fs);

% Peak analysis
% Define latency in seconds
latency = 1.25
[procMSNA, bandMSNA, pks, locs] = ...
	FindAllPeaks(rawMSNA, t, fs, rr, latency); 

% Calculate summary data
x = length(t)/fs/60;
y = length(locs)/x;
z = length(locs)/length(rr)*100;

% Data output
freq = round(y);
incidence = round(z);
freq
incidence
