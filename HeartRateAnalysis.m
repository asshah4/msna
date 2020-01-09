function [rr] = HeartRateAnalysis(fileName, ecg, bp, t, fs);
%
%
% Looking at both ECG and BP curves to determine RR intervals
% This RR series will not be adequate for HRV analysis

%% Extract RR intervals from ECG

% Initialize HRV params from toolbox
HRVparams = InitializeHRVparams(fileName);

% Detect RR intervals
% rr is in seconds
% t_RR is time stamp
[t_RR, rr, jqrs_ann, SQIjw, StartIdxSQIwindows_jw] = ... 
    ConvertRawDataToRRIntervals(ecg, HRVparams, fileName);

% Get R peaks 
tm = 0:1/fs:(length(ecg)-1)/fs;
r_peaks = jqrs(ecg, HRVparams);


%% Analyze BP data

% Add low pass filter to remove heart rate effects ~ 6Hz
fnorm = 6/(fs/2);
d = designfilt('lowpassfir', 'FilterOrder', 100, ...
'CutOffFrequency', fnorm);

% However filters add a delay to the signal, which can be corrected for
% grpdelay(d, 10*fs, fs); % used to visualize the delay
D = mean(grpdelay(d));

% Add filter
tmp = filter(d, [bp; zeros(D,1)]); % append D zeroes
bp = tmp(D+1:end); % shift the data 

% Identify all the peaks to generate BP intervals
% Rough draft of peaks
[pks, locs] = findpeaks(bp, tm, ...
    'MinPeakDistance', 0.3, ... % distance between peaks
    'MinPeakHeight', median(bp) + iqr(bp)/2, ... % SBP higher than mean only
    'MinPeakProminence', mean(bp)/2 ...
	);


%% Compare ECG and BP data to find the correct RR intervals

% Don't need RR intervals, just the time points
% SBP is delayed by ~ 0.6 seconds after R peak
% Also, BP has "calibration" periods that have no peaks
% BP = locs, ECG = t_RR


% Create a vector that sees if every ECG point has a BP point
lower_bound = t_RR + 0.3;
upper_bound = t_RR + 1.3;

%x = any(locs >= lower_bound(:) & locs <= upper_bound(:));
%y = t_RR(x);
%length(y);

% Check a stacked plot of ECG and BP
%ts = 1:10000;
%T = array2table([tm(ts)', bp(ts), ecg(ts)], ...
%    'VariableNames', {'Time', 'BP', 'ECG'});
%stackedplot(T, 'XVariable', 'Time');

%% Add missing data to both BP and ECG signal

% Find positions where there are time gaps
int = diff(locs); % intervals
idx = find(int > 2); % positions of gaps TOTAL
sp = round(int(idx)/median(int)) - 1; % bordered by peaks
newInt = int;
newLocs = locs;

% Insert this into interval spacing
for i = 1:length(idx)
    
    for j = 1:sp(i)
        k = j - 1;
        newLocs = [newLocs(1:idx(i)+k) ...
            newLocs(idx(i)) + j*nanmean(int) ...
            newLocs((idx(i)+j):end)];
    end
  
    % Shift position of next NaN by inserted amount
    idx = idx + sp(i);
end

% ECG rendition

% Find positions where there are time gaps
int = diff(t_RR); % intervals
idx = find(int > 2); % positions of gaps TOTAL
sp = round(int(idx)/median(int)) - 1; % bordered by peaks
newRR = int;
newTRR = t_RR;

% Insert this into interval spacing
for i = 1:length(idx)
    
    for j = 1:sp(i)
        k = j - 1;
        newTRR = [newTRR(1:idx(i)+k) ...
            newTRR(idx(i)) + j*nanmean(int) ...
            newTRR((idx(i)+j):end)];
    end
  
    % Shift position of next NaN by inserted amount
    idx = idx + sp(i);
end

%% Return a variable

% Use bp this time, and remove the 0.6 sec delay
% Get rid of any beats that occur before the ECG signal
x = newLocs - 0.6;
y = x(find(x >= 0));
rr = y;

rr = t_RR;


end
