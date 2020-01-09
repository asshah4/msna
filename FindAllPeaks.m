function [procMSNA, bandMSNA, pks, locs] = FindAllPeaks(rawMSNA, t, fs, rr)
%
% 	OVERVIEW
% 		Load raw MSNA data, time series, and sampling frequency
% 		Will identify all the peaks, with data clean up
%
%   INPUT:
%       rawMSNA 	: Raw vector signal, same length as time series
%       t 		    : Time series in seconds
%       fs 			: Frequency of sampling rate
%
%   OUTPUT:
%       pks 		: Raw vector signal, same length as time series
%       locs 		: Time series in seconds

%% Cleaning of MSNA data
N = length(rawMSNA);

% Before processing, identify missing samples
missing = isnan(rawMSNA);
fprintf('Missing %d samples of %d\n', sum(missing), max(N))

% Resample the data, simple interpolation
rawMSNA = resample(rawMSNA, t); 

% Low pass filter to adjust for maximum heart rate
fnorm = 6/(fs/2);
d = designfilt('lowpassfir', 'FilterOrder', 100, ...
	'CutOffFrequency', fnorm);

% However filters add a delay to the signal, which can be corrected for
D = mean(grpdelay(d));

% Compensate for the delay
lowfilterMSNA = filter(d, [rawMSNA; zeros(D,1)]); % append D zeroes
lowfilterMSNA = lowfilterMSNA(D+1:end); % shift the data 

% Detrend using median filter
x = lowfilterMSNA;

% Apply median filter
y = medfilt1(x, 200);
z = x - y;

% Much more detrended
medfiltMSNA = z;

%% Scaling the data and suggestive signal/noise patterns

% Which MSNA data to use for scaling?
temp = medfiltMSNA;

% Make peak of 100
x = maxk(temp, 100);
y = 100/mean(x); % scaling factor
scaleMSNA = temp * y;

% Find baseline signal voltage during non-bursting segment (NBS)
baseline = median(scaleMSNA);

% Noise...use IQR/2 of NBS (arbitrary choice)
noiseLevel = iqr(scaleMSNA)/2;

% Threshold for finding peaks 
% literature 3:1 ratio for signal to noise
threshold = baseline + noiseLevel*3;

% Establish a moving mean value
movingBaseline = movmean(scaleMSNA, fs*3);
movingThreshold = movstd(scaleMSNA, fs*3);

% Final form of the MSNA
procMSNA = scaleMSNA;

%% Peak detection

% Rough draft of peaks
[pks, locs_pks] = findpeaks(procMSNA, t, ...
    'MinPeakDistance', .2, ...
    'MinPeakProminence', threshold/2 ... % temp: divided by 2
    );

%% Remove peaks that are above a moving SD

x = any(t == locs_pks(:));
y = movingThreshold(x)';
locs_pks = locs_pks(pks >= y(:));
pks = pks(pks >= y(:));

%% Option: scaled and median filtered MSNA

pks = pks;
locs = locs_pks;

% Need to use a vector that accounts for latency from RR beat
% Centered around 1.25 s (Hamner and Taylor 2001)
latency = 1.25;
maxerr = 0.300;
lower_bounds = rr + latency - maxerr;
upper_bounds = rr + latency + maxerr;


% Limit to locations of interest and peaks
new_locs = any(locs >= lower_bounds(:) & locs <= upper_bounds(:));

% Redefine peaks with ECG limitations
pks_c = pks(new_locs);
locs_c = locs(new_locs);
	
%% Bandpass filter

% Do a 0.5 to 4k hz bandpass
bandMSNA = bandpass(rawMSNA, [0.5 4], fs);

% Cut off values
movingBaseline = movmean(bandMSNA, fs*3);
movingThreshold = movstd(bandMSNA, fs*3);
interval = min(diff(rr));
baseline = mean(bandMSNA);
noiseLevel = iqr(bandMSNA);
threshold = baseline + noiseLevel;

% Find baseline signal voltage during non-bursting segment (NBS)
baseline = median(scaleMSNA);

% Noise...use IQR/2 of NBS (arbitrary choice)
noiseLevel = iqr(scaleMSNA)/2;

% Threshold for finding peaks 
% literature 3:1 ratio for signal to no

% Locate peaks
[pks locs] = findpeaks(bandMSNA, t, ...
	'MinPeakDistance', interval, ...
	'MinPeakHeight', threshold/2 ...
	);


x = any(t == locs(:));
y = movingThreshold(x)';
locs = locs(pks >= y(:));
pks = pks(pks >= y(:));

% Find surrounding minima
[vals spots] = findpeaks(bandMSNA.*-1, t, ...
	'MinPeakDistance', 0.25 ...
	);

vals = vals.*-1;
spots_left = locs - 0.4;
spots_right = locs + 0.4;

x = any(spots >= spots_left(:) & spots <= spots_right(:));
spots = spots(x);
vals = vals(x);

% Visualize the plot PRN
plot(t, bandMSNA, ...
	locs, pks, 'o')

%% Option: bandpass filtered data
% Need to use a vector that accounts for latency from RR beat
% Centered around 1.25 s (Hamner and Taylor 2001)
latency = 1.25;
maxerr = 0.500;
lower_bounds = rr + latency - maxerr;
upper_bounds = rr + latency + maxerr;

% Limit to locations of interest and peaks
new_locs = any(locs >= lower_bounds(:) & locs <= upper_bounds(:));

% Redefine peaks with ECG limitations
pks = pks(new_locs);
locs = locs(new_locs);

end
