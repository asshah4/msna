%% Overview of project
% Purpose of code
    % First time looking at MSNA / signal processing on MATLAB
    % Need to learn how to extract data
    % Sample NVX_CKD file has channel 1 = MSNA, channel 2 = ECG
    % Sampling frequency is 1000 Hz, there is 10 minutes of data
    % Data already been amplified, band pass filtered (700-2000 Hz)
    % integrated at 100 ms
% Sections of code
    % Initial data definition
    % Initial data cleaning
    % Prepare MSNA data 

%% Initial data definitions
    %% Data definitions
    
    % Filename
    filename = 'NVX_CKD_018_073117.txt';
    patid = 'NVX_CKD_018_073117';

    % Sampling frequency in Hz
    fs = 1000;
    ts = 1/fs;

    % Create a table of the data, assuming column 5 is BP
    T = readtable(filename);
    T.Properties.VariableNames([1 2 3 5]) = {'timestamp' 'msna' 'ecg' 'bp'};
    
    %% Check to see if MSNA/ECG data is visible, first X seconds of data
    X = 10;
    S = stackedplot(T, {'msna','ecg', 'bp'})
    xlim([0, X*fs])
    
    X = 5;
    S = stackedplot(T, {'ecg', 'bp'})
    xlim([30000 50000])



    %% Identify MSNA data

    % For further processing, need a row vector of the signal
    vecMSNA = T.msna;
    N = length(vecMSNA);

    % Time vector is needed to correspond of real time duration of signal
    % Calculated by the steps from 0 to end of time using sampling Hz
    t = 0:ts:(N*ts)-ts;

    % Sampling frequency
    f0 = fs/N;

%% Inital data cleaning
    
    %% Interpolate missing data
    % Before processing, identify missing samples
    missing = isnan(vecMSNA);
    fprintf('Missing %d samples of %d\n', sum(missing), max(N))

    % Resample the data, simple interpolation
    vecMSNA = resample(vecMSNA, t); 

    %% Low pass filter to remove heart rate effect ~ 6 Hz

    % Normalized frequency
    fnorm = 6/(fs/2);
    d = designfilt('lowpassfir', 'FilterOrder', 100, ...
        'CutOffFrequency', fnorm);

    % However filters add a delay to the signal, which can be corrected for
    grpdelay(d, 10*fs, fs);
    D = mean(grpdelay(d));

    % Compensate for the delay
    lowfilterMSNA = filter(d, [vecMSNA; zeros(D,1)]); % append D zeroes
    lowfilterMSNA = lowfilterMSNA(D+1:end); % shift the data 

    % Plot out filtered data
    figure
    plot(t, vecMSNA, t, lowfilterMSNA)
    title('MSNA with Low pass filter of 6 Hz')
    legend('Original','Filtered')
    xlabel('Time (s)'), ylabel('Voltage (Units)'), grid on
    xlim([0 t(end)])

%% Prepare MSNA data for burst detection
    
    %% Median filter instead for detrending
    
    % Which MSNA
    x = lowfilterMSNA;
    
    % Apply median filter
    y = medfilt1(x, fs);
    z = x - y;
    
    % Visualize new data
    figure
    plot(t, x, t, y, t, z)
    legend('Original','Filtered', 'Subtracted')
    
    % Much more detrended, like this better
    medfiltMSNA = z;

    %% Simple approach to scale data to size
    
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

    % Scaled data to show which peaks are being identified
    figure
    plot(t, scaleMSNA)
    title('Scaled MSNA data')
    xlabel('Time (s)'), ylabel('Voltage (AU)')
    xlim([0 t(end)]), ylim([-10, 20])
    yline(baseline, '--r', 'LineWidth', 2)
    yline(noiseLevel, ':y', 'LineWidth', 2)
    yline(threshold, '-g', 'LineWidth', 2)

    %% Sliding/local scaling with threshold values of signal
    
    % Chosen signal
    temp = medfiltMSNA;
    
    % Scale mean/baseline data to 0
    temp = normalize(temp, 'center', 'mean');
    
    % Scale the max values to 100
    x = maxk(temp, 20);
    y = 100/mean(x); % scaling factor
    temp = temp * y;
    
    % Establish a moving mean value
    x = movmean(temp, fs*3);
    y = movstd(temp, fs*3);
    
    % Visualize this
    figure
    plot(t, temp, t, x, t, y)
    title('Scaled MSNA data')
    xlabel('Time (s)'), ylabel('Voltage (AU)')
    xlim([0 t(end)]), ylim([-20, 20])
    legend('MSNA','Baseline', 'Threshold')
    
    % Save scaled data
    normMSNA = temp;
    movingBaseline = x;
    movingThreshold = y;

%% Peak detection

    %% Load ECG data from prior
    
    % Simple vector
    vecECG = T.ecg;
    
    % Initialize HRV params from toolbox
    HRVparams = InitializeHRVparams(patid);
    HRVparams.sqi.LowQualityThreshold = 0.75; % Default: 0.9
    
    % Detect RR intervals
    [t_RR, rr, jqrs_ann, SQIjw, StartIdxSQIwindows_jw] = ... 
        ConvertRawDataToRRIntervals(vecECG, HRVparams, patid);
    
    % Get R peaks 
    tm = 0:1/fs:(length(vecECG)-1)/fs;
    r_peaks = jqrs(vecECG, HRVparams);
    
    % Plot
    figure
    plot(tm, vecECG)
    hold on;
    plot(r_peaks./fs, vecECG(r_peaks),'o')
    legend('ecg signal', 'detected R peaks')
    
    %% Create a way to limit peaks to within range
    
    % Need to use a vector that accounts for latency from RR beat
    lower_bounds = t_RR + 0.3;
    upper_bounds = t_RR + 1.3;
    
    % find peaks
    [pks, locs] = findpeaks(normMSNA, t, ...
        'MinPeakHeight', median(movingThreshold), ...
        'MinPeakDistance', 0.25, ...
        'MinPeakWidth', 0.25 ...
        );
    
    plot(t, normMSNA, locs, pks, 'o')
    xlim([0 t(end)]), ylim([-20, 50])
    
    % Limit to locations of interest and peaks
    new_locs = any(locs >= lower_bounds(:) & locs <= upper_bounds(:));
    
    % R edefine peaks with ECG limitations
    pks_corrected = pks(new_locs);
    locs_corrected = locs(new_locs);
    
    % Plot corrected peaks
    figure
    plot(t, normMSNA, locs_corrected, pks_corrected, 'o')
    xlim([0 t(end)]), ylim([-20, 50])
    legend('Signal', 'ECG-limited peaks')
    
    % Number of bursts over entire time series
    x = length(t)/fs/60;
    y = length(locs_corrected)/x;
    fprintf(['Burst rate (burst/minute) = %d bursts/minute, ' ...
        'total time = %d minutes\n'], round(y), round(x))
    
    % Bursts incidence = bursts / 100 beats
    x = length(locs_corrected)/length(rr)*100;
    fprintf(['Burst incidence (bursts/100beats) = %d bursts/100beats, ' ...
        'total beats = %d, average HR = %d\n'], ...
        round(x), length(rr), round(60/mean(rr)))