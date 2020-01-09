function [timestamp, rawMSNA, ecg, bp, t, N] = ExtractRawSignal(fileName, fs)
% function [timestamp, rawMSNA, ecg, bp] = ExtractRawSignal(rawdata, fs)
%
% Simple function to extract HRV from tables, and create basic variables to describe the time series

% Read in data table
T = readtable(fileName);

% Extract variables
T.Properties.VariableNames([1 2 3 5]) = {'timestamp' 'msna' 'ecg' 'bp'};

% Time stamp
timestamp = T.timestamp;

% Raw MSNA data
rawMSNA = T.msna;

% ECG data
ecg = T.ecg;

% Continuous blood pressure
bp = T.bp;

% Overall sample length
N = length(rawMSNA);

% Time steps 
ts = 1/fs;

% Sampling frequency
f0 = fs/N;

% Time vector in steps by frequency
t = 0:ts:(N*ts)-ts;
