% EEG preprocessing + source localization (frequency & time domain)
% This script processes multiple EEG recordings:
%   (1) Preprocessing per epoch:
%       - DC removal per chanel
%       - Common-average reference (CAR) per time sample
%   (2) Frequency-domain FPVS response:
%       - FFT per epoch
%       - Complex-average across epochs (preserve phase-locking)
%       - Extract harmonics using: small FFT-bin neighborhood complex 
%         average (leakage robustness), neighboorhodd-bin baseline
%         substraction (FPVS-style local noise)
%       - Build complex harmonics stacks for SBR-R1 localization
%   (3) Time-domain FPVS response:
%       - Average epochs in time
%       - Band-pass filtered signals for (T)RAP MUSIC localization
%   (4) Group scalp maps:
%       - Summed baseline-corrected harmonics amplitudes for F1 (F/n)
%         and F0 (F)
%
% External functions:
%   - mySBRR1_AltOpt:sparse block source localization(frequency domain)
%   - allmusic:MUSIC/TRAP localization(time domain)
%   - scalpmap_mesh_spline: scalp interpolation/plotting
%
% Data format:
%   - Each recording with .lw5 file(header) and a .mat containing "data"
%   - Data format as [samples x channels x epochs]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
clear; close all; clc;

% Forward model and mesh
modeldir = ('C:\Users\vivianahc\Documents\MATLAB\EEGData-Understanding Human Individuation\AdultMNI152model\');
load([modeldir 'GMSurfaceMesh.mat']);
load([modeldir 'HeadVolumeMesh.mat']);
load([modeldir 'ScalpSurfaceMesh.mat']);

% Lead field (Biosemi128)
load ('LFBiosemi128.mat');
LF = LFfin; % N sensors x (3*Ndip)
elnames = cellstr(labelsn(:));   % converts string to cell array of char

% EEG data repository
datadir = 'C:\Users\vivianahc\Documents\MATLAB\EEGData-Understanding Human Individuation\Data-review-FPVS-oddball\';
lw5Files = dir(fullfile(datadir, '*.lw5'));
matFiles = dir(fullfile(datadir, '*.mat'));

% File pairing (.lw5 -> .mat)
matMap = containers.Map('KeyType','char','ValueType','char');
for k = 1:numel(matFiles)
    [~,stem] = fileparts(matFiles(k).name);
    matMap(stem) = fullfile(matFiles(k).folder, matFiles(k).name);
end

% FPVS design (Hz): general visual(F0) / face recognition(F1) responses
F0_Base = 6.0;
N_Odd = 5; % oddball ratio
F1_odd = F0_Base / N_Odd;

% Harmonic rules
N_Harm = 4; % keep first 4 harmonics
F1_Max_Hz = 5.0; % oddball harmonics <= 5 Hz

% Time-domain bands for TRAP (adjust if needed)
Band_F1 = [0.8  5.0]; % oddball response 
Band_F0 = [5.0 25.0]; % base response

% TRAP stopping (avoid expGOF=1.0; it is a bad target)
expGOF_trap = 0.95;

% FPVS spectral extraction setting
sigHalfWidthBins = 1; % +/- bin around harmonic for complex averaging
noiseOffsetBins = 3; % how far away noise bands start
noiseHalfWidthBins = 2; % +/- bins per noise sideband
useHalfSpectrum = true; % only use [0...Fs/2]

% Outputs
nSub = numel(lw5Files);
all_freq = cell(nSub,1); % freq localization 
all_time = cell(nSub,1); % time localization
all_Yf1h = cell(nSub,1); % complex F1 harm stacks(channels x H1)
all_Yf0h = cell(nSub,1); % complex F0 harm stacks(channels x H0)
all_f1_chan_bc = cell(nSub,1); % baseline-corrected harm-sum/ch
all_f0_chan_bc = cell(nSub,1);
lw5_info = struct(); 
lw5_list = cell(nSub,1);

% Main loop
for i = 1:nSub
    lw5path = fullfile(lw5Files(i).folder, lw5Files(i).name);
    [~,stem] = fileparts(lw5Files(i).name);
    
    % Match numeric data
    if ~isKey(matMap, stem)
        warning('No matching .mat for %s. Skipping.', lw5Files(i).name);
        continue;
    end
    matpath = matMap(stem);
    
    % Acquisition info
    [lw5_info.Fs, lw5_info.nSamples, lw5_info.nCh, lw5_info.nEpochs] = readEEGinfo(lw5path);
    lw5_list{i} = lw5_info;
    Fs = lw5_info.Fs;

    % Load EEG data
    X = loadMatData(matpath); % [samples x channels x epochs]
    X = permute(X, [2 3 1]); % [channels x samples x epochs]
    X = X(1:128,:,:);  % keep 128 channels (Biosemi128 order)
    
    % Preprocessing per epoch
    Xss = X;
    for e = 1:size(Xss,3)
        Xe = Xss(:,:,e);
        Xe = Xe - mean(Xe, 2); % remove DC per channel
        Xe = Xe - mean(Xe, 1); % CAR per time sample
        Xss(:,:,e) = Xe;
    end

    % FFT per epoch then complex-average across epochs
    Ns = size(Xss, 2);
    T = Ns/Fs;
    nCycles = floor(T * F1_odd);  
    Ns_new  = round(nCycles * Fs / F1_odd);
    Xss = Xss(:, 1:Ns_new, :);  
    Ns = Ns_new;
    Fr = Fs / Ns; 
    fvec = (0:Ns-1) * Fr;
    
    FFT_ep = fft(Xss, [], 2); % ch x freq x epoch
    FFT_avg = mean(FFT_ep, 3); % complex mean across epochs
    %EEGviewer(abs(FFT_avg), 1/Fr, elnames) % optional visualization

    % Harmonic lists (Hz)
    harm_f0 = F0_Base * (1:floor((Fs/2)/F0_Base));
    harm_f0 = harm_f0(1:N_Harm);

    harm_f1 = F1_odd * (1:ceil(F1_Max_Hz / F1_odd));
    harm_f1 = harm_f1(1:N_Harm);

    % Convert to 1-based FFT bins
    hz2bin  = @(hz) max(2, round(hz / Fr) + 1); % avoid DC bin=1
    Id_f0 = unique(arrayfun(hz2bin, harm_f0));
    Id_f1 = unique(arrayfun(hz2bin, harm_f1));
    
     % Extract harmonic stacks (complex + baseline correction)
    [Yf1_cmp, A1_bc] = extractHarmonicsBC(FFT_avg, Id_f1, ...
        sigHalfWidthBins, noiseOffsetBins, noiseHalfWidthBins, useHalfSpectrum);
    [Yf0_cmp, A0_bc] = extractHarmonicsBC(FFT_avg, Id_f0, ...
        sigHalfWidthBins, noiseOffsetBins, noiseHalfWidthBins, useHalfSpectrum);

    % Scalp maps: sum across harmonics of baseline-corrected amplitude
    % Clamp negatives values(optional but common)
    f1_chan_bc = sum(max(A1_bc, 0), 2); 
    f0_chan_bc = sum(max(A0_bc, 0), 2);
    
    % Store these for group maps
    all_f1_chan_bc{i} = f1_chan_bc;
    all_f0_chan_bc{i} = f0_chan_bc;

    % Frequency-domain source localization (SBR-R1)
    EEGloc_freq = struct();

    [EEGloc_freq.inddipsbr_f1, EEGloc_freq.orrsbr_f1, EEGloc_freq.ampsbr_f1, EEGloc_freq.GOFsbr_f1] = ...
        mySBRR1_AltOpt(Yf1_cmp, LF, 0.01);

    [EEGloc_freq.inddipsbr_f0, EEGloc_freq.orrsbr_f0, EEGloc_freq.ampsbr_f0, EEGloc_freq.GOFsbr_f0] = ...
        mySBRR1_AltOpt(Yf0_cmp, LF, 0.01);

    % Time-domain source localization, (T)RAP/MUSIC
    % Average over epochs in time first, then filter
    Xe_mean = mean(Xss, 3); % [ch x samples]
    EEG_F1 = bandpass_filter(Xe_mean, Fs, Band_F1); % face recognition 
    EEG_F0 = bandpass_filter(Xe_mean, Fs, Band_F0); % base band

    ndip_f1 = max(1, numel(EEGloc_freq.inddipsbr_f1));
    ndip_f0 = max(1, numel(EEGloc_freq.inddipsbr_f0));

    EEGloc_time = struct(); 
    
    [EEGloc_time.inddiptrap_f1, EEGloc_time.orrtrap_f1, EEGloc_time.amptrap_f1, EEGloc_time.GOFtrap_f1] = ...
    allmusic(EEG_F1, LF, expGOF_trap, 'trap', 1, ndip_f1);

    [EEGloc_time.inddiptrap_f0, EEGloc_time.orrtrap_f0, EEGloc_time.amptrap_f0, EEGloc_time.GOFtrap_f0] = ...
    allmusic(EEG_F0, LF, expGOF_trap, 'trap', 1, ndip_f0);

    % Store per subject
    all_Yf1h{i} = Yf1_cmp;
    all_Yf0h{i} = Yf0_cmp;
    all_freq{i} = EEGloc_freq;
    all_time{i} = EEGloc_time;

    fprintf('[%s] Fs=%.3f | Fr=%.6f\n | Ns=%d', stem, Fs, Fr, Ns);
end
    
% Group scalp maps 
valid_f1 = ~cellfun('isempty', all_f1_chan_bc);
valid_f0 = ~cellfun('isempty', all_f0_chan_bc);

if any(valid_f1)
    f1_grp = mean(cat(2, all_f1_chan_bc{valid_f1}), 2);
    figure;
    scalpmap_mesh_spline(f1_grp, nnewXYZ, ScalpSurfaceMesh.node, ScalpSurfaceMesh.face);
    axis off; colorbar;
    title('Group FI (F/n): summed baseline-corrected harmonic amplitude');
end

if any(valid_f0)
    f0_grp = mean(cat(2, all_f0_chan_bc{valid_f0}), 2);
    figure;
    scalpmap_mesh_spline(f0_grp, nnewXYZ, ScalpSurfaceMesh.node, ScalpSurfaceMesh.face);
    axis off; colorbar;
    title('Group Base (F): summed baseline-corrected harmonic amplitude');
end

% ------------------------- Helper functions ------------------------------
function [Fs, nSamples, nChannels, nEpochs] = readEEGinfo(filename)
    % Read LetsWave .lw5 header fields needed for sampling and dimensions
    S = load(filename, '-mat');
    hdr = S.header;
    nEpochs = double(hdr.datasize(1));
    nChannels = double(hdr.datasize(2));
    nSamples = double(hdr.datasize(end));
    Fs = 1 / double(hdr.xstep);
end

function X = loadMatData(matpath)
    % Load LetsWave .mat "data" and ensure [samples x channels x epochs]
    S = load(matpath,'data');
    X = squeeze(S.data);
    if ndims(X) ~= 3
        error('Expected data to be 3D [samples x channels x epochs] in %s', matpath);
    end
end

function y = bandpass_filter(x, Fs, Fband)
    % Zero-phase Butterworth band-pass on channels (rows) over time (cols)
    Wn = Fband/(Fs/2);
    [b,a] = butter(4, Wn, 'bandpass');
    y = filtfilt(b,a,double(x.')).';
end

function [Ycmp, A_bc] = extractHarmonicsBC(FFT_avg, Ik, sigHalfWidthBins, noiseOffsetBins, noiseHalfWidthBins, useHalfSpectrum)
    % Extract FPVS harmonics as:
    %   - complex mean over a small bin neighborhood around each harmonic (leakage robustness)
    %   - baseline correction using magnitude in two sidebands around that harmonic
    %
    % Inputs:
    %   FFT_avg: [channels x freqBins] complex spectrum (already epoch-averaged)
    %   Ik: target harmonic bin indices (1-based)
    %
    % Outputs:
    %   Ycmp: [channels x H] complex harmonic stack
    %   A_bc: [channels x H] baseline-corrected amplitude |Ycmp| - localNoise

    [M,N] = size(FFT_avg);
    H     = numel(Ik);

    Ycmp = zeros(M,H);
    A_bc = zeros(M,H);

    if useHalfSpectrum
        kMax = floor(N/2) + 1;  
    else
        kMax = N;
    end

    for h = 1:H
        k0 = Ik(h);

        % Signal neighborhood (complex average)
        sigBins = (k0-sigHalfWidthBins):(k0+sigHalfWidthBins);
        sigBins = sigBins(sigBins >= 2 & sigBins <= kMax); % avoid DC bin=1
        Yh = mean(FFT_avg(:,sigBins),2);
        Ycmp(:,h) = Yh;

        % Noise sidebands (magnitude average)
        leftBins = (k0-noiseOffsetBins-noiseHalfWidthBins):(k0-noiseOffsetBins+noiseHalfWidthBins);
        rightBins = (k0+noiseOffsetBins-noiseHalfWidthBins):(k0+noiseOffsetBins+noiseHalfWidthBins);
        noiseBins = [leftBins rightBins];
        noiseBins = noiseBins(noiseBins >= 2 & noiseBins <= kMax);
        noiseBins = setdiff(noiseBins, sigBins);  % against overlap

        if isempty(noiseBins)
            localNoise = 0;
        else
            localNoise = mean(abs(FFT_avg(:,noiseBins)),2);
        end

        A_bc(:,h) = abs(Yh) - localNoise;
    end
end
