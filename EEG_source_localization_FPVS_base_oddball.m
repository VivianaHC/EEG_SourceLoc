% Simple EEG preprocessing + source localization
%
% For each recording:
% 1) load EEG
% 2) preprocess each epoch:
%      - remove DC per channel
%      - common average reference
% 3) compute epoch-averaged FFT
% 4) extract harmonics for F0 (base response) and F1 (face individuation)
% 5) compute baseline-corrected harmonic amplitudes for scalp maps
% 6) run frequency-domain localization (SBR-R1)
% 7) run time-domain localization (TRAP/MUSIC)

clear; close all; clc;

% Paths for forward model and EEG data 
modeldir = 'C:\Users\vivianahc\Documents\MATLAB\EEGData-Understanding Human Individuation\AdultMNI152model\';
datadir  = 'C:\Users\vivianahc\Documents\MATLAB\EEGData-Understanding Human Individuation\Data-review-FPVS-oddball\';

load([modeldir 'GMSurfaceMesh.mat']);
load([modeldir 'HeadVolumeMesh.mat']);
load([modeldir 'ScalpSurfaceMesh.mat']);
load('LFBiosemi128.mat'); % Lead-field, ch-labels

LF = LFfin;
elnames = cellstr(labelsn(:)); 

% FPVS parameters 
F0 = 6.0; % face presentation frequency
Nodd = 5; % oddball every 5 stimuli
F1 = F0 / Nodd; % face individuation frequency

Nharm = 4; % keep first 4 harmonics
F1max = 5.0; % oddball harmonics only up to 5 Hz

BandF1 = [0.5 5.0];
BandF0 = [5.0 25.0];
expGOF = 0.95;

% File list
lw5Files = dir(fullfile(datadir, '*.lw5'));
matFiles = dir(fullfile(datadir, '*.mat'));

matMap = containers.Map;
for k = 1:numel(matFiles)
    [~, name] = fileparts(matFiles(k).name);
    matMap(name) = fullfile(matFiles(k).folder, matFiles(k).name);
end

nSub = numel(lw5Files);

all_freq = cell(nSub,1);
all_time = cell(nSub,1);
all_f1_map = cell(nSub,1);
all_f0_map = cell(nSub,1);

for i = 1:nSub
    [~, stem] = fileparts(lw5Files(i).name);
    lw5path = fullfile(lw5Files(i).folder, lw5Files(i).name);

    if ~isKey(matMap, stem)
        warning('No matching .mat file for %s', stem);
        continue;
    end

    matpath = matMap(stem);

    % Read header and data
    EEGinfo = load(lw5path, '-mat');
    hdr = EEGinfo.header;
    Fs = 1 / double(hdr.xstep);

    S = load(matpath, 'data');
    X = squeeze(S.data); % expected [epochs x channels x samples]
    X = permute(X, [2 3 1]); % -> [channels x samples x epochs]
    X = X(1:128,:,:); % keep Biosemi128

    % Preprocess
    for e = 1:size(X,3)
        Xe = X(:,:,e);
        Xe = Xe - mean(Xe, 2); % remove DC per channel
        Xe = Xe - mean(Xe, 1); % common average reference
        X(:,:,e) = Xe;
    end

    % Frequency domain
    Ns = size(X,2);
    FFTep = fft(X, [], 2);
    FFTavg = mean(FFTep, 3);
    Fr = Fs / Ns;
    %EEGviewer(abs(FFTavg), 1/Fr, elnames) % optional visualization
    % Harmonic bins
    idxF0 = harmonic_bins(F0, Nharm, Fs, Fr);
    idxF1 = harmonic_bins(F1, Nharm, Fs, Fr, F1max);

    % Amplitudes baseline correction for scalp maps
    A1 = harmonic_amplitude_bc(FFTavg, idxF1);
    A0 = harmonic_amplitude_bc(FFTavg, idxF0);

    all_f1_map{i} = sum(max(A1,0), 2);
    all_f0_map{i} = sum(max(A0,0), 2);

    % Frequency-domain localization
    Yf1 = FFTavg(:, idxF1);
    Yf0 = FFTavg(:, idxF0);

    freq = struct();
    [freq.inddipsbr_f1, freq.orrsbr_f1, freq.ampsbr_f1, freq.GOFsbr_f1] = ...
        mySBRR1_AltOpt(Yf1, LF, 0.01);

    [freq.inddipsbr_f0, freq.orrsbr_f0, freq.ampsbr_f0, freq.GOFsbr_f0] = ...
        mySBRR1_AltOpt(Yf0, LF, 0.01);

    % Time-domain localization
    Xe = mean(X, 3); % epoch average in time
    EEG_F1 = bandpass_eeg(Xe, Fs, BandF1);
    EEG_F0 = bandpass_eeg(Xe, Fs, BandF0);

    ndip_f1 = max(1, numel(freq.inddipsbr_f1));
    ndip_f0 = max(1, numel(freq.inddipsbr_f0));

    time = struct();
    [time.inddiptrap_f1, time.orrtrap_f1, time.amptrap_f1, time.GOFtrap_f1] = ...
        allmusic(EEG_F1, LF, expGOF, 'trap', 1, ndip_f1);

    [time.inddiptrap_f0, time.orrtrap_f0, time.amptrap_f0, time.GOFtrap_f0] = ...
        allmusic(EEG_F0, LF, expGOF, 'trap', 1, ndip_f0);

    all_freq{i} = freq;
    all_time{i} = time;

    fprintf('[%s] Fs = %.3f Hz | Fr = %.6f Hz\n', stem, Fs, Fr);
end

% ------------------ Group scalp maps -------------------------------------
plot_group_scalpmap(all_f1_map, nnewXYZ, ScalpSurfaceMesh, ...
    'Group F1 (F/n): summed baseline-corrected harmonic amplitude');

plot_group_scalpmap(all_f0_map, nnewXYZ, ScalpSurfaceMesh, ...
    'Group F0 (F): summed baseline-corrected harmonic amplitude');

% --------------------- Support functions ---------------------------------
function idx = harmonic_bins(F, Nharm, Fs, Fr, Fmax)
    if nargin < 5
        harm = F * (1:floor((Fs/2)/F));
    else
        harm = F * (1:ceil(Fmax/F));
        harm = harm(harm <= Fmax);
    end
    harm = harm(1:min(Nharm, numel(harm)));
    idx = round(harm / Fr) + 1;
    idx(idx < 2) = 2; 
    idx = unique(idx);
end

function A = harmonic_amplitude_bc(FFTavg, idx)
    % Local baseline correction:
    % amplitude at harmonic minus mean of nearby bins
    nCh = size(FFTavg,1);
    nH = numel(idx);
    nF = size(FFTavg,2);
    kMax = floor(nF/2) + 1;

    A = zeros(nCh, nH);

    for h = 1:nH
        k = idx(h);
        signalAmp = abs(FFTavg(:,k));

        noiseBins = [k-5:k-3, k+3:k+5];
        noiseBins = noiseBins(noiseBins >= 2 & noiseBins <= kMax);

        if isempty(noiseBins)
            baseline = zeros(nCh,1);
        else
            baseline = mean(abs(FFTavg(:,noiseBins)), 2);
        end
        A(:,h) = signalAmp - baseline;
    end
end

function Y = bandpass_eeg(X, Fs, band)
    Wn = band / (Fs/2);
    [b,a] = butter(4, Wn, 'bandpass');
    Y = filtfilt(b, a, double(X.')).'; 
end

function plot_group_scalpmap(allMaps, nnewXYZ, ScalpSurfaceMesh, ttl)
    valid = ~cellfun('isempty', allMaps);
    if ~any(valid)
        return;
    end

    grp = mean(cat(2, allMaps{valid}), 2);

    figure;
    scalpmap_mesh_spline(grp, nnewXYZ, ScalpSurfaceMesh.node, ScalpSurfaceMesh.face);
    axis off;
    colorbar;
    title(ttl);
end
