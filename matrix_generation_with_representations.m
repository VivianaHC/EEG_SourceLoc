clear; clc;

load('C:\Users\vivianahc\Downloads\all_freq.mat'); % all_freq
load('C:\Users\vivianahc\Downloads\all_time.mat'); % all_time

% Cortex and source space for representations
modeldir = 'C:\Users\vivianahc\Documents\MATLAB\EEGData-Understanding Human Individuation\AdultMNI152model\';
load(fullfile(modeldir,'GMSurfaceMesh.mat'));
load('LFBiosemi128.mat'); % GMVoxels_mm

GM = double(GMVoxels_mm);
nPatients = numel(all_freq);
nVoxels = size(GM, 1);

% Output matrices
SBR_F1_bin = zeros(nVoxels, nPatients);
SBR_F0_bin = zeros(nVoxels, nPatients);
TRAP_F1_bin = zeros(nVoxels, nPatients);
TRAP_F0_bin = zeros(nVoxels, nPatients);

SBR_F1_amp = zeros(nVoxels, nPatients);
SBR_F0_amp = zeros(nVoxels, nPatients);
TRAP_F1_amp = zeros(nVoxels, nPatients);
TRAP_F0_amp = zeros(nVoxels, nPatients);

% Raw amplitudes (optional)
SBR_F1_amp_raw = cell(nPatients,1);
SBR_F0_amp_raw = cell(nPatients,1);
TRAP_F1_amp_raw = cell(nPatients,1);
TRAP_F0_amp_raw = cell(nPatients,1);

% SBR results
for p = 1:nPatients

    subj = all_freq{p};
    if isempty(subj), continue; end

    % SBR F1 
    if isfield(subj,'inddipsbr_f1') && ~isempty(subj.inddipsbr_f1)

        idx = double(subj.inddipsbr_f1(:));
        good = idx >= 1 & idx <= nVoxels;
        idx = idx(good);

        SBR_F1_bin(unique(idx), p) = 1;

        A = double(subj.ampsbr_f1);
        SBR_F1_amp_raw{p} = A;

        if ~isempty(A)
            A = A(good,:);

            nH = size(A,2)/2;
            Are = A(:,1:nH);
            Aim = A(:,nH+1:end);

            % magnitude per harmonic, then sum harmonics
            aScalar = sum(sqrt(Are.^2 + Aim.^2), 2);

            for k = 1:numel(idx)
                v = idx(k);
                SBR_F1_amp(v,p) = max(SBR_F1_amp(v,p), aScalar(k));
            end
        end
    end

    % SBR F0
    if isfield(subj,'inddipsbr_f0') && ~isempty(subj.inddipsbr_f0)

        idx = double(subj.inddipsbr_f0(:));
        good = idx >= 1 & idx <= nVoxels;
        idx = idx(good);

        SBR_F0_bin(unique(idx), p) = 1;

        A = double(subj.ampsbr_f0);
        SBR_F0_amp_raw{p} = A;

        if ~isempty(A)
            A = A(good,:);

            nH = size(A,2)/2;
            Are = A(:,1:nH);
            Aim = A(:,nH+1:end);

            aScalar = sum(sqrt(Are.^2 + Aim.^2), 2);

            for k = 1:numel(idx)
                v = idx(k);
                SBR_F0_amp(v,p) = max(SBR_F0_amp(v,p), aScalar(k));
            end
        end
    end
end

% TRAP MUSIC results
for p = 1:nPatients

    subj = all_time{p};
    if isempty(subj), continue; end

    % TRAP F1 
    if isfield(subj,'inddiptrap_f1') && ~isempty(subj.inddiptrap_f1)

        idx = double(subj.inddiptrap_f1(:));
        good = idx >= 1 & idx <= nVoxels;
        idx = idx(good);

        TRAP_F1_bin(unique(idx), p) = 1;

        A = double(subj.amptrap_f1);
        TRAP_F1_amp_raw{p} = A;

        if ~isempty(A)
            A = A(good,:);

            % RMS amplitude across time, one scalar per source
            aScalar = sqrt(mean(A.^2, 2));

            for k = 1:numel(idx)
                v = idx(k);
                TRAP_F1_amp(v,p) = max(TRAP_F1_amp(v,p), aScalar(k));
            end
        end
    end

    % TRAP F0
    if isfield(subj,'inddiptrap_f0') && ~isempty(subj.inddiptrap_f0)

        idx = double(subj.inddiptrap_f0(:));
        good = idx >= 1 & idx <= nVoxels;
        idx = idx(good);

        TRAP_F0_bin(unique(idx), p) = 1;

        A = double(subj.amptrap_f0);
        TRAP_F0_amp_raw{p} = A;

        if ~isempty(A)
            A = A(good,:);

            aScalar = sqrt(mean(A.^2, 2));

            for k = 1:numel(idx)
                v = idx(k);
                TRAP_F0_amp(v,p) = max(TRAP_F0_amp(v,p), aScalar(k));
            end
        end
    end
end

% Group summaries
% Number of participants reconstructing each voxel
SBR_F1_count = sum(SBR_F1_bin, 2);
SBR_F0_count = sum(SBR_F0_bin, 2);
TRAP_F1_count = sum(TRAP_F1_bin, 2);
TRAP_F0_count = sum(TRAP_F0_bin, 2);

% Visualization
figure('Color','w','Name','Histograms');
tiledlayout(2,2);
nexttile; title(sprintf('SBR F1')); hold on;
histogram(SBR_F1_count);

nexttile; title(sprintf('SBR F0')); hold on;
histogram(SBR_F0_count);

nexttile; title(sprintf('TRAP F1')); hold on;
histogram(TRAP_F1_count);

nexttile; title(sprintf('TRAP F0')); hold on;
histogram(TRAP_F0_count);

% Mean amplitude per voxel across all participants
SBR_F1_amp_mean = mean(SBR_F1_amp, 2);
SBR_F0_amp_mean = mean(SBR_F0_amp, 2);
TRAP_F1_amp_mean = mean(TRAP_F1_amp, 2);
TRAP_F0_amp_mean = mean(TRAP_F0_amp, 2);

% Mean amplitude only among participants where voxel was reconstructed
SBR_F1_amp_detected = mean_nonzero(SBR_F1_amp);
SBR_F0_amp_detected = mean_nonzero(SBR_F0_amp);
TRAP_F1_amp_detected = mean_nonzero(TRAP_F1_amp);
TRAP_F0_amp_detected = mean_nonzero(TRAP_F0_amp);

% Thresholded maps for cortex representations
% countThr = 10 means the voxel was reconstructed in at least 10/130 patients
% This keeps the representation conservative and avoids displaying very rare sources
countThr = 10;

SBR_F1_count_thr = threshold_map(SBR_F1_count, countThr);
SBR_F0_count_thr = threshold_map(SBR_F0_count, countThr);
TRAP_F1_count_thr = threshold_map(TRAP_F1_count, countThr);
TRAP_F0_count_thr = threshold_map(TRAP_F0_count, countThr);

% Amplitude maps are displayed only on voxels that pass the count threshold.
% This avoids interpreting large amplitudes from rarely reconstructed sources.
SBR_F1_amp_thr = SBR_F1_amp_detected;
SBR_F0_amp_thr = SBR_F0_amp_detected;
TRAP_F1_amp_thr = TRAP_F1_amp_detected;
TRAP_F0_amp_thr = TRAP_F0_amp_detected;

SBR_F1_amp_thr(SBR_F1_count < countThr) = 0;
SBR_F0_amp_thr(SBR_F0_count < countThr) = 0;
TRAP_F1_amp_thr(TRAP_F1_count < countThr) = 0;
TRAP_F0_amp_thr(TRAP_F0_count < countThr) = 0;

% Cortex representations using thresholded count values
ThresholdedVoxel(GMSurfaceMesh, GM, ...
    SBR_F1_count_thr, SBR_F0_count_thr, ...
    TRAP_F1_count_thr, TRAP_F0_count_thr, ...
    'Thresholded voxel count maps');

% Cortex representations using amplitude magnitudes on the same thresholded support
ThresholdedVoxel(GMSurfaceMesh, GM, ...
    SBR_F1_amp_thr, SBR_F0_amp_thr, ...
    TRAP_F1_amp_thr, TRAP_F0_amp_thr, ...
    'Thresholded voxel amplitude maps');

% Slices using count values
SlicesMiddleCross(GMSurfaceMesh, GM, ...
    SBR_F1_count_thr, SBR_F0_count_thr, ...
    TRAP_F1_count_thr, TRAP_F0_count_thr, ...
    'Counts');

% Slices using amplitude magnitudes
SlicesMiddleCross(GMSurfaceMesh, GM, ...
    SBR_F1_amp_thr, SBR_F0_amp_thr, ...
    TRAP_F1_amp_thr, TRAP_F0_amp_thr, ...
    'Amplitudes');

% Save output
save('source_localization_summary.mat', ...
    'SBR_F1_bin','SBR_F0_bin','TRAP_F1_bin','TRAP_F0_bin', ...
    'SBR_F1_amp','SBR_F0_amp','TRAP_F1_amp','TRAP_F0_amp', ...
    'SBR_F1_count','SBR_F0_count','TRAP_F1_count','TRAP_F0_count', ...
    'SBR_F1_count_thr','SBR_F0_count_thr','TRAP_F1_count_thr','TRAP_F0_count_thr', ...
    'SBR_F1_amp_mean','SBR_F0_amp_mean','TRAP_F1_amp_mean','TRAP_F0_amp_mean', ...
    'SBR_F1_amp_detected','SBR_F0_amp_detected', ...
    'TRAP_F1_amp_detected','TRAP_F0_amp_detected', ...
    'SBR_F1_amp_thr','SBR_F0_amp_thr','TRAP_F1_amp_thr','TRAP_F0_amp_thr', ...
    'SBR_F1_amp_raw','SBR_F0_amp_raw','TRAP_F1_amp_raw','TRAP_F0_amp_raw');

% Helper functions
function m = mean_nonzero(A)
% Mean across patients only where amplitude is nonzero
m = zeros(size(A,1),1);
for i = 1:size(A,1)
    x = A(i,:);
    x = x(x > 0);
    if ~isempty(x)
        m(i) = mean(x);
    end
end
end

function v_thr = threshold_map(v, thr)
% Keep only voxels reconstructed in at least thr patients
v_thr = zeros(size(v));
m = v >= thr;
v_thr(m) = v(m);
end

function ThresholdedVoxel(mesh, GM, v_sbr_f1, v_sbr_f0, v_trap_f1, v_trap_f0, figName)
figure('Color','w','Name',figName);
tiledlayout(2,2);

nexttile; title('SBR F1'); hold on;
CortexShell(mesh);
m = v_sbr_f1 > 0;
scatter3(GM(m,1), GM(m,2), GM(m,3), 30, v_sbr_f1(m), 'filled');
axis off equal; colorbar;

nexttile; title('SBR F0'); hold on;
CortexShell(mesh);
m = v_sbr_f0 > 0;
scatter3(GM(m,1), GM(m,2), GM(m,3), 30, v_sbr_f0(m), 'filled');
axis off equal; colorbar;

nexttile; title('TRAP F1'); hold on;
CortexShell(mesh);
m = v_trap_f1 > 0;
scatter3(GM(m,1), GM(m,2), GM(m,3), 30, v_trap_f1(m), 'filled');
axis off equal; colorbar;

nexttile; title('TRAP F0'); hold on;
CortexShell(mesh);
m = v_trap_f0 > 0;
scatter3(GM(m,1), GM(m,2), GM(m,3), 30, v_trap_f0(m), 'filled');
axis off equal; colorbar;
end

function SlicesMiddleCross(mesh, GM, v1, v0, t1, t0, modeLabel)
% Smooth slice rendering of scalar voxel values
sigmaMM = 7;
maxDistMM = 15;

s1 = smoothSourceValues(GM, v1, sigmaMM, maxDistMM);
s0 = smoothSourceValues(GM, v0, sigmaMM, maxDistMM);
u1 = smoothSourceValues(GM, t1, sigmaMM, maxDistMM);
u0 = smoothSourceValues(GM, t0, sigmaMM, maxDistMM);

figure('Color','w','Name',['Middle slabs - SBR F1 - ' modeLabel]);
SliceOnAxesMiddleCross(mesh, GM, s1);

figure('Color','w','Name',['Middle slabs - SBR F0 - ' modeLabel]);
SliceOnAxesMiddleCross(mesh, GM, s0);

figure('Color','w','Name',['Middle slabs - TRAP F1 - ' modeLabel]);
SliceOnAxesMiddleCross(mesh, GM, u1);

figure('Color','w','Name',['Middle slabs - TRAP F0 - ' modeLabel]);
SliceOnAxesMiddleCross(mesh, GM, u0);
end

function SliceOnAxesMiddleCross(mesh, GM, vals)
% Representation of one coronal, horizontal, and sagittal slab in the middle

V = mesh.node;
F = mesh.face;
cortex.faces = F;
cortex.vertices = V;

% Grids 
step = 0.5;
x = min(V(:,1)):step:max(V(:,1));
y = min(V(:,2)):step:max(V(:,2));
z = min(V(:,3)):step:max(V(:,3));
[Xc,Zc] = meshgrid(x,z);
[Xh,Yh] = meshgrid(x,y);
[Ys,Zs] = meshgrid(y,z);

% Interpolant from all values
Fi = scatteredInterpolant( ...
    GM(:,1), GM(:,2), GM(:,3), double(vals), ...
    'natural', 'nearest');

hold on;

% Middle coronal slab
yMin = min(V(:,2));
yMax = max(V(:,2));
yMid = 0.5 * (yMin + yMax);

% Define slab tickness and layers
coronalHalfThickness = 2;
coronalStep = 1;
coronalOffsets = -coronalHalfThickness:coronalStep:coronalHalfThickness;

for dy = coronalOffsets
    Yc = (yMid + dy) * ones(size(Xc));
    Qc = [Xc(:), Yc(:), Zc(:)];
    INc = inpolyhedron(cortex, Qc);
    INc = reshape(INc, size(Xc));
    if ~any(INc(:)), continue; end

    Gc = Fi(Xc,Yc,Zc);
    Gc(~INc) = NaN;

    surf(Xc, Yc, Zc, Gc, ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.35);
end

% Middle horizontal slab
zMin = min(V(:,3));
zMax = max(V(:,3));
zMid = 0.5 * (zMin + zMax);

% Define slab tickness and layers
horizontalHalfThickness = 2;
horizontalStep = 1;
horizontalOffsets = -horizontalHalfThickness:horizontalStep:horizontalHalfThickness;

for dz = horizontalOffsets
    Zh = (zMid + dz) * ones(size(Xh));
    Qh = [Xh(:), Yh(:), Zh(:)];
    INh = inpolyhedron(cortex, Qh);
    INh = reshape(INh, size(Xh));
    if ~any(INh(:)), continue; end

    Gh = Fi(Xh, Yh, Zh);
    Gh(~INh) = NaN;

    surf(Xh, Yh, Zh, Gh, ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.35);
end

% Middle sagittal slab
xMin = min(V(:,1));
xMax = max(V(:,1));
xMid = 0.5 * (xMin + xMax);

% Define slab tickness and layers
sagHalfThickness = 2;
sagStep = 1;
sagOffsets = -sagHalfThickness:sagStep:sagHalfThickness;

for dx = sagOffsets
    Xs = (xMid + dx) * ones(size(Ys));
    Qs = [Xs(:), Ys(:), Zs(:)];
    INs = inpolyhedron(cortex, Qs);
    INs = reshape(INs, size(Ys));
    if ~any(INs(:)), continue; end

    Gs = Fi(Xs, Ys, Zs);
    Gs(~INs) = NaN;

    surf(Xs, Ys, Zs, Gs, ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.35);
end
CortexShell(mesh)
end

function CortexShell(mesh)
% Cortex representation 
V = mesh.node;
F = mesh.face;
patch('Faces', F, 'Vertices', V, ...
    'FaceColor', [0.7 0.7 0.7], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.1);
axis('equal');
axis('off');
view(0,90) %axial
end

function valsSm = smoothSourceValues(GM, vals, sigmaMM, radiusMM)
% Smooth values defined on source points using a Gaussian kernel

idxNZ = vals > 0;
XYZ = GM(idxNZ,:);
V = vals(idxNZ);

valsSm = zeros(size(vals));

if isempty(V)
    return;
end

for i = 1:size(GM,1)
    d2 = sum((XYZ - GM(i,:)).^2, 2);
    m = d2 <= radiusMM^2;
    if ~any(m), continue; end

    w = exp(-d2(m) / (2*sigmaMM^2));
    valsSm(i) = sum(w .* V(m)) / sum(w);
end
end
