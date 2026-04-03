clear; close all; clc;

% -------------------------------------------------------------------------
% Visualization at the group level for source localization results 
% - For each source (voxel), count how many subjects reconstructed it
%   histogram representing #subjects that selected i source at least once
%
% - Keep only voxels with sufficient reproducibility: histo(i) >= countThr
%
% - Cluster voxels in source space using DBSCAN: constraint with a spatial 
%   distance criterion (epsMM = 10 mm)
%
% - Compute a cluster support: sum normalized of voxel counts over cluster 
%   members
%
% - Generate visualization maps where voxels inside retained clusters keep 
%   their normalized sum. All other voxels are reset to 0 
% -------------------------------------------------------------------------

% Load data
modeldir = 'C:\Users\vivianahc\Documents\MATLAB\EEGData-Understanding Human Individuation\AdultMNI152model\';
load(fullfile(modeldir,'HeadVolumeMesh.mat'));
load(fullfile(modeldir,'GMSurfaceMesh.mat'));
load(fullfile(modeldir,'ScalpSurfaceMesh.mat'));
load('LFBiosemi128.mat'); % GMVoxels_mm

direc = 'C:\Users\vivianahc\Documents\MATLAB\EEGData-Understanding Human Individuation';
load(fullfile(direc,'all_freq.mat')); % SBR results
load(fullfile(direc,'all_time.mat')); % TRAP results

GM = double(GMVoxels_mm);
Nsources = size(GM,1);
nSubjects = numel(all_freq);

% Parameters
useSubjectUnique = true; % count each subject max once per voxel
computeNonUnique = false; % redundancy diagnostics (optional)

countThr = 10; % voxel (source) reconstructed by at least n subjects
epsMM = 10; % neighborhood radius in mm
minPts = 2; % minimum number of candidate voxels to form a cluster
thrnbs = 5; % Threshold in raw histogram (only for visualization in figure)

% Histograms 
% Each subject contributes max 1 per voxel
histU_sbr_f1 = SubjectHistogram(all_freq, 'inddipsbr_f1',  Nsources);
histU_sbr_f0 = SubjectHistogram(all_freq, 'inddipsbr_f0',  Nsources);
histU_trap_f1 = SubjectHistogram(all_time, 'inddiptrap_f1', Nsources);
histU_trap_f0 = SubjectHistogram(all_time, 'inddiptrap_f0', Nsources);

% Count all detections including duplicates (optional)
if computeNonUnique
    histD_sbr_f1 = DetectionHistogram(all_freq, 'inddipsbr_f1',  Nsources);
    histD_sbr_f0 = DetectionHistogram(all_freq, 'inddipsbr_f0',  Nsources);
    histD_trap_f1 = DetectionHistogram(all_time, 'inddiptrap_f1', Nsources);
    histD_trap_f0 = DetectionHistogram(all_time, 'inddiptrap_f0', Nsources);

    eff_sbr_f1 = DetectionRedundancySubject(all_freq, 'inddipsbr_f1');
    eff_sbr_f0 = DetectionRedundancySubject(all_freq, 'inddipsbr_f0');
    eff_trap_f1 = DetectionRedundancySubject(all_time, 'inddiptrap_f1');
    eff_trap_f0 = DetectionRedundancySubject(all_time, 'inddiptrap_f0');
end

% Select histogram 
if useSubjectUnique
    hist_sbr_f1 = histU_sbr_f1;
    hist_sbr_f0 = histU_sbr_f0;
    hist_trap_f1 = histU_trap_f1;
    hist_trap_f0 = histU_trap_f0;
else
    hist_sbr_f1 = histD_sbr_f1;
    hist_sbr_f0 = histD_sbr_f0;
    hist_trap_f1 = histD_trap_f1;
    hist_trap_f0 = histD_trap_f0;
end

% Threshold and clustering 
% Cluster only voxels above the count threshold
C_SBR_F1 = clusterThresholdedHistogram(GM, hist_sbr_f1, countThr, epsMM, minPts);
C_SBR_F0 = clusterThresholdedHistogram(GM, hist_sbr_f0, countThr, epsMM, minPts);
C_TRAP_F1 = clusterThresholdedHistogram(GM, hist_trap_f1, countThr, epsMM, minPts);
C_TRAP_F0 = clusterThresholdedHistogram(GM, hist_trap_f0, countThr, epsMM, minPts);

% Compute cluster support from counts
C_SBR_F1 = ClusterCountStrength(C_SBR_F1);
C_SBR_F0 = ClusterCountStrength(C_SBR_F0);
C_TRAP_F1 = ClusterCountStrength(C_TRAP_F1);
C_TRAP_F0 = ClusterCountStrength(C_TRAP_F0);

% ============ Figures ==============
% Raw histograms
RawHistogram(GMSurfaceMesh, GM, hist_sbr_f1, hist_sbr_f0, ...
    hist_trap_f1, hist_trap_f0, thrnbs);

% Thresholded voxel maps (no clustering)
ThresholdedVoxel(GMSurfaceMesh, GM, ...
    C_SBR_F1.hist_thr, C_SBR_F0.hist_thr, ...
    C_TRAP_F1.hist_thr, C_TRAP_F0.hist_thr);

% Thresholded voxels and clustered 
ThresholdedVoxel(GMSurfaceMesh, GM, ...
    C_SBR_F1.hist_thr_clust, C_SBR_F0.hist_thr_clust, ...
    C_TRAP_F1.hist_thr_clust, C_TRAP_F0.hist_thr_clust);

% Planes of cortex, voxels clustered
SlicesMiddleCross(GMSurfaceMesh, GM, ...
    C_SBR_F1.hist_thr_clust, C_SBR_F0.hist_thr_clust, ...
    C_TRAP_F1.hist_thr_clust, C_TRAP_F0.hist_thr_clust);

% ========== Support functions ============= 
function histv = SubjectHistogram(allStruct, fieldname, Nsources)
% Count how many subjects reconstructed voxel i at least once
histv = zeros(Nsources,1);
for s = 1:numel(allStruct)
    subj = allStruct{s};
    if isempty(subj) || ~isfield(subj, fieldname) || isempty(subj.(fieldname))
        continue;
    end
    idx = unique(subj.(fieldname)(:));
    idx = idx(idx >= 1 & idx <= Nsources);
    histv(idx) = histv(idx) + 1;
end
end

function histv = DetectionHistogram(allStruct, fieldname, Nsources)
% Count all detections including duplicates
histv = zeros(Nsources,1);
for s = 1:numel(allStruct)
    subj = allStruct{s};
    if isempty(subj) || ~isfield(subj, fieldname) || isempty(subj.(fieldname))
        continue;
    end
    idx = subj.(fieldname)(:);
    idx = idx(idx >= 1 & idx <= Nsources);
    histv(idx) = histv(idx) + 1;
end
end

function dup = DetectionRedundancySubject(allStruct, fieldname)
% Number of repeated detections per subject
dup = zeros(numel(allStruct),1);
for s = 1:numel(allStruct)
    subj = allStruct{s};
    if isempty(subj) || ~isfield(subj, fieldname) || isempty(subj.(fieldname))
        dup(s) = 0;
        continue;
    end
    idx = subj.(fieldname)(:);
    dup(s) = numel(idx) - numel(unique(idx));
end
end

function C = clusterThresholdedHistogram(GM, histv, countThr, epsMM, minPts)
% Cluster voxels above the threshold
hist_thr = zeros(size(histv));
mask = histv >= countThr;
hist_thr(mask) = histv(mask);

idxCand = find(hist_thr > 0);
pts = GM(idxCand,:);

if isempty(idxCand)
    C = struct('idxCand',[],'labels',[],'K',0,'centroid',[], ...
               'members',{{}},'strength_count',[],'peak_count',[], ...
               'hist_thr', hist_thr);
    return;
end

labels = dbscan(pts, epsMM, minPts);
goodLabels = unique(labels(labels > 0));
K = numel(goodLabels);

if K == 0
    C = struct('idxCand',idxCand,'labels',labels,'K',0,'centroid',[], ...
               'members',{{}},'strength_count',[],'peak_count', [],...
               'hist_thr', hist_thr);
    return;
end

centroid = zeros(K,3);
members = cell(K,1);

for j = 1:K
    thisLabel = goodLabels(j);
    m = (labels == thisLabel);
    idxMembers = idxCand(m);
    members{j} = idxMembers;
    centroid(j,:) = mean(GM(idxMembers,:), 1);
end

C = struct('idxCand',idxCand,'labels',labels,'K',K,'centroid',centroid, ...
           'members',{members},'strength_count',[],'peak_count',[], ...
           'hist_thr', hist_thr);
end

function C = ClusterCountStrength(C)
% Compute cluster counts 
C.hist_thr_clust = C.hist_thr; 

if C.K < 1
    C.strength_count = [];
    C.peak_count = [];
    return;
end

for k = 1:C.K
    vals = C.hist_thr(C.members{k});
    n = numel(vals);
    if n == 0
        continue;
    end
    % Different metrics
    sum_val = sum(vals);
    mean_val = mean(vals);
    peak_val = max(vals);
    % Balanced measure, sum normalized  
    strenght_balanced = sum_val / sqrt(n); 
    
    % Store
    C.strength_count(k) = strenght_balanced; 
    C.sum_count(k) = sum_val;
    C.mean_count(k) = mean_val;
    C.peak_count(k) = peak_val; % strongest voxel count in cluster  
    % Assign to members of cluster
    C.hist_thr_clust(C.members{k}) = strenght_balanced; 
end
C.K = numel(C.strength_count);
end

function RawHistogram(mesh, GM, h_sbr_f1, h_sbr_f0, h_trap_f1, h_trap_f0, thrnbs)
figure('Color','w','Name','Raw histograms (thresholded only for visualization)');
tiledlayout(2,2);

nexttile; title(sprintf('SBR F1 (count > %d)', thrnbs)); hold on;
CortexShell(mesh);
scatter3(GM(h_sbr_f1 > thrnbs,1), GM(h_sbr_f1 > thrnbs,2), ...
    GM(h_sbr_f1 > thrnbs,3), 30, h_sbr_f1(h_sbr_f1 > thrnbs), 'filled');
axis off equal; colorbar;

nexttile; title(sprintf('SBR F0 (count > %d)', thrnbs)); hold on;
CortexShell(mesh);
scatter3(GM(h_sbr_f0 > thrnbs,1), GM(h_sbr_f0 > thrnbs,2), ...
    GM(h_sbr_f0 > thrnbs,3), 30, h_sbr_f0(h_sbr_f0 > thrnbs), 'filled');
axis off equal; colorbar;

nexttile; title(sprintf('TRAP F1 (count > %d)', thrnbs)); hold on;
CortexShell(mesh);
scatter3(GM(h_trap_f1 > thrnbs,1), GM(h_trap_f1 > thrnbs,2), ...
    GM(h_trap_f1 > thrnbs,3), 30, h_trap_f1(h_trap_f1 > thrnbs), 'filled');
axis off equal; colorbar;

nexttile; title(sprintf('TRAP F0 (count > %d)', thrnbs)); hold on;
CortexShell(mesh);
scatter3(GM(h_trap_f0 > thrnbs,1), GM(h_trap_f0 > thrnbs,2), ...
    GM(h_trap_f0 > thrnbs,3), 30, h_trap_f0(h_trap_f0 > thrnbs), 'filled');
axis off equal; colorbar;
end

function ThresholdedVoxel(mesh, GM, v_sbr_f1, v_sbr_f0, v_trap_f1, v_trap_f0)
figure('Color','w','Name','Thresholded voxel count maps');
tiledlayout(2,2);

nexttile; title('SBR F1 thresholded'); hold on;
CortexShell(mesh);
m = v_sbr_f1 > 0;
scatter3(GM(m,1), GM(m,2), GM(m,3), 30, v_sbr_f1(m), 'filled');
axis off equal; colorbar;

nexttile; title('SBR F0 thresholded'); hold on;
CortexShell(mesh);
m = v_sbr_f0 > 0;
scatter3(GM(m,1), GM(m,2), GM(m,3), 30, v_sbr_f0(m), 'filled');
axis off equal; colorbar;

nexttile; title('TRAP F1 thresholded'); hold on;
CortexShell(mesh);
m = v_trap_f1 > 0;
scatter3(GM(m,1), GM(m,2), GM(m,3), 30, v_trap_f1(m), 'filled');
axis off equal; colorbar;

nexttile; title('TRAP F0 thresholded'); hold on;
CortexShell(mesh);
m = v_trap_f0 > 0;
scatter3(GM(m,1), GM(m,2), GM(m,3), 30, v_trap_f0(m), 'filled');
axis off equal; colorbar;
end

function SlicesMiddleCross(mesh, GM, v1, v0, t1, t0)
% Smooth slices
sigmaMM = 7; % smooth projection from source to cortex
maxDistMM = 15; % neighbors in this radius contribute to smothness 

s1 = smoothSourceValues(GM, v1, sigmaMM, maxDistMM);
s0 = smoothSourceValues(GM, v0, sigmaMM, maxDistMM);
u1 = smoothSourceValues(GM, t1, sigmaMM, maxDistMM);
u0 = smoothSourceValues(GM, t0, sigmaMM, maxDistMM);

figure('Color','w','Name','Middle slabs - SBR F1');
SliceOnAxesMiddleCross(mesh, GM, s1);

figure('Color','w','Name','Middle slabs - SBR F0');
SliceOnAxesMiddleCross(mesh, GM, s0);

figure('Color','w','Name','Middle slabs - TRAP F1');
SliceOnAxesMiddleCross(mesh, GM, u1);

figure('Color','w','Name','Middle slabs - TRAP F0');
SliceOnAxesMiddleCross(mesh, GM, u0);

end

function SliceOnAxesMiddleCross(mesh, GM, vals)
% Representation of one coronal, horizontal, and sagittal slab in the middle

V = mesh.node;
F = mesh.face;
cortex.faces = F;
cortex.vertices = V;

% Grids 
step = 0.3;
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

% ---- Middle coronal slab
yMin = min(V(:,2));
yMax = max(V(:,2));
yMid = 0.5 * (yMin + yMax);

% Define slab tickness and layers
coronalHalfThickness = 3;
coronalStep = 0.5;
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

% ----- Middle horizontal slab
zMin = min(V(:,3));
zMax = max(V(:,3));
zMid = 0.5 * (zMin + zMax);

% Define slab tickness and layers
horizontalHalfThickness = 3;
horizontalStep = 0.5;
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

% ---- Middle sagittal slab
xMin = min(V(:,1));
xMax = max(V(:,1));
xMid = 0.5 * (xMin + xMax);

% Define slab tickness and layers
sagHalfThickness = 3;
sagStep = 0.5;
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
view(102,0) %lateral representation
% view(90,0) %sagittal
% view(0,0) %coronal
% view(0,90) %axial
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
