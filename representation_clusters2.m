clear; close all; clc;

% -------------------------------------------------------------------------
% This script produces 3D visualization of source localization results in 
% cortex space:
% 1) Histogram per source =  #subjects that reconstructed it. It shows where
%    reconstructions happen and how often.
% 2) Spatial summary using DBSCAN clustering of nonzero sources (<10 mm).
%    It provides a spatially summary of hotspots.
% 3) Cluster strength where each cluster summarized by a single scalar:
%    strength = sum(histogram values over cluster members). It representes
%    how much subject-level support a spatial hotspot has
% 4) Cortex slices  
% -------------------------------------------------------------------------

% Load data
modeldir = 'C:\Users\vivianahc\Documents\MATLAB\EEGData-Understanding Human Individuation\AdultMNI152model\';
load(fullfile(modeldir,'HeadVolumeMesh.mat'));
load(fullfile(modeldir,'GMSurfaceMesh.mat')); % GMSurfaceMesh.node, .face 
load(fullfile(modeldir,'ScalpSurfaceMesh.mat'));
load LFBiosemi128.mat % GMVoxels_mm(Nsources x 3)

direc = 'C:\Users\vivianahc\Documents\MATLAB\EEGData-Understanding Human Individuation\Outputs';
load(fullfile(direc,'all_freq.mat')); % frequency domain results (SBR)
load(fullfile(direc,'all_time.mat')); % time domain results (TRAP/MUSIC)

GM = double(GMVoxels_mm);
Nsources = size(GM,1);
nSubjects = numel(all_freq);

% Choose histogram mode
useSubjectUnique = true; % group suport interpretation 
computeNonUnique = false; % quantifying redundancy/efficency on algs 

% Choose what visualize('count' subject count, 'amp' based on source amps
strengthMode = 'count';%'count';

% ------------------------ Histograms -------------------------------------
% Each subject contributes max 1 per source 
histU_sbr_f1 = buildSubjectHistogram(all_freq, 'inddipsbr_f1', Nsources);
histU_sbr_f0 = buildSubjectHistogram(all_freq, 'inddipsbr_f0', Nsources);
histU_trap_f1 = buildSubjectHistogram(all_time, 'inddiptrap_f1', Nsources);
histU_trap_f0 = buildSubjectHistogram(all_time, 'inddiptrap_f0', Nsources);

% Count all detections (duplicates included)
if computeNonUnique
    histD_sbr_f1 = buildDetectionHistogram(all_freq, 'inddipsbr_f1', Nsources);
    histD_sbr_f0 = buildDetectionHistogram(all_freq, 'inddipsbr_f0', Nsources);
    histD_trap_f1 = buildDetectionHistogram(all_time, 'inddiptrap_f1', Nsources);
    histD_trap_f0 = buildDetectionHistogram(all_time, 'inddiptrap_f0', Nsources);

    % How often the same source repeats
    % redundancy = total_detections - unique_detections
    eff_sbr_f1 = detectionRedundancyPerSubject(all_freq, 'inddipsbr_f1');
    eff_sbr_f0 = detectionRedundancyPerSubject(all_freq, 'inddipsbr_f0');
    eff_trap_f1 = detectionRedundancyPerSubject(all_time, 'inddiptrap_f1');
    eff_trap_f0 = detectionRedundancyPerSubject(all_time, 'inddiptrap_f0');

    fprintf('SBR F1: mean=%g duplicates/subject\n', mean(eff_sbr_f1));
    fprintf('SBR F0: mean=%g duplicates/subject\n', mean(eff_sbr_f0));
    fprintf('TRAP F1: mean=%g duplicates/subject\n', mean(eff_trap_f1));
    fprintf('TRAP F0: mean=%g duplicates/subject\n\n', mean(eff_trap_f0));
end

% Select which histogram is used for clustering + slice maps
% subject-unique: cluster strength means how many subjects support region
% detection: cluster strength means how many detections support region
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

% Mean normalized amplitude maps: mean over subjects of "normalized amplitude
% mass" assigned to voxel i
meanAmpNorm_sbr_f1 = buildMeanNormalizedAmplitudeMap(all_freq, 'inddipsbr_f1', 'ampsbr_f1', Nsources);
meanAmpNorm_sbr_f0 = buildMeanNormalizedAmplitudeMap(all_freq, 'inddipsbr_f0', 'ampsbr_f0', Nsources);
meanAmpNorm_trap_f1 = buildMeanNormalizedAmplitudeMap(all_time, 'inddiptrap_f1', 'amptrap_f1', Nsources);
meanAmpNorm_trap_f0 = buildMeanNormalizedAmplitudeMap(all_time, 'inddiptrap_f0', 'amptrap_f0', Nsources);

% Clustering - Spatial Summary
% Represent spatial group of voxels close each other
epsMM = 10; % <= 10mm for same cluster 
minPts = 3; % cluster must have at least 3 supporting sources (adjustable)

C_SBR_F1 = clusterFromHistogram(GM, hist_sbr_f1, epsMM, minPts);
C_SBR_F0 = clusterFromHistogram(GM, hist_sbr_f0, epsMM, minPts);
C_TRAP_F1 = clusterFromHistogram(GM, hist_trap_f1, epsMM, minPts);
C_TRAP_F0 = clusterFromHistogram(GM, hist_trap_f0, epsMM, minPts);

% Compute strength per cluster
% Strength based on count 
C_SBR_F1.strength_count = clusterStrengthFromMap(C_SBR_F1, hist_sbr_f1);
C_SBR_F0.strength_count = clusterStrengthFromMap(C_SBR_F0, hist_sbr_f0);
C_TRAP_F1.strength_count = clusterStrengthFromMap(C_TRAP_F1, hist_trap_f1);
C_TRAP_F0.strength_count = clusterStrengthFromMap(C_TRAP_F0, hist_trap_f0);

% Strength based on amplitude 
C_SBR_F1.strength_amp = clusterStrengthFromMap(C_SBR_F1, meanAmpNorm_sbr_f1);
C_SBR_F0.strength_amp = clusterStrengthFromMap(C_SBR_F0, meanAmpNorm_sbr_f0);
C_TRAP_F1.strength_amp = clusterStrengthFromMap(C_TRAP_F1, meanAmpNorm_trap_f1);
C_TRAP_F0.strength_amp =clusterStrengthFromMap(C_TRAP_F0, meanAmpNorm_trap_f0);

% Maps for Visualization
% Count-based voxel maps 
map_count_SBR_F1 = clusterStrengthMap(Nsources, C_SBR_F1, 'count');
map_count_SBR_F0 = clusterStrengthMap(Nsources, C_SBR_F0, 'count');
map_count_TRAP_F1 = clusterStrengthMap(Nsources, C_TRAP_F1, 'count');
map_count_TRAP_F0 = clusterStrengthMap(Nsources, C_TRAP_F0, 'count');

% Amplitude-based voxel maps
map_amp_SBR_F1 = clusterStrengthMap(Nsources, C_SBR_F1, 'amp');
map_amp_SBR_F0 = clusterStrengthMap(Nsources, C_SBR_F0, 'amp');
map_amp_TRAP_F1 = clusterStrengthMap(Nsources, C_TRAP_F1, 'amp'); 
map_amp_TRAP_F0 = clusterStrengthMap(Nsources, C_TRAP_F0, 'amp'); 

% Select which one to plot
switch lower(strengthMode)
    case 'count'
        map_strength_SBR_F1 = map_count_SBR_F1;
        map_strength_SBR_F0 = map_count_SBR_F0;
        map_strength_TRAP_F1 = map_count_TRAP_F1;
        map_strength_TRAP_F0 = map_count_TRAP_F0;
        centroidStrength = @(C) C.strength_count;
    case 'amp'
        map_strength_SBR_F1  = map_amp_SBR_F1;
        map_strength_SBR_F0  = map_amp_SBR_F0;
        map_strength_TRAP_F1 = map_amp_TRAP_F1;
        map_strength_TRAP_F0 = map_amp_TRAP_F0;
        centroidStrength = @(C) C.strength_amp;
    otherwise
        error('strengthMode must be ''count'' or ''amp''.');
end

% ------------------- Figure 1: Raw histograms -------------------------
thrnbs = 5; % threshold only for readability
plotRawHistogramFigure(GMSurfaceMesh, GM, ...
    hist_sbr_f1, hist_sbr_f0, hist_trap_f1, hist_trap_f0, thrnbs);

% -------------------- Figure 2: Cluster centroids --------------------- 
plotClusterCentroidsFigure(GMSurfaceMesh, ...
    C_SBR_F1, C_SBR_F0, C_TRAP_F1, C_TRAP_F0, strengthMode);

% ---------------- Figure 3: Inside cortex cut planes -------------- 
plotSlices(GMSurfaceMesh, GM, ...
    map_strength_SBR_F1, map_strength_SBR_F0, ...
    map_strength_TRAP_F1, map_strength_TRAP_F0);

% -------------------------------------------------------------------------
% FUNCTIONS
% -------------------------------------------------------------------------
function histv = buildSubjectHistogram(allStruct, fieldname, Nsources)
% #subjects that reconstructed source i at least once (unique per subject)
histv = zeros(Nsources,1);
for s = 1:numel(allStruct)
    subj = allStruct{s};
    if isempty(subj) || ~isfield(subj, fieldname) || isempty(subj.(fieldname)), continue; end
    idx = unique(subj.(fieldname)(:));
    idx = idx(idx>=1 & idx<=Nsources);
    histv(idx) = histv(idx) + 1;
end
end

function histv = buildDetectionHistogram(allStruct, fieldname, Nsources)
% total detections across all subjects (duplicates included)
histv = zeros(Nsources,1);
for s = 1:numel(allStruct)
    subj = allStruct{s};
    if isempty(subj) || ~isfield(subj, fieldname) || isempty(subj.(fieldname)), continue; end
    idx = subj.(fieldname)(:);
    idx = idx(idx>=1 & idx<=Nsources);
    histv(idx) = histv(idx) + 1;
end
end

function dup = detectionRedundancyPerSubject(allStruct, fieldname)
dup = zeros(numel(allStruct),1);
for s = 1:numel(allStruct)
    subj = allStruct{s};
    if isempty(subj) || ~isfield(subj, fieldname) || isempty(subj.(fieldname))
        dup(s) = 0; continue;
    end
    idx = subj.(fieldname)(:);
    dup(s) = numel(idx) - numel(unique(idx));
end
end

function meanAmpNorm = buildMeanNormalizedAmplitudeMap(allStruct, idxField, ampField, Nsources)
acc = zeros(Nsources,1);
nValid = 0;

for s = 1:numel(allStruct)
    subj = allStruct{s};
    if isempty(subj) || ~isfield(subj, idxField) || isempty(subj.(idxField)), continue; end

    idx = subj.(idxField)(:);
    idx = idx(idx>=1 & idx<=Nsources);

    A = subj.(ampField);

    % Make amplitude a scalar per detection:
    Nd = numel(idx);
    m = min(size(A,1), Nd);
    idx = idx(1:m);
    a = sqrt(sum(abs(double(A(1:m,:))).^2, 2));

    % Accumulate duplicates for this subject (voxel-wise sum amplitude)
    subjMap = accumarray(double(idx), a, [Nsources 1], @sum, 0);

    % Normalize within subject
    den = sum(subjMap);
    subjMap = subjMap ./ den;

    acc = acc + subjMap;
    nValid = nValid + 1;
end
    meanAmpNorm = acc ./ nValid;
end

function C = clusterFromHistogram(GM, histv, epsMM, minPts)
% DBSCAN over coordinates of sources with histv>0
idxNZ = find(histv > 0);
pts   = GM(idxNZ,:);

labels = dbscan(pts, epsMM, minPts);
K = max(labels);

centroid = zeros(K,3);
members  = cell(K,1);

for k = 1:K
    m = (labels == k);
    idxMembers = idxNZ(m);
    members{k} = idxMembers;
    centroid(k,:) = mean(GM(idxMembers,:),1);
end

C = struct('idxNZ',idxNZ,'labels',labels,'K',K,'centroid',centroid,...
    'members',{members}); 
end

function strength = clusterStrengthFromMap(C, voxelMap)
% Sum(voxelMap(member voxels of cluster k))
strength = zeros(C.K,1);
voxelMap = double(voxelMap(:));

for k = 1:C.K
    strength(k) = sum(voxelMap(C.members{k}));
end
end

function vals = clusterStrengthMap(Nsources, C, whichStrength)
% Each voxel in cluster k receives the scalar cluster strength(count or amp)
vals = zeros(Nsources,1);

switch lower(whichStrength)
    case 'count'
        if ~isfield(C,'strength_count'), return; end
        svec = C.strength_count;
    case 'amp'
        if ~isfield(C,'strength_amp'), return; end
        svec = C.strength_amp;
    otherwise
        error('whichStrength must be ''count'' or ''amp''.');
end

for k = 1:C.K
    vals(C.members{k}) = svec(k);
end
end

function plotRawHistogramFigure(mesh, GM, h_sbr_f1, h_sbr_f0, h_trap_f1, h_trap_f0, thrnbs)
figure('Color','w','Name','Raw histograms (thresholded for readability)');
tiledlayout(2,2);

nexttile; title('SBR F1'); hold on;
plotCortexShell(mesh);
scatter3(GM(h_sbr_f1>thrnbs,1),GM(h_sbr_f1>thrnbs,2),GM(h_sbr_f1>thrnbs,3), ...
    30, h_sbr_f1(h_sbr_f1>thrnbs),'filled');
axis off equal; colorbar;

nexttile; title('SBR F0'); hold on;
plotCortexShell(mesh);
scatter3(GM(h_sbr_f0>thrnbs,1),GM(h_sbr_f0>thrnbs,2),GM(h_sbr_f0>thrnbs,3), ...
    30, h_sbr_f0(h_sbr_f0>thrnbs),'filled');
axis off equal; colorbar;

nexttile; title('TRAP F1'); hold on;
plotCortexShell(mesh);
scatter3(GM(h_trap_f1>thrnbs,1),GM(h_trap_f1>thrnbs,2),GM(h_trap_f1>thrnbs,3), ...
    30, h_trap_f1(h_trap_f1>thrnbs),'filled');
axis off equal; colorbar;

nexttile; title('TRAP F0'); hold on;
plotCortexShell(mesh);
scatter3(GM(h_trap_f0>thrnbs,1),GM(h_trap_f0>thrnbs,2),GM(h_trap_f0>thrnbs,3), ...
    30, h_trap_f0(h_trap_f0>thrnbs),'filled');
axis off equal; colorbar;
end

function plotCortexShell(mesh)
V = mesh.node; F = mesh.face;
patch('Faces',F,'Vertices',V, ...
    'FaceColor',[0.82 0.82 0.82], 'EdgeColor','none', 'FaceAlpha',0.08);
axis equal;
camlight headlight; camlight right;
end

function plotClusterCentroidsFigure(mesh, C1, C2, C3, C4, strengthMode)
figure('Color','w','Name','Cluster centroids (colored by chosen strength)');
tiledlayout(2,2);

nexttile; title('SBR F1'); hold on; plotCortexShell(mesh);
plotCentroidsOnAxes(gca, C1, strengthMode); axis off equal; colorbar;

nexttile; title('SBR F0'); hold on; plotCortexShell(mesh);
plotCentroidsOnAxes(gca, C2, strengthMode); axis off equal; colorbar;

nexttile; title('TRAP F1'); hold on; plotCortexShell(mesh);
plotCentroidsOnAxes(gca, C3, strengthMode); axis off equal; colorbar;

nexttile; title('TRAP F0'); hold on; plotCortexShell(mesh);
plotCentroidsOnAxes(gca, C4, strengthMode); axis off equal; colorbar;
end

function plotCentroidsOnAxes(ax, C, strengthMode)
% Choose which strength to color with
switch lower(strengthMode)
    case 'count'
        if ~isfield(C,'strength_count') || isempty(C.strength_count), return; end
        cval = C.strength_count(:);
    case 'amp'
        if ~isfield(C,'strength_amp') || isempty(C.strength_amp), return; end
        cval = C.strength_amp(:);
    otherwise
        error('strengthModeToPlot must be ''count'' or ''amp''.');
end
markerSize = 40; % constant size (only color encodes strength)
scatter3(ax, C.centroid(:,1), C.centroid(:,2), C.centroid(:,3), ...
    markerSize, cval, 'filled');
end

function plotSlices(mesh, GM, v1, v0, t1, t0)
figure('Color','w','Name','Inside-cortex cut planes');
tiledlayout(2,2);

ax1 = nexttile; title(ax1,'SBR F1');  plotSliceOnAxes(ax1, mesh, GM, v1);
ax2 = nexttile; title(ax2,'SBR F0');  plotSliceOnAxes(ax2, mesh, GM, v0);
ax3 = nexttile; title(ax3,'TRAP F1'); plotSliceOnAxes(ax3, mesh, GM, t1);
ax4 = nexttile; title(ax4,'TRAP F0'); plotSliceOnAxes(ax4, mesh, GM, t0);
end

function plotSliceOnAxes(ax, GMSurfaceMesh, GMVoxels_mm, vals)
% Coronal slices 

ySlices = -150:10:150; % coronal slice positions along Y (mm)
stepXZ  = 0.3; % grid step in X/Z

V = GMSurfaceMesh.node;
F = GMSurfaceMesh.face;
cortex.faces = F;
cortex.vertices = V;

x = min(V(:,1)):stepXZ:max(V(:,1));
z = min(V(:,3)):stepXZ:max(V(:,3));
[X,Z] = meshgrid(x,z);

Fi = scatteredInterpolant( ...
    GMVoxels_mm(:,1), GMVoxels_mm(:,2), GMVoxels_mm(:,3), double(vals), ...
    'natural','nearest');

hold(ax,'on');
axis(ax,'equal');
view(ax,[115 18]);

for y0 = ySlices
    Y = y0 * ones(size(X));
    Q  = [X(:), Y(:), Z(:)];
    IN = inpolyhedron(cortex, Q);
    IN = reshape(IN, size(X));
    if ~any(IN(:)), continue; end

    G = Fi(X,Y,Z);
    G(~IN) = NaN;

    surf(ax, X, Y, Z, G, 'EdgeColor','none');
end

axis(ax,'off');
colormap(ax, parula);
colorbar(ax);
end