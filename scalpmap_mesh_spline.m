function scalpmap_mesh_spline(velec, re, meshpoints, faces, m)
% Spherical spline interpolation of EEG scalp but plotted on the original 
% head mesh using plotmesh.
%   Inputs:
%     velec = potentials at electrodes (N x 1)
%     re = electrode positions (N x 3)
%     meshpoints = scalp surface vertices (M x 3)
%     faces = scalp surface faces (K x 3)
%     m = spline parameter (default = 4)

if nargin < 5, m = 4; end

% Keep original mesh for plotting
meshpoints_orig = meshpoints;

% Normalize electrode and mesh positions to unit sphere ---
re_unit = re ./ vecnorm(re,2,2);
mesh_unit = meshpoints ./ vecnorm(meshpoints,2,2);

% Electrode-electrode kernel
cosang = re_unit * re_unit';
G = spline_kernel(cosang, m);

% Regularized solve for spline coefficients
lambda = 1e-6;
c = (G + lambda*eye(size(G))) \ velec;

% Mesh-electrode kernel (on unit sphere)
cosang_mesh = mesh_unit * re_unit';
Gmesh = spline_kernel(cosang_mesh, m);
qc = Gmesh * c;

% Ensure non-negativity
qc(qc < 0) = 0;

% Plot using original scalp geometry with plotmesh
plotmesh([meshpoints_orig qc], faces, 'Edgealpha', 0.05);
hold on;
plot3(re(:,1), re(:,2), re(:,3), 'k*'); % electrode markers
axis equal; axis off;
camlight; lighting gouraud; material dull;
colorbar; caxis([0 max(qc)]);
view(90,0) %sagittal
%view(0,0) %coronal
%view(0,90) %axial
title('Scalp map');
end

function G = spline_kernel(cosang, m)
% Compute spherical spline kernel matrix
% cosang = dot products between unit vectors (cos(theta))
% m = flexibility parameter

nmax = 50; % number of Legendre terms
G = zeros(size(cosang));

% Clip values to [-1,1] for numerical stability
cosang = min(max(cosang,-1), 1);

for n = 1:nmax
    % Evaluate Legendre polynomial P_n at cos(theta)
    Pn = legendre(n, cosang(:)'); 
    Pn = squeeze(Pn(1,:,:)); % take the zonal (m=0) term
    coef = (2*n+1)/((n^m)*(n+1)^m);
    G = G + coef * reshape(Pn, size(cosang));
end
end
