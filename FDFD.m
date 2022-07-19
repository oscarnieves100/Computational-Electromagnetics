function [R,T,A] = FDFD(ER,UR,grid_width,grid_height,NPML,wavelength,...
    angle_incidence,polarization,scaling_factor,plotsOn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finite-Difference Frequency Domain solver in 2D using a PML in both ends
% of the y-direction (top and bottom) and periodic boundary conditions
% along the x-direction (left and right). Assuming a plane wave of
% wavelength "wavelength" (in meters), angle of incidence (in degrees),
% scaling factor in terms of units (e.g. scaling_factor = 1e-9 for
% nanometers). The other inputs into the program include two double-grid
% matrices ER (relative permittivity values) and UR (relative permeability 
% values) across the entire grid.
%
% Note that this program assumes a linear isotropic material where each
% cell of the grid (or each pixel in the image plot) is made up of a single
% isotropic material property, for instance epsilon_r = 5.8. For
% anisotropic materials, this particular program won't work because each
% grid cell will need to be made of a 3x3 matrix with all the tensor
% components of the permittivity and permeability. For such applications, a
% different implementation (slightly more complicated) will be required.
%
% The basic procedure for implementing this program is as follows: you
% create the matrices ER and UR by assigning relative permittivity and
% permeability values to the grid based on the materials/structures
% involved, and they must have size [Nx2, Ny2] (because we are using a
% double-grid approach). The size is selected by the user, but there are
% certain checks in this program which will determine whether the
% resolution you have selected is appropriate or not. For an example of
% this, see the script FDFD_cylinder.m in my GitHub repository. Generally
% speaking, you want to create a structure surrounded by free space (or
% vacuum) as the Perfectly Matched Layer (PML) is built around free space.
% If by any chance material or structure of your device happens to overlap
% with the PML, then this will lead to incorrect results for the
% reflectance and transmittance propertes, so special caution must be taken
% with this.
%
% Another input is polarization, which is either 'E' or 'H' (for TE and TM
% polarizations respectively). The program calculates reflectance
% properties based on the polarization.
%
% The outputs to this program are the Reflectance R, Transmittance T and
% Absorbance A of the structure in a single vector [R, T, A] and all in
% terms of single units rather than percentages (e.g. R=0.55 and T=0.45).
% You may also choose to set plotsOn = 1 or any number other than zero if
% you want to display plots, but the default is 0 so you can use this
% function in a loop (e.g. to compute the spectral response across multiple
% wavelengths).
%
% Written by: Oscar A. Nieves
% Last updated: 19/07/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Handle Input Arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Nx2, Ny2] = size(ER);
Nx = ceil(Nx2/2);
Ny = ceil(Ny2/2);

if ~exist('NPML','var')||isempty(NPML)
    NPML = ceil(Ny/10);
end

if ~exist('wavelength','var')||isempty(wavelength)
    wavelength = 600e-9;
end

if ~exist('angle_incidence','var')||isempty(angle_incidence)
    angle_incidence = 0;
end

if ~exist('polarization','var')||isempty(polarization)
    polarization = 'E';
end

if ~exist('scaling_factor','var')||isempty(scaling_factor)
    scaling_factor = 1e-9;
end

if ~exist('plotsOn','var')||isempty(plotsOn)
    plotsOn = 0;
end

% Normalize everything to wavelength
normal     = 1/wavelength;
wavelength = normal*wavelength;
width      = normal*grid_width;
maxsize    = normal*grid_height;

% Grid size (Nx must be an odd number)
aspect_ratio = maxsize/width;
height = ceil(aspect_ratio);

% Step size
dx = maxsize/Nx;
dy = maxsize*height/Ny;

% Assert that grid sampling is fine enough for wavelength
err_msg = ['Error: poor grid sampling. Either increase wavelength or ' ...
    'increase grid resolution.'];
assert(max([dx,dy]) <= wavelength/20, err_msg);
% Note: this check is important because poor sampling of the grid with
% respect to the incident wavelength can lead to inaccurate results. You
% may, if you wish, relax the condition to < wavelength/10 if it becomes
% impractical to simulate your problem based on available hardware, but it
% is best to keep dx and dy at least below 10% of the incident wavelength

% Axes
xa = [0:Nx-1]*dx;
ya = [0:Ny-1]*dy;
Nx2 = 2*Nx;
Ny2 = 2*Ny;

% Create wavevector
k0    = 2*pi/wavelength;
theta = angle_incidence*pi/180;
k     = k0*[sin(theta); cos(theta)];

% Calculate refractive index in reflection and transmission regions
er_ref = mean(ER(:,1));
ur_ref = mean(UR(:,1));
er_trn = mean(ER(:,Ny2));
ur_trn = mean(UR(:,Ny2)); 
nref = sqrt(er_ref*ur_ref);
ntrn = sqrt(er_trn*ur_trn);
plotER = real(sqrt(ER.*UR));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Absorbing Boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PML parameters
N0   = 376.73032165;    %free space impedance
amax = 5;
cmax = 1;
p    = 5;

% Initialize PML to problem space
sx = ones(Nx2,Ny2);
sy = ones(Nx2,Ny2);

% Construct PML in y-direction only
N = 2*NPML;
for n = 1 : N
    % compute PML value
    ay = 1 + amax*(n/N)^p;
    cy = cmax*sin(0.5*pi*n/N)^2;
    s  = ay*(1-1i*cy*N0/k0);
    
    % incorporate value into PML
    sy(:,N-n+1)   = s;
    sy(:,Ny2-N+n) = s;
end

% Incorporate PML into tensor components of material
ER2xx = ER ./ sx .* sy;
ER2yy = ER .* sx ./ sy;
ER2zz = ER .* sx .* sy;
UR2xx = UR ./ sx .* sy;
UR2yy = UR .* sx ./ sy;
UR2zz = UR .* sx .* sy;

% Overlay materials back into 1X grid
ERxx = ER2xx(2:2:Nx2,1:2:Ny2);
ERyy = ER2yy(1:2:Nx2,2:2:Ny2);
ERzz = ER2zz(1:2:Nx2,1:2:Ny2);
URxx = UR2xx(1:2:Nx2,2:2:Ny2);
URyy = UR2yy(2:2:Nx2,1:2:Ny2);
URzz = UR2zz(2:2:Nx2,2:2:Ny2);

clear N0 amax cmax p sx sy n N ay cy s;
clear ER2xx ER2yy ER2zz UR2xx UR2yy UR2zz;

% Form diagonal matrices
ERxx = diag(sparse(ERxx(:)));
ERyy = diag(sparse(ERyy(:)));
ERzz = diag(sparse(ERzz(:)));
URxx = diag(sparse(URxx(:)));
URyy = diag(sparse(URyy(:)));
URzz = diag(sparse(URzz(:)));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative matrices dimensions
M = Nx*Ny;

% DEX
DEX = sparse(M,M);
DEX = spdiags(-ones(M,1),0,DEX);
DEX = spdiags(+ones(M,1),+1,DEX);
for ny = 1 : Ny-1
    neq = Nx*(ny-1) + Nx;
    DEX(neq,neq+1) = 0;
end

% periodic BCs
for ny = 1 : Ny
    neq = Nx*(ny-1) + Nx;
    nv  = Nx*(ny-1) + 1;
    DEX(neq,nv) = +exp(-1i*k(1)*Nx*dx);
end
DEX = DEX / (k0*dx);

% DEY
DEY = sparse(M,M);
DEY = spdiags(-ones(M,1),0,DEY);
DEY = spdiags(+ones(M,1),+Nx,DEY);
DEY = DEY / (k0*dy);

% DHX
DHX = sparse(M,M);
DHX = spdiags(+ones(M,1),0,DHX);
DHX = spdiags(-ones(M,1),-1,DHX);
for ny = 2 : Ny
    neq = Nx*(ny-1) + 1;
    DHX(neq,neq-1) = 0;
end
for ny = 1 : Ny
   neq = Nx*(ny-1) + 1;
   nv  = Nx*(ny-1) + Nx;
   DHX(neq,nv) = -exp(+1i*k(1)*Nx*dx);
end
DHX = DHX / (dx*k0);

% DHY
DHY = sparse(M,M);
DHY = spdiags(+ones(M,1),0,DHY);
DHY = spdiags(-ones(M,1),-Nx,DHY);
DHY = DHY / (dy*k0);

% Field matrix A
switch polarization
    case 'E'
        A = DHX/URyy*DEX + DHY/URxx*DEY + ERzz;
    case 'H'
        A = DEX/ERyy*DHX + DEY/ERxx*DHY + URzz;
    otherwise
        error('Unrecognized polarization.');
end
clear DHX DHY DEX DEY;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total Field / Scattered Field Interface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TF/SF Masking Matrix Q
Q = zeros(Nx,Ny);
Q(:,1:NPML+2) = 1;
Q = diag(sparse(Q(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implement source
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create source
[Y,X] = meshgrid(ya,xa);
source = exp(-1i*(k(1)*X + k(2)*Y));
plotsource = source;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FDFD Main Solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Incorporate TF/SF into source
source = source(:);
b      = (Q*A - A*Q)*source;
clear Q source;

% Calculate resulting field
Field = A\b;
clear A b;
Field = full(Field);
Field = reshape(Field,Nx,Ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate plots
switch polarization
    case 'E'
        fieldtype = 'Ez';
    case 'H'
        fieldtype = 'Hz';
end
title1 = ['Source ', fieldtype, ' amplitude'];
title2 = ['Scattered ', fieldtype, ' amplitude'];
title3 = ['Scattered ', fieldtype, ' phase'];
plots = {'plotER', 'real(plotsource)', 'real(Field)', 'angle(Field)'};
titles = {'Refractive index', title1, title2, title3};

% Rescaling factor for plots
rescaling = 1/scaling_factor/normal;
if scaling_factor == 1
    units = 'm';
elseif scaling_factor == 1e-3
    units = 'mm';
elseif scaling_factor == 1e-6
    units = '\mum';
elseif scaling_factor == 1e-9
    units = 'nm';
elseif scaling_factor == 1e-12
    units = 'pm';
elseif scaling_factor == 1e-15
    units = 'fm';
else
    error('Scaling factor not recognized');
end

% Generate plots
if plotsOn ~= 0
    figure('units','normalized','outerposition',[0 0 1 1]);
    for loop = 1:length(plots)
      subplot(1,4,loop);
      imagesc(xa*rescaling,ya*rescaling,eval(plots{loop})'); 
      xlabel(['x (' units ')']); ylabel(['y (' units ')']);
      colormap jet; colorbar; axis equal tight;
      title(titles{loop});
    end
end
clear plotER plotsource plots titles plotfield;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate transmittance and reflectance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create reflected and transmitted field vectors
reflected   = Field(:,NPML+1);
transmitted = Field(:,Ny-NPML);
clear Field;

% Adjust phase tilt
phase       = exp(+1i*k(1)*xa');
reflected   = reflected.*phase;
transmitted = transmitted.*phase;

% Compute spatial harmonics using FFT
reflected   = fftshift(fft(reflected))/Nx;
transmitted = fftshift(fft(transmitted))/Nx;

% Calculate wavevector components
m = -floor(Nx/2):floor(Nx/2);
m = m';
kx_m = k(1) - 2*pi*m/(Nx*dx);
kref = conj(sqrt((k0*nref)^2 - kx_m.^2));
ktrn = conj(sqrt((k0*ntrn)^2 - kx_m.^2));

% Calculate reflectance and transmittance
switch polarization
    case 'E'
        R =  abs(reflected).^2 .* real(kref/k(2));
        T =  abs(transmitted).^2 .* real(ktrn*ur_ref/k(2)/ur_trn);
    case 'H'
        R =  abs(reflected).^2 .* real(kref/k(2));
        T =  abs(transmitted).^2 .* real(ktrn*er_ref/k(2)/er_trn);
end
clear kref ktrn kx_m m reflected transmitted;

% Calculate total transmittance and reflectance
T = sum(T);
R = sum(R);
A = 1 - (R+T);
end