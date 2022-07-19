%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finite Difference Frequency Domain 2D Solver for a periodic multilayer
% stack consisting of two alternating materials with corresponding
% permittivities epsilon1 and epsilon2, surrounded by free-space. The stack
% is made such that epsilon1 is the 1st layer and the last, so it has N+1
% layers, where N = 2*M and M is the number of periods or repetition of the
% two materials.
%
% The program also admits the inclusion of surface roughness via the RMS
% roughness parameters rms1 and rms2, and then generates Gaussian smoothed
% random interfaces between each set of layers using the random1Dsurface.m
% script from my GitHub page.
% 
% Program by: Oscar A. Nieves
% Updated: 19/07/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clc; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Units and constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Length units
meters      = 1;
centimeters = meters*1e-2;
millimeters = meters*1e-3;
micrometers = meters*1e-6;
nanometers  = meters*1e-9;
picometers  = meters*1e-12;

% time and frequency
seconds   = 1;
hertz     = 1/seconds;
kilohertz = hertz*1e+3;
megahertz = hertz*1e+6;
gigahertz = hertz*1e+9;
terahertz = hertz*1e+12;

% constants
c0    = 299792458*meters/seconds;      % speed of light in vacuum
e0    = 8.8541878176e-12*1/meters;     % permittivity of free space
u0    = 1.2566370614e-6*1/meters;      % permeability of free space
e_air = 1.0;
u_air = 1.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source parameters
wavelengths     = nanometers*linspace(450,840,10); % Wavelength of planar wave source
grid_width      = 500*nanometers;       % Maximum physical size of grid (x-axis)
rms1            = 00*nanometers;        % RMS roughness of interface 1-2
rms2            = rms1;                 % RMS roughness of interface 2-1
thickness1      = 100*nanometers;       % thickness of layer 1
thickness2      = 100*nanometers;       % thickness of layer 2
angle_incidence = 00;                   % Angle of incidence to device (degrees)
stack_periods   = 5;                    % number of layer-pairs in the stack

% Calculate maximum size of grid based on thickness of layers
period_length = thickness1 + thickness2;
stack_height  = stack_periods*period_length + thickness1;
grid_height   = 1.3*stack_height; % Maximum physical size of grid (y-axis)

% Note: maxsize should be at least 20% larger than the size of the device,
% which in this case is stack_length; otherwise you will get part of the 
% device overlapping with the PML and this will lead to wrong results like 
% very high absorbance when using dielectric materials only, which is not 
% a physical effect. You can do an easy check by using a purely dielectric 
% material like glass and seeing if your final absorbance value is 
% significant. In general, due to numerical error, you would expect 
% absorbance to be somewhere around 0.1% or less for dielectrics, but if 
% you end up with absorbance of +-1%, then most likely your device is 
% leaking into the PML and causing trouble.

% Material properties 
epsilon_medium = 1.0;      % Relative permittivity surrounding medium
mu_medium      = 1.0;      % Relative permeability surrounding medium
epsilon1       = 12.4;     % Relative permittivity of layer 1
epsilon2       = 6.8;      % Relative permittivity of layer 2 
mu_device      = 1.0;      % Relative permeability of device

% Grid size (Nx must be odd)
aspect_ratio = grid_height/grid_width;
height = ceil(aspect_ratio);
Nx     = 301;
Ny     = height*Nx;

% Polarization (type either 'E' or 'H' to display Electric or Magnetic field)
polarization = 'E';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid construction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust epsilons
epsilon1 = real(epsilon1) - abs(imag(epsilon1))*1i;
epsilon2 = real(epsilon2) - abs(imag(epsilon2))*1i;
% Note: for some reason, when the imaginary part of the permittivity is 
% positive, we get nonsense results for reflectance, transmittance, etc.
% for example transmittance 290% and absorbance -300%. This adjustment is
% done to prevent that from happening.

% Size of PML (Absorbing boundary)
empty_grid = ceil(Ny*(1 - stack_height/grid_height));
NPML = ceil(empty_grid/4);

% step size
dx = grid_width/Nx;
dy = grid_width*height/Ny;

% Axes
xa = [0:Nx-1]*dx;
ya = [0:Ny-1]*dy;

% 2X grid
Nx2 = 2*Nx;
Ny2 = 2*Ny;
dx2 = dx/2;
dy2 = dy/2;
xa2 = [0:Nx2-1]*dx2;
ya2 = [0:Ny2-1]*dy2;
[Y2,X2] = meshgrid(ya2,xa2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Device
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Device parameters
x0 = ceil(Nx2/2);
y0 = ceil(Ny2/4);
centering = max(ya2)/2;
rms_vec = repmat([rms1, rms2], [1,stack_periods+1]);
eps_vec = repmat([epsilon1, epsilon2], [1,stack_periods + 1]);
thicknesses = repmat([thickness1, thickness2], [1,stack_periods + 1]);
for nn = 1:stack_periods*2+2
    if rms_vec(nn) == 0.0 % in case random1surface returns NaNs
        rough{nn} = zeros(1,Nx2);
    else
        [~, rough{nn}, ~] = random1Dsurface(grid_width, grid_width/20, ...
            rms_vec(nn), Nx2);
    end
end
total_length = stack_periods * (thickness1 + thickness2) + thickness1;
top_layer = centering + total_length/2 + rough{1};
prev_layer = top_layer;
top_line = centering + total_length/2;
bottom_line = centering - total_length/2;
device = epsilon_medium * ones(Nx2, Ny2);
for nn = 1:stack_periods*2 + 1
    top_curve = prev_layer;
    bottom_curve = mean(top_curve) - thicknesses(nn) + rough{nn+1};
    prev_layer = bottom_curve;
    refractive_curves{nn} = top_curve;
    device(Y2>=bottom_curve.' & Y2<=top_curve.') = eps_vec(nn);
end
refractive_curves{stack_periods*2 + 2} = bottom_curve;
device0 = sqrt(device(1:2:end,1:2:end)); % convert to refractive index

% Overlay material properties
ER = device;
UR = ones(Nx2, Ny2);
clear device;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call FDFD Solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_TE = zeros(1,length(wavelengths));
T_TE = zeros(1,length(wavelengths));
A_TE = zeros(1,length(wavelengths));
R_TM = zeros(1,length(wavelengths));
T_TM = zeros(1,length(wavelengths));
A_TM = zeros(1,length(wavelengths));

tic;
for nn = 1:length(wavelengths)
    wavelength = wavelengths(nn);
    [R_TE(nn),T_TE(nn),A_TE(nn)] = FDFD(ER,UR,grid_width,grid_height,NPML,...
        wavelength,angle_incidence,'E',nanometers);
    [R_TM(nn),T_TM(nn),A_TM(nn)] = FDFD(ER,UR,grid_width,grid_height,NPML,...
        wavelength,angle_incidence,'H',nanometers);
    disp(['loop = ' num2str(nn) '/' num2str(length(wavelengths))]);
end
time0 = toc;
disp(['Time elapsed = ' num2str(time0/60) ' minutes']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Spectral response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color','white');
LW = 2;
wavelengths = wavelengths*1e9; 

subplot(2,2,1);
plot(wavelengths,R_TE,'color','blue','linewidth',LW);
xlabel('Wavelength (nm)');
ylabel('Reflectance TE');
xlim([min(wavelengths), max(wavelengths)]);
ylim([0, 1]);

subplot(2,2,3);
plot(wavelengths,T_TE,'color','black','linewidth',LW);
xlabel('Wavelength (nm)');
ylabel('Transmittance TE');
xlim([min(wavelengths), max(wavelengths)]);
ylim([0, 1]);

subplot(2,2,2);
plot(wavelengths,R_TM,'color','blue','linewidth',LW);
xlabel('Wavelength (nm)');
ylabel('Reflectance TM');
xlim([min(wavelengths), max(wavelengths)]);
ylim([0, 1]);

subplot(2,2,4);
plot(wavelengths,T_TM,'color','black','linewidth',LW);
xlabel('Wavelength (nm)');
ylabel('Transmittance TM');
xlim([min(wavelengths), max(wavelengths)]);
ylim([0, 1]);