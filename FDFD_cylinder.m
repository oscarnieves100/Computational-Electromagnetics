%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finite Difference Frequency Domain 2D Cylinder solver
%
% Program by: Oscar A. Nieves
% Updated: 19/07/2022
%
% This script calls the function FDFD.m to solve for the reflectance and
% transmittance properties of a cylinder surrounded by free-space.
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
wavelength      = 500*nanometers;   % Wavelength of planar wave source
grid_width      = 600*nanometers;   % Maximum physical size of grid (x-axis)
grid_height     = 1000*nanometers;  % Maximum physical size of grid (y-axis)
diameter        = 400*nanometers;   % diameter of cylinder
angle_incidence = 0;                % Angle of incidence to device (degrees)

% Material properties 
epsilon_medium = 1.0;      % Relative permittivity surrounding medium
mu_medium      = 1.0;      % Relative permeability surrounding medium
epsilon_device = 6.0;      % Relative permittivity of device
mu_device      = 1.0;      % Relative permeability of device

% Grid size (Nx must be an odd number)
aspect_ratio = grid_height/grid_width;
height = ceil(aspect_ratio);
Nx     = 301;
Ny     = height*Nx;

% Polarization (type either 'E' or 'H' to display Electric or Magnetic field)
polarization = 'E';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid construction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Size of PML (Absorbing boundary)
NPML = ceil(Ny/10);

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
x0 = mean(xa2);
y0 = mean(ya2);
device = ((X2-x0).^2 + (Y2-y0).^2) <= (diameter/2)^2;

% Overlay material properties
ER = epsilon_medium + (epsilon_device - epsilon_medium)*device;
UR = mu_medium + (mu_device - mu_medium)*device;
clear device;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call FDFD Solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[R,T,A] = FDFD(ER,UR,grid_width,grid_height,NPML,wavelength, ...
    angle_incidence,polarization,nanometers,1);

disp(['Reflectance = ' num2str(R)]);
disp(['Transmittance = ' num2str(T)]);
disp(['Absorbance = ' num2str(A)]);
disp(['Total energy = ' num2str(R+A+T)]);