function [xdomain, smoothgaussian, H] = random1Dsurface(length0, ...
                                corrlength, rms_height, Nx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1D spatially smooth Gaussian random surface generator.
% 
% Program by: Oscar A. Nieves
% Updated: 19/07/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate spatially smooth Gaussian surface
if nargin==0
    length0 = 100e-9;
    corrlength = 7e-9;
    rms_height = 10e-9;
    plots = 1;
    factor = 50;
    N0 = ceil( log10(factor*length0/corrlength)/log10(2) );
    N = 2^N0;
else
    N = Nx;
    plots = 0;
end

% Length and correlation length parameters
L = length0;
Lc = corrlength;

% Generate surface and normalize height to match RMS
xdomain = linspace(-L/2,L/2,N);
H = rms_height * randn(1,N);
G = exp(-2*(xdomain.^2)/Lc^2);
f = 2*L/sqrt(pi)/Lc/N;
smoothgaussian = f * real( ifft( fft(H).*fft(G) ));
smoothgaussian = smoothgaussian - mean(smoothgaussian);
rms_smooth = sqrt( mean(smoothgaussian.^2) );
amplifier = rms_height./rms_smooth;
smoothgaussian = amplifier*smoothgaussian;

if plots ~= 0
    close all;
    plot(xdomain, smoothgaussian);
end
end