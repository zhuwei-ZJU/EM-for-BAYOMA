% a 3 closely-spaced modes
% other synthetic, laboratory and field data can be found in:
% https://dataverse.harvard.edu/dataverse.xhtml?alias=harvard&q=Bayoma

%% synthetic model
f = [0.98 1 1.02];   % frequencies, Hz
z = [0.8 1 1.2]/100;  % damping ratios, 1
phi = [1 2 2;2 1 -2;1 -2 2]'/3;    % mode shapes
S = blkdiag([1 0.5*exp(1i*pi/4); 0.5*exp(-1i*pi/4) 1],1);  % modal force PSD, \mug^2/Hz
Se = 100; % channel noise PSD \mug^2/Hz

fs = 100;   % sampling frequency, Hz
t0 = 500;   % initial burn-in time for stationary state
T = t0 + 5000;  % sampling period
nt = T*fs;  % #samples
t = 0:1/fs:(T-1/fs);    % time instances

%% identification - bayoma
clear,clc
in = load('modes3.mat');
% in.m = 1; in.tdata = tdata; in.fs = fs;
in.f0 = {[0.98 1.0 1.02]};
in.f1f2 = [0.85 1.15];

in.tol_cvg = 1e-3;
in.alg = 'P-EM';
out1 = bayoma_main(in);

