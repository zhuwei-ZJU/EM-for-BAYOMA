# EM-for-BAYOMA
Bayesian operational modal analysis based on the expectation-maximization algorithm.
An expectation-maximization is proposed to speed up the original Bayesian FFT method for operational modal analysis.
Fast computations for both the most probable value (MPV) and the posterior covariance matrix (PCM) are invloved.
Before running the algorithm, the frequency bands should be manually selected by conbining the information in the spectra.
The significant peak in the power spectrum density (PSD) tells where the structural mode possibly locates.
The singular value (SV) spectrum consists of eigenvalues of the PSD matrix at each frequency point.
In a resonance band, the number of lines significantly above others reveals the number of modes, and other lines indicate the noise.

An example can be found as below:

clear all; close all; clc
addpath(genpath(pwd));  % setpath

%% true values
f = [0.98 1 1.02];   % frequencies, Hz
z = [0.8 1 1.2]/100;  % damping ratios, 1
phi = [1 2 2;2 1 -2;1 -2 2]'/3;    % mode shapes
S = blkdiag([1 0.5*exp(1i*pi/4); 0.5*exp(-1i*pi/4) 1],1);  % modal force PSD, \mug^2/Hz
Se = 100; % channel noise PSD \mug^2/Hz


%% identification - bayoma
in = load('modes3.mat');
% in.tdata contains measured data;
% in.fs is the sampling frequency, Hz;
in.f0 = {[0.98 1.0 1.02]}; % initial guess of frequencies
% if there are multiple bands, in.f0 = {[f01, f02]; [f03 f04]; ...}
in.f1f2 = [0.85 1.15]; % selected frequency bands
% if there are multiple bands, in.f1f2 =[lower1 upper1; lower2 upper2; ...]

% IN contains the following mandatory fields:
%   tdata = (nt,n) ambient data; nt = no. of time points; n = no. of dofs
%   fs = scalar, sampling rate (Hz)
%   f1f2 = (nb,2), f1f2(:,1) and f1f2(:,2) gives the lower and upper bound 
%          of frequency bands used for modal identification, nb = no. of
%          bands
%   f0 = {nb,1}, cell array of initial guesses of natural frequencies
%
% The following optional fields can be supplied in IN:
%   maxiter = max. no. of iterations in determining MPV
%   p = power, fft data is multiplied by (2*pi*f)^p.
%       This power converts the data to acceleration data assumed in the
%       theory of the program. 
%       E.g., p=0 for acceleration data (default)
%             p=1 converts velocity data to acceleration data
%             p=2 converts displacment data to acceleration data
%     The measurementn error of data is assumed to have constant PSD in each
%     selected band.


out = em4bayoma(in);
% OUT is a structure with the following fields:
% 1. MPV of modal parameters:
%  f = (1,m) natural frequencies (Hz)
%  z = (1,m) damping ratios
%  phi = (n,m) mode shape, each column normalized to unity
%  Se = scalar, prediction error variance
%  S = (m,m) PSD matrix of modal force
%
% Note: If tdata is given in g and fs in Hz, then Se and S are two-sided
% with a unit of g^2/Hz
% 
% 2. Posterior uncertainty
%  coefv = cell containing blocks of c.o.v. for 
%    f(:),z(:),S(:),Se,phi(:,1),phi(:,2),...,phi(:,m),
%  Note that off-diagonal terms of c.o.v. of S is the c.o.v. of coherence

