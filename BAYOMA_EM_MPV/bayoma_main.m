function out = bayoma_main(in)
% Fast Bayesian FFT modal identification using EM algorithm
%--------------------------------------------------------------------------
% out = bayoma_main(in)
% perform Bayesian modal identification using FFT on a specified frequency 
% band limited by f1f2.
%
% If you use this program, please cite
%   1. Li, Binbin; Au, Siu-Kui. An expectation-maximization algorithm for 
%       Bayesian operational modal analysis with multiple (possibly close)
%       modes. Mechanical Systems and Signal Processing, 2019,132:490-511.
%   2. Zhu, Wei and Li, Binbin, EM-Aided Fast Posterior Covariance Computation
%       in Bayesian FFT Method. Mechanical Systems and Signal Processing.
% IN contains the following mandatory fields:
%   tdata = (nt,n) ambient data; nt = no. of time points; n = no. of dofs
%   fs = scalar, sampling rate (Hz)
%   f1f2 = (nb,2), f1f2(:,1) and f1f2(:,2) gives the lower and upper bound 
%          of frequency bands used for modal identification
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
%     The measurementn error of data is assumed to have constant PSD in the
%     selected band.
%
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
%
%  cmat = cell of covariance matrix. Except for phi, all other entries are
%  wrt multipliers (w/ MPV=1). This implies that the diagonal entries of
%  cmat gives directly the sq. c.o.v.
%  The covariance matrices of different bands are arranged row-wise as
%   [f,z,diag of S, Re of off diag of S, Im of off diag of S, Se, phi] of 
%   band 1, [...] of band 2, [...] of band 3, etc
%
% Note: 
% The EM algorithm may be difficult to converge to true parameter values when
% the modal signal-to-noise ratio S/(4z^2*Se) is greater than 1e5. A warning 
% will be issued if it happens. 

% copyright statement
disp('============================================');
disp('BAYOMA_EM - BAYesian Operational Modal Analysis Using EM Algorithm');
disp('-- Binbin Li (bbl@zju.edu.cn) - Nov 1, 2019');
disp('-- Wei Zhu (weiz.20@intl.zju.edu.cn) - Feb 10, 2024');
disp('============================================');

%% defaults
in_default.maxiter = 500;  % max. no. of iterations in determining MPV
in_default.tol_cvg = 1e-5;  % convergence criterion
in_default.tol_hess = 1e-4;
in_default.z0 = 1e-2; % initial guess of damping
in_default.zlim = [1e-4 1]; % range of z outside which result will be nan
in_default.snmin = 1; % s/n ration below which result will be nan
in_default.options_fmins = optimset('TolX',in_default.tol_cvg/10,...
                   'TolFun',in_default.tol_cvg/10);
in_default.p = 0;   % acceleration data

%------------------------------------------------------------------------
in = optfield(in,in_default);

%% ==========================================================================
% check mandatory fields from IN
if ~isfield(in,'tdata')
  error('tdata is required field.');
end
if ~isfield(in,'fs')
  error('fs is required field');
end
if ~isfield(in,'f0')
  error('f0 is required field.');
end
if ~isfield(in,'f1f2')
  error('f1f2 is required field.');
end

% fs
%=====
if length(in.fs(:))>1
  warning('Only fs(1) will be used.');
end
in.fs = in.fs(1);

% f1f2
%======
if ~isnumeric(in.f1f2)
  error('f1f2 should be a numerical array.');
end
if size(in.f1f2,2)~=2
  error('f1f2 should have 2 columns')
end
nb = size(in.f1f2,1); % no. of freq. bands

% f0
%======
if ~iscell(in.f0)
  error('f0 should be a cell.');
end
if ~isequal(size(in.f0),[nb 1])
  error('The 1st dim. of f0 and f1f2 should be equal.');
end


%% ==========================================================================  
% prepare the data
in.tdata = detrend(in.tdata); % detrend
nt = size(in.tdata,1);

nf0 = floor((nt+1)/2); % no. of ffts up to Nyquist, will change later
ff0 = (0:nf0-1).'*in.fs/nt; % (nf0,1), ordinates of freq.

F = sqrt(1/in.fs/nt)*fft(in.tdata);    % scaled FFT - two sided 

in.nt = nt;
in.ff = ff0(2:nf0); % (:,1)
in.F = F(2:nf0,:); % (:,n) skip zero freq., up to Nyquist

in0 = in;
in = rmfield(in,'tdata'); % clear tdata to save space

% run through different bands
o = cell(1,nb); % store results for each band
% in0 = in;
for ii = 1:nb
  in = in0;
  in.f0 = in0.f0(ii);
  in.f1f2 = in0.f1f2(ii,:);
  o{ii} = bayoma_em(in);
end

% combine fields in different bands
out = o{1};
varx = fieldnames(o{1}).'; % cell array of field names
for ii=2:nb
  for xx = varx
    if isstruct(out.(xx{:}))
      vary = fieldnames(out.(xx{:})).';
      for yy=vary
        out.(xx{:}).(yy{:}) = [out.(xx{:}).(yy{:}),o{ii}.(xx{:}).(yy{:})];
      end
    else
      out.(xx{:}) = [out.(xx{:}),o{ii}.(xx{:})];
    end
  end
end

% Posterior covariance
fprintf('Computing Hessian for posterior covariance matrix ... ');
fprintf('by analytical formula ... ');
fprintf('\n');
outPCM = PostCOV_scalarSe(out.f.',out.z.',out.phi,out.S,out.Se,in0);
out.time_req(2,:) = outPCM.time_req;
out.coefv = outPCM.coefv;


% arrange field orders
% varx = {'f','z','Sii','Se','phi','rms','sn','coefv',...
%         'S','time_req','cmat','loglik'};
varx = {'f','z','Sii','Se','phi','rms','sn','coefv',...
        'S','time_req','loglik'};
out = orderfields(out,varx);

for m = 1:length(out.f)
    if out.sn(m) > 4e4
        warnmsg = ['Identification results for the frequency = ', num2str(out.f(m)), ...
            ' Hz is not reliable, consider other methods for double-check'];
        warning(warnmsg);
        fprintf('\n');
    end
end