%% Calculate posterior cov of modal parameters by matrix calculus (Louis identity)
% Single setup and scalar Se
% Assemble all the functions before, can be applied to multiple bands
%========================================================================
% 2209014-Firstly written by Wei at ZJUI
% 221111-Consistent with the page-wise operation
% 221121-Generalized to multi-bands
% 230816-Fixed bug for in.iband
% 231031-Detrend tdata
%========================================================================
% Input
%========================================================================
% fnb: MPV of natural frequency [m,1] or [1,m]
% znb: MPV of damping ratio [m,1] or [1,m]
% PHInb: MPV of mode shape [n,m]
% Snb: MPV of modal force PSD [m,m]
% Senb: MPV of prediction error [m,1] or [1,m]
% in: Struct contains:
% # Mandatory fields: tdata, fs, f1f2, f0
% # Optional fields: nseg, iband, hess_alg, cons_alg, sing_alg
%   ## nseg: number of segment (1 by default)
%   ## iband: bands to calculate
%   ## hess_alg: way of calculating hessian matrix, can be 'par' and
%      'page' (default)
%      ### 'par' => use parllel function of MATLAB
%      ### 'page' => use page-wise operations
%   ## cons_alg: way of considering the mode shape norm constraints,
%      can be 'comp' and 'lag' (default)
%      ### 'lag' => use lagrange multiplier and pseudo inverse, only
%          applying to the optimum or the uncertainty can't be quantified in this way
%      ### 'comp' => use the concept of the composite function and pseudo inverse,
%          applying to all parameter values, not just the optimum
%      ### 'nullproj' => use nullspace projection method
%========================================================================
% Output
%========================================================================
% Struct contains 2 fields: coefv and time_req (computational time)
% coefv is a struct contains the following fields:
% f: Posterior cov of natural frequency [1,m]
% z: Posterior cov of damping ratio [1,m]
% PHI: Posterior cov of mode shape [1,m]
% S: Posterior cov of diagonal entries of modal force PSD [1,m]
% Se: Posterior cov of prediction error scalar [1,m]
%========================================================================
function oo = PostCOV_scalarSe(fnb,znb,PHInb,Snb,Senb,in)
%========================================================
% Defaults setting
%========================================================
if ~isfield(in,'nseg') || isempty(in.nseg)
    in.nseg = 1;
end
if isfield(in,'hess_alg') && isempty(in.hess_alg)
    in.hess_alg = 'page';
elseif ~isfield(in,'hess_alg') || ~contains(in.hess_alg,'par')
    in.hess_alg = 'page';
end
if isfield(in,'cons_alg') && isempty(in.cons_alg)
    in.cons_alg = 'nullproj';
elseif ~isfield(in,'cons_alg') || ~contains(in.cons_alg,'comp') && ~contains(in.cons_alg,'lag')
    in.cons_alg = 'nullproj';
end
%========================================================
% Check input (Refer to the 'commonproc' function)
%========================================================
% check modal parameters
m = length(fnb);
[r,c] = cellfun(@size,Snb);
rr = sum(r);   cc = sum(c);
if length(znb)~=m || rr~=m || cc~=m || size(PHInb,2)~=m || length(Senb)~=m
    error('The dimension of modal parameters must be consistent!');
end
% check mandatory fields from in
if ~isfield(in,'tdata')
    error('tdata is a required field!');
end
if ~isfield(in,'fs')
    error('fs is a required field!');
end
if ~isfield(in,'f0')
    error('f0 is a required field!');
end
if ~isfield(in,'f1f2')
    error('f1f2 is a required field!');
end
% check legitimacy
% fs
if length(in.fs(:))>1
    warning('Only fs(1) will be used.');
end
in.fs = in.fs(1);
% f1f2
if ~isnumeric(in.f1f2)
    error('f1f2 should be a numerical array!');
end
if size(in.f1f2,2)~=2
    error('f1f2 should have 2 columns!')
end
% f0
if ~iscell(in.f0)
    error('f0 should be a cell!');
end
nb = size(in.f1f2,1);
if ~isequal(size(in.f0),[nb 1])
    error('The 1st dim. of f0 and f1f2 should be equal!');
end
% nseg
if rem(in.nseg,1) || (in.nseg<1) || (in.nseg>size(in.tdata,1))
    error('nseg must be an integer betwen 1 and first dim. of tdata');
end
M = cellfun(@length,in.f0);   % number of modes in every band [nb,1]
modind(:,2) = cumsum(M);
modind(:,1) = modind(:,2)-M+1;   % mode index
% iband
if isfield(in,'iband')
    if ~isempty(in.iband)
        in.iband = round(in.iband); % force to integer
        if any(in.iband(:)<1 | in.iband(:)>size(in.f1f2,1))
            tmpstr = 'Each entry in iband(:) must be an integer ';
            tmpstr = [tmpstr,'between 1 and 1st dim. of f1f2.'];
            error(tmpstr);
        end
        in.f0 = in.f0(in.iband);
        in.f1f2 = in.f1f2(in.iband,:);
        M = M(in.iband);
        modind = modind(in.iband,:);
        nb = length(in.iband);   % no. of freq. bands
    end
end
%========================================================
% Initialization
%========================================================
n = size(PHInb,1);   % number of dofs
Ntheta = (M+1).^2 + M*n;   % number of parameters in every band [nb,1]
%========================================================
% Calculate posterior cov band-by-band
%========================================================
in.tdata = detrend(in.tdata);
nt = size(in.tdata,1);
nf0 = floor((nt+1)/2); % no. of ffts up to Nyquist, will change later
ff0 = (0:nf0-1).'*in.fs/nt; % (nf0,1), ordinates of freq.
Fnb = sqrt(1/in.fs/nt)*fft(in.tdata);    % scaled FFT - two sided
Fnb = Fnb(2:nf0,:);
ffnb = ff0(2:nf0); % (:,1)
for ib = 1:nb
%     [F,ff,~] = proc_fft(in.tdata,in.f1f2(ib,:),in.fs,in.nseg);
    I1I2 = round(in.f1f2(ib,:)*nt/in.fs);
    II = I1I2(1):I1I2(2);
    F = Fnb(II,:); % (nf0,n,nseg)
    ff = ffnb(II); % (nf0,1)

    m = M(ib);   ntheta = Ntheta(ib);
    modindib = modind(ib,1):modind(ib,2);
    f = fnb(modindib);   z = znb(modindib);   PHI = PHInb(:,modindib);
    Se = Senb(modindib(1));
    if isfield(in,'iband')
        S = Snb{in.iband(ib)};
    else
        S = Snb{ib};
    end
    if contains(in.hess_alg,'par')
        if ib == 1
            flg = input('Parllel pool is ready (0/1)?\n');
        end
        %============================================================
        % Use codes to start the parpool, proving to be low-efficient
        %============================================================
        %         if ib == 1
        %             delete(gcp("nocreate"));   % stop the process before starting a new run
        %             numCore = feature('numcores');   % get the maximum core num of PC
        %             parpool(numCore);   % start parpool
        %         end
        if ~flg
            tmpstr1 = 'Please press the ''Start Parallel Pool'' button';
            tmpstr = [tmpstr1,' in the lower left quater!\n'];
            error(tmpstr);
        else
            tic
            [g,H] = NLLFHess_scalarSe3(f,z,PHI,S,Se,ff,F);
        end
    else   % page-wise operation
        tic
        [g,H] = NLLFHess_scalarSe(f,z,PHI,S,Se,ff,F);
    end
    %==========================================================
    % Bypassing singularity due to mode shape norms
    %===========================================================
    if contains(in.cons_alg,'nullproj')   % null space projection
        % Eliminate the effect of the modal parameters units
        T = diag([f;z;ones(m*n,1);real(diag(S));real(vec(S,'hh'));imag(vec(S,'hh'));Se]);
        H = T*H*T;
        %         % Nullspace of the jacobian matrix of the constraints, [ntheta,ntheta-m]
        %         IPHI = 2*m+1:2*m+m*n;
        %         Jmat = zeros(m,ntheta);
        %         %         Jmat = spalloc(m,ntheta,m*n);
        %         Jmat(:,IPHI) = vec(PHI).'.*repelem(eye(m),1,n);
        %         U = null(Jmat);
        U = cstgradnull(PHI);
        PCM = U/(U.'*H*U)*U.';
    elseif contains(in.cons_alg,'lag')   % Lagrange multilpier and pseudoinverse
        % Mapping the parameters to a set of free variables
        %==========================================================
        % Full form, ntheta won't be very large in single setup
        %===========================================================
        % Gradient of the mapping function
        IPHI = 2*m+1:2*m+m*n;
        gu = eye(ntheta);
        gu(IPHI,IPHI) = speye(m*n)-PHI(:)*PHI(:).'.*repelem(speye(m),n,n);
        sumcons = ...
            sparse(IPHI,IPHI,repelem(diag(reshape(g(IPHI),n,m).'*PHI),n),ntheta,ntheta);
        % Hessian of the composite function
        Hc = gu.'*(H+sumcons)*gu;
        %==========================================================
        % Sparse form, maybe can save space a lot in multiple setup
        %===========================================================
        %     % Gradient of the mapping function
        %     Ixx = [1:2*m,2*m+m*n:ntheta];
        %     gu = speye(ntheta);
        %     gu(IPHI,IPHI) = speye(m*n)-PHI(:)*PHI(:).'.*repelem(speye(m),n,n);
        %     ii = [Ixx,repmat(IPHI,1,m*n)];
        %     jj = [Ixx,repelem(IPHI,m*n)];
        %     vv = [ones(ntheta-m*n,1);vec(speye(m*n)-PHI(:)*PHI(:).'.*repelem(speye(m),n,n))];
        %     gu = sparse(ii,jj,vv,ntheta,ntheta);
        %     % Sum of parts related to the constraints
        %     sumcons = sparse(IPHI,IPHI,repelem(diag(reshape(g(IPHI),n,m).'*PHI),n),ntheta,ntheta);
        % Eliminate the effect of the modal parameters units
        T = diag([f;z;ones(m*n,1);diag(S);real(vec(S,'hh'));imag(vec(S,'hh'));Se]);
        Hc = T*Hc*T;
        % pseudoinverse
        [V,D] = eig(Hc);
        D = diag(D);
        [D,I] = sort(D);
        V = V(:,I);
        Nnpos = length(find(D<=0));
        I = max(Nnpos,m)+1:ntheta;
        PCM = gu*V(:,I)*diag(1./D(I))*V(:,I).'*gu.';
    else   % Composite function and pseudoinverse
        gu = eye(ntheta);
        Hu = zeros(ntheta^2,ntheta);
        PHIHu = zeros(n^2,ntheta);
        for ii = 1:m
            IPHIii = 2*m+1+(ii-1)*n:2*m+ii*n;
            phi = PHI(:,ii);
            gu(IPHIii,IPHIii) = eye(n)-PHI(:,ii)*PHI(:,ii).';
            [Er1,Ec1] = pick_vec(ntheta,ntheta,IPHIii,IPHIii);
            PHIHu(:,IPHIii) = kron(3*(phi*phi.')-eye(n),phi)...
                - kron(phi,eye(n)) - vec(eye(n))*phi.';
            Hu = Hu + kron(Ec1,Er1.')*PHIHu;
        end
        % Hessian of the composite function
        Hc = gu.'*H*gu + kron(eye(ntheta),g)*Hu;
        % Eliminate the effect of the modal parameters units
        T = diag([f;z;ones(m*n,1);diag(S);real(vec(S,'hh'));imag(vec(S,'hh'));Se]);
        Hc = T*Hc*T;
        % pseudoinverse
        [V,D] = eig(Hc);
        D = diag(D);
        [D,I] = sort(D);
        V = V(:,I);
        Nnpos = length(find(D<=0));
        I = max(Nnpos,m)+1:ntheta;
        PCM = gu*V(:,I)*diag(1./D(I))*V(:,I).'*gu.';
    end
    oo.time_req(modindib) = toc;
    %========================================================
    % Form the posterior COV of the modal parameters
    %========================================================
    diagPCM = real(diag(PCM));
    oo.coefv.f(modindib) = sqrt(diagPCM(1:m))';
    oo.coefv.z(modindib) = sqrt(diagPCM(m+1:2*m))';
    oo.coefv.Sii(modindib) = sqrt(diagPCM(2*m+m*n+1:2*m+m*n+m)');
    oo.coefv.Se(modindib) = sqrt(diagPCM(end));
    oo.coefv.phi(modindib) = sqrt(sum(reshape(diagPCM(2*m+1:2*m+m*n)',n,m)));
end
delzeros = @(a) a(a~=0);
oo.coefv = structfun(delzeros,oo.coefv,'UniformOutput',false);
end