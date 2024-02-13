function out = bayoma_em(in)
% Fast Bayesian FFT modal identification using EM algorithm
%
% 9/29/2018 - written by Binbin Li
%

tic; % start counting computational time

% extract fields for convenience
fs = in.fs; f0 = in.f0; f1f2 = in.f1f2;
z0 = in.z0; maxiter = in.maxiter; %tol_svd = in.tol_svd;
tol_cvg = in.tol_cvg; %tol_hess = in.tol_hess; 
zlim = in.zlim; snmin = in.snmin; nt = in.nt;

% for single band, convert f0 to a vector
f0 = cell2mat(f0); % (no. of modes in the band,1)
m = length(f0(:)); % # modes in selected band

fprintf([repmat('=',[1 60]),'\n']);

%% Preparing data
% calc freq. indices for selected band
I1I2 = round(f1f2*nt/fs);
II = I1I2(1):I1I2(2);

% confine to selected band
ff = in.ff(II); % (nf0,1)
F = in.F(II,:); % (nf0,n,nseg)
[nf,n] = size(F);

if m > n
    errmsg = 'This program can only be used to identify the case with more';
    errmsg = [errmsg newline 'sensors than the number of modes within the selected band'];
    error(errmsg);
end
    
% multiply by power
if in.p~=0
  disp(['FFT data is multiplied by (2*pi*f)^',num2str(in.p)]);
  F = F.*repmat((2*pi*ff).^in.p,[1,n]); % (nf,n)
end

%% Initial setting
% initial guess of freq. & damping
f = f0;
z = z0*ones(1,m);

time_req = nan(2,m);

% singular values at peak
sv0 = zeros(1,m);
phi = zeros(n,m);
for ii = 1:m
  I0 = round((f0(ii)-ff(1))*nt/fs)+1; % index corresponding to initial guess freq.
  dn0 = round(0.005*nt*f0(ii)/fs); % half bandwidth round f0 for taking average
  II = find([1:nf]<=I0+dn0 & [1:nf]>=I0-dn0);
  DD = real(F(II,:).'*conj(F(II,:))); % (n,n)  
  DD = DD/length(II);
  [vtmp,stmp] = svd(DD);
  phi(:,ii) = vtmp(:,1);
  sv0(ii) = stmp(ii,ii);
end

% arrange f in descending order of magnitude, to be consistent w/ BA1
[~,I] = sort(sv0,'descend');
f = f(I);

% calculate A0 & determine the modeshape space
%======================================================================
A0 = real(F.'*conj(F)); % (n,n)

[~,dA] = svd(A0);
dA = diag(dA); % (n,1), singular values in descending order
fprintf('Eigenvalues of A0:');
fprintf(['\n',repmat('%1.1e ',[1 10])],dA);
fprintf('\n');

mp = min(m,n); % modeshape space dim. mp
fprintf('Assume mode shape space dimension = %2i\n',mp);  
dA1 = dA(1:mp); 

% initialize S
ffT = ff';
betak = f'./ffT;
ihik = 1-betak.^2 - 1i*2*z'.*betak; % hik^-1
mf0i = F*phi/(phi'*phi).*ihik.';    % (hik.^-1*Qk).'
S = (mf0i'*mf0i).'/nf;
sqr_iS = chol(S,'lower')^-1;    % inverse of modal force PSD, S^-1/2, lower triangular matrix
 
    
% asymp. optimal value for predictio error variance = Se
traceA0 = trace(A0);
ResA0 = (traceA0 - sum(dA1)); % residual sum of eigenvaleus of A0

if n == mp % note that n>=mp always
    Se1 = eigs(real(F(1:10,:).'*conj(F(1:10,:)))/10,m);
    Se = min(Se1);
else
    Se = ResA0/nf/(n-mp);
end
sqr_Se = sqrt(Se);

%% ========================================================================
% optimize modal parameters sequentially by EM iteration

% array used for deciding convergence
aalast = [f(:).',z(:).',phi(:).',Se,diag(S).'];
fprintf('Convergence based on f(:),z(:),phi(:),Se,diag(S))\n');

iter = 0;
fprintf(['Iter.\tfreq (Hz) ',repmat(' ',1,10*(m-1))]);
fprintf(['damping   ',repmat(' ',1,10*(m-1))]);
fprintf('Se        ');
fprintf(['Sii       ',repmat(' ',1,10*(m-1))]);
fprintf('\n');
fprintf('%-3i\t\t',iter);
fprintf('%-6.2e ',[f(:).',z(:).',Se,diag(S).']);
fprintf('\n');


% calculate modal fft, used frequently later
d = trace(F*F');   % sum(Fk'*Fk)
exitflag = 1; % will turn to 0 if any fminsearch not converge
fSe_last = [f Se];

loglik = zeros(maxiter,1);

miniter = 10;
while iter < maxiter
    iter = iter + 1;
    
    % commonly used quantities
    betak = f'./ffT;
    ihk = 1-betak.^2 - 1i*2*z'.*betak; % hik^-1
    sqr_phi = chol(phi'*phi);
    rk = (F*phi/sqr_phi).';    % sqr_phi'^-1*phi'*Fk
    
    % E-step
    EQk = zeros(m,nf); % E[Qk]
    EQk2 = zeros(m,m,nf); % E[QkQk']
    S0 = zeros(m,m,nf);  % used in computing Si in M step
    logliki = zeros(1,nf);
    
    for k0 = 1:nf
        R = qr([sqr_phi rk(:,k0);sqr_Se*sqr_iS*diag(ihk(:,k0)) zeros(m,1)]);
        r12 = R(1:m,end); R11 = triu(R(1:m,1:m));
        Rr = R11\[sqr_Se*eye(m) r12];
        hRr = diag(ihk(:,k0))*Rr;
        EQk(:,k0) = Rr(:,end);  % iR11*r12;
        EQk2(:,:,k0) = Rr*Rr';
        S0(:,:,k0) = hRr*hRr';
        logliki(k0) = -2*sum(log(abs(diag(R11))))+r12'*r12/Se;
    end
    loglik(iter) = -n*nf*log(pi)-(n-m)*log(Se)-d/Se ...
        + 2*sum(sum(log(abs(ihk))))+2*nf*sum(log(abs(diag(sqr_iS)))) + sum(logliki);
    
    
    % M-step
    FkEQk = real(EQk*conj(F))'/nf;  % sum(Fk*EQk')
    EQk2mean = mean(real(EQk2),3);
    phi = FkEQk/EQk2mean; 
    Se = trace(A0/nf-phi*EQk2mean*phi')/n;
    
%     L = chol([EQk2mean FkEQk'; FkEQk A0/nf],'lower');
%     L11 = tril(L(1:m,1:m)); L21 = L(m+1:end,1:m);
%     L22 = tril(L(m+1:end,m+1:end));
%     phi = L21/L11;
%     V = trace(L22*L22')/n;
    sqr_Se = sqrt(Se);
    
    S = mean(S0,3); % modal force PSD
    if m == 1; S = real(S); end
    sqr_iS = chol(S,'lower')^-1;    % inverse of modal force PSD, S^-1/2, lower triangular matrix
    
    q2 = @(x)Q2(x,[f z],ffT,sqr_iS'*sqr_iS,EQk2);
    [x,~,tmp] = fminsearch(q2,ones(1,m*2));
    %         x = fminunc(q2,[fi zi],options);
    x = x.*[f z];
    f = x(1:m);   % freqncies
    z = x(m+1:end);   % damping ratios
    
    if tmp~=1
      exitflag = 0;
    end 
    
    % return to original parameters
    rnorm = sqrt(sum(phi.^2));
    phi = phi./rnorm;
    S = diag(rnorm)*S*diag(rnorm);
    sqr_iS = diag(1./rnorm)*sqr_iS;    % S^-1
    
    % display
    fprintf('%-3i\t\t',iter);
    fprintf('%-6.2e ',[f(:).',z(:).',Se,diag(S).']);
    fprintf('\n');

    
    % decide convergence    
    aa = [f(:).',z(:).',phi(:).',Se,diag(S).'];

    I = find(aalast);
    err = aa(I)./aalast(I) - 1; % relative error
    if all(abs(err)<tol_cvg) && (exitflag==1) && (iter >= miniter)
        break;
    else
        aalast = aa;
    end
    
    fSe = [f Se]; % preliminary convergence test
    if all(abs(fSe./fSe_last-1) < tol_cvg*10) && iter >= miniter && m > 1        
        break;
    else
        fSe_last = fSe;
    end
    
    
end

%% P-EM acceleration - only used for multi-mode
if m > 1    % multi-mode acceleration
    
    disp('P-EM algorithm begins ... ...');
    % initialization
    P0 = [f; z; phi; sqr_iS; sqr_Se*ones(1,m)];
    P1 = oneEMstep(F,ffT,P0,d,A0);
    P2 = oneEMstep(F,ffT,P1,d,A0);
    
    iter2 = 0;
    while iter < maxiter
        iter = iter + 1;
        iter2 = iter2 + 1;
        P_best = P2;
        L_best = oneloglik(F,ffT,P2,d);
        kk = 0; tt = 1+1.5^kk*0.1;
        P_new = (1-tt)^2*P0+2*tt*(1-tt)*P1+tt^2*P2;
        L_new = oneloglik(F,ffT,P_new,d);
        while L_new > L_best
            P_best = P_new; L_best = L_new;
            kk = kk + 1; tt = 1+1.5^kk*0.1;
            P_new = (1-tt)^2*P0+2*tt*(1-tt)*P1+tt^2*P2;
            L_new = oneloglik(F,ffT,P_new,d);
        end
        
        
        % convergence testing
        P_best = P2P(P_best,m,n);
        loglik(iter) = L_best;
        f = P_best(1,:); z = P_best(2,:); phi = P_best(3:2+n,:);    % frequency, damping, mode shapes
        sqr_iS = P_best(3+n:2+n+m,:);    % inverse of modal force PSD, S^-1/2
        sqr_S = sqr_iS^-1; S = sqr_S*sqr_S';     % modal force PSD, S
        sqr_Se = P_best(end,1);    % channel noise, Se^1/2
        Se = sqr_Se^2; % channel noise, Se
        
        % display
        fprintf('%-3i\t\t',iter);
        fprintf('%-6.2e ',[f(:).',z(:).',Se,diag(S).']);
        fprintf('\n');

        aa = [f(:).',z(:).',phi(:).',sqr_Se,diag(sqr_iS).'];
        I = find(aalast);
        err = aa(I)./aalast(I) - 1; % relative error
        if all(abs(err)<tol_cvg) && (exitflag==1) && (iter2 > 10)
            break;
        else
            aalast = aa;
        end
        
        % updating parameters
        P0 = P1;  P1 = P2;
        [P2,loglik1] = oneEMstep(F,ffT,P_best,d,A0);
        [P2,loglik2] = oneEMstep(F,ffT,P2,d,A0);
        
        %         if loglik2 < loglik1; disp('Likelihood decreases!!!'); break; end
    end
    
end

loglik(iter:end) = [];   

% rearrange in ascending order of frequency
%=========================================================================
% f, z, S
%--------
[f,I]=sort(f);
z = z(I);
phi = phi(:,I);
S = S(I,I);

[BA1,alph1,alph2] = svd(phi);    % mode shape space representation
BA1 = BA1(:,1:m);
alph = alph1(1:m,:)*alph2';

toc
time_req(1,:) = toc;

if all(abs(err)<tol_cvg)
    fprintf('Results have converged to within a relative tolerance of %5.3e\n',tol_cvg);
else
    fprintf('Results have not converged to within a relative tolerance of %5.3e\n',tol_cvg);
end


% calc. s/n ratios
sn = (diag(S).*(1/4./z(:).^2)./Se(:)).'; % (1,m)

% modal rms
rms = sqrt(pi/4*diag(S).'.*f./z); % g, rms, S is two-sided g/sqrt(Hz)

%% =========================================================================
% output

% expand V to all modes
Se = repmat(Se,[1,m]); % (1,m)

% calc. s/n ratios
sn = (diag(S).*(1/4./z(:).^2)./Se(:)).'; % (1,m)


%========================================================================
% output
%========================================================================
out = [];
out.f = f;
out.z = z;
out.S = S;
out.Sii = diag(S).';
out.Se = Se;
out.phi = phi;
out.rms = rms;
out.sn = sn;
out.time_req = time_req;
% out.cmat = cmat;
out.loglik = loglik;

% new outputs
% varx = {'f','z','Sii','Se','phi'};
% for i=1:5
%   out.coefv.(varx{i}) = coefv{i};
% end
% out.coefv.rms = rmscov;

% set to nan for either 
% 1. freq. outside band
% 2. damping outside limits
% 3. sn limit too small
%-------------------------------------
I = out.f(:)<f1f2(1)|out.f(:)>f1f2(2);
I = I|out.z(:)<zlim(1)|out.z(:)>zlim(2);
I = I|out.sn(:)<snmin;
if any(I)
  fieldname = {'f','z','S','Sii','Se','phi','rms','sn'};
  for x2 = fieldname
    out.(x2{:})(:) = nan;
  end

%   varx = fieldnames(out.coefv).';
%   for x2 = varx
%     out.coefv.(x2{:})(:)=nan;
%     out.coefv0.(x2{:})(:)=nan;
%     out.coefv1.(x2{:})(:)=nan;
%   end  
end

out.S = {out.S}; % cell
out.loglik = {out.loglik}; % cell

end
