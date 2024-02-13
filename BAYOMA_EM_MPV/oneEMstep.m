function [P,loglik,Sei] = oneEMstep(Fki,fki,P0,d,A0)
% Frequency-domain Bayesian Identification: EM step

[Nf,No] = size(Fki);
mmodes = size(P0,2);   % #modes in the selected band

fi = P0(1,:); zi = P0(2,:); phii = P0(3:2+No,:);    % frequency, damping, mode shapes
sqr_iSi = P0(3+No:2+No+mmodes,:);    % inverse of modal force PSD, S^-1/2
sqr_Sei = P0(end,1);    % channel noise, Se^1/2
Sei = sqr_Sei^2;

betaik = fi'./fki;
ihik = 1-betaik.^2 - 1i*2*zi'.*betaik; % hik^-1
sqr_phii = chol(phii'*phii);
rki = (Fki*phii/sqr_phii).';    % sqr_phi'^-1*phi'*Fk

%% E step
% E-step
EQk = zeros(mmodes,Nf); % E[Qk]
EQk2 = zeros(mmodes,mmodes,Nf); % E[QkQk']
S0 = zeros(mmodes,mmodes,Nf);  % used in computing Si in M step
logliki = zeros(1,Nf);

for k0 = 1:Nf
    R = qr([sqr_phii rki(:,k0);sqr_Sei*sqr_iSi*diag(ihik(:,k0)) zeros(mmodes,1)]);
    r12 = R(1:mmodes,end); R11 = triu(R(1:mmodes,1:mmodes));
    Rr = R11\[sqr_Sei*eye(mmodes) r12];
    hRr = diag(ihik(:,k0))*Rr;
    EQk(:,k0) = Rr(:,end);  % iR11*r12;
    EQk2(:,:,k0) = Rr*Rr';
    S0(:,:,k0) = hRr*hRr';
    %             EQk2(:,:,k0) = iR11*(Sei*eye(mmodes)+r12*r12')*iR11';
    logliki(k0) = -2*sum(log(abs(diag(R11))))+r12'*r12/Sei;
end
loglik = -No*Nf*log(pi)-(No-mmodes)*log(Sei)- d/Sei ...
    + 2*sum(sum(log(abs(ihik))))+2*Nf*sum(log(abs(diag(sqr_iSi)))) + sum(logliki);


%% M-step
FkEQk = real(EQk*conj(Fki))'/Nf;  % sum(Fk*EQk')
EQk2mean = mean(real(EQk2),3);
% L = chol([EQk2mean FkEQk'; FkEQk A0/Nf],'lower');
% L11 = tril(L(1:mmodes,1:mmodes)); L21 = L(mmodes+1:end,1:mmodes);
% L22 = tril(L(mmodes+1:end,mmodes+1:end));
% phii = L21/L11;     % mode shape
% Sei = trace(L22*L22')/No;   % channel noise PSD
phii = FkEQk/EQk2mean;
Sei = trace(A0/Nf-phii*EQk2mean*phii')/No;
% sqr_Sei = sqrt(Sei);

Si = mean(S0,3); % modal force PSD
if length(Si) == 1; Si = real(Si); end
sqr_iSi = chol(Si,'lower')^-1;    % inverse of modal force PSD, S^-1/2, lower triangular matrix

q2 = @(x)Q2(x,[fi zi],fki,sqr_iSi'*sqr_iSi,EQk2);
x = fminsearch(q2,ones(1,mmodes*2));
x = x.*[fi zi];
fi = x(1:mmodes);   % freqncies
zi = x(mmodes+1:end);   % damping ratios

% return to original parameters
rnorm = sqrt(sum(phii.^2));
phii = phii./rnorm;
sqr_iSi = diag(1./rnorm)*sqr_iSi;    % S^-1

P = [fi; zi; phii; sqr_iSi; sqr_Sei*ones(1,mmodes)];     % output

end