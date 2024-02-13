function L = oneloglik(Fki,fki,P,d)
% Calculate log likelihood for Frequency-domain Bayesian Identification
% robust implementation

[Nf,No] = size(Fki);
mmodes = size(P,2);   % #modes in the selected band

fi = P(1,:); zi = P(2,:); phii = P(3:2+No,:);    % frequency, damping, mode shapes
sqr_iSi = P(3+No:2+No+mmodes,:);    % inverse of modal force PSD, S^-1/2
sqr_Sei = P(end,1);    % channel noise, Se^1/2
Sei = sqr_Sei^2;

betaik = fi'./fki;
ihik = 1-betaik.^2 - 1i*2*zi'.*betaik; % hik^-1
sqr_phii = chol(phii'*phii);
rki = (Fki*phii/sqr_phii).';    % sqr_phi'^-1*phi'*Fk

logliki = zeros(1,Nf);
for k0 = 1:Nf
    R = qr([sqr_phii rki(:,k0);sqr_iSi*diag(sqr_Sei*ihik(:,k0)) zeros(mmodes,1)]);
    r12 = R(1:mmodes,end); R11 = triu(R(1:mmodes,1:mmodes));
    logliki(k0) = -2*sum(log(abs(diag(R11))))+r12'*r12/Sei;
end
L = -No*Nf*log(pi)-(No-mmodes)*log(Sei)-d/Sei ...
    + 2*sum(sum(log(abs(ihik))))+2*Nf*sum(log(abs(diag(sqr_iSi)))) + sum(logliki);