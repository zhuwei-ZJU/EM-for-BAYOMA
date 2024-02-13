function P = P2P(P,m,n)
% Transform to original parameters

phi = P(3:2+n,:);
sqr_iSi = P(3+n:2+n+m,:);

rnorm = sqrt(sum(phi.^2));
phi = phi./rnorm;
sqr_iSi = diag(1./rnorm)*sqr_iSi;    % S^-1

P = [P(1:2,:); phi; sqr_iSi; P(end,:)];     % output