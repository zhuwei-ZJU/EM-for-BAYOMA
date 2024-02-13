%% A function to analytically obtain null space of the constraints gradient
%==========================================================
% 221206-Firstly written by Wei at ZJUI
%===========================================================
% Input
%===========================================================
% PHI: mode shape matrix within a selected band [n,m]
%===========================================================
% Output
%===========================================================
% U: null space of the constaints gradient [ntheta,ntheta-m]
%===========================================================
function U = cstgradnull(PHI)
[n,m] = size(PHI);
PHIbar = -PHI./PHI(1,:);   % [n-1,m]
%==========================================================
% Sparse form
%===========================================================
% ii = ones(m*n,1)*(1:n:m*n);
% jj = repmat((1:m*n),m,1);
% vv = repelem(speye(m),1,n).*vec(PHIbar).';   % [m,m*n]
% KsiPHI = sparse(ii,jj,vv,m*n,m*n);
% Ksi = speye(m*n) + KsiPHI;   % [m*n,m*n]
% Ksi(:,1:n:end) = [];   % [m*n,m*n-m]
% U = blkdiag(speye(2*m),Ksi,speye(m^2+1));
%==========================================================
% Full form
%===========================================================
KsiPHI = zeros(m*n,m*n);
KsiPHI(1:n:m*n,:) = repelem(speye(m),1,n).*vec(PHIbar).';   % [m*n,m*n]
Ksi = speye(m*n) + KsiPHI;   % [m*n,m*n]
Ksi(:,1:n:end) = [];   % [m*n,m*n-m]
U = blkdiag(eye(2*m),Ksi,eye(m^2+1));
U = orth(U);
end





