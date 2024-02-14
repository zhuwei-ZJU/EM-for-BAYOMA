%% Elimilation matrix transforming vec(A) into vech(A) or vechh(A)
function L = elim_mat(n,varargin)
%===========================================================
% Ref: https://ww2.mathworks.cn/matlabcentral/answers/274256
%===========================================================
% 220629-Firstly written by Wei at ZJUI
% 220712-Add hollow-half-vectorization
%===========================================================
% Input
%===========================================================
% n: dimension of the square matrix
%===========================================================
% Output
%===========================================================
% L: [n^2,n(n+1)/2] L*vec(A) = vech(A)
% Note I = find(tril(ones(n)))
% L = pick_vec(n^2,1,I,1) also return the same L
% L: [n^2,n(n-1)/2] L*vec(A) = vechh(A)
% Note I = find(tril(ones(n)) - diag(ones(n,1)))
% L = pick_vec(n^2,1,I,1) also return the same L
%===========================================================
if nargin < 2
    T = tril(ones(n));   % Lower triangle of 1's
    I = find(T);   % Get linear indexes of 1's
    L = zeros(n^2,n*(n+1)/2);   % Start with L'
    I = I + n^2*(0:n*(n+1)/2-1)';   % Linear indexes of the 1's within L'
elseif varargin{1} == 'hh'
    T = tril(ones(n)) - diag(ones(n,1));
    I = find(T);
    L = zeros(n^2,n*(n-1)/2);
    I = I + n^2*(0:n*(n-1)/2-1)';
else
    error('Wrong input!');
end
L(I) = 1;   % Put the 1's in place
L = L';   % Now transpose to actual L
end