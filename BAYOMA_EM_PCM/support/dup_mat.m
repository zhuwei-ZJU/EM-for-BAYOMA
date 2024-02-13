%% Duplication matrix connecting vec and vech, vechh
%=============================================================
% ref: https://ww2.mathworks.cn/matlabcentral/answers/473737
%=============================================================
% 220629-Firstly written by Wei at ZJUI
% 220712-Add hollow-half-vectorization
%=============================================================
% Input
%=============================================================
% n: dimension of the square matrix
%=============================================================
% Output
%=============================================================
% D: [n^2,n(n+1)/2] D*vech(A)=>the same dimension of vec(A)
% i.e.: vertorization of a symmetric matrix A
% D: [n^2,n*(n-1)/2] D*vechh(a)=>the same dimension of vec(A)
% i.e.: vertorization of a hollow symmetric matrix A
%=============================================================
function D = dup_mat(n,varargin)
if nargin < 2
    m = n * (n + 1) / 2;
    nsq = n^2;
    r = 1;
    a = 1;
    v = zeros(1, nsq);
    cn = cumsum(n:-1:2);
    for i = 1:n
        v(r:r + i - 2) = i - n + cn(1:i - 1);
        r = r + i - 1;
        v(r:r + n - i) = a:a + n - i;
        r = r + n - i + 1;
        a = a + n - i + 1;
    end
    D = sparse(1:nsq, v, 1, nsq, m);
elseif contains(varargin,'hh')
    D = zeros(n^2,n*(n-1)/2);
    [r,c] = find(tril(ones(n)) - diag(ones(n,1)));
    I1 = sub2ind([n,n],r,c);
    I2 = sub2ind([n,n],c,r);
    D(I1,:) = eye(length(I1));
    D(I2,:) = eye(length(I2));
elseif contains(varargin,'sk')
    D = zeros(n^2,n*(n-1)/2);
    [r,c] = find(tril(ones(n)) - diag(ones(n,1)));
    I1 = sub2ind([n,n],r,c);
    I2 = sub2ind([n,n],c,r);
    D(I1,:) = eye(length(I1));
    D(I2,:) = -eye(length(I2));
else
    error('Wrong input!');
end
end