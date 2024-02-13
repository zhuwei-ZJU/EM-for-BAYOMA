%% A matrix transforing the vectorization into diag
function R = vec2diag(n)
%===========================================
% 220511-Firstly written by Wei at ZJUI
%===========================================
% Input
%===========================================
% n: size of the squared matrix
%===========================================
% Output:
%===========================================
% R: [n,n^2] R*vec(X) = diag(X)
%===========================================
R = zeros(n,n^2);
for ii = 1:n
    R(ii,n*(ii-1)+ii) = 1;
end
R = sparse(R);
end