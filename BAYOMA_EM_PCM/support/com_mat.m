%% Commutation matrix tranforming vec(A) into vec(A')
function K = com_mat(n,m)
%======================================================
% ref: https://en.wikipedia.org/wiki/Commutation_matrix
%======================================================
% 220511-Firstly written by Wei at ZJUI
%======================================================
% Input
%======================================================
% n: number of rows
% m: number of columns
%======================================================
% Output
%======================================================
% K: [mn,mn] K*vec(A) = vec(A')
%======================================================
A = reshape(1:n*m,n,m);
v = reshape(A',1,[]);
K = eye(n*m);
K = K(v,:);
end