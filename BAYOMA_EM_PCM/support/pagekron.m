%% Page-wise kronecker product
%============================================================================
% 221107-Firstly written by Wei at ZJUI
%============================================================================
% Input
%============================================================================
% A: Muti-dimensional array 
% B: Muti-dimensional array with the same dimension as A except the first
% two dimensions
%============================================================================
% Output
%============================================================================
% K: put the kronecker product of the first two dimension of AB
% into the first two dimensions
%============================================================================
function K = pagekron(A,B)
IA = size(A);
IB = size(B);
K = repelem(A,IB(1),IB(2)).*repmat(B,IA(1),IA(2));
end