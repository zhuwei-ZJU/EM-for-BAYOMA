%% Stack the matrix (or its lower triangle parts) along the column direction
%===========================================================
% 220629-Firstly written by Wei at ZJUI
% 220712-Add vechh
% 221109-Add the basic vec
%===========================================================
% Input
%===========================================================
% A: a matrix
% Optional input: 'h' (half-vertorization) or 
%                 'hh' (half-hallow-vertorization) 
%===========================================================
% Output
%===========================================================
% vecA: [n^2,1]   vectorize A
% vechA: [n(n+1)/2,1]   vectorize the lower triangle parts
% vechhA: [n(n-1)/2,1]   vectorize the lower triangle parts
% (not including the diagonal entries)  n: dimension of A
%===========================================================
% Syntax
%===========================================================
% vecA = vec(A)   <=> A(:)
% vechA = vec(A,'h')
% vechhA = vec(A,'hh')
%===========================================================
function a = vec(A,varargin)
[nr,nc] = size(A);
if nargin > 2
    error('Too much input!');
elseif nargin < 2
    a = A(:);
elseif nr ~= nc || length(size(A))~=2
    error('Input matrix must be square in this input form!');
elseif contains('h',varargin)
    II = tril(ones(nr));
    a = A(II~=0);
elseif contains('hh',varargin)
    II = tril(ones(nr)) - eye(nr);
    a = A(II~=0);
else
    error('Unknown input!');
end
end



