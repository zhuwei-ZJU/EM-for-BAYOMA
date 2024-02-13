%% Page-wise stacking the matrix (or its lower triangle parts) along the column direction
%===============================================================
% 221110-Firstly written by Wei at ZJUI
%===============================================================
% Input
%===============================================================
% A: Multi-dimension array
% Optional input: 'h' (half-vertorization) or 
%                 'hh' (half-hallow-vertorization) 
%===============================================================
% Output
%===============================================================
% vecA: [n^2,1,...]   vectorize A
% vechA: [n(n+1)/2,1,...]   vectorize the lower triangle parts
% vechhA: [n(n-1)/2,1,...]   vectorize the lower triangle parts
% (not including the diagonal entries)  n: dimension of A
%===============================================================
% Syntax
%===============================================================
% vecA = pagevec(A) 
% vechA = pagevec(A,'h')
% vechhA = pagevec(A,'hh')
%===============================================================
function a = pagevec(A,varargin)
% High-dimension version of vec
IA = size(A);
% if length(IA) < 3
%     error('Please use the ''vec'' function!');
% end
nr = IA(1);   nc = IA(2);
if nargin > 2
    error('Too much input!');
elseif nargin < 2
    a = reshape(A,[nr*nc,1,IA(3:end)]);
%     a = zeros([nr*nc,1,IA(3:end)]);
%     a(:) = A(:);
elseif nr ~= nc
    error('The first two dimensions of the input array must be the same in this input form!');
elseif contains('h',varargin)
    IA3end = IA(3:end);
    II = repmat(tril(ones(nr)),[1,1,IA3end]);
    a = zeros([nr*(nr+1)/2,1,IA3end]);
    a(:) = A(II~=0);
elseif contains('hh',varargin)
    IA3end = IA(3:end);
    II = repmat(tril(ones(nr)) - eye(nr),[1,1,IA3end]);
    a = zeros([nr*(nr-1)/2,1,IA3end]);
    a(:) = A(II~=0);
else
    error('Unknown input!');
end
end
 