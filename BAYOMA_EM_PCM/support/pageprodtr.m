%% Page-wise matrix product trace
%============================================================================
% 221107-Firstly written by Wei at ZJUI
% 221110-Allow the input array to be a matrix
%============================================================================
% Input
%============================================================================
% A: Muti-dimensional array with a square matrix in the first two dimensions 
% B: Muti-dimensional array with the same first two dimensions as A
%============================================================================
% Output
%============================================================================
% traceAB: put the trace of the first two dimension of AB into the first two
% dimensions
%============================================================================
function traceAB = pageprodtr(A,B)
% Set the sencond input to the identity array with appropriate size, one can
% get the page-wise matrix trace of the first input

% Check input
IA = size(A);
IB = size(B);
if IA(1) ~= IA(2)
    error('The first two dimensions of the first input array must be the same!');
end
if IB(1) ~= IB(2)
    error('The first two dimensions of the second input array must be the same!');
end
if IB(1) ~= IA(1) || IB(2) ~= IA(2)
    error('The first two dimesions of A and B must be the same!');
end
if length(IA) < 3
    vecAdp = vec(A.');
else
    vecAdp = pagevec(pagetranspose(A));
end
if length(IB) < 3
    vecB = vec(B);
else
    vecB = pagevec(B);
end
traceAB = pagemtimes(vecAdp,'transpose',vecB,'none');
end