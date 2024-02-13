%% Calculate the mode assurance criterion (MAC) between columns of 2 matrices
%===========================================================
% 230522-Firstly written by Wei at ZJUI
%===========================================================
% Input
%===========================================================
% A: [m,n]
% B: [m,p]
%===========================================================
% Output
%===========================================================
% M: [n,p]
%===========================================================
% Syntax
%===========================================================
% M = calMAC(A)
% M = calMAC(A,b)
%===========================================================
function M = calMAC(A,B)
if nargin == 1
    M = vec(A.'*A./(vecnorm(A).'*vecnorm(A)),'hh');
elseif size(A,1)~=size(B,1)
    error('The dimensions of input matrices must be consistent!');
else
    M = A.'*B./(vecnorm(A).'*vecnorm(B));
end
end