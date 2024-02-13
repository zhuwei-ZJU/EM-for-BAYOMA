%% Creat a diagnoal matrix or extract the diagonal elements in a page-wise manner
%============================================================================
% 221108-Firstly written by Wei at ZJUI
% 221109-Remove the ndSparse since it can't apply to array multiplications
%============================================================================
% Input
%============================================================================
% A: Muti-dimensional array with a square matrix or a vector
% in the first two dimensions
%============================================================================
% Output
%============================================================================
% diagA: Extract the digonal elements or use the vertor components as the
% diagnonal elements to creat a matrix in the first two dimensions of diagA
%============================================================================
function diagA = pagediag(A)
% Check input
IA = size(A);
if length(A) < 3
    error('The dimension of the input must be at least 3!')
end
if IA(1)~=IA(2) && IA(1)~=1 && IA(2)~=1
    error('The first 2 dimensions of the input must be a square matrix or a vector!');
end

if IA(1)==1
    diagA = repelem(A,IA(2),1).*eye(IA(2));
end
if IA(2)==1
    diagA = repelem(A,1,IA(1)).*eye(IA(1));
end

if IA(1)==IA(2)
    IA3end = IA(3:end);
    diagA = zeros([IA(1),1,IA3end]);
    IIA = repmat(eye(IA(1)),[1,1,IA3end]);
    diagA(1:end) = A(IIA~=0);
end
end

