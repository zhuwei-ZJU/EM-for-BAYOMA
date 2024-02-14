%% A trivial vector that can pick some specific elements in a matrix
%===============================================================
% 220525-Firstly written by Wei at ZJUI
% 220701-Exptend to the general matrix form
% 220902-Change the output into sparse form
%===============================================================
% Input
%===============================================================
% Two type of input
% #1: pick_vec(m,r) e.g.: pick_vec(3,2)=>[0 1 0]'
% #2: pick_vec(n,m,a,b): pick the [a,b] part of the [n,m] matrix
%===============================================================
% Output:
%===============================================================
% #1: e: picking vector
% #2: [e1,e2]: e1: picking matrix for rows, e2: for columns
% Satisfy the equality: e1*A*e2 = A(a,b)
%===============================================================
function varargout = pick_vec(varargin)
if nargin < 2
    error('Too few input arguments!');
elseif nargin == 2
    if length(varargin{2}) > 1
        error('A vector can only pick up one element at a time!');
    else
        e = zeros(varargin{1},1);
        e(varargin{2}) = 1;
        varargout{1} = sparse(e);
    end
elseif nargin == 4
    e1 = zeros(length(varargin{3}),varargin{1});
    e2 = zeros(varargin{2},length(varargin{4}));
    e1(:,varargin{3}) = eye(length(varargin{3}));
    e2(varargin{4},:) = eye(length(varargin{4}));
    varargout{1} = sparse(e1);
    varargout{2} = sparse(e2);
elseif nargin > 4
    error('Too many input arguments!')
else
    error('Wrong input!');
end
end