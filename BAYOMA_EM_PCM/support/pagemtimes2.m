%% Page-wise matrix multiplication: allow arbitrary number of input matrices
%============================================================================
% 221109-Firstly written by Wei at ZJUI
%============================================================================
% Input
%============================================================================
% A: Muti-dimensional array with a square matrix or a vector
% in the first two dimensions
%============================================================================
% Output
%============================================================================
% Multiplication of matrices in a page-wise manner
%============================================================================
% Syntax
%============================================================================
% Z = pagemtime2(X1,X2,X3)
% Z = pagemtime2(X1,'transpX1',X2,'transpX2',X3,'transpX3')
% transpX must be 'transpose', 'ctranspose', or 'none'
%============================================================================
function Z = pagemtimes2(varargin)
if ischar(varargin{2})
    Z = pagemtimes(varargin{1},varargin{2},varargin{3},varargin{4});
    ii = 5;
    while ii <= nargin-1
        if ischar(varargin{ii+1})
            Z = pagemtimes(Z,'none',varargin{ii},varargin{ii+1});
            ii = ii+2;
        else
            Z = pagemtimes(Z,varargin{ii});
            ii = ii+1;
        end
    end
    if ii == nargin
        Z = pagemtimes(Z,varargin{ii});
    end

else
    Z = pagemtimes(varargin{1},varargin{2});
    ii = 3;
    while ii <= nargin-1
        if ischar(varargin{ii+1})
            Z = pagemtimes(Z,'none',varargin{ii},varargin{ii+1});
            ii = ii+2;
        else
            Z = pagemtimes(Z,varargin{ii});
            ii = ii+1;
        end
    end
    if ii == nargin
        Z = pagemtimes(Z,varargin{ii});
    end
end
end

