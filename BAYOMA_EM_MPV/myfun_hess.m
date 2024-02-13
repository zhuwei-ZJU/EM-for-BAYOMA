function [hh,gg]=myfun_hess(f,z,S,V,BA1,alph,ff,F,dimx)
% calculate hessian using analytical formalae
%
% THIS FUNCTION MAKES USE OF fft02_funLcchess or fft02_funLhess, 
% BUT THERE IS SOME DIFFERENCE
% IN DEFINITION, SO THIS FUNCTION SERVES TO BRIDGE THE GAP
%
% differences:
% fft02_funL
% fft02_funLcchess             bayoma
%--------------------------------------------
% deriv wrt original var       wrt ratio
%
% Notes:
% fft02_funLhess - original formulation, for n<=m
% fft02_funLcchess - fast formulation, for n>=m
 

% [n,m]=size(phi);
% 
% if n<=m % original formulation
%   [hh,gg]=fft02_funLhess(f,z,S,V,phi,ff,F,dimx); % cell  
% else % n>=m, fast formulation
%   [hh,gg]=fft02_funLcchess(f,z,S,V,phi,ff,F,dimx); % cell
% end

n = size(BA1,1);
m = size(alph,2);
if n <= m % original form
    [hh,gg]=fft02_funLhess(f,z,S,V,BA1,alph,ff,F,dimx); % cell
else % n>m, fast form
    [hh,gg]=fft02_funLcchess(f,z,S,V,BA1*alph,ff,F,dimx); % cell
end

% non-dimensionalize gg and hh
%-------------------------------------------------------------------------
[diagS,ReS,ImS]=herm2vec(S,1);
xx = {f(:).',z(:).',diagS(:).',ReS(:).',ImS(:).',V,ones(1,n*m)}; % {1,nblk} cell

if dimx(5)==0 % exclude ImS
  xx{5}=[];
end
if dimx(7)==0 % exclude phi
  xx{7}=[];
end
xx = cell2mat(xx); % (1,nx)
% xx = abs(xx); % important! otherwise inconsistent for Tij that can be -ve
  
gg = gg.*xx; % (1,nx)
hh = diag(xx)*hh*diag(xx); % (nx,nx)
end