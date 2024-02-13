%% Calculate the 1st and 2nd moments of the latent variable (Modal response)
%========================================================================
% 220725-Firstly written by Wei at ZJUI
% 221105-Use 'pagemtimes' function to avoid the 'for loop'
% 221121-Fixed bug in function 'pagemldivide'
%========================================================================
% Input
%========================================================================
% f: Natural frequency [m,1]
% z: Damping ratio [m,1]
% PHI: Mode shape [n,m]
% S: Modal force PSD [m,m]
% Se: Prediction error [n,n] diagonal
% ff: Frequency, [nf,1]
% F: FFT of measured data, [nf,n]
%========================================================================
% Output
%========================================================================
% The 1st and 2nd original and central moment of the latent variable
% Mu: Mean of Yitak [m,1,nf]
% Sigma: Covariance matrix of Yitak [m,m,nf]
% Mu2: Second original moment of Yitak, Mu2(:,:,kk)= E[Yitak*Yitak']
% Nu2: Vectorization of the sencond original moment [m^2,1,nf]
%========================================================================
function [Mu,Sigma,Mu2,Nu2] = EYitak12(f,z,PHI,S,Se,ff,F)
% Used for calculating the gradient and hessian matrix of the Q function:
% Not only will be used for the posterior COV calculation process, but also
% for the MPV calculation process if some gradient-based optimization 
% methods are used!

beta(:,1,:) = f./ff.';   % [m,1,nf]
h = 1./(1 - beta.^2 - 1i*2*beta.*z);   % [m,1,nf]
H = h.*S.*pagectranspose(h);   % [m,m,nf]
iSe = diag(1./diag(Se));
P = pageinv(H) + PHI.'*iSe*PHI;   % [m,m,nf]
F3(:,1,:) = F.';   % [n,1,nf]
% Mu = pagemtimes(pagemldivide(P,PHI.')*iSe,F3);   % [m,1,nf]
Mu = pagemtimes2(pageinv(P),PHI.',iSe,F3);   % [m,1,nf]
Sigma = pageinv(P);   % [m,m,nf]
Mu2 = pagemtimes(Mu,'none',Mu,'ctranspose') + Sigma;   % [m,m,nf]
Nu2 = pagevec(Mu2);   % [m^2,1,nf]
% Nu2temp = shiftdim(Mu2,2);   % [nf,m,m]
% Nu2(:,1,:) = Nu2temp(:,:).';   % [m^2,1,nf]
end