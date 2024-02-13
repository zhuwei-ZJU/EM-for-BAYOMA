%% Calculate the 3nd and 4th moments of the latent variable (Modal response)
%========================================================================
% 220725-Firstly written by Wei at ZJUI
% 221105-Use 'pagemtimes' function to avoid the 'for loop'
%========================================================================
% Input
%========================================================================
% Mu: Mean of Yitak [m,1,nf]
% Sigma: Covariance matrix of Yitak [m,m,nf]
% Mu2: Second original moment of Yitak [m,m,nf]
%========================================================================
% Output
%========================================================================
% The 3-4th original and central moment of the latent variable
% Mu_tilde2(:,:,kk) = E[Yitak*Yitak.']
% Nu3(:,:,kk) = E[vec(Yitak*Yitak')*Yitak']
% Nu_tilde3(:,:,kk) = E[vec(Yitak*Yitak')*Yitak.']
% Nu4(:,:,kk) = E[vec(Yitak*Yitak')*vec(Yitak*Yitak').']
%========================================================================
function [Mu_tilde2,Nu3,Nu_tilde3,Nu4] = EYitak34(Mu,Sigma,Mu2)
% Used for calculating the posterior COV of modal parameters:
% The product of the first derivative and the expected product of the first
% derivative

[m,~,nf] = size(Mu);
Mu_tilde2 = pagemtimes(Mu,'none',Mu,'transpose');
Im = eye(m);
Im111 = repelem(eye(m),1,m);
er1(:,1,1,:) = Im;   % [m,1,1,m]
er1 = repmat(er1,[1,1,m^2,1,nf]);   % [m,1,m^2,m,nf]
er2(:,1,:,1) = Im;   % [m,1,m,1]
er2 = repmat(er2,[1,1,m,m,nf]);   % [m,1,m^2,m,nf]
er3(:,1,:,1) = Im111;   % [m,1,m^2,1]
er3 = repmat(er3,[1,1,1,m,nf]);   % [m,1,m^2,m,nf]
MU(:,:,1,1,:) = Mu;   % [m,1,1,1,nf]
MU = repmat(MU,[1,1,m^2,m,1]);   % [m,1,m^2,m,nf]
SIGMA(:,:,1,1,:) = Sigma;   % [m,m,1,1,nf]
SIGMA = repmat(SIGMA,[1,1,m^2,m,1]);   % [m,m,m^2,m,nf]
MU2(:,:,1,1,:) = Mu2;   % [m,m,1,1,nf]
MU2 = repmat(MU2,[1,1,m^2,m,1]);   % [m,m,m^2,m,nf]
er1er2 = pagemtimes(er1,'none',er2,'transpose');   % [m,m,m^2,m,nf]
MUer1er2 = pagemtimes(MU,'ctranspose',er1er2,'none');   % [1,m,m^2,m,nf]
SIGMAer1er2 = pageprodtr(SIGMA,er1er2);   % [1,1,m^2,m,nf]
Nu3 = pagemtimes(pagemtimes(MUer1er2,MU2)+pagemtimes(SIGMAer1er2,'none',MU,'ctranspose'),er3);
Nu3 = shiftdim(Nu3,2);   % [m^2,m,nf]

er4 = er2;   er5 = er3;   er6 = er1;
er5er6 = pagemtimes(er5,'none',er6,'transpose');   % [m,m,m^2,m,nf]
er5er6MU = pagemtimes(er5er6,MU);
er5er6SIGMA = pageprodtr(er5er6,SIGMA);
Nu_tilde3 = pagemtimes(er4,'transpose',pagemtimes(MU2,er5er6MU)+pagemtimes(MU,er5er6SIGMA),'none');
Nu_tilde3 = shiftdim(Nu_tilde3,2);   % [m^2,m,nf]

er7 = repmat(er2,[1,1,1,m,1]);   % [m,1,m^2,m^2,nf]
er8 = repmat(er3,[1,1,1,m,1]);   % [m,1,m^2,m^2,nf]
er9(:,1,1,:) = Im;   % [m,1,1,m]
er9 = repmat(er9,[1,1,m^2,m,nf]);   % [m,1,m^2,m^2,nf]
er10(:,1,1,:) = Im111;   % [m,1,1,m^2]
er10 = repmat(er10,[1,1,m^2,1,nf]);   % [m,1,m^2,m^2,nf]
MU = repmat(MU,[1,1,1,m,1]);   % [m,1,m^2,m^2,nf]
SIGMA = repmat(SIGMA,[1,1,1,m,1]);   % [m,m,m^2,m^2,nf]
MU2 = repmat(MU2,[1,1,1,m,1]);   % [m,m,m^2,m^2,nf]
er8er9 = pagemtimes(er8,'none',er9,'transpose');   % [m,m,m^2,m^2,nf]
er8er9SIGMA = pageprodtr(er8er9,SIGMA);   % [1,1,m^2,m^2,nf]
MU2MU2 = pagemtimes2(MU2,er8er9,MU2);   % [m,m,m^2,m^2,nf]
MUMUSIGMA = pagemtimes2(MU,'ctranspose',er8er9,'none',MU,SIGMA);
MU2SIGMA = pagemtimes(MU2,er8er9SIGMA);
Nu4 = pagemtimes2(er7,'transpose',MU2MU2+MUMUSIGMA+MU2SIGMA,'none',er10);
Nu4 = shiftdim(Nu4,2);   % [m^2,m^2,nf]
end