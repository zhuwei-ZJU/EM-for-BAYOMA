%% Calculate the gradient and hessian matrix of the NLLF (negative log-liklihood function)
% An indirect way by Louis Identity
%========================================================================
% 220904-Firstly written by Wei at ZJUI
% 221111-Consistent with the page-wise operation
%========================================================================
% Input
%========================================================================
% f: MPV of natural frequency [m,1]
% z: MPV of damping ratio [m,1]
% PHI: MPV of mode shape [n,m]
% S: MPV of modal force PSD [m,m]
% Se: MPV of prediction error scalar
% ff: Frequency, [nf,1]
% F: FFT of measured data, [nf,n]
%========================================================================
% Output
%========================================================================
% grad: Gradient of the NLLF
% Hess: Hessian matrix of the NLLF
%========================================================================
function [grad,Hess] = NLLFHess_scalarSe(f,z,PHI,S,Se,ff,F)
[n,m] = size(PHI);   % [number of dof, number of modes]
ntheta = m + m + m*n + m + m*(m-1) + 1;   % number of parameters
Rm = full(vec2diag(m));  Dh = dup_mat(m,'hh');   Ds = dup_mat(m,'sk');
[Mu,Sigma,Mu2,Nu2] = EYitak12(f,z,PHI,S,Se,ff,F);
[Mu_tilde2,Nu3,Nu_tilde3,Nu4] = EYitak34(Mu,Sigma,Mu2);
oo = ABterms_ScalarSe(f,z,PHI,S,Se,ff,F,Mu,Mu2);
A1 = oo.Alpha1; A2 = oo.Alpha2; A3 = oo.Alpha3; A4 = oo.Alpha4;
A5 = oo.Alpha5; A6 = oo.Alpha6; A7 = oo.Alpha7; A8 = oo.Alpha8;
A9 = oo.Alpha9; A10 = oo.Alpha10; A11 = oo.Alpha11; A12 = oo.Alpha12;
A13 = oo.Alpha13; A14 = oo.Alpha14; A15 = oo.Alpha15; A16 = oo.Alpha16;
B1 = oo.Beta1; B2 = oo.Beta2; B3 = oo.Beta3; B4 = oo.Beta4;
B5 = oo.Beta5; B6 = oo.Beta6; B7 = oo.Beta7; B8 = oo.Beta8;
B9 = oo.Beta9; B10 = oo.Beta10; B11 = oo.Beta11; B12 = oo.Beta12;
B13 = oo.Beta13; B14 = oo.Beta14; B15 = oo.Beta15; B16 = oo.Beta16;
B17 = oo.Beta17; B18 = oo.Beta18; B19 = oo.Beta19; B20 = oo.Beta20;
B21 = oo.Beta21; B22 = oo.Beta22; B23 = oo.Beta23;  B24 = oo.Beta24;
B25 = oo.Beta25;
%========================================================================
% The first derivative of Q function = the frist derivative of LLF
%========================================================================
dLLF_df = pagemtimes(A1,Nu2) + pagemtimes(A2,Nu2) + A3;
dLLF_dz = pagemtimes(A4,Nu2) + pagemtimes(A5,Nu2) + A6;
dLLF_dPHI = pagemtimes(A7,conj(Mu)) + pagemtimes(A8,Mu) + pagemtimes(A9,Nu2) + pagemtimes(A10,Nu2);
dLLF_dS = pagemtimes(A11,Nu2) + A12;
dLLF_dSe = pagemtimes(A13,conj(Mu)) + pagemtimes(A14,Mu) + pagemtimes(A15,Nu2) + A16;
%========================================================================
% Fisher's identity
%========================================================================
grad = -real([dLLF_df;dLLF_dz;dLLF_dPHI;pagemtimes(Rm,dLLF_dS);...
    pagemtimes(Dh,'transpose',dLLF_dS,'none');1i*pagemtimes(Ds,'transpose',dLLF_dS,'none');dLLF_dSe]);
grad = sum(grad,3).';
if nargout > 1
    %========================================================
    % Initialization
    %========================================================
    Hess = zeros(ntheta,ntheta);
    If = 1:m;   Iz = m+1:2*m;   IPHI = 2*m+1:2*m+m*n;
    IdiagS = 2*m+m*n+1:2*m+m*n+m;   IReSij = 2*m+m*n+m+1:2*m+m*n+m+m*(m-1)/2;
    IImSij = 2*m+m*n+m+m*(m-1)/2+1:2*m+m*n+m+m*(m-1); ISe = ntheta;
    %========================================================================
    % The expected product of 2 first derivatives
    %========================================================================
    % ff
    E1_ff = pagemtimes2(A1,Nu4,A1,'transpose') + pagemtimes2(A1,Nu4,A2,'transpose') + pagemtimes2(A1,Nu2,A3,'transpose') ...
        + pagemtimes2(A2,Nu4,A1,'transpose') + pagemtimes2(A2,Nu4,A2,'transpose') + pagemtimes2(A2,Nu2,A3,'transpose') ...
        + pagemtimes2(A3,'none',Nu2,'transpose',A1,'transpose') + pagemtimes2(A3,'none',Nu2,'transpose',A2,'transpose') + pagemtimes(A3,'none',A3,'transpose');
    % fz
    E1_fz = pagemtimes2(A1,Nu4,A4,'transpose') + pagemtimes2(A1,Nu4,A5,'transpose') + pagemtimes2(A1,Nu2,A6,'transpose') ...
        + pagemtimes2(A2,Nu4,A4,'transpose') + pagemtimes2(A2,Nu4,A5,'transpose') + pagemtimes2(A2,Nu2,A6,'transpose') ...
        + pagemtimes2(A3,'none',Nu2,'transpose',A4,'transpose') + pagemtimes2(A3,'none',Nu2,'transpose',A5,'transpose') + pagemtimes(A3,'none',A6,'transpose');
    % fPHI
    E1_fPHI = pagemtimes2(A1,Nu3,A7,'transpose') + pagemtimes2(A1,Nu_tilde3,A8,'transpose') + pagemtimes2(A1,Nu4,A9,'transpose') + pagemtimes2(A1,Nu4,A10,'transpose') ...
        + pagemtimes2(A2,Nu3,A7,'transpose') + pagemtimes2(A2,Nu_tilde3,A8,'transpose') + pagemtimes2(A2,Nu4,A9,'transpose') + pagemtimes2(A2,Nu4,A10,'transpose') ...
        + pagemtimes2(A3,'none',Mu,'ctranspose',A7,'transpose') + pagemtimes2(A3,'none',Mu,'transpose',A8,'transpose') + pagemtimes2(A3,'none',Nu2,'transpose',A9,'transpose') + pagemtimes2(A3,'none',Nu2,'transpose',A10,'transpose');
    % fS
    E1_fS = pagemtimes2(A1,Nu4,A11,'transpose') + pagemtimes2(A1,Nu2,A12,'transpose')  ...
        + pagemtimes2(A2,Nu4,A11,'transpose') + pagemtimes2(A2,Nu2,A12,'transpose')  ...
        + pagemtimes2(A3,'none',Nu2,'transpose',A11,'transpose') + pagemtimes(A3,'none',A12,'transpose');
    % fSe
    E1_fSe = pagemtimes2(A1,Nu3,A13,'transpose') + pagemtimes2(A1,Nu_tilde3,A14,'transpose') + pagemtimes2(A1,Nu4,A15,'transpose') + pagemtimes2(A1,Nu2,A16,'transpose') ...
        + pagemtimes2(A2,Nu3,A13,'transpose') + pagemtimes2(A2,Nu_tilde3,A14,'transpose') + pagemtimes2(A2,Nu4,A15,'transpose') + pagemtimes2(A2,Nu2,A16,'transpose')...
        + pagemtimes2(A3,'none',Mu,'ctranspose',A13,'transpose') + pagemtimes2(A3,'none',Mu,'transpose',A14,'transpose') + pagemtimes2(A3,'none',Nu2,'transpose',A15,'transpose') + pagemtimes(A3,'none',A16,'transpose');
    % zz
    E1_zz = pagemtimes2(A4,Nu4,A4,'transpose') + pagemtimes2(A4,Nu4,A5,'transpose') + pagemtimes2(A4,Nu2,A6,'transpose') ...
        + pagemtimes2(A5,Nu4,A4,'transpose') + pagemtimes2(A5,Nu4,A5,'transpose') + pagemtimes2(A5,Nu2,A6,'transpose') ...
        + pagemtimes2(A6,'none',Nu2,'transpose',A4,'transpose') + pagemtimes2(A6,'none',Nu2,'transpose',A5,'transpose') + pagemtimes(A6,'none',A6,'transpose');
    % zPHI
    E1_zPHI = pagemtimes2(A4,Nu3,A7,'transpose') + pagemtimes2(A4,Nu_tilde3,A8,'transpose') + pagemtimes2(A4,Nu4,A9,'transpose') + pagemtimes2(A4,Nu4,A10,'transpose') ...
        + pagemtimes2(A5,Nu3,A7,'transpose') + pagemtimes2(A5,Nu_tilde3,A8,'transpose') + pagemtimes2(A5,Nu4,A9,'transpose') + pagemtimes2(A5,Nu4,A10,'transpose') ...
        + pagemtimes2(A6,'none',Mu,'ctranspose',A7,'transpose') + pagemtimes2(A6,'none',Mu,'transpose',A8,'transpose') + pagemtimes2(A6,'none',Nu2,'transpose',A9,'transpose') + pagemtimes2(A6,'none',Nu2,'transpose',A10,'transpose');
    % zS
    E1_zS = pagemtimes2(A4,Nu4,A11,'transpose') + pagemtimes2(A4,Nu2,A12,'transpose')  ...
        + pagemtimes2(A5,Nu4,A11,'transpose') + pagemtimes2(A5,Nu2,A12,'transpose')  ...
        + pagemtimes2(A6,'none',Nu2,'transpose',A11,'transpose') + pagemtimes(A6,'none',A12,'transpose');
    % zSe
    E1_zSe = pagemtimes2(A4,Nu3,A13,'transpose') + pagemtimes2(A4,Nu_tilde3,A14,'transpose') + pagemtimes2(A4,Nu4,A15,'transpose') + pagemtimes2(A4,Nu2,A16,'transpose') ...
        + pagemtimes2(A5,Nu3,A13,'transpose') + pagemtimes2(A5,Nu_tilde3,A14,'transpose') + pagemtimes2(A5,Nu4,A15,'transpose') + pagemtimes2(A5,Nu2,A16,'transpose')...
        + pagemtimes2(A6,'none',Mu,'ctranspose',A13,'transpose') + pagemtimes2(A6,'none',Mu,'transpose',A14,'transpose') + pagemtimes2(A6,'none',Nu2,'transpose',A15,'transpose') + pagemtimes(A6,'none',A16,'transpose');
    % PHIPHI
    E1_PHIPHI = pagemtimes2(A7,conj(Mu_tilde2),A7,'transpose') + pagemtimes2(A7,'none',Mu2,'transpose',A8,'transpose') + pagemtimes2(A7,'none',Nu3,'transpose',A9,'transpose') + pagemtimes2(A7,'none',Nu3,'transpose',A10,'transpose') ...
        + pagemtimes2(A8,Mu2,A7,'transpose') + pagemtimes2(A8,Mu_tilde2,A8,'transpose') + pagemtimes2(A8,'none',Nu_tilde3,'transpose',A9,'transpose') + pagemtimes2(A8,'none',Nu_tilde3,'transpose',A10,'transpose') ...
        + pagemtimes2(A9,Nu3,A7,'transpose') + pagemtimes2(A9,Nu_tilde3,A8,'transpose') + pagemtimes2(A9,Nu4,A9,'transpose') + pagemtimes2(A9,Nu4,A10,'transpose') ...
        + pagemtimes2(A10,Nu3,A7,'transpose') + pagemtimes2(A10,Nu_tilde3,A8,'transpose') + pagemtimes2(A10,Nu4,A9,'transpose') + pagemtimes2(A10,Nu4,A10,'transpose');
    % PHIS
    E1_PHIS = pagemtimes2(A7,'none',Nu3,'transpose',A11,'transpose') + pagemtimes2(A7,conj(Mu),A12,'transpose')  ...
        + pagemtimes2(A8,'none',Nu_tilde3,'transpose',A11,'transpose') + pagemtimes2(A8,Mu,A12,'transpose')  ...
        + pagemtimes2(A9,Nu4,A11,'transpose') + pagemtimes2(A9,Nu2,A12,'transpose') ...
        + pagemtimes2(A10,Nu4,A11,'transpose') + pagemtimes2(A10,Nu2,A12,'transpose');
    % PHISe
    E1_PHISe = pagemtimes2(A7,conj(Mu_tilde2),A13,'transpose') + pagemtimes2(A7,'none',Mu2,'transpose',A14,'transpose') + pagemtimes2(A7,'none',Nu3,'transpose',A15,'transpose') + pagemtimes2(A7,conj(Mu),A16,'transpose') ...
        + pagemtimes2(A8,Mu2,A13,'transpose') + pagemtimes2(A8,Mu_tilde2,A14,'transpose') + pagemtimes2(A8,'none',Nu_tilde3,'transpose',A15,'transpose') + pagemtimes2(A8,Mu,A16,'transpose')...
        + pagemtimes2(A9,Nu3,A13,'transpose') + pagemtimes2(A9,Nu_tilde3,A14,'transpose') + pagemtimes2(A9,Nu4,A15,'transpose') + pagemtimes2(A9,Nu2,A16,'transpose') ...
        + pagemtimes2(A10,Nu3,A13,'transpose') + pagemtimes2(A10,Nu_tilde3,A14,'transpose') + pagemtimes2(A10,Nu4,A15,'transpose') + pagemtimes2(A10,Nu2,A16,'transpose');
    % SS
    E1_SS = pagemtimes2(A11,Nu4,A11,'transpose') + pagemtimes2(A11,Nu2,A12,'transpose')  ...
        + pagemtimes2(A12,'none',Nu2,'transpose',A11,'transpose') + pagemtimes(A12,'none',A12,'transpose');
    % SSe
    E1_SSe = pagemtimes2(A11,Nu3,A13,'transpose') + pagemtimes2(A11,Nu_tilde3,A14,'transpose') + pagemtimes2(A11,Nu4,A15,'transpose') + pagemtimes2(A11,Nu2,A16,'transpose') ...
        + pagemtimes2(A12,'none',Mu,'ctranspose',A13,'transpose') + pagemtimes2(A12,'none',Mu,'transpose',A14,'transpose') + pagemtimes2(A12,'none',Nu2,'transpose',A15,'transpose') + pagemtimes(A12,'none',A16,'transpose');
    % SeSe
    E1_SeSe = pagemtimes2(A13,conj(Mu_tilde2),A13,'transpose') + pagemtimes2(A13,'none',Mu2,'transpose',A14,'transpose') + pagemtimes2(A13,'none',Nu3,'transpose',A15,'transpose') + pagemtimes2(A13,conj(Mu),A16,'transpose') ...
        + pagemtimes2(A14,Mu2,A13,'transpose') + pagemtimes2(A14,Mu_tilde2,A14,'transpose') + pagemtimes2(A14,'none',Nu_tilde3,'transpose',A15,'transpose') + pagemtimes2(A14,Mu,A16,'transpose')...
        + pagemtimes2(A15,Nu3,A13,'transpose') + pagemtimes2(A15,Nu_tilde3,A14,'transpose') + pagemtimes2(A15,Nu4,A15,'transpose') + pagemtimes2(A15,Nu2,A16,'transpose') ...
        + pagemtimes2(A16,'none',Mu,'ctranspose',A13,'transpose') + pagemtimes2(A16,'none',Mu,'transpose',A14,'transpose') + pagemtimes2(A16,'none',Nu2,'transpose',A15,'transpose') + pagemtimes(A16,'none',A16,'transpose');
    %========================================================================
    % The expectation of the second derivatives
    %========================================================================
    E2_ff = B1 + B2 + B3 + B4 + B5;
    E2_fz = B6 + B7 + B8 + B9 + B10;
    E2_fS = B11 + B12;
    E2_zz = B13 + B14 + B15;
    E2_zS = B16 + B17;
    E2_PHIPHI = B18 + B19;
    E2_PHISe = B20;
    E2_SS = B21 + B22 + B23;
    E2_SeSe = B24 + B25;
    %========================================================================
    % The product of the first derivative of Q function
    %========================================================================
    % ff
    E3_ff = pagemtimes(dLLF_df,'none',dLLF_df,'transpose');
    % fz
    E3_fz = pagemtimes(dLLF_df,'none',dLLF_dz,'transpose');
    % fPHI
    E3_fPHI = pagemtimes(dLLF_df,'none',dLLF_dPHI,'transpose');
    % fS
    E3_fS = pagemtimes(dLLF_df,'none',dLLF_dS,'transpose');
    % fSe
    E3_fSe = pagemtimes(dLLF_df,'none',dLLF_dSe,'transpose');
    % zz
    E3_zz = pagemtimes(dLLF_dz,'none',dLLF_dz,'transpose');
    % zPHI
    E3_zPHI = pagemtimes(dLLF_dz,'none',dLLF_dPHI,'transpose');
    % zS
    E3_zS = pagemtimes(dLLF_dz,'none',dLLF_dS,'transpose');
    % zSe
    E3_zSe = pagemtimes(dLLF_dz,'none',dLLF_dSe,'transpose');
    % PHIPHI
    E3_PHIPHI = pagemtimes(dLLF_dPHI,'none',dLLF_dPHI,'transpose');
    % PHIS
    E3_PHIS = pagemtimes(dLLF_dPHI,'none',dLLF_dS,'transpose');
    % PHISe
    E3_PHISe = pagemtimes(dLLF_dPHI,'none',dLLF_dSe,'transpose');
    % SS
    E3_SS = pagemtimes(dLLF_dS,'none',dLLF_dS,'transpose');
    % SSe
    E3_SSe = pagemtimes(dLLF_dS,'none',dLLF_dSe,'transpose');
    % SeSe
    E3_SeSe = pagemtimes(dLLF_dSe,dLLF_dSe);
    %========================================================================
    % Louis identity
    %========================================================================
    E_ff = - E1_ff - E2_ff + E3_ff;
    E_fz= - E1_fz - E2_fz + E3_fz;
    E_fPHI = - E1_fPHI + E3_fPHI;
    E_fS = - E1_fS - E2_fS + E3_fS;
    E_fSe = - E1_fSe + E3_fSe;
    E_zz = - E1_zz - E2_zz + E3_zz;
    E_zPHI = - E1_zPHI + E3_zPHI;
    E_zS = - E1_zS - E2_zS + E3_zS;
    E_zSe = - E1_zSe + E3_zSe;
    E_PHIPHI = - E1_PHIPHI - E2_PHIPHI + E3_PHIPHI;
    E_PHIS = - E1_PHIS + E3_PHIS;
    E_PHISe = - E1_PHISe - E2_PHISe + E3_PHISe;
    E_SS= - E1_SS - E2_SS + E3_SS;
    E_SSe = - E1_SSe + E3_SSe;
    E_SeSe = - E1_SeSe - E2_SeSe + E3_SeSe;
    %========================================================
    % Form the hessian matrix
    %========================================================
    E_ff = sum(E_ff,3); E_fz = sum(E_fz,3); E_fPHI = sum(E_fPHI,3);
    E_fS = sum(E_fS,3); E_fSe = sum(E_fSe,3); E_zz = sum(E_zz,3);
    E_zPHI = sum(E_zPHI,3); E_zS = sum(E_zS,3); E_zSe = sum(E_zSe,3);
    E_PHIPHI = sum(E_PHIPHI,3); E_PHIS = sum(E_PHIS,3); E_PHISe = sum(E_PHISe,3);
    E_SS = sum(E_SS,3); E_SSe = sum(E_SSe,3); E_SeSe = sum(E_SeSe);
    Hess(If,If) = E_ff;
    Hess(If,Iz) = E_fz;
    Hess(If,IPHI) = E_fPHI;
    Hess(If,IdiagS) = E_fS*Rm.';
    Hess(If,IReSij) = E_fS*Dh;
    Hess(If,IImSij) = 1i*E_fS*Ds;
    Hess(If,ISe) = E_fSe;
    Hess(Iz,Iz) = E_zz;
    Hess(Iz,IPHI) = E_zPHI;
    Hess(Iz,IdiagS) = E_zS*Rm.';
    Hess(Iz,IReSij) = E_zS*Dh;
    Hess(Iz,IImSij) = 1i*E_zS*Ds;
    Hess(Iz,ISe) = E_zSe;
    Hess(IPHI,IPHI) = E_PHIPHI;
    Hess(IPHI,IdiagS) = E_PHIS*Rm.';
    Hess(IPHI,IReSij) = E_PHIS*Dh;
    Hess(IPHI,IImSij) = 1i*E_PHIS*Ds;
    Hess(IPHI,ISe) = E_PHISe;
    Hess(IdiagS,IdiagS) = Rm*E_SS*Rm.';
    Hess(IdiagS,IReSij) = Rm*E_SS*Dh;
    Hess(IdiagS,IImSij) = 1i*Rm*E_SS*Ds;
    Hess(IdiagS,ISe) = Rm*E_SSe;
    Hess(IReSij,IReSij) = Dh.'*E_SS*Dh;
    Hess(IReSij,IImSij) = 1i*Dh.'*E_SS*Ds;
    Hess(IReSij,ISe) = Dh.'*E_SSe;
    Hess(IImSij,IImSij) = -Ds.'*E_SS*Ds;
    Hess(IImSij,ISe) = 1i*Ds.'*E_SSe;
    Hess(ISe,ISe) = E_SeSe;
    Hess = triu(Hess);
    Hess = real(Hess + Hess.' - Hess.*eye(ntheta));
end
end