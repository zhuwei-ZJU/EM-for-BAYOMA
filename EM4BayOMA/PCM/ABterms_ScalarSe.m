%% Calculate some frequently used 1st and 2nd derivative terms of Q function
% Note only the A13,14,15,16 and B18,19,20,24,25 are different from ABterms
%========================================================================
% 220924-Firstly written by Wei at ZJUI
% 221108-Use 'pagemtimes' function to get ABterms at all frequency points
%========================================================================
% Input
%========================================================================
% f: Natural frequency [m,1]
% z: Damping ratio [m,1]
% PHI: Mode shape [n,m]
% S: Modal force PSD [m,m]
% Se: Prediction error: scalar
% ff: Frequency, [nf,1]
% F: FFT of measured data, [nf,n]
% Mu,Mu2: The 1st and 2nd central mement of the latent variable
%========================================================================
% Output
%========================================================================
% A structure contains the derivative terms used for calculate the
% gradient and hessian matrix of the Q function
%========================================================================
function oo = ABterms_ScalarSe(f,z,PHI,S,Se,ff,F,Mu,Mu2)
[n,m] = size(PHI);   % [number of dof, number of modes]
beta = f./ff.';   % [m,nf]
h = 1./(1 - beta.^2 - 1i*2*beta.*z);  ih = 1./h;    % [m,nf]
D = 1./((1-beta.^2).^2+4*z.^2.*beta.^2);   D2 = D.^2;   % [m,nf]
ff2 = ff.^2; ff4 = ff.^4; ff6 = ff.^6;   % [nf,1]
f2 = f.^2; f3 = f.^3; f4 = f.^4;   % [m,1]
z2 = z.^2;

f_ff = f./ff.';   f_ff2 = f./ff2.';   z_ff = z./ff.';   % [m,nf]
dhdf(:,1,:) = 2*f_ff2 + 2i*z_ff;   % -d(diag(inv(hk)))/df   [m,1,nf]
dhdf = pagediag(dhdf);   % [m,m,nf]
dhdz(:,1,:) = 2i*f_ff;   % -d(diag(inv(hk)))/dz   [m,1,nf]
dhdz = pagediag(dhdz);   % [m,m,nf]
dhpdf = conj(dhdf);   % -d(diag(inv(hk')))/df   [m,m,nf]
iS = inv(S);   veciS = vec(iS.');   iSiS = kron(iS.',iS);
ih3(:,1,:) = ih;   % [m,1,nf]
iSdpihp = iS.'.*pagectranspose(ih3);   iSih = iS.*pagetranspose(ih3);    % [m,m,nf]
Rm = full(vec2diag(m));   Kmm = com_mat(m,m);  Im = eye(m);  In = eye(n);   RmKmm = Rm*Kmm;
iSe = 1/Se;   iSe2 = 1/Se^2;   iSe3 = 1/Se^3;
ImiSePHI = kron(Im,iSe*PHI);
vecPHIPHI = vec(PHI.'*PHI).';   % for scalar Se
%========================================================================
% First derivative terms that does not depend on the latent variable
%========================================================================
% f
oo.Alpha1 = pagemtimes2(Rm,pagekron(dhdf,iSdpihp),Kmm);   % [m,m^2,nf]
oo.Alpha2 = pagemtimes2(Rm,pagekron(pagemtimes(dhpdf,iSih),Im),Kmm);   % [m,m^2,nf]
oo.Alpha3(:,1,:) = -4*D.*f./ff2.' + 4*D.*f3./ff4.' + 8*D.*f.*z2./ff2.';   % [m,1,nf]
% z
oo.Alpha4 = pagemtimes2(Rm,pagekron(dhdz,iSdpihp),Kmm);   % [m,m^2,nf]
oo.Alpha5 = -pagemtimes2(Rm,pagekron(pagemtimes(dhdz,iSih),Im),Kmm);   % [m,m^2,nf]
oo.Alpha6(:,1,:) = 8*D.*f2.*z./ff2.';   % [m,1,nf]
% PHI
F3(:,1,:) = F.';   % [n,1,nf]
oo.Alpha7(:,1:m,:) = iSe*pagekron(Im,F3);   % [m*n,m,nf]
oo.Alpha7 = reshape(oo.Alpha7,m*n,m,[]);   % 231221 add for single dof
oo.Alpha8 = conj(oo.Alpha7);   % [m*n,1,nf]
oo.Alpha9 = -ImiSePHI;   % [m*n,m^2]
oo.Alpha10 = -ImiSePHI*Kmm;   % [m*n,m^2]
% S
oo.Alpha11 = pagemtimes(Kmm,pagekron(iSdpihp,iSih));   % [m^2,m^2,nf]
oo.Alpha12 = -veciS;   % [m^2,1]
% Se
oo.Alpha13(1,:,:) = -iSe2*PHI.'*F.';   % [1,m,nf]
oo.Alpha14 = conj(oo.Alpha13);   % [1,m,nf]
oo.Alpha15 = iSe2*vecPHIPHI;   % [m^2,1]
diagFF = diag(F*F');
oo.Alpha16(1,1,:) = -n*iSe + iSe2*diagFF;   % [1,1,nf]
%========================================================================
% Second derivative terms
%========================================================================
% ff
Mu2dp = pagetranspose(Mu2);
oo.Beta1 = -pagemtimes(pagemtimes(dhdf,Mu2).*iS.',dhpdf);   % [m,m,nf]
oo.Beta2 = -pagemtimes(pagemtimes(dhpdf,iS).*Mu2dp,dhdf);   % [m,m,nf]
fff2(1,:,:) = ff2.'; 
oo.Beta3 = 2*pagemtimes(iSdpihp,Mu2dp).*Im./fff2;   % [m,m,nf]
oo.Beta4 = 2*pagemtimes(Mu2,'transpose',iSdpihp,'ctranspose').*Im./fff2;   % [m,m,nf]
oo.Beta5(:,1,:) = D2.*(-16*f2./ff4.'+16*f4./ff6.'+32*z2.*f2./ff4.')-4*D./ff2.'  ...
    + D2.*(16*f4./ff6.'-16*f.^6./ff.^8.'-32*z2.*f4./ff6.')+12*D.*f2./ff4.'   ...
    + D2.*(32*z2.*f2./ff4.'-32*z2.*f4./ff6.'-64*z.^4.*f2./ff4.')+8*D.*z2./ff2.';   % [m,1,nf]
oo.Beta5 = pagediag(oo.Beta5);   % [m,m,nf]
% fz
oo.Beta6 = pagemtimes(pagemtimes(dhdf,Mu2).*iS.',dhdz);   % [m,m,nf]
fff(1,:,:) = ff.';
oo.Beta7 = 2i*pagemtimes(iSdpihp,Mu2dp)./fff.*Im;   % [m,m,nf]
oo.Beta8 = -pagemtimes(pagemtimes2(dhpdf,iS).*Mu2dp,dhdz);   % [m,m,nf]
oo.Beta9 = -2i*pagemtimes(Mu2dp,'none',iSih,'transpose').*Im./fff;   % [m,m,nf]
oo.Beta10(:,1,:) = 32*D2.*f3.*z./ff4.' - 32*D2.*f.^5.*z./ff6.'   ...
-64*D2.*f3.*z.^3./ff4.' + 16*D.*f.*z./ff2.';  % [m,1,nf]
oo.Beta10 = pagediag(oo.Beta10);   % [m,m,nf]
% fs
Mu2ihp = Mu2.*pagectranspose(ih3);
Mu2dpih = Mu2dp.*pagetranspose(ih3);
oo.Beta11 = -pagemtimes2(RmKmm,pagekron(Im,pagemtimes(dhdf,Mu2ihp)),iSiS); % [m,m^2,nf]
oo.Beta12 = -pagemtimes2(RmKmm,pagekron(Mu2dpih,dhpdf),iSiS); % [m,m^2,nf]
% zz
oo.Beta13(:,1,:) = -64*D2.*f4.*z2./ff4.' + 8*D.*f2./ff2.';   % [m,1,nf]
oo.Beta13 = pagediag(oo.Beta13);   % [m,m,nf]
oo.Beta14 = pagemtimes(pagemtimes(dhdz,Mu2).*iS.',dhdz);   % [m,m,nf]
oo.Beta15 = pagemtimes(pagemtimes(dhdz,iS).*Mu2dp,dhdz);   % [m,m,nf]
% zS
oo.Beta16 = -pagemtimes2(RmKmm,pagekron(Im,pagemtimes(dhdz,Mu2ihp)),iSiS);   % [m,m,nf]
oo.Beta17 = pagemtimes2(RmKmm,pagekron(Mu2dpih,dhdz),iSiS);   % [m,m,nf] 
% PHIPHI
oo.Beta18 = -iSe*pagekron(Mu2dp,In);   % [mn,mn,nf]
oo.Beta19 = conj(oo.Beta18);   % [mn,mn,nf]
% PHISe
oo.Beta20 = iSe2*pagevec(-pagemtimes(F3,'none',Mu,'ctranspose')- ...
    pagemtimes(conj(F3),'none',Mu,'transpose')+pagemtimes(PHI,2*real(Mu2)));   % [mn,1,nf]
% SS
iSihMu2ihpiS = pagemtimes2(iSih,Mu2,iSih,'ctranspose');
oo.Beta21 = -pagemtimes(Kmm,pagekron(pagetranspose(iSihMu2ihpiS),iS));   % [m^2,m^2,nf]
oo.Beta22 = -pagemtimes(Kmm,pagekron(iS.',iSihMu2ihpiS));   % [m^2,m^2,nf]
oo.Beta23 = Kmm*iSiS;   % [m^2,m^2]
% SeSe
oo.Beta24(1,1,:) = n*iSe2 - 2*iSe3*diagFF;   % [1,1,nf]
MupPHIpF3 = pagemtimes2(Mu,'ctranspose',PHI,'transpose',F3);
oo.Beta25 =  2*iSe3*(2*real(MupPHIpF3)-pageprodtr(PHI.'*PHI,Mu2));   % [1,1,nf]
end
