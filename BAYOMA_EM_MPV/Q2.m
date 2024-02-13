function q2 = Q2(x,fz,fki,iSi,EQk2)
%% Q2 function to optimize f & z

mmodes = length(x)/2;    %#modes
x = x.*fz;
f = x(1:mmodes);    % frequencies
z = x(mmodes+1:end); % damping ratios

Nf = length(fki);   % # frequency points in selected band

betak = f'./fki;
ihk = 1-betak.^2 - 1i*2*z'.*betak; % hik^-1
iDk = (1-betak.^2).^2 + 4*z'.^2.*betak.^2; % Dk^-1


S0i = zeros(mmodes,mmodes,Nf);
for m = 1:mmodes
    for mm = m:mmodes
        S0i(m,mm,:) = ihk(m,:).*conj(ihk(mm,:)).*reshape(EQk2(m,mm,:),1,[]);
    end
end
S0i = sum(S0i,3)/Nf;
S0i = S0i+S0i';
S0i = S0i-diag(real(diag(S0i)))/2;


q2 = -sum(log(prod(iDk,1))) + Nf*real(trace(iSi*S0i));


end