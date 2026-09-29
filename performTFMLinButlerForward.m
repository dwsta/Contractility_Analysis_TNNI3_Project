

function [Ox,Oy] = performTFMLinButlerForward(NxD,NyD,dxD,dyD,Ix,Iy)
    
    E=8e3;nu=0.48;
    h=150/0.325; 

    k = [-NxD/2:NxD/2-1]*((2*pi)/(dxD*NxD));
    l = [-NyD/2:NyD/2-1]*((2*pi)/(dyD*NyD));
    [K,L] = meshgrid(k,l);
    
    Ks = fftshift(K); % alpha
    Ls = fftshift(L); % beta
    K2s = (Ks.^2 + Ls.^2);
    q = sqrt(K2s); % q is the scalar, modulus of 2D wave vector (k,l)
    
    % Defintions
    c = cosh(q.*h);
    s = sinh(q.*h);
    
    % Constructing M matrix for vector u_hat = [u_k ; u_l];
    f1 = E.* (c.*q)./(2.*(1+nu).*s);
    f2 = ( E./(2.*(1-nu.^2).*q.*s) ) .* ...
            ( (3-4*nu)*nu.*s.*(c.^2) - (1-nu).*c.*q.*h + (1-2*nu).^2.*s + s.*(q.*h).^2 ) ./ ...
                    ( (3-4*nu).*s.*c + q.*h );

    % FFT of data
    Ixf = fft2(Ix);
    Iyf = fft2(Iy);
    
    Oxf = zeros(size(Ixf));
    Oyf = zeros(size(Iyf));
    g1 = 1./f1;
    g2 = -f2./(f1.*(f1+q.^2 .*f2));

    D1 = g1 + g2 .* Ks.^2;
    D2 = g1 + g2 .* Ls.^2;
    offDiag =  g2 .* Ks .* Ls;

    Oxf = D1 .* Ixf + offDiag .* Iyf;
    Oyf = offDiag .* Ixf + D2 .* Iyf;
    
    inds = isnan(Oxf) | isnan(Oyf);
    
    Oxf(inds) = 0;
    Oyf(inds) = 0;
    
    Ox = real(ifft2(Oxf));
    Oy = real(ifft2(Oyf));
end
