function [Ufilt,Vfilt,dUdx,dVdy] = FourierFiltAndDiff(X,Y,U,V,cutoff)
    
    % Hanning window from Ricardo's codes
    alpha = 0.2;

    HWx = zeros(size(X)) ;
    HWy=HWx;
    Lx = range(X(:));
    Ly = range(Y(:));

    HWx(X<=alpha*(Lx)/2) = (1 + cos(pi*(2*X(X<=alpha*(Lx)/2)/(alpha*(Lx))-1)))/2 ;
    HWx((X>alpha*(Lx)/2).*(X<=(Lx)*(1-alpha/2))>0) = 1 ;
    HWx(X>(Lx)*(1-alpha/2)) = (1 + cos(pi*(2*X(X>(Lx)*(1-alpha/2))/(alpha*(Lx))-2/alpha+1)))/2 ;

    HWy(Y<=alpha*(Ly)/2) = (1 + cos(pi*(2*Y(Y<=alpha*(Ly)/2)/(alpha*(Ly))-1)))/2 ;
    HWy((Y>alpha*(Ly)/2).*(Y<=(Ly)*(1-alpha/2))>0) = 1 ;
    HWy(Y>(Ly)*(1-alpha/2)) = (1 + cos(pi*(2*Y(Y>(Ly)*(1-alpha/2))/(alpha*(Ly))-2/alpha+1)))/2 ;

    HW=HWx.*HWy;

    U = U .* HW;
    V = V .* HW;

    % Zero-padding
    padSize = 30;
    U = padarray(U,[padSize padSize]);
    V = padarray(V,[padSize padSize]);

    [ny,nx] = size(U);
    
    % Take out the max cutoff% of the frequencies (Arbitrary)   
    if nargin <= 4
        cutoff = 0.7;
    end
    Uf = fft2(U);
    Vf = fft2(V);
    k = [-nx/2:nx/2-1]*((2*pi)/(Lx));
    l = [-ny/2:ny/2-1]*((2*pi)/(Ly));
    [K,L] = meshgrid(k,l);
    
    K = fftshift(K); % alpha
    L = fftshift(L); % beta
    
    Kmax = max(abs(K(:))); Lmax= max(abs(L(:))); 
    inds = abs(K) > (cutoff*Kmax) | abs(L) > (cutoff*Lmax);
    
    Uf(inds) = 0; Vf(inds) = 0;
    K(inds) = 0; L(inds) = 0;
    
    Ufilt = real(ifft2(Uf));
    Vfilt = real(ifft2(Vf));
    
    % Take the first derivative
    dUfdx = 1i.*K.*Uf;
    dVfdy = 1i.*L.*Vf;
    
    dUfdx(inds) = 0;dVfdy(inds) = 0;
    
    dUdx = real(ifft2(dUfdx));
    dVdy = real(ifft2(dVfdy));
    
    Ufilt = Ufilt(padSize+1:end-padSize,padSize+1:end-padSize);
    Vfilt = Vfilt(padSize+1:end-padSize,padSize+1:end-padSize);
    
    dUdx = dUdx(padSize+1:end-padSize,padSize+1:end-padSize);
    dVdy = dVdy(padSize+1:end-padSize,padSize+1:end-padSize);
    
end