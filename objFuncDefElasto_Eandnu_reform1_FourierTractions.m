function [costVal,uu] = objFuncDefElasto_Eandnu_reform1_FourierTractions(E,nu,M,deltaS,deps,uBeads,delInds)
    
    % (13.2) eqn from L.L
    % Assumes that E and nu are constant 
    
    dsxxdX1 = (E./(1-nu.^2)) .* (deps.dExxdX1 + nu .* deps.dEyydX1);
    dsyydX2 = (E./(1-nu.^2)) .* (deps.dEyydX2 + nu .* deps.dExxdX2);
    dsxydX1 = (E./(1+nu)) .* (deps.dExydX1);
    dsxydX2 = (E./(1+nu)) .* (deps.dExydX2);
    
    % -Tx = dsxx/dx + dsigxy/dy
    ttx = ( -(dsxxdX1 + dsxydX2) );
    tty = ( -(dsxydX1 + dsyydX2) );
    
%     uu = M * [ttx(:); tty(:)]*(10/0.325);
ttx = ttx * (10/0.325);
tty = tty * (10/0.325);

% Try with Fourier method (Forward)
% Zero-pad with HW (Ricardo's code)
alpha = 0.2;
X = delInds{1};
Y = delInds{2};

Xind = zeros(size(X)); Yind = zeros(size(Y));
Xind(5:end-4,5:end-4) = X(5:end-4,5:end-4);
Yind(5:end-4,5:end-4) = Y(5:end-4,5:end-4);
inds = (X == Xind) & (Y == Yind);

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

ttx = ttx .* HW;
tty = tty .* HW;

padSize = 20;
ttx = padarray(ttx,[padSize padSize]);
tty = padarray(tty,[padSize padSize]);

[ny,nx] = size(ttx);

[uux,uuy] = performTFMLinButlerForward(nx,ny,deltaS,deltaS,ttx,tty);

uux = uux(padSize+1:end-padSize,padSize+1:end-padSize);
uuy = uuy(padSize+1:end-padSize,padSize+1:end-padSize);

[uux,uuy] = FourierFiltAndDiff(X,Y,uux,uuy,0.5);
% uux = imgaussfilt(uux,0.8);
% uuy = imgaussfilt(uuy,0.8);

uu = [uux(:); uuy(:)];
%     uu(delInds) = nan; % To remove contribution from close to 0 cell deformation
    
    costVal =  sqrt(mean( (uBeads(inds) - uu(inds)).^2 ) );
%       costVal = sum( (uBeads - uu).^ 2 .*abs(uBeads)) ./sum(abs(uBeads)) ;
    
end

