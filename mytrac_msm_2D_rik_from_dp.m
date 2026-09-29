function [tx,ty,nxx,nxy,nyy,ut0,vt0,X0,Y0,Lx,Ly]...
    = mytrac_msm_2D_rik_from_dp(E,h,h0,nx0,ny0,x,y,u,v,fcalx,Nyq)

%   ^ x (also remember that when doing imagesc, the vertical axis is upside down)
%   |
%   |
%   |
%   |
%   --------> Y
%   
% E Young's modulus of the gel;
% h = thickness of the gel
% h0 = height at which beads are (if beads on surface h = h0)
% [nx0,ny0] = output resolution (in case you want over-sampled output)
% fcalx =  Calibration factor of the microscope (in microns / px)
% Nyq = Nyquist limit for cutoff filter since output can't contain frequencies higher than those measured.

x = x'; y = y'; u = u'; v = v'; % Legacy PIV had all transposed

n1= size(x,1);
n2= size(y,2);
w = zeros(size(x));
% median
u(isnan(u))=0;
u=u-median(u(:));
v(isnan(v))=0;
v=v-median(v(:));
% %%%%%%%%%%%%

xr = sort(unique(x))*fcalx;
yr = sort(unique(y))*fcalx;

ix0 =  1;
ixe = n1;
jy0 =  1;
jye = n2;

vecx=ix0:length(xr);
vecy=jy0:length(yr);

%theory behind it
KXM = pi/(xr(2)-xr(1))*2/Nyq;
KYM = pi/(yr(2)-yr(1))*2/Nyq;

%new xy has real length, x has the values 16,32,64... -1/2 to put in center
x  = xr(vecx) - xr(1)/2;
y  = yr(vecy) - yr(1)/2;

u = u(vecx,vecy)*fcalx;
v = v(vecx,vecy)*fcalx;
w = w(vecx,vecy)*1;


Lx0 = (x(end) - x(1)); %bego
Ly0 = (y(end) - y(1)); %bego

xgap = Lx0/(nx0-1);
ygap = Ly0/(ny0-1);

%X0 Y0 are in um

% Old (would cause problems with Lx != Ly
% [Y0,X0] = meshgrid(x(1):xgap:x(end),y(1):ygap:y(end));
% [y0,x0] = meshgrid(x,y);

% Correct, counterintuitive but consequence of how X and Y are defined from PIV output 
[Y0,X0] = meshgrid(y(1):ygap:y(end),x(1):xgap:x(end));
[y0,x0] = meshgrid(y,x);


% Corrected
ut = interp2(x0',y0',u',X0',Y0')';
vt = interp2(x0',y0',v',X0',Y0')';
wt = interp2(x0',y0',w',X0',Y0')';

ut0=ut;
vt0=vt;
wt0=wt;

%       fudges for extrap

ut(ut~=ut)=0;
vt(vt~=vt)=0;
wt(wt~=wt)=0;


[Y,X] = meshgrid(0:ygap:Ly0,0:xgap:Lx0);


% filtering
alpha = 0.2;

HWx = zeros(size(X)) ;
HWy=HWx;

HWx(X<=alpha*(Lx0)/2) = (1 + cos(pi*(2*X(X<=alpha*(Lx0)/2)/(alpha*(Lx0))-1)))/2 ;
HWx((X>alpha*(Lx0)/2).*(X<=(Lx0)*(1-alpha/2))>0) = 1 ;
HWx(X>(Lx0)*(1-alpha/2)) = (1 + cos(pi*(2*X(X>(Lx0)*(1-alpha/2))/(alpha*(Lx0))-2/alpha+1)))/2 ;

HWy(Y<=alpha*(Ly0)/2) = (1 + cos(pi*(2*Y(Y<=alpha*(Ly0)/2)/(alpha*(Ly0))-1)))/2 ;
HWy((Y>alpha*(Ly0)/2).*(Y<=(Ly0)*(1-alpha/2))>0) = 1 ;
HWy(Y>(Ly0)*(1-alpha/2)) = (1 + cos(pi*(2*Y(Y>(Ly0)*(1-alpha/2))/(alpha*(Ly0))-2/alpha+1)))/2 ;

HW=HWx.*HWy;

ut = ut.*HW;
vt = vt.*HW;
wt = wt.*HW;
%
nx = nx0;
ny = ny0;
Lx = Lx0;
Ly = Ly0;



Theta_x = 2*pi*X/Lx0;
Theta_y = 2*pi*Y/Ly0;

%	Transforms displacements to Fourier space
%
ut =     fft2(ut);  %the zero-frequency component is in the upper-left corner of the two-dimensional FFT
% ut = fftshift(ut);
%
vt =     fft2(vt);  %the zero-frequency component is in the upper-left corner of the two-dimensional FFT
% vt = fftshift(vt);
%
wt =     fft2(wt);  %the zero-frequency component is in the upper-left corner of the two-dimensional FFT


ut = fftshift(ut);  %the zero-frequency component is near the center of the matrix
vt = fftshift(vt);  %the zero-frequency component is near the center of the matrix
wt = fftshift(wt);  %the zero-frequency component is near the center of the matrix

kx = [-nx/2:nx/2-1]*2*pi/Lx;
ky = [-ny/2:ny/2-1]*2*pi/Ly;
k2max = max(kx)^2 + max(ky)^2;

KXM = pi/(xr(2)-xr(1))*2/Nyq;
KYM = pi/(yr(2)-yr(1))*2/Nyq;


nfil = 2; %width of the Gaussian filter in multiples of (dx^2 + dy^2)^0.5
%	Compute the tractions from the displacements

sig = 0.49;

I = sqrt(-1);
tx = zeros(nx,ny);
ty = zeros(nx,ny);
tz = zeros(nx,ny);
nxx = zeros(nx,ny);
nxy = zeros(nx,ny);
nyy = zeros(nx,ny);


for ix=1:nx
    if (abs(kx(ix))<=KXM)
        for iy=1:ny
            if (abs(ky(iy))<=KYM)
                k2 = kx(ix)^2+ky(iy).^2;
                filter = exp(nfil*k2/k2max/4*pi^2);
                %            [txz,tyz,tzz,uu,uv,uw] = forcesonlytzz(kx(ix),ky(iy),h,sig,squeeze(tz(ix,iy)),E);
                %            [txz,tyz,tzz] = forcestzz(kx(ix),ky(iy),h,h,sig,ut(ix,iy),vt(ix,iy));
                % [txz,tyz,tzz] = forces3D(kx(ix),ky(iy),h,h,sig,ut(ix,iy),vt(ix,iy),mic*wt(ix,iy));
                [txz,tyz,tzz,u,v,w] = forces3D(kx(ix),ky(iy),h,h0,sig,ut(ix,iy),vt(ix,iy),wt(ix,iy));
                tx(ix,iy) = E*txz/filter;
                ty(ix,iy) = E*tyz/filter;
                tz(ix,iy) = E*tzz/filter;
            end
        end
    end
end

%
% zero modes
%
tx(nx/2+1,ny/2+1)=ut(nx/2+1,ny/2+1)*E/h/(1+sig)/2;
ty(nx/2+1,ny/2+1)=vt(nx/2+1,ny/2+1)*E/h/(1+sig)/2;
tz(nx/2+1,ny/2+1)=wt(nx/2+1,ny/2+1)*E*(1-sig)/h/(1+sig)/(1-2*sig);

for ix=1:nx
    if (abs(kx(ix))<=KXM)
        for iy=1:ny
            if (abs(ky(iy))<=KYM)
                k2 = kx(ix)^2+ky(iy).^2;
               
                % Lateral Monolayer Stresses
                nxx(ix,iy) = -I*kx(ix)*tx(ix,iy)/k2 + I*ky(iy)*ty(ix,iy)/k2 - I*(ky(iy)^2*(1+sig)*(kx(ix)*tx(ix,iy)+ky(iy)*ty(ix,iy)))/k2^2; 
                nxy(ix,iy) = -I*ky(iy)*tx(ix,iy)/k2 - I*kx(ix)*ty(ix,iy)/k2 + I*(kx(ix)*ky(iy)*(1+sig)*(kx(ix)*tx(ix,iy)+ky(iy)*ty(ix,iy)))/k2^2;
                nyy(ix,iy) = I*kx(ix)*tx(ix,iy)/k2 - I*ky(iy)*ty(ix,iy)/k2 - I*(kx(ix)^2*(1+sig)*(kx(ix)*tx(ix,iy)+ky(iy)*ty(ix,iy)))/k2^2; 
            end
        end
    end
end
nxx(nx/2+1,ny/2+1)=0;
nxy(nx/2+1,ny/2+1)=0;
nyy(nx/2+1,ny/2+1)=0;


%


tx = ifftshift(tx);
ty = ifftshift(ty);
tz = ifftshift(tz);


nxx = ifftshift(nxx);
nxy = ifftshift(nxy);
nyy = ifftshift(nyy);

ut = ifftshift(ut);
vt = ifftshift(vt);
wt = ifftshift(wt);

% elaser=sum(sum(tx*ut+ty*vt+tz*wt));
% elaser=real(elaser);


tx = real(ifft2(tx))';
ty = real(ifft2(ty))';
tz = real(ifft2(tz))';


nxx = real(ifft2(nxx))';
nxy = real(ifft2(nxy))';
nyy = real(ifft2(nyy))';

X0 = X0';
Y0 = Y0';
ut0 = ut0';
vt0 = vt0';