%
%       baldo and juanc 07/08/09
%	3D displacements prescribed on h
%	Condition tauzz = 0 on h is NOT IMPOSED NOW!!!
%
% juanca/baldo July 2009
% Last time modified by: baldo 11/28/09
%

function [tauxz,tauyz,tauzz,u,v,w] = forces3D(alpha,beta,h,h0,sigma,uph0,vph0,wph0);

I = sqrt(-1);
k = sqrt(alpha^2+beta^2);
%
if (k~=0) & (k*h < 100)
   %
   solup0(1) = -1/4*alpha^2.*h/k^2/(-1+sigma).*cosh(k*h) - ...
                1/4*(3*alpha^2-4*alpha^2*sigma+4*beta^2-4*beta^2*sigma)/k^3/(-1+sigma).*sinh(k*h);
   solup0(2) = -1/4*alpha*beta/k^2.*h/(-1+sigma).*cosh(k*h)+1/4*alpha*beta.*sinh(k*h)/k^3/(-1+sigma);
   solup0(3) =  1/4*I.*sinh(k*h).*h*alpha/(-1+sigma)/k;
   solup0(4) = -1/4*(4*k-4*k*sigma)/k/(-1+sigma).*cosh(k*h)-1/4.*h*alpha^2.*sinh(k*h)/k/(-1+sigma);
   solup0(5) = -1/4.*sinh(k*h).*h*beta/(-1+sigma)/k*alpha;
   solup0(6) =  1/4*I.*h*alpha/(-1+sigma).*cosh(k*h)+1/4*I.*sinh(k*h)*alpha/k/(-1+sigma);
   %
   solvp0(1) = -1/4*alpha*beta/k^2.*h/(-1+sigma).*cosh(k*h)+1/4*alpha*beta.*sinh(k*h)/k^3/(-1+sigma);
   solvp0(2) = -1/4*beta^2/k^2.*h/(-1+sigma).*cosh(k*h) - ...
                1/4*(3*beta^2-4*beta^2*sigma+4*alpha^2-4*alpha^2*sigma)/k^3/(-1+sigma).*sinh(k*h);
   solvp0(3) =  1/4*I.*sinh(k*h).*h*beta/k/(-1+sigma);
   solvp0(4) = -1/4.*sinh(k*h).*h*beta/(-1+sigma)/k*alpha;
   solvp0(5) = -1/4*(4*k-4*k*sigma)/k/(-1+sigma).*cosh(k*h)-1/4.*h*beta^2.*sinh(k*h)/k/(-1+sigma);
   solvp0(6) =  1/4*I.*h*beta/(-1+sigma).*cosh(k*h)+1/4*I.*sinh(k*h)*beta/k/(-1+sigma);
   %
   solwp0(1) =  1/2*I*alpha/k.*h.*sinh(k*h)/(-1+2*sigma);
   solwp0(2) =  1/2*I*beta/k.*h.*sinh(k*h)/(-1+2*sigma);
   solwp0(3) =  1/2.*h/(-1+2*sigma).*cosh(k*h)+1/2*(-3+4*sigma)/k/(-1+2*sigma).*sinh(k*h);
   solwp0(4) =  1/2*I*alpha.*h/(-1+2*sigma).*cosh(k*h)+1/2*I*alpha.*sinh(k*h)/k/(-1+2*sigma);
   solwp0(5) =  1/2*I.*h*beta/(-1+2*sigma).*cosh(k*h)+1/2*I.*sinh(k*h)*beta/k/(-1+2*sigma);
   solwp0(6) =  1/2*(-2+4*sigma)/(-1+2*sigma).*cosh(k*h)+1/2*k.*h.*sinh(k*h)/(-1+2*sigma);
   %
   C1N = (-8*I*(-1+sigma)*(alpha-I*beta)*(alpha+I*beta)*alpha*h0*k^3*cosh(2*k*h0)+8*I*(-1+sigma)*(alpha-I*beta)*(alpha+I*beta)*alpha*h0*k^3)*wph0+(8*alpha*h0*k^4*beta*(-1+sigma)*vph0+8*alpha^2*h0*k^4*(-1+sigma)*uph0)*sinh(2*k*h0)+(-2*alpha*beta*(-3+4*sigma)*k^3*vph0+2*(-3+4*sigma)*(4*alpha^2*sigma+4*beta^2*sigma-4*alpha^2-3*beta^2)*k^3*uph0)*cosh(2*k*h0)+2*alpha*beta*(2*h0^2*alpha^2+2*h0^2*beta^2-3+4*sigma)*k^3*vph0-2*(9*beta^2-28*alpha^2*sigma+2*h0^2*beta^2*alpha^2+2*h0^2*beta^4-24*beta^2*sigma+12*alpha^2+16*alpha^2*sigma^2+16*beta^2*sigma^2)*k^3*uph0;
   %
   C2N = (-8*I*(-1+sigma)*(alpha-I*beta)*(alpha+I*beta)*beta*h0*k^3*cosh(2*k*h0)+8*I*(-1+sigma)*(alpha-I*beta)*(alpha+I*beta)*beta*h0*k^3)*wph0+(8*k^4*h0*beta^2*(-1+sigma)*vph0+8*alpha*h0*k^4*beta*(-1+sigma)*uph0)*sinh(2*k*h0)+(2*k^3*(-3+4*sigma)*(4*alpha^2*sigma+4*beta^2*sigma-4*beta^2-3*alpha^2)*vph0-2*alpha*beta*(-3+4*sigma)*k^3*uph0)*cosh(2*k*h0)-2*k^3*(16*beta^2*sigma^2+16*alpha^2*sigma^2+12*beta^2+9*alpha^2+2*h0^2*alpha^4+2*h0^2*beta^2*alpha^2-24*alpha^2*sigma-28*beta^2*sigma)*vph0+2*alpha*beta*(2*h0^2*alpha^2+2*h0^2*beta^2-3+4*sigma)*k^3*uph0;
   %
   C3N = -4*k^2*h0*(alpha^2+beta^2)*(-1+2*sigma)*wph0*cosh(k*h0)+4*k*(alpha^2+beta^2)*(-1+2*sigma)*(-3+4*sigma)*sinh(k*h0)*wph0-4*I*(-1+2*sigma)*(alpha-I*beta)*(alpha+I*beta)*beta*k*h0*sinh(k*h0)*vph0-4*I*(-1+2*sigma)*(alpha-I*beta)*(alpha+I*beta)*alpha*k*h0*sinh(k*h0)*uph0;
   %
   C12D = -(alpha^2+beta^2)^2*(4*h0^2*alpha^2+48*sigma^2+27+4*h0^2*beta^2-72*sigma)*sinh(k*h0)+sinh(3*k*h0)*(alpha^2+beta^2)^2*(-3+4*sigma)^2;
   %
   C3D = (-24*alpha^2*sigma+9*alpha^2+16*alpha^2*sigma^2+16*beta^2*sigma^2+9*beta^2-24*beta^2*sigma)*cosh(2*k*h0)-16*alpha^2*sigma^2+24*alpha^2*sigma-9*beta^2-9*alpha^2-2*alpha^2*(alpha^2+beta^2)*h0^2+24*beta^2*sigma-16*beta^2*sigma^2-2*h0^2*(alpha^2+beta^2)*beta^2;
   %
   %
   %
   C1 = C1N/C12D;
   C2 = C2N/C12D;
   C3 = C3N/C3D;
   %
   u  = C1*solup0(1) + C2*solvp0(1) + C3*solwp0(1); 
   v  = C1*solup0(2) + C2*solvp0(2) + C3*solwp0(2);
   w  = C1*solup0(3) + C2*solvp0(3) + C3*solwp0(3);
%   u = uph0;
%   v = vph0;
%   w = wph0;
   uz = C1*solup0(4) + C2*solvp0(4) + C3*solwp0(4);
   vz = C1*solup0(5) + C2*solvp0(5) + C3*solwp0(5);
   wz = C1*solup0(6) + C2*solvp0(6) + C3*solwp0(6);
   %
   tauxz = (uz + I*alpha*w)/(1+sigma)/2;
   tauyz = (vz + I*beta *w)/(1+sigma)/2;
   tauzz = ((1-sigma)*wz + sigma*(I*alpha*u + I*beta*v))/((1+sigma)*(1-2*sigma));
   %
elseif k==0;
   u = uph0;
   v = vph0;
   w = wph0;
   tauxz = uph0/(1+sigma)/h0/2;
   tauyz = vph0/(1+sigma)/h0/2;
   tauzz = wph0*(1-sigma)/(1+sigma)/(1-2*sigma)/h0;
else
   %
   % we apply the limit kh>>1
   %
   %
   u = uph0;
   v = vph0;
   w = wph0;
   % % % THIS LIMIT IS NOT CORRECT % % %
%   tauxz = (1/2)*(-4*alpha^2-3*beta^2+4*alpha^2*sigma+4*beta^2*sigma)*uph0/((-3+4*sigma)*(1+sigma)*sqrt(alpha^2+beta^2))+(1/2*I)*(4*wph0*sigma*sqrt(alpha^2+beta^2)+I*vph0*beta-2*wph0*sqrt(alpha^2+beta^2))*alpha/((-3+4*sigma)*(1+sigma)*sqrt(alpha^2+beta^2));
%   tauyz = -(1/2)*uph0*alpha*beta/((-3+4*sigma)*(1+sigma)*sqrt(alpha^2+beta^2))+(1/2)*(-3*vph0*alpha^2+4*vph0*alpha^2*sigma-(2*I)*wph0*beta*sqrt(alpha^2+beta^2)-4*vph0*beta^2+4*vph0*beta^2*sigma+(4*I)*wph0*beta*sigma*sqrt(alpha^2+beta^2))/((-3+4*sigma)*(1+sigma)*sqrt(alpha^2+beta^2));
%   tauzz = 0;
   uz = ((-3+4*sigma)*k^2-alpha^2)*uph0/(k*(-3+4*sigma))-alpha*beta*vph0/(k*(-3+4*sigma))+I*alpha*wph0/(-3+4*sigma);
   vz = -alpha*beta*uph0/(k*(-3+4*sigma))+((-3+4*sigma)*k^2-beta^2)*vph0/(k*(-3+4*sigma))+I*beta*wph0/(-3+4*sigma);
   wz = I*uph0*alpha/(-3+4*sigma)+I*vph0*beta/(-3+4*sigma)+2*k*(-1+2*sigma)*wph0/(-3+4*sigma);
   tauxz = (uz + I*alpha*w)/(1+sigma)/2;
   tauyz = (vz + I*beta *w)/(1+sigma)/2;
   tauzz = ((1-sigma)*wz + sigma*(I*alpha*u + I*beta*v))/((1+sigma)*(1-2*sigma));
end

