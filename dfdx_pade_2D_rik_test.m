Nx = 50;
Ny = 30;
x = linspace(0,2*pi,Nx);
y = linspace(0,2*pi,Ny);

[X,Y] = meshgrid(x,y);
dx = x(2)-x(1);
dy = y(2)-y(1);

u = sin(4*X).*cos(2*Y);

dudx = dfdx_pade_2D_rik(u,dx);
dudy = dfdy_pade_2D_rik(u,dy);

figure(1);
subplot(3,1,1)
imagesc(u,'XData',x,'YData',y); colorbar
title('f(x,y)')
subplot(3,1,2)
imagesc(dudx,'XData',x,'YData',y); colorbar
title('df(x,y)/dx')
subplot(3,1,3)
imagesc(dudy,'XData',x,'YData',y); colorbar
title('df(x,y)/dy')