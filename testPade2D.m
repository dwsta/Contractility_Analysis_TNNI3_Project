% Use Scheme from Parvif Moin (Ch 2.4)
% Let fp_j denote the first derivative, of function f evaluated at point j
%
% Pade's 4th order approximation:
% fp_j+1 + 4*fp_j + fp_j-1  = (3/h)*(f_j+1-f_j-1)
%
% Boundaries 3rd order
% fp_0 + 2 fp_1 = (1/h) * (-5/2 f0 + 2f1 + 1/2 f2)
% fp_n + 2fp_n-1 = (1/h)* (5/2 fn - 2 f_n-1 - 1/2 f_n-2)

% Calculate df/d(column index)
% [1 2 0  ...  0]
% [1 4 1  ...  0]
% [0 1 4 1 ... 0]
% [0 0 ...  ...0]
% [0 0 ... 1 4 1]
% [0 0 ... 0 2 1] for each row or Y

function dfdx = dfdx_1st_derivative_rik_2D(f,dx)
x = linspace(0,2*pi,100);
y = linspace(0,2*pi,8);

[X,Y] = meshgrid(x,y);
f = sin(8*X);
dx = X(1,2)-X(1,1);
dy = Y(2,1)-Y(1,1);

Nx = length(x);
Ny = length(y);

% This would differentiate a row. Need to block repeat to whole matrix
lhs =  diag(4*ones(1,Nx)) + diag(ones(1,Nx-1),1) + diag(ones(1,Nx-1),-1);
lhs(1,1) = 1;
lhs(1,2) = 2;
lhs(end,end) = 1;
lhs(end,end-1) = 2;

% Block diagonal for the left hand side
LHS = kron(eye(Ny),sparse(lhs));

% Right hand side can be computed, per "row"... 
rhs(:,1) = 1/dx * (-5/2 * f(:,1) + 2 * f(:,2) + 1/2 * f(:,3));
rhs(:,Nx) = 1/dx * (5/2 * f(:,Nx) - 2 * f(:,Nx-1) - 1/2 * f(:,Nx-2));
rhs(:,2:Nx-1) = 3/dx * (f(:,3:Nx)-f(:,1:Nx-2));
% ... and then re-ordered as one single column
RHS = reshape(rhs',[],1);

% The solution needs some un-tangling
dfdx = reshape(LHS\RHS,Nx,Ny)';

% Plot to check
figure(1)
subplot(2,1,1)
imagesc(f,'XData',x,'YData',y); colorbar
title('f(x)')
subplot(2,1,2)
imagesc(dfdx,'XData',x,'YData',y); colorbar
title('df(x)/dx')

