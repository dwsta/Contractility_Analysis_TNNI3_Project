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

x = linspace(0,2*pi);
f = sin(8*x);
dx = x(2)-x(1);
Nx = length(x);

LHS =  diag(4*ones(1,Nx)) + diag(ones(1,Nx-1),1) + diag(ones(1,Nx-1),-1);
LHS(1,1) = 1;
LHS(1,2) = 2;
LHS(end,end) = 1;
LHS(end,end-1) = 2;

RHS(1)  = 1/dx * (-5/2 * f(1) + 2 * f(2) + 1/2 * f(3));
RHS(Nx) = 1/dx * (5/2 * f(Nx) - 2 * f(Nx-1) - 1/2 * f(Nx-2));
RHS(2:Nx-1) = 3/dx * (f(3:Nx)-f(1:Nx-2));

dfdx = inv(LHS)*RHS';
figure(1)
plot(x,f,'b');
hold on
plot(x,dfdx,'r')
legend({'f(x)','df(x)/dx'})