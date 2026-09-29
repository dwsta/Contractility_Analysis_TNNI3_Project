% Use Scheme from Parvif Moin (Ch 2.4)
% Pade's 4th order approximations:
% fp_j+1 + fp_j-1 + fp_j = (3/h)*(f_j+1-f_j-1)
%
% Boundaries 3rd order
% fp_0 + 2 fp_1 = (1/h) * (-5/2 f0 + 2f1 + 1/2 f2)
% fp_n + 2fp_n-1 = (1/h)* (5/2 fn - 2 fn - 1/2 f_n-2)

function [dfdx,dfdy] = Func_1stDer_PadesScheme(f,dx,dy)

    Ny = size(f,2); Nx = size(f,1);

    % Calculate df/d(column index)
    % [1 2 0  ...  0]
    % [1 4 1  ...  0]
    % [0 1 4 1 ... 0]
    % [0 0 ...  ...0]
    % [0 0 ... 1 4 1]
    % [0 0 ... 0 2 1] for each row or Y

    DBlock1 = sparse(1:Nx,1:Nx, [1 4*ones(1,Nx-2) 1],Nx,Nx) + sparse(1:Nx-1,2:Nx, ...
        [2 ones(1,Nx-2)],Nx,Nx)  + ...
                sparse(2:Nx,1:Nx-1, [ones(1,Nx-2) 2],Nx,Nx);
    D = [];

    fX = reshape(f',numel(f),1);
    fY = f(:);
    bx = nan(numel(f),1);
    by = nan(numel(f),1);

    % This loop works only if Nx = Ny
    for i = 1 : Ny 
        D = blkdiag(D,DBlock1);
        fT = fX([1:Nx]+(i-1)*Ny);
        bx([1:Nx]+(i-1)*Ny,1) = [ -5/2*fT(1) + 2 * fT(2) + 0.5*fT(3); 3 * ( fT(3:end) - fT(1:end-2) ); ...
                                        (5/2)*fT(Ny)-2*fT(Ny-1)-0.5*fT(Ny-2) ] ;
        fT1 = fY([1:Nx]+(i-1)*Ny);
        by([1:Nx]+(i-1)*Ny,1) = [ -5/2*fT1(1) + 2 * fT1(2) + 0.5*fT1(3); 3 * ( fT1(3:end) - fT1(1:end-2) ); ...
            (5/2)*fT1(Ny)-2*fT1(Ny-1)-0.5*fT1(Ny-2) ] ;
    end
    bx = bx ./ dx;
    by = by ./ dy;

    dfdx = reshape( D \ bx , Ny, Nx)' ;
    dfdy = reshape( D \ by , Ny, Nx);

end



