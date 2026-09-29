function [costVal,uu] = objFuncDefNH(K,G,M,deltaX,deltaY,shTerm,bulkTerm,uBeads,delInds)
    
    % Compressible N-H formulation in terms of First Piola-Kirchoff tensor
    % P = mu * (F - F^-T) + lam * (J-1) * J * F^-T 
    % P = mu * shTerm + lam * bulkTerm;
    
    N = numel(uBeads)./2;
    
    
    P = G * shTerm + K * bulkTerm;
    
    P11 = squeeze( P(:,:,1,1) ) ;
    P12 = squeeze( P(:,:,1,2) ) ;
    P21 = squeeze( P(:,:,2,1) ) ;
    P22 = squeeze( P(:,:,2,2) ) ;
    
    [dP11dx,~] = Func_1stDer_PadesScheme(P11,deltaX,deltaY); % sigma_xx
    [~,dP22dy] = Func_1stDer_PadesScheme(P22,deltaX,deltaY); % sigma_yy
    [~,dP12dy] = Func_1stDer_PadesScheme(P12,deltaX,deltaY); % sigma_xy
    [dP21dx,~] = Func_1stDer_PadesScheme(P21,deltaX,deltaY); % sigma_yx
    
    dP11dx = imgaussfilt(dP11dx,0.5); dP22dy = imgaussfilt(dP22dy,0.5);
    dP21dx = imgaussfilt(dP21dx,0.5); dP12dy = imgaussfilt(dP12dy,0.5);
    
    ttx = -(dP11dx+dP12dy);
    tty = -(dP22dy+dP21dx);

    uu =  M * [ttx(:); tty(:)] *(10/0.35);
    %     uu = [ttx(:); tty(:)] * 8e3; 
    
    uu(delInds) = nan; % To remove contribution from close to 0 cell deformation
    
    %     eps = 1e-3;
    
    %     costVal = nansum( abs( (uBeads - uu).^ 2 ./ (uBeads + uu + eps).^2 ));
    costVal = nanmean( ( uBeads - uu ).^2 );
    
end
