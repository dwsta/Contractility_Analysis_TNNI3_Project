function [LinEsol, NHsol] = runoptimNHLinE_Pades(UF,VF,UBF,VBF,dV,X,Y,inds,runLinE,runNH)
% Setup and run linear Elasticity (E,nu) model and/or Neo-Hookean model
% (G,K) for monolayer
% INPUTS
% UF,VF - Cell displacements
% UBF, VBF - Bead displacements
% dV - Divergence field
% X,Y - Positions for UF/UBF
% inds = Time indices to run the fitting (3rd index) in the UF array

    LinEsol = struct;  NHsol = struct;
    
    for i = 1 : size(UF,3)
        
        UU = UF(:,:,i); VV = VF(:,:,i);
        UUB = UBF(:,:,i); VVB = VBF(:,:,i);
        dVV = dV(:,:,i);
        XX = X; YY = Y;

        dx = XX(1,2) - XX(1,1);
        dy = YY(2,1) - YY(1,1);
        
        % 2D strain matrix (epsilon) 
%         [dUdX1,dUdX2] = Func_1stDer_PadesScheme(squeeze(UBF(:,:,i)),deltaS,deltaS);
%         [dVdX1,dVdX2] = Func_1stDer_PadesScheme(squeeze(VBF(:,:,i)),deltaS,deltaS);
        dUdX1 = dfdx_pade_2D_rik(squeeze(UBF(:,:,i)),dx);
        dUdX2 = dfdy_pade_2D_rik(squeeze(UBF(:,:,i)),dy);
        dVdX1 = dfdx_pade_2D_rik(squeeze(VBF(:,:,i)),dx);
        dVdX2 = dfdy_pade_2D_rik(squeeze(VBF(:,:,i)),dy);
        

        % F - deformation gradient tensor
        % J - Jacobi matrix (|F|)
        % shTerma and bulkTerm are multiples of G and K in th energy term
        % for NH model
        F = nan([ size(UU,[1 2]) 2 2]);
        shTerm = nan([ size(UU,[1 2]) 2 2]);
        bulkTerm = nan([ size(UU,[1 2]) 2 2]);
        J = nan([ size(UU,[1 2])]);
        F11N = 1 + dUdX1; F12N = dUdX2; 
        F21N = dVdX1; F22N = 1 + dVdX2;
        F(:,:,1,1) = F11N;
        F(:,:,1,2) = F12N;
        F(:,:,2,1) = F21N;
        F(:,:,2,2) = F22N;

        for ii = 1:size(F11N,1)
            for jj = 1 : size(F11N,2)
                FT = [F11N(ii,jj) F12N(ii,jj); F21N(ii,jj) F22N(ii,jj)];
                JT = det(FT);
                J(ii,jj) = JT;
                FinvT = inv(FT).';
                shTerm(ii,jj,:,:) = FT-FinvT;
                bulkTerm(ii,jj,:,:) = (JT-1)*JT*FinvT;
            end
        end

        % Derivtive of Strain matrix
        dE = struct();
        Exx = dUdX1; Eyy = dVdX2; Exy = (dUdX2+dVdX1)./2;
%         [dExxdX1,dExxdX2] = Func_1stDer_PadesScheme(Exx,deltaS,deltaS);
%         [dEyydX1,dEyydX2] = Func_1stDer_PadesScheme(Eyy,deltaS,deltaS);
%         [dExydX1,dExydX2] = Func_1stDer_PadesScheme(Exy,deltaS,deltaS);
        
        dExxdX1 = dfdx_pade_2D_rik(Exx,dx);
        dExxdX2 = dfdy_pade_2D_rik(Exx,dy);
        
        dEyydX1 = dfdx_pade_2D_rik(Eyy,dx);
        dEyydX2 = dfdy_pade_2D_rik(Eyy,dy);
        
        dExydX1 = dfdx_pade_2D_rik(Exy,dx);
        dExydX2 = dfdy_pade_2D_rik(Exy,dy);
        
        dE.dExxdX1 = dExxdX1;
        dE.dExxdX2 = dExxdX2;
        dE.dEyydX1 = dEyydX1;
        dE.dEyydX2 = dEyydX2;
        dE.dExydX1 = dExydX1;
        dE.dExydX2 = dExydX2;

        uBeads = [reshape(UUB,[],1); reshape(VVB,[],1)];
        
        if runLinE 
            [solE,exitFlagE,fvalE,uPredictE] = runoptimsolveElasto_Eandnu([],...
                                        dx,dE,uBeads,{XX,YY},@objFuncDefElasto_Eandnu_reform1_FourierTractions);
            uuXE = reshape(uPredictE(1:end/2),size(XX)); uuYE = reshape(uPredictE(1+end/2:end),size(XX));
            LinEsol(i).sol = solE;
            LinEsol(i).fval = fvalE;
            LinEsol(i).uuX = uuXE;
            LinEsol(i).uuY = uuYE;
            LinEsol(i).exitflag = exitFlagE;

            % Calculate the stress tensor
            
            sxx = (solE.E./(1-solE.nu.^2)) .* (Exx + solE.nu .* Eyy);
            syy = (solE.E./(1-solE.nu.^2)) .* (Eyy + solE.nu .* Exx);
            sxy = (solE.E./(1+solE.nu)) .* (Exy);
            
            LinEsol(i).sxx = sxx;
            LinEsol(i).syy = syy;
            LinEsol(i).sxy = sxy;
        end 
        
        if runNH
            [sol,exitflag,fval,uPredict] = runoptimsolveNH_Pades([],...
                                    dx,dx,shTerm,bulkTerm,uBeads,{XX,YY});

            uuX = reshape(uPredict(1:end/2),size(XX)); uuY = reshape(uPredict(1+end/2:end),size(XX));
        
            NHsol(i).sol = sol;
            NHsol(i).fval = fval;
            NHsol(i).uuX = uuX;
            NHsol(i).uuY = uuY;
            NHsol(i).exitflag = exitflag;
            
            % Calculate the stress tensor (First Piola-Kirchoff tensor)
            
            Pxx = sol.G * squeeze(shTerm(:,:,1,1)) + sol.K * squeeze( shTerm(:,:,1,1) ) ;
            Pyy = sol.G * squeeze(shTerm(:,:,2,2)) + sol.K * squeeze( shTerm(:,:,2,2) ) ;
            Pxy = sol.G * squeeze(shTerm(:,:,1,2)) + sol.K * squeeze( shTerm(:,:,1,2) ) ;
            Pyx = sol.G * squeeze(shTerm(:,:,2,1)) + sol.K * squeeze( shTerm(:,:,2,1) ) ;
            NHsol(i).Pxx = Pxx;
            NHsol(i).Pyy = Pyy;
            NHsol(i).Pxy = 0.5*(Pxy + Pyx);
        end
        
    end
end