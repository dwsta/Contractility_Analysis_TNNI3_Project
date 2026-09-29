function [LinEsol, NHsol] = runoptimNHLinE(UF,VF,UBF,VBF,dV,X,Y,inds,runLinE,runNH)
% Setup and run linear Elasticity (E,nu) model and/or Neo-Hookean model
% (G,K) for monolayer
% INPUTS
% UF,VF - Cell displacements
% UBF, VBF - Bead displacements
% dV - Divergence field
% X,Y - Positions for UF/UBF
% inds = Time indices to run the fitting (3rd index) in the UF array

    k=1;
    LinEsol = struct;  NHsol = struct;
    
    for ii = inds

        UU = UF(:,:,ii); VV = VF(:,:,ii);
        UUB = UBF(:,:,ii); VVB = VBF(:,:,ii);
        dVV = dV(:,:,ii);
        XX = X; YY = Y;

        deltaS = XX(1,2) - XX(1,1);
        % 2D strain matrix (epsilon)
        [~,~,dUdX1,dVdX2] = FourierFiltAndDiff(XX,YY,UU,VV);
        [~,~,dVdX1,dUdX2] = FourierFiltAndDiff(XX,YY,VV,UU);
        
        
        
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
        [~,~,dExxdX1,dExxdX2] = FourierFiltAndDiff(XX,YY,Exx,Exx);
        [~,~,dEyydX1,dEyydX2] = FourierFiltAndDiff(XX,YY,Eyy,Eyy);
        [~,~,dExydX1,dExydX2] = FourierFiltAndDiff(XX,YY,Exy,Exy);
        dE.dExxdX1 = dExxdX1;
        dE.dExxdX2 = dExxdX2;
        dE.dEyydX1 = dEyydX1;
        dE.dEyydX2 = dEyydX2;
        dE.dExydX1 = dExydX1;
        dE.dExydX2 = dExydX2;

        uBeads = [reshape(UUB,[],1); reshape(VVB,[],1)];
        
        if runLinE 
            [solE,exitFlagE,fvalE,uPredictE] = runoptimsolveElasto_Eandnu([],...
                                        deltaS,dE,uBeads,{XX,YY},@objFuncDefElasto_Eandnu_reform1_FourierTractions);
            uuXE = reshape(uPredictE(1:end/2),size(XX)); uuYE = reshape(uPredictE(1+end/2:end),size(XX));
            LinEsol(k).sol = solE;
            LinEsol(k).fval = fvalE;
            LinEsol(k).uuX = uuXE;
            LinEsol(k).uuY = uuYE;
            LinEsol(k).exitflag = exitFlagE;
        end 
        
        if runNH
            [sol,exitflag,fval,uPredict] = runoptimsolveNH([],...
                                    deltaS,deltaS,shTerm,bulkTerm,uBeads,{XX,YY});

            uuX = reshape(uPredict(1:end/2),size(XX)); uuY = reshape(uPredict(1+end/2:end),size(XX));
        
            NHsol(k).sol = sol;
            NHsol(k).fval = fval;
            NHsol(k).uuX = uuX;
            NHsol(k).uuY = uuY;
            NHsol(k).exitflag = exitflag;
        end
        
        k=k+1;
    end
end