function [sol,exitflag,fval,uPredict] = runoptimsolveNH(M,deltaX,deltaY,shTerm,bulkTerm,uBeads,inds)
%     Kfix=-10^4;
%     K = optimvar('K',1,'UpperBound',Kfix+1,'LowerBound',Kfix-1);
    K = optimvar('K',1);
    G = optimvar('G',1);%,'UpperBound',0.1,'LowerBound',0);

    if ~iscell(inds)
        objFunc = @(K,G) objFuncDefNH(K,G,M,deltaX,deltaY,shTerm,bulkTerm,uBeads,inds);
    else
        objFunc = @(K,G) objFuncDefNH_FourierTraction(K,G,M,deltaX,deltaY,shTerm,bulkTerm,uBeads,inds);
    end

    obj = fcn2optimexpr( objFunc, K ,G );
    prob = optimproblem('Objective',obj);

    x0 = struct();
    x0.K = 1e6;
    x0.G = 1e6;

    options = prob.optimoptions;
    options.Display = 'none';
    options.OptimalityTolerance = 1e-8;
    
    
%     tic;
    [sol,fval,exitflag,output] = solve(prob,x0,'Options',options);
%     toc;
    
    if ~iscell(inds)
        [~,uPredict] = objFuncDefNH(sol.K,sol.G,M,deltaX,deltaY,shTerm,bulkTerm,uBeads,inds);
    else
        [~,uPredict] = objFuncDefNH_FourierTraction(sol.K,sol.G,M,deltaX,deltaY,shTerm,bulkTerm,uBeads,inds);
    end
end
