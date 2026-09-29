function [sol,exitflag,fval,uPredict] = runoptimsolveElasto_Eandnu(M,deltaS,deps,uBeads,inds,objFuncDef)
    E = optimvar('E',1);
    nu = optimvar('nu',1,'UpperBound',0.5,'LowerBound',-1);

    objFunc = @(E,nu) objFuncDef(E,nu,M,deltaS,deps,uBeads,inds);

    obj = fcn2optimexpr( objFunc, E, nu);
    prob = optimproblem('Objective',obj);

    x0 = struct();
    x0.E = 1e8;
    x0.nu = 0.45;

    options = prob.optimoptions;
    options.Display = 'none';
    options.OptimalityTolerance = 1e-8;

%     tic;
    [sol,fval,exitflag,output] = solve(prob,x0,'Options',options);
%     toc;
        
    [~,uPredict] = objFuncDef(sol.E,sol.nu,M,deltaS,deps,uBeads,inds);
    
end
