function solveContestBinary
    p.Aeq = computeAeq;
    p.beq = ones(345,1);
    n    = 345;
    myOptions = optimoptions(@fmincon,'algorithm','active-set','Display','iter','PlotFcns',@optimplotfval);
    p.lb = 0*ones(1,n);
    p.ub = 1*ones(1,n);
    p.nonlcon = @(x) constraint(x);
    p.objective = @(x) cost(x);
    p.options = myOptions;
    p.solver = 'fmincon';
    p.x0 = rand(37*345,1);
    x = fmincon(p);
end

function [c,ceq] = constraint(x)
    
    x = obtainIntengerVariable(x);
    [w,v1,v2] = ISCSO_2021(x',0);
    c = max([v1,v2])
    ceq = [];
end

function ix = obtainIntengerVariable(x)
    x = reshape(x,37,345);
    [x,ix] = max(x);
    ix = ix';
end

function c = cost(x)
    x = obtainIntengerVariable(x);
    [w,v1,v2] = ISCSO_2021(x',0);
    c = w
end

function Aeq = computeAeq
    nSections = 37;
    nel = 345;
    Aeq = sparse(nel,nSections*nel);
    nValues = ones(1,nSections);
    
    for iElem = 1:nel
        n1 = (iElem-1)*nSections + 1;
        n2 = (n1-1) + nSections;
        Aeq(iElem,n1:n2) = nValues;
    end

end