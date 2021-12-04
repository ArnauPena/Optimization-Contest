function solveContest
%     plotHomogeneousSolution();
    solveProblem();
end

function solveProblem()
    n    = 345;
    myOptions = optimoptions(@fmincon, 'algorithm', 'active-set', ... 
                'SpecifyConstraintGradient', true, ...
                'SpecifyObjectiveGradient', true, ...
                'Display','iter','PlotFcns', @optimplotfval);
    p.lb = 0*ones(1,n);
    p.ub = 1*ones(1,n);
    p.nonlcon = @(x) constraint(x);
    p.objective = @(x) cost(x);
    p.options = myOptions;
    p.solver = 'fmincon';
    p.x0 = 0.99*ones(1,345);
    x = fmincon(p);
    p.objective(x)
    p.nonlcon(x)
end

function plotHomogeneousSolution()
    sects = 1:37;
    for i = sects
        x0 = i*ones(1,345);
        [w(i), v1(i), v2(i)] = ISCSO_2021(x0,0);
    end
    figure();
    plot (sects, w);
    figure();
    plot (sects, [v1; v2]);
end

function [Si] = computeData()
    run("DataSections.m")
end

function [sections] = computeSectionData()
    Si = computeData();
    sections = Si(:,1);
end

function A = computeSection(s)
    sections = computeSectionData();
    Amin = min(sections);
    Amax = max(sections);
    p = 3;
    A = Amin * (1-s.^p) + s.^p * Amax;
end

function gradA = computeSectionGradient(s)
    sections = computeSectionData();
    Amin = min(sections);
    Amax = max(sections);
    p = 3;
    gradA = - Amin * p * s.^(p-1) + Amax * p*s.^(p-1);
end

function [c, dC] = cost(s)
    x = s*(37-1) + 1;
    x = round(x);
    [w,v1,v2] = ISCSO_2021(x,0);
    c = w;
    dA = computeSectionGradient(s);
    dC = dA;
end

function [c,ceq,dC,dceq] = constraint(s)
%[c,ceq,DC,DCeq]
    x = s*(37-1) + 1;
    x = round(x);
    [w,v1,v2] = ISCSO_2021(x,0);
    c = [v1,v2];
    ceq = [];
    dA = computeSectionGradient(s);
    dC = [1./dA; 1./dA]';
    dceq = [];
end