function solveContest
%     plotHomogeneousSolution();
    solveProblem();
end

function solveProblem()
    n    = 345;
    s0 = 0.99*ones(1,345);
    sec = computeSectionData();
    myOptions = optimoptions(@fmincon, 'algorithm', 'interior-point', ... 
                'SpecifyConstraintGradient', true, ...
                'SpecifyObjectiveGradient', true, ...
                'Display','iter','PlotFcns', @optimplotfval);
    p.lb = 0*ones(1,n);
    p.ub = 1*ones(1,n);
    p.nonlcon = @(x) constraint(x, s0,sec);
    p.objective = @(x) cost(x, s0,sec);
    p.options = myOptions;
    p.solver = 'fmincon';
    p.x0 = s0;
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

function A = computeSection(s,sec)
    sections = sec;
    Amax     = max(sections);
    Amin     = min(sections)/Amax;
    Amax     = 1;
    p        = 3;
    A        = Amin * (1-s.^p) + s.^p * Amax;
end

function gradA = computeSectionGradient(s,sec)
    sections = sec;
    Amax     = max(sections);
    Amin     = min(sections)/Amax;
    Amax     = 1;
    p        = 3;
    gradA    = - Amin * p * s.^(p-1) + Amax * p*s.^(p-1);
end

function [c, dC] = cost(s,s0,sec)
    x = s*(37-1) + 1;
    x = round(x);
    [w,v1,v2] = ISCSO_2021(x,0);
    c = w;
%     dA = computeSectionGradient(s);
    [dC, ~, ~] = computeFiniteDiffDerivatives(s,s0,sec);
%     dC = dA;
end

function [c,ceq,dC,dceq] = constraint(s,s0,sec)
%[c,ceq,DC,DCeq]
    x = s*(37-1) + 1;
    x = round(x);
    [w,v1,v2] = ISCSO_2021(x,0);
    c = [v1,v2];
    ceq = [];
%     dA = computeSectionGradient(s);
    [~, dc_u, dc_sig] = computeFiniteDiffDerivatives(s,s0,sec);
    dC = [dc_u; dc_sig]';
%     dC = [1./dA; 1./dA]';

    dceq = [];
end

function [dc, dc_u, dc_sig] = computeFiniteDiffDerivatives(s,s0,sec)
    A      = computeSection(s,sec);
    A0     = computeSection(s0,sec);
    dA     = computeSectionGradient(s,sec);
    dA0    = computeSectionGradient(s0,sec);
    A = A./A0;
    dA = dA./dA0;
    dc     = load("dC.mat", "dC");
    dc = dc.dC;
    dc_sig = load("dSig.mat", 'dV1');
    dc_sig = dc_sig.dV1;
    dc_u   = load("dU.mat", 'dV2');
    dc_u = dc_u.dV2;
    dc     = dc.*dA;
    dc_u   = dc_u.*dA./A.^2;
    dc_sig = dc_sig.*dA./A.^2;
end