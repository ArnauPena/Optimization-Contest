%% Augmented Langrangian method %%

function solveAugLagr
    [sec]   = computeSectionData;
    tol     = 1e-3;
    l_u     = 30;
    l_sig   = 30;
    rho_u   = 1;
    rho_sig = 1;
    s0      = 0.99*ones(1,345);
    dS      = 5;
    [c,c_u,c_sig,La0]    = computeParameters(s0,l_u,l_sig,rho_u,rho_sig);
    [dc,dc_u,dc_sig,dLa] = computeParametersGradient(s0,sec,l_u,l_sig,rho_u,rho_sig,c_u,c_sig);
    iter = 0;
    while dS > tol || c_u > 0 || c_sig > 0
        disp(num2str(c));
        disp(num2str(c_u));
        disp(num2str(c_sig));
        if iter == 0
            tau = norm(s0)/norm(dLa);
        else
            tau = 10*tau;
        end
        rightStepLength = 0;
        while rightStepLength == 0
            smean = s0 - tau*dLa;
            for isec = 1:length(smean)
                s(isec)   = max(0,min(1,smean(isec)));
            end
            [c,c_u,c_sig,La] = computeParameters(s,l_u,l_sig,rho_u,rho_sig);
            if La < La0
                rightStepLength = 1;
            else
%                 La0 = La;
                tau = tau/2;
                if tau < 1e-20
                    rightStepLength = 1;
                end
            end
        end
        l_u0   = l_u;
        l_sig0 = l_sig;
        l_u    = l_u0 + rho_u*max(-l_u0/rho_u,c_u);
        l_sig  = l_sig0 + rho_sig*max(-l_sig0/rho_sig,c_sig);
        [c,c_u,c_sig,La] = computeParameters(s,l_u,l_sig,rho_u,rho_sig);
        [~,~,~,dLa] = computeParametersGradient(s,sec,l_u,l_sig,rho_u,rho_sig,c_u,c_sig);
        dS = norm(s-s0)/norm(s0);
        s0 = s;
        La0 = La;
    end
end

function [c,c_u,c_sig,La] = computeParameters(s,l_u,l_sig,rho_u,rho_sig)
    x = s*(37-1) + 1;
    x = round(x);
    [w,v1,v2] = ISCSO_2021(x,0);
    c     = w;
    c_sig = v1;
    c_u   = v2;
    La    = c + l_u*c_u + 0.5*rho_u*c_u^2 + l_sig*c_sig + 0.5*rho_sig*c_sig^2;        
end

function [dc,dc_u,dc_sig,dLa] = computeParametersGradient(s,sec,l_u,l_sig,rho_u,rho_sig,c_u,c_sig)
    dA     = computeSectionGradient(s,sec);
    A      = computeSection(s,sec);
    dc     = dA;
    dc_u   = -dA./A.^2;
    dc_sig = -dA./A.^2;
    dLa    = dc + (l_u + rho_u*c_u).*dc_u + (l_sig + rho_sig*c_sig).*dc_sig;
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
    Amin     = min(sections);
    Amax     = max(sections);
    p        = 3;
    A        = Amin * (1-s.^p) + s.^p * Amax;
end

function gradA = computeSectionGradient(s,sec)
    sections = sec;
    Amin     = min(sections);
    Amax     = max(sections);
    p        = 3;
    gradA    = - Amin * p * s.^(p-1) + Amax * p*s.^(p-1);
end
