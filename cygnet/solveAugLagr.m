%% Augmented Langrangian method %%

function solveAugLagr
    [sec]   = computeSectionData;
    run("auglagrData.m");
    [w,~,~] = ISCSO_2021(x0,0);
    w0      = w;
    wNew    = w0;
    monitor = figure(1);
    iter = 0;
    [c,c_u,c_sig,LaNew] = computeParameters(s0,l_u,l_sig,rho_u,rho_sig);
    [~,~,~,dLa]         = computeParametersGradient(s0,sec,l_u,l_sig,rho_u,rho_sig,c_u,c_sig);
    [cost, const_u, const_sig, lu, lsig, tauV] = computeMonitoring(c,c_u,c_sig,l_u,l_sig,tau,w0,iter,s0,...
            cost,const_u,const_sig,lu,lsig,tauV,monitor);
    iter = iter + 1;
    while dS > tol || c_u > 0 || c_sig > 0
        displayInfo(iter,c,c_u,c_sig);
        [tau, c_u, c_sig, s]  = determineStepLength(s0,l_u,l_sig,rho_u,rho_sig,tau,LaNew,dLa);
        [l_u, l_sig]          = updateMultipliers(l_u, l_sig, rho_u, rho_sig, c_u, c_sig);
        [c,c_u,c_sig,LaTrial] = computeParameters(s,l_u,l_sig,rho_u,rho_sig);
        [~,~,~,dLa]           = computeParametersGradient(s,sec,l_u,l_sig,rho_u,rho_sig,c_u,c_sig);
        [cost, const_u, const_sig, lu, lsig, tauV] = computeMonitoring(c,c_u,c_sig,l_u,l_sig,tau,w0,iter,s,...
            cost, const_u, const_sig, lu, lsig, tauV, monitor);
        cost_gradient = (c - wNew)/dS;
        wNew          = c;
        dS            = norm(s-s0)/norm(s0);
        s0            = s;
        LaNew         = LaTrial;
        iter          = iter + 1;
    end
end

function displayInfo(iter,c,c_u,c_sig)
     disp('Iteration: ');
     disp(num2str(iter));
     disp('Weight: ');
     disp(num2str(c));
     disp('Displacement constraint: ');
     disp(num2str(c_u));
     disp('Stress constraint: ');
     disp(num2str(c_sig));
end

function [tau, c_u, c_sig, s] = determineStepLength(s0, l_u,l_sig,rho_u,rho_sig,tau,LaNew,dLa)
        if tau < 1e-8
            tau = 0.1*norm(s0)/norm(dLa);
            else
            tau = 10*tau;
        end
        rightStepLength = 0;
        while rightStepLength == 0
            smean = s0 - tau*dLa;
            s     = max(0,min(1,smean));
            [c,c_u,c_sig,LaTrial] = computeParameters(s,l_u,l_sig,rho_u,rho_sig);
            if LaTrial < LaNew
                rightStepLength = 1;
            else
                tau = tau/2;
                if tau < 1e-8
                    rightStepLength = 1;
                end
            end
        end
    
end

function x = calculateSectionID(s)
    x = s*(37-1) + 1;
    x = round(x);
end

function [c,c_u,c_sig,La] = computeParameters(s,l_u,l_sig,rho_u,rho_sig)
    x = calculateSectionID(s);
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


function [l_u, l_sig] = updateMultipliers(l_u, l_sig, rho_u, rho_sig, c_u, c_sig)
    l_u0   = l_u;
    l_sig0 = l_sig;
    l_u    = l_u0   + rho_u*max(-l_u0/rho_u,c_u);
    l_sig  = l_sig0 + rho_sig*max(-l_sig0/rho_sig,c_sig);
end

function [cost, const_u, const_sig, lu, lsig, tauV] = computeMonitoring(c, c_u, c_sig, l_u, l_sig, tau, w0, iter, s,...
    cost, const_u, const_sig, lu, lsig, tauV, monitor)
        cost(iter+1)      = c;
        const_u(iter+1)   = c_u;
        const_sig(iter+1) = c_sig;
        lu(iter+1)        = l_u;
        lsig(iter+1)      = l_sig;
        tauV(iter+1)      = tau;
        x                 = calculateSectionID(s);
        auglagr_monitoring(monitor, 1:(iter+1), x, cost, const_sig, const_u, lu, lsig ,tauV, w0);
end

function [dC, dV1, dV2] = calculateFiniteGradient(s)
    x0 = calculateSectionID(s);
    [c0, v10, v20] = ISCSO_2021(x0,0);
    for i = 1:numel(s)
        s(i) = s(i)+1/37;
        x = calculateSectionID(s);
        [c, v1, v2] = ISCSO_2021(x,0);
        ds = 1/37;
        dC(i)  = (c-c0) / (ds);
        dV1(i) = (v1-v10) / (ds);
        dV2(i) = (v2-v20) / (ds);
    end
end

