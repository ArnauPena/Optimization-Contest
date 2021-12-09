%% Augmented Langrangian method %%

function solveAugLagr
clear
clc
close all
    [sec]   = computeSectionData;
    method = 'finite differences';
    run("auglagrData.m");
    [w,~,~] = ISCSO_2021(x0,0);
    w0      = w;
    wNew    = w0;
    monitor = figure(1);
    iter = 0;
    s_ini = s0;
    [c,c_u,c_sig,LaNew] = computeParameters(s0,l_u,l_sig,rho_u,rho_sig);
    dC = 0;
    dU = 0;
    dSig = 0;
    update = 1;
    [dc,dc_u,dc_sig,dLa,dC,dU,dSig]         = computeParametersGradient(s0,s_ini,sec,l_u,l_sig,rho_u,rho_sig,c_u,c_sig,method,update,dC,dU,dSig);
    [cost, const_u, const_sig, lu, lsig, tauV] = computeMonitoring(c,c_u,c_sig,l_u,l_sig,tau,w0,iter,s0,...
            cost,const_u,const_sig,lu,lsig,tauV,monitor,dc, dc_u, dc_sig);
    iter = iter + 1;
    cNew = c;
    while dS > tol || c_u > 0 || c_sig > 0
        displayInfo(iter,c,c_u,c_sig);
        multiple      = (iter/10);
        multiple      = multiple/fix(multiple);
        update        = multiple;
        [tau, c_u, c_sig, s]   = determineStepLength(s0,l_u,l_sig,rho_u,rho_sig,tau,LaNew,dLa,cNew);
        [l_u, l_sig]           = updateMultipliers(l_u, l_sig, rho_u, rho_sig, c_u, c_sig);
        [c,c_u,c_sig,LaTrial]  = computeParameters(s,l_u,l_sig,rho_u,rho_sig);
        [dc,dc_u,dc_sig,dLa,dC,dU,dSig] = computeParametersGradient(s,s_ini,sec,l_u,l_sig,rho_u,rho_sig,c_u,c_sig,method,update,dC,dU,dSig);
        [cost, const_u, const_sig, lu, lsig, tauV] = computeMonitoring(c,c_u,c_sig,l_u,l_sig,tau,w0,iter,s,...
            cost, const_u, const_sig, lu, lsig, tauV, monitor,dc,dc_u,dc_sig);
        dS            = norm(s-s0)/norm(s0);
        s0            = s;
        cNew          = c;
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

function [tau, c_u, c_sig, s] = determineStepLength(s0, l_u,l_sig,rho_u,rho_sig,tau,LaNew,dLa,cNew)
%         if tau < 1e-8
%             tau = 0.05*norm(s0)/norm(dLa);
            tau = 1/max(abs(dLa))*4/37;
%         else
%             tau = 10*tau;
%         end
        rightStepLength = 0;
        while rightStepLength == 0
            smean = s0 - tau*dLa;
            s     = max(0,min(1,smean));
            [cTrial,c_u,c_sig,LaTrial] = computeParameters(s,l_u,l_sig,rho_u,rho_sig);
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

function [dc,dc_u,dc_sig,dLa,dC,dU,dSig] = computeParametersGradient(s,s0,sec,l_u,l_sig,rho_u,rho_sig,c_u,c_sig,method,update,dC,dU,dSig)
    switch method
        case 'approach'
            dA     = computeSectionGradient(s,sec);
            A      = computeSection(s,sec);
            dc     = dA;
            dc_u   = -dA./A.^2;
            dc_sig = -dA./A.^2;
        case 'finite differences'
            if update == 1
                [dc, dc_u, dc_sig,dC,dU,dSig] = updateFiniteDiffDerivatives(s,s0,sec);
            else
                [dc, dc_u, dc_sig] = computeFiniteDiffDerivatives(s,s0,sec,dC,dU,dSig);
            end
%             [dc, dc_sig, dc_u] = calculateFiniteGradient(s);
        otherwise
            error('Method not implemented')
    end
    dLa    = dc + (l_u + rho_u*c_u).*dc_u + (l_sig + rho_sig*c_sig).*dc_sig;
end

function [dc, dc_u, dc_sig] = computeFiniteDiffDerivatives(s,s0,sec,dC,dU,dSig)
    A      = computeSection(s,sec);
    A0     = computeSection(s0,sec);
    dA     = computeSectionGradient(s,sec);
    dA0    = computeSectionGradient(s0,sec);
    A = A./A0;
    dA = dA./dA0;
    dc     = dC.*dA;
    dc_u   = dU.*dA./A.^2;
    dc_sig = dSig.*dA./A.^2;
end

function [dc, dc_u, dc_sig,dC,dU,dSig] = updateFiniteDiffDerivatives(s,s0,sec)
    A      = computeSection(s,sec);
    A0     = computeSection(s0,sec);
    dA     = computeSectionGradient(s,sec);
    dA0    = computeSectionGradient(s0,sec);
    A      = A./A0;
    dA     = dA./dA0;
    [dC, dSig, dU] = calculateFiniteGradient(s);
    dc     = dC.*dA;
    dc_u   = dU.*dA./A.^2;
    dc_sig = dSig.*dA./A.^2;
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


function [l_u, l_sig] = updateMultipliers(l_u, l_sig, rho_u, rho_sig, c_u, c_sig)
    l_u0   = l_u;
    l_sig0 = l_sig;
    l_u    = l_u0   + rho_u*max(-l_u0/rho_u,c_u);
    l_sig  = l_sig0 + rho_sig*max(-l_sig0/rho_sig,c_sig);
end

function [cost, const_u, const_sig, lu, lsig, tauV] = computeMonitoring(c, c_u, c_sig, l_u, l_sig, tau, w0, iter, s,...
    cost, const_u, const_sig, lu, lsig, tauV, monitor,dc,dc_u,dc_sig)
        cost(iter+1)      = c;
        const_u(iter+1)   = c_u;
        const_sig(iter+1) = c_sig;
        lu(iter+1)        = l_u;
        lsig(iter+1)      = l_sig;
        tauV(iter+1)      = tau;
        x                 = calculateSectionID(s);
        auglagr_monitoring(monitor, 1:(iter+1), x, cost, const_sig, const_u, lu, lsig ,tauV, w0, dc,dc_u,dc_sig);
end

function [dC, dV1, dV2] = calculateFiniteGradient(s)
    x0 = calculateSectionID(s);
    [c0, v10, v20] = ISCSO_2021(x0,0);
    for i = 1:numel(s)
        s1    = s;
        s2    = s;
        s1(i) = s(i)+1/37;
        s2(i) = s(i)-1/37;

        if s1(i) > 1
            s1(i) = 1;
        elseif s2(i) < 0
            s2(i) = 0;
        end
            x1 = calculateSectionID(s1);
            x2 = calculateSectionID(s2);
            [c1, v1_1, v2_1] = ISCSO_2021(x1,0);
            [c2, v1_2, v2_2] = ISCSO_2021(x1,0);
            ds = 1/37;
            dC1(i)   = (c1-c0) / (ds);
            dV1_1(i) = (v1_1-v10) / (ds);
            dV2_1(i) = (v2_1-v20) / (ds);
            dC2(i)   = (c2-c0) / (ds);
            dV1_2(i) = (v1_2-v10) / (ds);
            dV2_2(i) = (v2_2-v20) / (ds);
            dC(i)    = (dC1(i) + dC2(i)) / 2;
            dV1(i)   = (dV1_1(i) + dV1_2(i)) / 2;
            dV2(i)   = (dV2_1(i) + dV2_2(i)) / 2;
    end
end



