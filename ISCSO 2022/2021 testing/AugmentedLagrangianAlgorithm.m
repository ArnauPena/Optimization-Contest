classdef AugmentedLagrangianAlgorithm

    properties (Access = public)
        xFinal
    end

    methods (Access = public)

        function obj = AugmentedLagrangianAlgorithm()
            dataSectionsFile = "DataSections.m";
            [sec]   = obj.computeSectionData(dataSectionsFile);
            method  = 'finite differences';
            run("auglagrData.m");
            monitor = figure(1);
            s_ini   = s0;
            iter    = 0;
            dC      = 0;
            dU      = 0;
            dSig    = 0;
            update  = 1;
            [c,c_u,c_sig,LaNew]             = obj.computeParameters(s0,l_u,l_sig,rho_u,rho_sig);
            w0                              = c;
            [dc,dc_u,dc_sig,dLa,dC,dU,dSig] = obj.computeParametersGradient(s0,s_ini,...
                sec,l_u,l_sig,rho_u,rho_sig,c_u,c_sig,method,update,dC,dU,dSig);
            [cost, const_u, const_sig, lu, lsig, tauV] = obj.computeMonitoring(c,c_u,...
                c_sig,l_u,l_sig,tau,w0,iter,s0,...
                cost,const_u,const_sig,lu,lsig,tauV,monitor,dc, dc_u, dc_sig);
            iter = iter + 1;          
            while dS > tol || c_u > 0 || c_sig > 0
                obj.displayInfo(iter,c,c_u,c_sig);
                if iter == 1
                    [s0,cNew,c_u,c_sig] = obj.makeFeasible(s0, dc_u, dc_sig);
                    [dc,dc_u,dc_sig,dLa,dC,dU,dSig] = obj.computeParametersGradient(s0,...
                        s_ini,sec,l_u,l_sig,rho_u,rho_sig,c_u,c_sig,method,update,dC,dU,dSig);
                end
                update                 = obj.checkUpdate(iter);
                [tau, c_u, c_sig, s]   = obj.determineStepLength(s0,l_u,l_sig,...
                    rho_u,rho_sig,tau,LaNew,dLa,cNew);
                [l_u, l_sig]           = obj.updateMultipliers(l_u, l_sig, rho_u, ...
                    rho_sig, c_u, c_sig);
                [c,c_u,c_sig,LaTrial]  = obj.computeParameters(s,l_u,l_sig,rho_u,...
                    rho_sig);
                [dc,dc_u,dc_sig,dLa,dC,dU,dSig] = obj.computeParametersGradient(s,s_ini,...
                    sec,l_u,l_sig,rho_u,rho_sig,c_u,c_sig,method,update,dC,dU,dSig);
                [cost, const_u, const_sig, lu, lsig, tauV] = obj.computeMonitoring(c,c_u,...
                    c_sig,l_u,l_sig,tau,w0,iter,s,...
                    cost, const_u, const_sig, lu, lsig, tauV, monitor,dc,dc_u,dc_sig);
                [iter,dS,s0,cNew,LaNew] = obj.updateParameters(s,s0,c,LaTrial,iter);
            end
            obj.xFinal = calculateSectionID(s);
        end
    end

    methods (Access = private)

        function [tau, c_u, c_sig, s] = determineStepLength(obj,s0, l_u,l_sig,rho_u,...
                rho_sig,tau,LaNew,dLa,cNew)
            tau1 = 1/mean(abs(dLa))*1/37;
            tau2 = 1/max(abs(dLa))*2/37;
            tau  = max(tau1,tau2);
            rightStepLength = 0;
            while rightStepLength == 0
                smean = s0 - tau*dLa;
                s     = max(0,min(1,smean));
                [~,c_u,c_sig,LaTrial] = obj.computeParameters(s,l_u,l_sig,rho_u,rho_sig);
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

        function [c,c_u,c_sig,La] = computeParameters(obj,s,l_u,l_sig,rho_u,rho_sig)
            x = obj.calculateSectionID(s);
            [w,v1,v2] = ISCSO_2021(x,0);
            c     = w;
            c_sig = v1;
            c_u   = v2;
            La    = c + l_u*c_u + 0.5*rho_u*c_u^2 + l_sig*c_sig + 0.5*rho_sig*c_sig^2;
        end

        function [dc,dc_u,dc_sig,dLa,dC,dU,dSig] = computeParametersGradient(obj,s,s0,...
                sec,l_u,l_sig,rho_u,rho_sig,c_u,c_sig,method,update,dC,dU,dSig)
            switch method
                case 'approach'
                    dA     = obj.computeSectionGradient(s,sec);
                    A      = obj.computeSection(s,sec);
                    dc     = dA;
                    dc_u   = -dA./A.^2;
                    dc_sig = -dA./A.^2;
                case 'finite differences'
                    if update == 1
                        [dc, dc_u, dc_sig,dC,dU,dSig] = obj.updateFiniteDiffDerivatives(s,...
                            s0,sec);
                    else
                        [dc, dc_u, dc_sig] = obj.computeFiniteDiffDerivatives(s,s0,sec,...
                            dC,dU,dSig);
                    end
                otherwise
                    error('Method not implemented')
            end
            dLa    = dc + (l_u + rho_u*c_u).*dc_u + (l_sig + rho_sig*c_sig).*dc_sig;
        end

        function [dc, dc_u, dc_sig] = computeFiniteDiffDerivatives(obj,s,s0,sec,dC,dU,dSig)
            A      = obj.computeSection(s,sec);
            A0     = obj.computeSection(s0,sec);
            dA     = obj.computeSectionGradient(s,sec);
            dA0    = obj.computeSectionGradient(s0,sec);
            A      = A./A0;
            dA     = dA./dA0;
            dc     = dC.*dA;
            dc_u   = dU.*dA./A.^2;
            dc_sig = dSig.*dA./A.^2;
        end

        function [dc, dc_u, dc_sig,dC,dU,dSig] = updateFiniteDiffDerivatives(obj,s,s0,sec)
            A      = obj.computeSection(s,sec);
            A0     = obj.computeSection(s0,sec);
            dA     = obj.computeSectionGradient(s,sec);
            dA0    = obj.computeSectionGradient(s0,sec);
            A      = A./A0;
            dA     = dA./dA0;
            [dC, dSig, dU] = obj.calculateFiniteGradient(s);
            dc     = dC.*dA;
            dc_u   = dU.*dA./A.^2;
            dc_sig = dSig.*dA./A.^2;
        end

        function [cost, const_u, const_sig, lu, lsig, tauV] = computeMonitoring(obj,...
                c, c_u, c_sig, l_u, l_sig, tau, w0, iter, s,...
                cost, const_u, const_sig, lu, lsig, tauV, monitor,dc,dc_u,dc_sig)
            cost(iter+1)      = c;
            const_u(iter+1)   = c_u;
            const_sig(iter+1) = c_sig;
            lu(iter+1)        = l_u;
            lsig(iter+1)      = l_sig;
            tauV(iter+1)      = tau;
            x                 = obj.calculateSectionID(s);
            auglagr_monitoring(monitor, 1:(iter+1), x, cost, const_sig, ...
                const_u, lu, lsig ,tauV, w0, dc,dc_u,dc_sig);
        end

        function [dC, dV1, dV2] = calculateFiniteGradient(obj,s)
            x0 = obj.calculateSectionID(s);
            [c0, v10, v20] = ISCSO_2021(x0,0);
            for i = 1:numel(s)
                su    = s;
                sl    = s;
                su(i) = s(i)+1/37;
                sl(i) = s(i)-1/37;

                if su(i) > 1
                    su(i) = 1;
                elseif sl(i) < 0
                    sl(i) = 0;
                end
                x1 = obj.calculateSectionID(su);
                x2 = obj.calculateSectionID(sl);
                [cu, v1_u, v2_u] = ISCSO_2021(x1,0);
                [cl, v1_l, v2_l] = ISCSO_2021(x2,0);
                ds = 1/37;
                dCu(i)   = (cu-c0) / (ds);
                dV1_u(i) = (v1_u-v10) / (ds);
                dV2_u(i) = (v2_u-v20) / (ds);
                dCl(i)   = (c0-cl) / (ds);
                dV1_l(i) = (v10-v1_l) / (ds);
                dV2_l(i) = (v20-v2_l) / (ds);
                dC(i)    = (dCu(i) + dCl(i)) / 2;
                dV1(i)   = (dV1_u(i) + dV1_l(i)) / 2;
                dV2(i)   = (dV2_u(i) + dV2_l(i)) / 2;
            end
        end

        function [s,c,c_u,c_sig] = makeFeasible(obj, s, dc_u, dc_sig)
            tau        = 1/max( max(abs(dc_u), abs(dc_sig)) )*3/37;
            isFeasible = false;
            while ~isFeasible
                smean         = s + tau*max(-dc_u,-dc_sig);
                s             = max(0,min(1,smean));
                x             = obj.calculateSectionID(s);
                [c,c_u,c_sig] = ISCSO_2021(x,0);
                isFeasible    = max(c_u, c_sig) <= 0;
                tau           = 2*tau;
            end
        end

    end

    methods (Static, Access = private)

        function x = calculateSectionID(s)
            x = s*(37-1) + 1;
            x = round(x);
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
            Amax     = max(sec);
            Amin     = min(sec)/Amax;
            Amax     = 1;
            p        = 3;
            gradA    = - Amin * p * s.^(p-1) + Amax * p*s.^(p-1);
        end

        function update = checkUpdate(iter)
            multiple      = (iter/10);
            multiple      = multiple/fix(multiple);
            update        = multiple;
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

        function [l_u, l_sig] = updateMultipliers(l_u, l_sig, rho_u, rho_sig, c_u, c_sig)
            l_u0   = l_u;
            l_sig0 = l_sig;
            l_u    = l_u0   + rho_u*max(-l_u0/rho_u,c_u);
            l_sig  = l_sig0 + rho_sig*max(-l_sig0/rho_sig,c_sig);
        end

        function [sections] = computeSectionData(dataSectionsFile)
            run(dataSectionsFile);
            sections = Si(:,1);
        end

        function [iter,dS,s0,cNew,LaNew] = updateParameters(s,s0,c,LaTrial,iter)
            dS     = norm(s-s0)/norm(s0);
            s0     = s;
            cNew   = c;
            LaNew  = LaTrial;
            iter   = iter + 1;
        end

    end
end