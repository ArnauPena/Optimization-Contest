classdef FminconAlgorithm < handle
    
    properties (Access = private)
        historial
        searchdir
        preiter
        currentiter
        dc, dc_sig, dc_u
        dCApprox, dUApprox, dSigApprox, sApprox
    end
    
    methods (Access = public)

        function obj = FminconAlgorithm()
            n    = 345;
            obj.historial.fval = [];
            obj.historial.x = [];
            update = 1;
            obj.currentiter = 1;
            s0 = 0.2*ones(1,345);
            sec = obj.computeSectionData();
            [dc,dc_u,dc_sig] = computeParametersGradient(obj,s0,sec,update);
            obj.dc = dc;
            obj.dc_u = dc_u;
            obj.dc_sig = dc_sig;
            myOptions = optimoptions(@fmincon,...
                        'OutputFcn',@obj.outfun, ...
                        'algorithm', 'active-set', ... 
                        'SpecifyConstraintGradient', true, ...
                        'SpecifyObjectiveGradient', true, ...
                        'Display','iter','PlotFcns', @optimplotfval);
            pr.lb = 0*ones(1,n);
            pr.ub = 1*ones(1,n);
            pr.nonlcon = @(x) obj.constraint(x, sec);
            pr.objective = @(x) obj.cost(x, sec);
            pr.options = myOptions;
            pr.solver = 'fmincon';
            pr.x0 = 0.2*ones(1,345);
            x = fmincon(pr);
            pr.objective(x)
            pr.nonlcon(x)
        end

    end

    methods (Access = private)

        function [c, dC] = cost(obj, s,sec)
            x = s*(37-1) + 1;
            x = round(x);
            [w,v1,v2] = ISCSO_2021(x,0);
            c = w;
            update = 0;
            if (mod(obj.currentiter,10) == 0)
                update = 1;
            else
                update = 0;
            end
            [dC, ~, ~] = obj.computeParametersGradient(s,sec,update);
            obj.currentiter = obj.currentiter + 1;
        end
        
        function [c,ceq,dC,dceq] = constraint(obj, s,sec)
            x = s*(37-1) + 1;
            x = round(x);
            [w,v1,v2] = ISCSO_2021(x,0);
            c = [v1,v2];
            ceq = [];
            update = 0;
            if (mod(obj.currentiter,10) == 0)
                update = 1;
            else
                update = 0;
            end
            [~, dc_u, dc_sig] = obj.computeParametersGradient(s,sec,update);
            dC = [dc_u; dc_sig]';
            dceq = [];
        end

        function stop = outfun(obj, x,optimValues,state)
            hist = obj.historial;
            stop = false;
            switch state
                 case 'init'
                     hold on
                 case 'iter'
                   hist.fval = [hist.fval; optimValues.fval];
                   hist.x = [hist.x; x];
                   obj.searchdir = [obj.searchdir;... 
                                optimValues.searchdirection'];
                   hist.iter = optimValues.iteration;
                   obj.historial = hist;
                 case 'done'
                     hold off
                 otherwise
             end
         end

        function [dc,dc_u,dc_sig] = computeParametersGradient(obj,s,sec,update)
            if update == 1
                [dC, dSig, dU] = obj.calculateFiniteGradient(s);
                obj.dCApprox = dC;
                obj.dSigApprox = dSig;
                obj.dUApprox = dU;
                obj.sApprox = s;
            end
            [dc, dc_u, dc_sig] = obj.computeFiniteDiffDerivatives(s,sec);
        end
        
        function [c,c_u,c_sig] = computeParameters(obj,s)
            x = obj.calculateSectionID(s);
            [w,v1,v2] = ISCSO_2021(x,0);
            c     = w;
            c_sig = v1;
            c_u   = v2;
        end
        
        function [sections] = computeSectionData(obj)
            Si = obj.computeData();
            sections = Si(:,1);
        end

        function [dC, dV1, dV2] = calculateFiniteGradient(obj, s)
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

        
        function [dc, dc_u, dc_sig] = computeFiniteDiffDerivatives(obj,s,sec)
            dC     = obj.dCApprox;
            dU     = obj.dUApprox;
            dSig   = obj.dSigApprox;
            s0     = obj.sApprox;
            A      = obj.computeSection(s,sec);
            A0     = obj.computeSection(s0,sec);
            dA     = obj.computeSectionGradient(s,sec);
            dA0    = obj.computeSectionGradient(s0,sec);
            A = A./A0;
            dA = dA./dA0;
            dc     = dC.*dA;
            dc_u   = dU.*dA./A.^2;
            dc_sig = dSig.*dA./A.^2;
        end
                    
    end
   

    methods (Static, Access = private)

        function x = calculateSectionID(s)
            x = s*(37-1) + 1;
            x = round(x);
        end

        function [Si] = computeData()
            run("DataSections.m")
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

    end

end

