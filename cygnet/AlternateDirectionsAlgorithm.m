classdef AlternateDirectionsAlgorithm < handle
    
    properties (Access = private)
        dCApprox, dSigApprox, dUApprox
        sApprox
    end
    
    methods (Access = public)
        


        function obj = AlternateDirectionsAlgorithm()
            
         obj.computeProjectedGradient()
        end

    end

    methods (Access = private)
        
        function computeProjectedGradient(obj)
            monitor = AlternateDirectionsMonitor();
            s0        = 0.1*ones(1,345);
            method = 'finite differences';
            update = 1;
            sec    = obj.computeSectionData();             
           
             dc = load('dc.mat', 'dc').dc;
             dc_u = load('dc_uAltDir.mat', 'dc_u').dc_u;
             dc_sig = load('dc_sigAltDir.mat', 'dc_sig').dc_sig;
             
             obj.dCApprox = dc;
             obj.dSigApprox = dc_sig;
             obj.dUApprox = dc_u;
             obj.sApprox = s0;
             
             s = obj.makeFeasible(s0, dc_u,dc_sig);            
             [c,c_u,c_sig] = obj.computeParameters(s);
             iter = 1;
             while iter < 10000
                 
               if (iter == 0 || mod(iter,100) == 0)
                    update = 1;
                    disp(iter)
                    disp('fem update')
                else
                    update = 0;
                end
                 
                 
                 [dc,dc_u,dc_sig] = obj.computeParametersGradient(s,sec,method,update);   

                 h =  1/max(dc)*3/37;
                 tau = h;
                 hasCostDescend = false;
                 i = 1;
                 while ~hasCostDescend
                     smean = s - tau*abs(dc);                
                     sT = max(0,min(1,smean));
                     [cNew,c_u,c_sig] = obj.computeParameters(sT);
                     if (c_u > 0 || c_sig >0)
                      [~,dc_u,dc_sig] = obj.computeParametersGradient(sT,sec,method,1);   
                      sT = obj.makeFeasible(sT, dc_u,dc_sig); 
                      [cNew,c_u,c_sig] = obj.computeParameters(sT);                     
                     end
                     if cNew < 1.1*c
                        hasCostDescend = true;
                        if (c_u ~= 0 && c_sig ~=0)
                          sT = obj.makeFeasible(sT, dc_u,dc_sig); 
                        end
                     end
                     cP(i) = cNew;
                     tauP(i) = tau;
                     i = i+1;
                     tau = tau/2;
                     
                 end
                 
                 s = sT;
                 c = cNew;
                 objective(iter) = c;
                 stressvio(iter) = c_sig;
                 dispvio(iter) = c_u;
                 z = zeros(size(s));
                 iter = iter+1;
                                 sections = obj.calculateSectionID(s);

                 monitor.update(1:iter-1, sections, objective, stressvio, dispvio,dc,dc_sig,dc_u,z)

                 
                 
             end
             
             
        end
        
        function computeAlternate(obj)
             monitor = AlternateDirectionsMonitor();
             s0        = 0.1*ones(1,345);
         %   s0        = rand(1,345);

            method = 'finite differences';
            update = 1;
            sec    = obj.computeSectionData();
%            [dc,dc_u,dc_sig] = obj.computeParametersGradient(s0,sec,method,update);

             dc = load('dc.mat', 'dc').dc;
             dc_u = load('dc_uAltDir.mat', 'dc_u').dc_u;
             dc_sig = load('dc_sigAltDir.mat', 'dc_sig').dc_sig;
             
             obj.dCApprox = dc;
             obj.dSigApprox = dc_sig;
             obj.dUApprox = dc_u;
             obj.sApprox = s0;
             
             s = obj.makeFeasible(s0, dc_u,dc_sig);
             %s = s0;
             
            [c,c_u,c_sig] = obj.computeParameters(s);
            iter = 1;
            RECORD_C = 1000000;
            RECORD_S = [];
            DCGLOBAL = [];
            
            
            
            z = zeros(size(s));
            tau = 1;

            while iter < 1000
                
                if (iter == 1 || mod(iter,10) == 0 || mod(iter-1,10)==0)
                    update = 1;
                    disp(iter)
                    disp('fem update')
                else
                    update = 0;
                end
                
                [dc,~,~] = obj.computeParametersGradient(s,sec,method,0);   
                s = obj.makeCostDecrease(s, c,dc+z,z);
                [c,c_u,c_sig] = obj.computeParameters(s);
                objective(iter) = c;
                stressvio(iter) = c_sig;
                dispvio(iter) = c_u;
                sections = obj.calculateSectionID(s);
                iterations = 1:1:length(objective);
                monitor.update(iterations, sections, objective, stressvio, dispvio,dc,dc_sig,dc_u,z)
                DCGLOBAL(:,iter) = dc;                
                iter = iter + 1;
                if (c<RECORD_C && c_u == 0 && c_sig == 0)
                    RECORD_C = c;
                    RECORD_S = s;
                end 
                
                
                
                if (c_u ~= 0 && c_sig ~=0)
                    [~,dc_u,dc_sig] = obj.computeParametersGradient(s,sec,method,update);              
                    sM = s + 1/tau*z;
                    sF = obj.makeFeasible(sM, dc_u,dc_sig);
                    [c,c_u,c_sig] = obj.computeParameters(sF);
                    objective(iter) = c;
                    stressvio(iter) = c_sig;
                    dispvio(iter) = c_u;
                    update = 1;
                    iter = iter + 1;
                else
                    sF = s;
                end

                if (c<RECORD_C && c_u == 0 && c_sig == 0)
                    RECORD_C = c;
                    RECORD_S = s;
                end
                
                z = z + tau*(s - sF);
                s = sF;
               
            end
            DCGLOBAL;  
            
            
            
        end        

        function [dc,dc_u,dc_sig] = computeParametersGradient(obj,s,sec,method,update)
            switch method
                case 'approach'
                    dA     = obj.computeSectionGradient(s,sec);
                    A      = obj.computeSection(s,sec);
                    dc     = dA;
                    dc_u   = -dA./A.^2;
                    dc_sig = -dA./A.^2;
                case 'finite differences'

                    if update == 1
                        [dC, dSig, dU] = obj.calculateFiniteGradient(s);
                        obj.dCApprox = dC;
                        obj.dSigApprox = dSig;
                        obj.dUApprox = dU;
                        obj.sApprox = s;
                    end
                    [dc, dc_u, dc_sig] = obj.computeFiniteDiffDerivatives(s,sec);

                otherwise
                    error('Method not implemented')
            end

        end

        
        function [c,c_u,c_sig] = computeParameters(obj,s)
            x = obj.calculateSectionID(s);
            [w,v1,v2] = ISCSO_2021(x,0);
            c     = w;
            c_sig = v1;
            c_u   = v2;
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
        
        function [sections] = computeSectionData(obj)
            Si = obj.computeData();
            sections = Si(:,1);
        end

        function [s] = makeFeasible(obj, s, dc_u, dc_sig)
            
                [a,b] = max(abs([-dc_u;-dc_sig]));
                is1 = b == 1;
                g(is1) = -dc_u(is1);
                g(~is1) = -dc_sig(~is1);            
            
            h =  1/max(g)*1/37;
            tau = h;

            isFeasible = false;
            i = 1;
            while ~isFeasible

                smean = s + tau*abs(g);                
                s = max(0,min(1,smean));
                [c,c_u,c_sig] = obj.computeParameters(s);
                isFeasible = max(c_u, c_sig) <= 0;
                cp(i) = max(c_u, c_sig);
                t(i) = tau;
                  i = i+1;
                tau = 1.5*tau;
            end
        end

        function [s] = makeCostDecrease(obj, s, c, dc,z)
            tau = 1/max(dc)*3/37;
            c = c+z*s';
            hasCostDecreased = false;
            while ~hasCostDecreased
                smean = s - tau*dc;
                s = max(0,min(1,smean));
                [cNew,~,~] = obj.computeParameters(s);
                hasCostDecreased = cNew + z*s' <= c;
                tau = tau/2;
            end
        end

%         function [dC, dV1, dV2] = calculateFiniteGradient(obj, s)
%             x0 = obj.calculateSectionID(s);
%             [c0, v10, v20] = ISCSO_2021(x0,0);
%             for i = 1:numel(s)
%                 s(i) = s(i)+1/37;
%                 if s(i) > 1
%                     s(i) = 1;
%                 end
%                 x = obj.calculateSectionID(s);
%                 [c, v1, v2] = ISCSO_2021(x,0);
%                 ds = 1/37;
%                 dC(i)  = (c-c0) / (ds);
%                 dV1(i) = (v1-v10) / (ds);
%                 dV2(i) = (v2-v20) / (ds);
%                 s(i) = s(i) - ds;
%             end
%             
%         end

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

