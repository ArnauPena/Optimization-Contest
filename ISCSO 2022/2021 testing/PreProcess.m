classdef PreProcess < handle

    properties (Access = public)
        information
        functionEvaluations
        stepLength
    end

    properties (Access = private)
        
    end

    methods (Access = public)

        function obj = PreProcess()
            obj.init()
        end

        function [f,c1,c2] = computeCostAndConstraints(obj,x)
            [f,c1,c2] = ISCSO_2021(x,0);
            obj.functionEvaluations = obj.functionEvaluations + 1;
        end

        function c = computeConstraint(obj,x)
            [~,c1,c2] = ISCSO_2021(x,0);
            c         = [c1,c2];
            obj.functionEvaluations = obj.functionEvaluations + 1;
        end

        function g = computeCostGradient(obj,x)
            info     = obj.information;
            l        = info.barsLength;
            gL       = info.lengthLessCostGradient;
            x        = obj.contToDiscrete(x);
            g        = zeros(1,length(x));
            x        = max(1,min(37,x));
            intenger = fix(x);
            incr     = x - intenger;
            for i = 1:length(x)
                if intenger(i) == 37
                    g(i) = gL(intenger(i))*l(i);
                else
                    k1   = gL(intenger(i));
                    k2   = gL(intenger(i) + 1);
                    g(i) = (k1*(1-incr(i)) + k2*incr(i))*l(i);
                end
            end
            g = g/norm(g,2);
        end

        function DC = computeConstraintGradient(obj,x)
            x  = round(obj.contToDiscrete(x));
            x  = max(1,min(37,x));
            DC = zeros(length(x),2);
            c0 = obj.computeConstraint(x);
            for i = 1:length(x)
                x1 = x;
                x2 = x;
                if x(i) == 1
                    x1(i)   = x(i) + 1;
                    c1      = obj.computeConstraint(x1);
                    DC(i,:) = 0.5*(c1 - c0);
                elseif x(i) == 37
                    x1(i)   = x(i) - 1;
                    c1      = obj.computeConstraint(x1);
                    DC(i,:) = 0.5*(c0 - c1);
                else
                    x1(i)   = x(i) - 1;
                    x2(i)   = x(i) + 1;
                    c1      = obj.computeConstraint(x1);
                    c2      = obj.computeConstraint(x2);
                    DC(i,:) = 0.5*(c2 - c1);
                end
            end
            DC(:,1) = DC(:,1)/norm(DC(:,1),2);
            DC(:,2) = DC(:,2)/norm(DC(:,2),2);
        end

        function x = updatePrimal(obj,x,g)
            t = obj.stepLength;
            x = x - t*g;
            x = min(1,max(0,x));
        end

        function is = stepLengthIsTooSmall(obj)
            is = obj.stepLength < 1e-5;
        end

        function decreaseStepLength(obj)
            obj.stepLength = obj.stepLength/2;
        end

        function increaseStepLength(obj,factor)
            obj.stepLength = obj.stepLength*factor;
        end

        function computeFirstStepLength(obj,g,x,factor)
            obj.stepLength = norm(g,2)/norm(x,2)*factor;
        end

    end

    methods (Access = private)

        function init(obj)
            info.numOfDesignVariables   = 345;
            info.numOfSections          = 37;
            info.density                = 7.85;
            obj.functionEvaluations     = 0;
            info.sections               = load("sectionsValue","Si").Si;
            info.barsLength             = obj.computeBarsLength(info);
            info.costGradient           = @(x) obj.computeCostGradient(x);
            info.lengthLessCostGradient = obj.computeLengthLessCostGradient(info);
            obj.information             = info;
        end
            
        function l = computeBarsLength(obj,info)
            n         = info.numOfDesignVariables;
            sections0 = ones(1,n);
            l         = zeros(1,n);
            [f0,~,~]  = obj.computeCostAndConstraints(sections0);
            for i = 1:n
                sections1      = sections0;
                sections1(1,i) = 2;
                [f1,~,~]       = obj.computeCostAndConstraints(sections1);
                l(i)           = obj.computeLength(f0,f1,info);
            end
        end

    end

    methods (Static, Access = public)

        function l = computeLength(f0,f1,info)
            Si  = info.sections;
            rho = info.density;
            l   = (f1 - f0)/((Si(2) - Si(1))*rho);
        end

        function x = contToDiscrete(x)
            x = x*36 + 1;
            x = max(1,min(x,37));
        end

        function x = discreteToCont(x)
            x = x*1/36 - 1/36;
        end

        function g = computeLengthLessCostGradient(info)
            Si     = info.sections;
            rho    = info.density;
            g      = zeros(1,length(Si));
            g(1)   = rho*(Si(2) - Si(1));
            g(end) = rho*(Si(end) - Si(end-1));
            for i = 2:length(Si)-1
                g(i) = 0.5*rho*(Si(i+1) - Si(i-1));
            end
        end

    end

end