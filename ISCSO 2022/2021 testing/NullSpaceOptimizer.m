classdef NullSpaceOptimizer < handle

    properties(Access = public)
        bestSolution
        functionEvaluations
    end

    properties(Access = private)
        s
        information
        cost
        constraint
    end

    methods(Access = public)

        function obj = NullSpaceOptimizer()
            obj.init();
            obj.solve();
        end

    end

    methods (Access = private)

        function init(obj)
            obj.functionEvaluations   = 0;
            info.numOfDesignVariables = 345;
            info.numOfSections        = 37;
            info.density              = 7.85;
            info.sections             = load("sectionsValue","Si").Si;
            info.barsLength           = obj.computeBarsLength(info);
            info.costGradient         = @(x) obj.computeCostGradient(x);
            info.lengthLessCostGradient = obj.computeLengthLessCostGradient(info);
            obj.information           = info;
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

        function [f,c1,c2] = computeCostAndConstraints(obj,x)
            [f,c1,c2] = ISCSO_2021(x,0);
            obj.functionEvaluations = obj.functionEvaluations + 1;
        end

        function g = computeCostGradient(obj,x)
            info     = obj.information;
            l        = info.barsLength;
            gL       = info.lengthLessCostGradient;
            x        = contToDiscrete(x);
            g        = zeros(1,length(x));
            intenger = fix(x);
            incr     = x - intenger;
            for i = 1:length(x)
                c1   = gL(intenger(i));
                c2   = gL(intenger(i) + 1);
                g(i) = (c1*(1-incr) + c2*incr)*l(i);
            end
        end

    end

    methods (Static, Access = private)

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