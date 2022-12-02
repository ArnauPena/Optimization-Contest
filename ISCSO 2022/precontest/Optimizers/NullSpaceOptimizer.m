classdef NullSpaceOptimizer < handle

    properties(Access = public)
        bestSolution
        functionEvaluations
        nIter
    end

    properties(Access = private)
        designVariable
        dualVariable
        data
        cost
        constraint
        hasConverged
        mOld
        meritNew
        acceptableStep
        lineSearchTrials
        meritGradient
        tol = 1e-8;
        monitor
        fEval
        stepLength
        dualVarVect
    end

    methods(Access = public)

        function obj = NullSpaceOptimizer()
            obj.init();
            obj.monitor = OptimizerMonitor();
            obj.solve();
        end

    end

    methods (Access = private)

        function init(obj)
            obj.data           = PreProcess;
            obj.designVariable = 0.1*ones(1,345);
            obj.dualVariable   = zeros(2,1);
            obj.nIter          = 0;
        end

        function solve(obj)
            obj.hasConverged = false;
            while ~obj.hasConverged
                obj.update();
                obj.increaseIter();
                obj.updateMonitoring();
                obj.checkConvergence();
            end
        end

        function update(obj)
            x0 = obj.designVariable;
            obj.mOld = obj.computeMeritFunction(x0);
            if obj.nIter == 0
                obj.calculateInitialStep();
            end
            obj.acceptableStep   = false;
            obj.lineSearchTrials = 0;
            DJ = obj.cost.gradient;
            Dg = obj.constraint.gradient;
            while ~obj.acceptableStep
                obj.computeMeritGradient(DJ,Dg);
                x = obj.data.updatePrimal(x0,obj.meritGradient');
                obj.computeFunctionAndGradient(x);
                obj.updateDual();
                obj.checkStep(x);
            end
            obj.updateOldValues(x);
        end

        function checkStep(obj,x)
            mNew = obj.computeMeritFunction(x);
            if obj.nIter == 0
                obj.acceptableStep = true;
                obj.meritNew       = mNew;
            end
            if mNew < obj.mOld
                obj.acceptableStep = true;
                obj.meritNew = mNew;
                factor = 1.2;
                obj.data.increaseStepLength(factor);
            elseif obj.data.stepLengthIsTooSmall()
                error('Convergence could not be achieved (step length too small)')
            else
                obj.data.decreaseStepLength();
                obj.lineSearchTrials = obj.lineSearchTrials + 1;
            end
        end

        function updateMonitoring(obj)
            it   = obj.nIter;
            iter = 1:it;
            obj.constraint.vect(it,:) = obj.constraint.value;
            obj.cost.vect(it)         = obj.cost.value;
            obj.fEval(it)             = obj.data.functionEvaluations;
            obj.stepLength(it)        = obj.data.stepLength;
            obj.dualVarVect(it,:)     = obj.dualVariable';
            s    = obj.data.contToDiscrete(obj.designVariable);
            obj.monitor.update(iter,s,obj.cost.vect,obj.constraint.vect,...
                obj.stepLength,obj.fEval,obj.dualVarVect)
        end

        function updateOldValues(obj,x)
            obj.designVariable = x;
        end

        function obj = checkConvergence(obj)
            if abs(obj.meritNew - obj.mOld) < obj.tol && obj.checkConstraint()
                obj.hasConverged = true;
            else

            end

        end

        function isZero = checkConstraint(obj)
            c      = obj.constraint.value;
            isZero = 0 == norm(c,2);
        end

        function calculateInitialStep(obj)
            x       = obj.designVariable;
            obj.computeGradient(x);
            l       = obj.dualVariable;
            DJ      = obj.cost.gradient';
            Dg      = obj.constraint.gradient;
            aJ      = 1;
            DmF     = aJ*(DJ + Dg*l);
            factor  = 1e-3;
            obj.data.computeFirstStepLength(DmF,x,factor);
        end

        function computeMeritGradient(obj,DJ,Dg)
            l       = obj.dualVariable;
            aJ      = 1;
            DmF     = aJ*(DJ' + Dg*l);
            obj.meritGradient = DmF;
        end

        function mF = computeMeritFunction(obj,x)
            obj.computeFunction(x);
            J  = obj.cost.value;
            h  = obj.constraint.value;
            l  = obj.dualVariable;
            aJ = 1;
            AJ = aJ*(J + l'*h);
            mF = AJ;
        end

        function computeFunctionAndGradient(obj,x)
            obj.computeFunction(x);
            obj.computeGradient(x);
        end

        function computeFunction(obj,x)
            x                       = obj.data.contToDiscrete(x);
            x                       = round(x);
            [f,c1,c2]               = obj.data.computeCostAndConstraints(x);
            obj.cost.value          = f;
            obj.constraint.value    = [c1,c2]';
        end

        function computeGradient(obj,x)
            obj.cost.gradient       = obj.data.computeCostGradient(x);
            if ~mod(obj.nIter,5)
                obj.constraint.gradient = obj.data.computeConstraintGradient(x);
            end
        end

        function updateDual(obj)
            DJ = obj.cost.gradient';
            Dh = obj.constraint.gradient;
            h  = obj.constraint.value;
            S  = (Dh'*Dh)^-1;
            aJ = 1;
            aC = 1;
            f2 = 1;
            AC = f2*aC/aJ*S*h;
            AJ = -aC/aJ*S*Dh'*DJ;
            l  = AC + AJ;
            obj.dualVariable = l;
        end

        function increaseIter(obj)
            obj.nIter = obj.nIter + 1;
        end

    end

    methods (Static, Access = private)

        function x = contToDiscrete(x)
            x = x*36 + 1;
            x = max(1,min(x,37));
        end

        function x = discreteToCont(x)
            x = x*1/36 - 1/36;
        end

    end
end