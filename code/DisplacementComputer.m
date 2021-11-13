classdef DisplacementComputer < handle

    properties(SetAccess = private, GetAccess = public)
        displacement
    end
    
    properties(Access = private)
        F, K
        solvertype
        DOFManager
    end
    
    methods(Access = public)

        function obj = DisplacementComputer(cParams)
            obj.init(cParams);
        end

        function obj = compute(obj)
            obj.calculateFreeDisplacement();
            obj.calculateFixedDisplacement();
        end

    end
    
    methods(Access = private)

        function init(obj, cParams)
            obj.K          = cParams.K;
            obj.F          = cParams.F;
            obj.solvertype = cParams.solvertype;
            obj.DOFManager = cParams.DOFManager;
        end

        function calculateFreeDisplacement(obj)
            freeDOFs = obj.DOFManager.freeDOFs;
            solver = Solver.create(obj.solvertype);
            solution = solver.solve(obj.K, obj.F);
            obj.displacement(freeDOFs,1) = solution;
        end

        function calculateFixedDisplacement(obj)
            fixedDOFs = obj.DOFManager.fixedDOFs;
            fixedDisp = obj.DOFManager.fixedDisp;
            obj.displacement(fixedDOFs,1) = fixedDisp;
        end

    end
end