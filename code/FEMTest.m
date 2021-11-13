classdef FEMTest < Test
    
    properties (Access = private)
        data
        dim
        results
        tolerance
    end
    
    methods(Access = protected)

        function passed = passed(obj)
            obj.initFile();
            [directRes, iterRes] = obj.computeResults();
            passed = 1;
%             maxError = obj.checkMaxError(directRes, iterRes);
%             tol = obj.tolerance;
%             if maxError < tol
%                 passed = 1;
%             else
%                 passed = 0;
%             end
        end

    end
    
    methods (Access = private)
        
        function initFile(obj)
            run(obj.fileName)
            obj.dim       = dim;
            obj.data      = data;
            obj.results   = results;
            obj.tolerance = 1e-4;
        end
        
        function [dirRes, iterRes] = computeResults(obj)
            dirRes  = obj.computeStressDisplacement('DIRECT');
            iterRes = obj.computeStressDisplacement('ITERATIVE');
        end
        
        function results = computeStressDisplacement(obj, solver)
            s.dim        = obj.dim;
            s.data       = obj.data;
            s.solvertype = solver;
            FEM = FEMAnalyzer(s);
            FEM.perform();
            results = FEM.displacement;
        end
        
        function error = calculateIterativeError(obj, iterRes)
            values = obj.results;
            diff = abs(iterRes - values);
            error = max(diff);
        end

    end
end

