classdef StiffnessTest < Test

    properties (Access = private)
        data
        dim
        connectivities
        KComp
        DOFManager
    end
    
    methods(Access = protected)

        function passed = passed(obj)
            obj.initFile();
            [Ke, Kg] = obj.computeResults();
            maxError = obj.checkMaxError(Ke, Kg);
            tolerance = 1e-10;
            if maxError < tolerance
                passed = 1;
            else
                passed = 0;
            end
        end

    end

    methods (Access = private)
        function initFile(obj)
            run(obj.fileName)
            obj.dim            = dim;
            obj.data           = data;
            obj.KComp          = KComp;
            obj.DOFManager     = DOFManager;
            obj.connectivities = connectivities;
        end

        function [Ke, Kg] = computeResults(obj)
            s.dim            = obj.dim;
            s.data           = obj.data;
            s.DOFManager     = obj.DOFManager;
            s.connectivities = obj.connectivities;
            SMC = StiffnessMatrixComputer(s);
            SMC.compute();
            Ke = SMC.KElem;
            Kg = SMC.KGlobal;
        end

        function maxError = checkMaxError(obj, Ke, Kg)
            KEError  = obj.calculateKElemError(Ke);
            KGError  = obj.calculateKGlobalError(Kg);
            maxError = max(KEError, KGError);
        end

        function error = calculateKElemError(obj, Ke)
            value = obj.KComp.KElem;
            diff = abs(value - Ke);
            error = max(diff, [], 'all');            
        end

        function error = calculateKGlobalError(obj, Kg)
            value = obj.KComp.KGlobal;
            diff = abs(value - Kg);
            error = max(diff, [], 'all');            
        end
        
    end
end