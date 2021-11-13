classdef ElementStiffnessComputer < handle

    properties(SetAccess = private, GetAccess = public)
        KBase
    end
    
    properties(Access = private)
        bar
    end
    
    methods(Access = public)
        function obj = ElementStiffnessComputer(cParams)
            obj.init(cParams);
        end
                
        function obj = compute(obj)
            n  = obj.bar;
            le = n.length;
            [E, A, ~] = n.getMaterialData();
            K = sparse(2,2);
            coef = A * E / le;
            K(1,1) = coef;
            K(1,2) = -coef;
            K(2,1) = -coef;
            K(2,2) = coef;
            obj.KBase = K;
        end
    end
    
    methods(Access = private)   
        function init(obj, cParams)
            obj.bar = cParams;
        end
    end
end