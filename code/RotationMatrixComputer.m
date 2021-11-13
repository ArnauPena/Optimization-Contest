classdef RotationMatrixComputer < handle
    properties(SetAccess = private, GetAccess = public)
        RotationMatrix
    end
    
    properties(Access = private)
        bar
    end
    
    methods(Access = public)
        function obj = RotationMatrixComputer(cParams)
            obj.init(cParams);
        end
                
        function obj = compute(obj)
            n = obj.bar;
            le = n.length;
            [x1, y1, z1, x2, y2, z2] = n.getNodeCoordinates();
            coef = 1/le;
            Re = sparse(2,6);
            Re(1,1) =   coef * (x2-x1);
            Re(1,2) =   coef * (y2-y1);
            Re(1,3) =   coef * (z2-z1);
            Re(2,4) =   coef * (x2-x1);
            Re(2,5) =   coef * (y2-y1);
            Re(2,6) =   coef * (z2-z1);
            obj.RotationMatrix = Re;
        end
    end
    
    methods(Access = private)   
        function init(obj, cParams)
            obj.bar = cParams;
        end            
    end
end