classdef NodeDOFConverter < handle

    properties(SetAccess = private, GetAccess = public)
        
    end
    
    properties(Access = private)
        dim
        data
        index
    end
    
    methods(Access = public)
        
        function obj = NodeDOFConverter(dim, data, i)
            obj.dim   = dim;
            obj.data  = data;
            obj.index = i;
        end
        
        function DOF = convert(obj)
            [nod, dir] = obj.processData();
            DOF = obj.calculateDOF(nod, dir);
        end
    end
    
    methods(Access = private)
        
        function [nod, dir] = processData(obj)
           M = obj.data;
           j = obj.index;
           nod = M(j,1);
           dir = M(j,2);
        end
        
        function DOF = calculateDOF(obj, nod, dir)
            ni = obj.dim.ni;
            DOF = ni*(nod-1)+dir;
        end

    end
end

