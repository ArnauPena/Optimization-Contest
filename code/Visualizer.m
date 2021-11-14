classdef Visualizer < handle
       
    properties (Access = private)
        data
        dim
        X
        Y
        Z
    end
    
    methods (Access = public)
        
        function obj = Visualizer(cParams)
            obj.init(cParams)           
        end
        
        function obj = compute(obj)
            obj.preCompute();
            obj.visualize();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.data = cParams.data;
            obj.dim  = cParams.dim;
        end
        
        function obj = preCompute(obj)
            nodes = obj.data.nodes;
            obj.X = nodes(:,1);
            obj.Y = nodes(:,2);
            obj.Z = nodes(:,3);
        end
       
        function obj = visualize(obj)
            Tn   = obj.data.nodalconnec;
            x    = obj.X;
            y    = obj.Y;
            z    = obj.Z;
            figure
            hold on
            axis equal;
            plot3(x(Tn)',y(Tn)',z(Tn)','-k','linewidth',0.5);
        end
        
    end
    
end