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
        
        function obj = initialize(obj)
            obj.getNodeData();
            obj.visualize();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.data = cParams.data;
            obj.dim  = cParams.dim;
        end
        
        function obj = getNodeData(obj)
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
            plot3(x(Tn)',y(Tn)',z(Tn)','-k','linewidth',0.5);
            view(5,5);
            xlim([0 18*4e3]);
            zlim([-25e3 5e3]);
            ylim([0 4000]);
            hold off

        end
        
    end
    
end