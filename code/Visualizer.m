classdef Visualizer < handle
       
    properties (Access = private)
        data
        dim
        FEM
        X, Y, Z
    end
    
    methods (Access = public)
        
        function obj = Visualizer(cParams)
            obj.init(cParams)
        end
        
        function obj = initialize(obj)
            obj.getNodeData();
        end
        
        function obj = plotStress(obj)
            Tn   = obj.data.nodalconnec;
            sig = obj.FEM.stress;
            x    = obj.X;
            y    = obj.Y;
            z    = obj.Z;
            obj.plotNodes();
            patch(x(Tn)',y(Tn)',z(Tn)',[sig';sig'],'edgecolor','flat','linewidth',2);
            cbar = colorbar('Ticks',linspace(min(sig),max(sig),5));
            title(cbar,{'Stress';'(MPa)'});
            hold off
        end
       
        function obj = plotNodes(obj)
            Tn   = obj.data.nodalconnec;
            x    = obj.X;
            y    = obj.Y;
            z    = obj.Z;
            figure
            hold on
            colormap jet;
            plot3(x(Tn)',y(Tn)',z(Tn)','-k','linewidth',0.5);
            view(5,5);
            xlim([0 18*4e3]);
            zlim([-25e3 5e3]);
            ylim([0 4000]);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.data = cParams.data;
            obj.dim  = cParams.dim;
            obj.FEM  = cParams.FEM;
        end
        
        function obj = getNodeData(obj)
            nodes = obj.data.nodes;
            obj.X = nodes(:,1);
            obj.Y = nodes(:,2);
            obj.Z = nodes(:,3);
        end
        
    end
    
end