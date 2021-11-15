classdef PostProcessor < handle
    
    properties (Access = public)
    end
    
    properties (Access = private)
        data
        dim
        FEM
        visualizer
    end
    
    methods (Access = public)
        
        function obj = PostProcessor(cParams)
            obj.init(cParams);
        end

    end
    
    methods (Access = private)
        
        function init(obj, cParams)
            obj.data = cParams.data;
            obj.dim  = cParams.dim;
            obj.FEM = cParams.FEM;
            obj.initVisualizer();
            buckle = obj.computeBuckling();
            mass   = obj.calculateMass();
        end
        
        function initVisualizer(obj)
            s.data = obj.data;
            s.dim  = obj.dim;
            obj.visualizer = Visualizer(s);
            obj.visualizer.initialize();
        end

        function mass = calculateMass(obj)
            barsMat = obj.FEM.mesh.barsMatrix;
            rho = barsMat(:,8);
            A   = barsMat(:,10);
            len = barsMat(:,7);
            massMat = rho.*A.*len;
            mass = sum(massMat);
        end
        
        function buckle = computeBuckling(obj)
            stress_cr = obj.FEM.mesh.barsMatrix(:,12);
            stress = obj.FEM.stress;
            nel = size(stress_cr,1);
            buckle = zeros(nel,1);
            for iElem = 1:size(stress_cr,1)
                sig = stress(iElem,1);
                sig_cr = stress_cr(iElem,1);
                diff = sig_cr - abs(sig);
                if (sig < 0 && diff < 0)
                    buckle(iElem,1) = 1;
                end
            end
        end

    end

    
end