classdef FEMAnalyzer < handle
    
    properties(SetAccess = private, GetAccess = public)
        displacement
        stress
    end
    
    properties(Access = private)
        data
        DOFManager
        solvertype
        connectivities
        K, KElem, KGlobal
        F, Fext
        dim
        mesh
    end
    
    methods(Access = public)
        
        function obj = FEMAnalyzer(cParams)
            obj.init(cParams);
        end

        function obj = perform(obj)
            obj.createMesh();
            obj.computeStiffnessMatrix();
            obj.computeForces();
            obj.computeDisplacements();
            obj.computeStress();
        end

    end
    
    methods(Access = public)
        function obj = compute(obj)
            obj.perform();
        end
    end
    
    methods(Access = private)
        
        function init(obj, cParams)
            obj.data       = cParams.data;
            obj.dim        = cParams.dim;
            obj.solvertype = cParams.solvertype;
            obj.createDOFManager();
        end

        function createDOFManager(obj)
            s.dim         = obj.dim;
            s.data.fixnod = obj.data.fixnod;
            DOFfixer = DOFFixer(s);
            DOFfixer.fix();
            obj.DOFManager = DOFfixer;
        end

        function createMesh(obj)
            s.dim = obj.dim;
            s.data = obj.data;
            obj.mesh = Mesh(s);
            obj.connectivities = obj.mesh.connectivities;
        end
        
        function computeStiffnessMatrix(obj)
            s.dim            = obj.dim;
            s.data           = obj.data;
            s.DOFManager     = obj.DOFManager;
            s.connectivities = obj.connectivities;
            s.mesh           = obj.mesh;
            KComp = StiffnessMatrixComputer(s);
            KComp.compute();
            obj.KElem   = KComp.KElem;
            obj.KGlobal = KComp.KGlobal;
            obj.K = KComp.calculateFreeStiffnessMatrix();
        end
        
        function computeForces(obj)
            s.dim        = obj.dim;
            s.data       = obj.data;
            s.KGlobal    = obj.KGlobal;
            s.DOFManager = obj.DOFManager;
            FC = ForcesComputer(s);
            FC.compute();
            obj.Fext = FC.Fext;
            obj.F = FC.F;
        end
        
        function computeDisplacements(obj)
            s.K          = obj.K;
            s.F          = obj.F;
            s.solvertype = obj.solvertype;
            s.DOFManager = obj.DOFManager;
            DC = DisplacementComputer(s);
            DC.compute();
            obj.displacement = DC.displacement;
        end
        
        function computeStress(obj)
            s.dim            = obj.dim;
            s.data           = obj.data;
            s.mesh           = obj.mesh;
            s.displacement   = obj.displacement;
            s.connectivities = obj.connectivities;
            s.KElem = obj.KElem;
            SC = StressComputer(s);
            SC.compute();
            obj.stress = SC.stress;
        end
        
    end
end