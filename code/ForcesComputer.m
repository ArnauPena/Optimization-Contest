classdef ForcesComputer < handle

    properties(SetAccess = private, GetAccess = public)
        Fext
        F
    end
    
    properties(Access = private)
        dim
        data
        fdata
        KGlobal
        DOFManager
    end
    
    methods(Access = public)
        
        function obj = ForcesComputer(cParams)
            obj.init(cParams);
        end

        function obj = compute(obj)
            obj.calculateExternalForcesMatrix();
            obj.calculateForceMatrix();
        end
        
    end
    
    methods(Access = private)

        function init(obj, cParams)
            obj.dim        = cParams.dim;
            obj.fdata      = cParams.data.fdata;
            obj.KGlobal    = cParams.KGlobal;
            obj.DOFManager = cParams.DOFManager;
        end
        
        function calculateExternalForcesMatrix(obj)
            Fdata = obj.fdata;
            d = obj.dim;
            ndof = d.ndof;
            forces = zeros(ndof,1);
            for i = 1:height(Fdata)
               val = Fdata(i,3);
               converter = NodeDOFConverter(d, Fdata, i);
               DOF = converter.convert();
               forces(DOF,1) = val;
            end
            obj.Fext = forces;
        end
        
        function calculateForceMatrix(obj) 
            DOFMgr = obj.DOFManager;
            freeDOFS  = DOFMgr.freeDOFs;
            fixedDOFS = DOFMgr.fixedDOFs;
            fixedDis  = DOFMgr.fixedDisp;
            freeFixedK = obj.KGlobal(freeDOFS, fixedDOFS);
            freeFext   = obj.Fext(freeDOFS,1);
            obj.F = freeFext - freeFixedK*fixedDis;
        end
        
    end
end