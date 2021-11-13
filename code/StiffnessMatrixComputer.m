classdef StiffnessMatrixComputer < handle

    properties(SetAccess = private, GetAccess = public)
        KGlobal
        KElem
    end
    
    properties(Access = private)
        dim
        data
        DOFManager
        connectivities
        mesh
    end

    methods(Access = public)

        function obj = StiffnessMatrixComputer(cParams)
            obj.init(cParams);
        end
        
    end
    
    methods(Access = public)   
        
        function obj = compute(obj)
            obj.createElementStiffness();
            obj.assembleGlobalMatrix();
        end
        
        function KLL = calculateFreeStiffnessMatrix(obj)
            freeDOFs = obj.DOFManager.freeDOFs;
            KLL   = obj.KGlobal(freeDOFs, freeDOFs);
        end
        
        function KRL = calculateFreeFixedStiffnessMatrix(obj)
            freeDOFs = obj.DOFManager.freeDOFs;
            KRL   = obj.KGlobal(freeDOFs, freeDOFs);
        end
        
    end
    
    methods(Access = private)
        
        function init(obj, cParams)
            obj.dim            = cParams.dim;
            obj.data           = cParams.data;
            obj.DOFManager     = cParams.DOFManager;
            obj.connectivities = cParams.connectivities;
            obj.mesh           = cParams.mesh;
        end

        function createDOFFixer(obj)
            s.dim         = obj.dim;
            s.data.fixnod = obj.data.fixnod;
            DOFfixer = DOFFixer(s);
            DOFfixer.fix();
            obj.DOFManager = DOFfixer;
        end

        function createElementStiffness(obj)
            nel = obj.dim.nel;
            nne = obj.dim.nne;
            ni  = obj.dim.ni;
            Kelem = zeros(nne*ni,nne*ni,nel);
            bars = obj.mesh.bars;
            for iElem = 1:nel
                bar   = bars(iElem);
                Ke    = bar.getElementStiffnessMatrix();
                Kelem(:,:,iElem) = Ke;
            end
            obj.KElem = Kelem;
        end

        function assembleGlobalMatrix(obj)
            ndof = obj.dim.ndof;
            nne  = obj.dim.nne;
            ni   = obj.dim.ni;
            nel  = obj.dim.nel;
            Kg = zeros(ndof,ndof);
            for elem = 1:nel
                for i = 1:nne*ni
                    I = obj.connectivities(elem,i);
                    for j = 1:nne*ni
                        J = obj.connectivities(elem,j);
                        Kel = obj.KElem(i, j, elem); 
                        Kg(I, J) = Kg(I, J) + Kel;
                    end
                end
            end
            obj.KGlobal = Kg;
        end

    end

end