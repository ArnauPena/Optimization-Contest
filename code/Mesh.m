classdef Mesh < handle

    properties (Access = public)
        bars
        connectivities
        mass
    end

    properties (Access = private)
        dim
        data
    end
    
    methods (Access = public)

        function obj = Mesh(cParams)
            obj.init(cParams);
        end

        function calculateMass(obj)
            nel = obj.dim.nel;
            m = 0;
            for iElem = 1:nel
                bar = obj.bars(iElem);
                len = bar.length;
                sec = bar.getSection();
                rho = 3;
                m = m + rho * len * sec;
            end
            obj.mass = m;
        end

    end
    
    methods (Access = private)
        
        function init(obj, cParams)
            obj.dim = cParams.dim;
            obj.data = cParams.data;
            obj.createBars();
            obj.computeConnectivities();
            obj.calculateMass();
        end
        
        function createBars(obj)
            nel = obj.dim.nel;
            for iElem = 1:nel
                s.data = obj.data;
                bar = Bar(s);
                bar.create(iElem);
                elems(iElem, 1) = bar;
            end
            obj.bars = elems;
        end
        
        function computeConnectivities(obj)
            nel = obj.dim.nel;
            nne = obj.dim.nne;
            ni  = obj.dim.ni;
            T = zeros(nel,nne*ni);
            for e = 1:nel
                for i = 1:nne
                    for j = 1:ni
                        I = ni*(i-1)+j;
                        Tn = obj.data.nodalconnec(e,i);
                        T(e,I) = ni*(Tn-1)+j;
                    end
                end
            end
            obj.connectivities = T;
        end

    end
end

