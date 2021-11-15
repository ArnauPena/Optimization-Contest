classdef Mesh < handle

    properties (Access = public)
        bars
        barsMatrix
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

    end
    
    methods (Access = private)
        
        function init(obj, cParams)
            obj.dim = cParams.dim;
            obj.data = cParams.data;
            obj.createBars();
            obj.computeConnectivities();
        end
        
        function createBars(obj)
            nel = obj.dim.nel;
            barsMat = zeros(nel, 13);
            for iElem = 1:nel
                s.data = obj.data;
                bar = Bar(s);
                bar.create(iElem);
                elems(iElem, 1) = bar;
                barsMat(iElem,:) = obj.fillBarsMatrix(bar);
            end
            obj.bars = elems;
            obj.barsMatrix = barsMat;
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

    methods (Static,Access = private)

        function a = fillBarsMatrix(bar)
            [x1, y1, z1, x2, y2, z2] = bar.getNodeCoordinates();
            [rho, E, A, Iz] = bar.getMaterialData();
            [buck, yield] = bar.getCriticalStresses();
            a(1) = x1;
            a(2) = x2;
            a(3) = y1;
            a(4) = y2;
            a(5) = z1;
            a(6) = z2;
            a(7) = bar.length;
            a(8) = rho;
            a(9) = E;
            a(10) = A;
            a(11) = Iz;
            a(12) = buck;
            a(13) = yield;
        end
    end
end

