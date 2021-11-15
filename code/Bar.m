classdef Bar < handle

    properties(SetAccess = private, GetAccess = public) % fer privat
        length
    end
    
    properties(Access = private)
        data
        x1, y1, z1
        x2, y2, z2
        E, A, Iz, rho
        RotationMatrix, KBase, KElem
        bucklingStress, yield
    end
    
    methods(Access = public)
        
        function obj = Bar(cParams)
            obj.init(cParams)
        end
        
        function obj = create(obj, e)
            X    = obj.data.nodes;
            nod1 = obj.data.nodalconnec(e,1);
            nod2 = obj.data.nodalconnec(e,2);
            mat  = obj.data.materials;
            obj.x1 = X(nod1, 1);
            obj.x2 = X(nod2, 1);
            obj.y1 = X(nod1, 2);
            obj.y2 = X(nod2, 2);
            obj.z1 = X(nod1, 3);
            obj.z2 = X(nod2, 3);
            obj.E   = mat(e,1);
            obj.A   = mat(e,2);
            obj.rho = mat(e,3);
            obj.Iz  = mat(e,4);
            obj.yield = mat(e,5);
            obj.length         = obj.calculateBarLength();
            obj.RotationMatrix = obj.calculateRotationMatrix();
            obj.KBase          = obj.calculateEuclideanStiffnessMatrix();
            obj.KElem          = obj.rotateStiffnessMatrix();
            obj.bucklingStress = obj.calculateCriticalBucklingStress();
        end
        
        function Re = getRotationMatrix(obj)
            Re = obj.RotationMatrix;
        end
        
        function Ke = getElementStiffnessMatrix(obj)
            Ke = obj.KElem;
        end
        
        function Re = calculateRotationMatrix(obj)
            RM = RotationMatrixComputer(obj);
            RM.compute();
            Re = RM.RotationMatrix;
            obj.RotationMatrix = Re;
        end
        
        function KB = calculateEuclideanStiffnessMatrix(obj)
            ESC = ElementStiffnessComputer(obj);
            ESC.compute();
            KB = ESC.KBase;
        end
        
        function [x1, y1, z1, x2, y2, z2] = getNodeCoordinates(obj)
            x1 = obj.x1;
            y1 = obj.y1;
            z1 = obj.z1;
            x2 = obj.x2;
            y2 = obj.y2;
            z2 = obj.z2;
        end
        
        function [rho, E, A, Iz] = getMaterialData(obj)
            rho = obj.rho;
            E = obj.E;
            A = obj.getSection();
            Iz = obj.Iz;
        end
        
        function A = getSection(obj)
            A = obj.A;
        end
        
        function [buck, yield] = getCriticalStresses(obj)
            buck = obj.bucklingStress;
            yield = obj.yield;
        end

    end

    methods(Access = private)
        
        function init(obj, cParams)
            obj.data = cParams.data;
        end
        
        function le = calculateBarLength(obj)
            X1 = obj.x1;
            Y1 = obj.y1;
            Z1 = obj.z1;
            X2 = obj.x2;
            Y2 = obj.y2;
            Z2 = obj.z2;
            le = sqrt((X2 - X1)^2 + (Y2 - Y1)^2 + (Z2 - Z1)^2);
        end
        
        function Ke = rotateStiffnessMatrix(obj)
            R = obj.RotationMatrix;
            K = obj.KBase;
            Ke = transpose(R) * K * R;
        end

        function sig = calculateCriticalBucklingStress(obj)
            sig = (pi^2 * obj.E * obj.Iz)/(obj.A * obj.length^2);
        end
        
    end
end