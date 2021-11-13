classdef Bar < handle

    properties(SetAccess = private, GetAccess = public) % fer privat
        length
    end
    
    properties(Access = private)
        data
        x1, y1, z1
        x2, y2, z2
        E, A, Iz
        RotationMatrix, KBase, KElem
    end
    
    methods(Access = public)
        
        function obj = Bar(cParams)
            obj.init(cParams)
        end
        
        function obj = create(obj, e)
            X    = obj.data.nodes;
            nod1 = obj.data.nodalconnec(e,1);
            nod2 = obj.data.nodalconnec(e,2);
            Mat  = obj.data.materials;
            Tm   = obj.data.matconnec(e);
            obj.x1 = X(nod1, 1);
            obj.x2 = X(nod2, 1);
            obj.y1 = X(nod1, 2);
            obj.y2 = X(nod2, 2);
            obj.z1 = X(nod1, 3);
            obj.z2 = X(nod2, 3);
            obj.E  = Mat(Tm,1);
            obj.A  = Mat(Tm,2);
            obj.Iz = Mat(Tm,3);
            obj.length  = obj.calculateBarLength();
            obj.RotationMatrix = obj.calculateRotationMatrix();
            obj.KBase = obj.calculateEuclideanStiffnessMatrix();
            obj.KElem = obj.rotateStiffnessMatrix();
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
        
        function [E, A, Iz] = getMaterialData(obj)
            E = obj.E;
            A = obj.getSection();
            Iz = obj.Iz;
        end
        
        function A = getSection(obj)
            A = obj.A;
        end

    end

    methods(Access = protected)
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
        
    end
end