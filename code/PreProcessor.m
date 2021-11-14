classdef PreProcessor < handle
    
    properties (Access = public)
        data
        dim
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = PreProcessor(fdatafile)
            obj.initFixedData(fdatafile);        
        end

        function obj = initVariableData(obj,vdatafile)
            run(vdatafile);
            obj.data.z_node = data.z_node;
            obj.data.Sel    = data.Sel;
        end

        function obj = updateVariableData(obj,cParams)
            obj.data.z_node = cParams.z_node;
            obj.data.Sel    = cParams.Sel;
        end

        function obj = computeForceCase(obj,opt)
            switch opt
                case 1
                    obj.data.fdata = obj.data.Fext(:,:,1);
                case 2
                    obj.data.fdata = obj.data.Fext(:,:,2);
                case 3
                    obj.data.fdata = obj.data.Fext(:,:,3);
                otherwise
                    error('Not a valid external force case');
            end
        end

        function obj = computeInitialData(obj,vdatafile)
            obj.initVariableData(vdatafile);
            obj.computeFixedCoor();
        end

        function obj = compute(obj)
            obj.computeVariableCoor();
            obj.connectNodes();
        end
        
    end
    
    methods (Access = private)

        function obj = initFixedData(obj,fdatafile)
            run(fdatafile);
            obj.data = data;
            obj.dim  = dim;
        end
        
        function obj = computeFixedCoor(obj)
            nNodes  = obj.dim.nnod;
            nDOF    = obj.dim.ni;
            nodes   = zeros(nNodes,nDOF);
            x       = obj.data.x;
            y12     = obj.data.y12;
            y34     = obj.data.y34;
            zsup    = obj.data.zsup;
            nodes   = obj.assemblyXandYandZsup(x,y12,y34,nodes,zsup);
            obj.data.nodes = nodes;
        end

        function obj = computeVariableCoor(obj)
            nNodes  = obj.dim.nnod;
            nodes   = obj.data.nodes;
            z_node  = obj.data.z_node;
            Nz_node = length(z_node);
            for k = 1:Nz_node
                m = 4*(k-1);
                n = nNodes-3-m;
                nodes(m+1,3) = z_node(k);
                nodes(n,3)   = z_node(k);
                nodes(m+4,3) = z_node(k);
                nodes(n+3,3) = z_node(k);
            end
            obj.data.nodes = nodes;
        end

        function obj = connectNodes(obj)
            x    = obj.data.x;
            k1   = length(x) - 1;
            k2   = k1 + 1;
            k3   = k1/2;
            iEl  = [1:k1];
            iEl2 = [1:k2];
            iEl3 = [1:k3];           
            Tn1(1:k1,1) = 4*(iEl-1) + 1;
            Tn1(1:k1,2) = 4*iEl + 1;
            Tn2(1:k1,1) = 4*iEl;
            Tn2(1:k1,2) = 4*iEl + 4;
            Tn3(1:k1,1) = 4*(iEl-1) + 2;
            Tn3(1:k1,2) = 4*iEl + 2;
            Tn4(1:k1,1) = 4*iEl - 1;
            Tn4(1:k1,2) = 4*iEl + 3;
            Tn5(1:k2,1) = 4*(iEl2-1) + 1;
            Tn5(1:k2,2) = 4*iEl2;
            Tn6(1:k2,1) = 4*(iEl2-1) + 2;
            Tn6(1:k2,2) = 4*iEl2 - 1;
            Tn7(1:k2,1) = 4*(iEl2-1) + 1;
            Tn7(1:k2,2) = 4*(iEl2-1) + 2;            
            Tn8(1:k2,1) = 4*iEl2 - 1;
            Tn8(1:k2,2) = 4*iEl2;
            Tn9(1:k3,1) = 4*(iEl3-1) + 2;
            Tn9(1:k3,2) = 4*iEl3 + 1;
            Tn10(1:k3,1) = 4*iEl3 - 1;
            Tn10(1:k3,2) = 4*iEl3 + 4;
            Tn11(1:k3,1) = 4*(iEl3-1) + 37;
            Tn11(1:k3,2) = 4*iEl3 + 38;
            Tn12(1:k3,1) = 4*(iEl3-1) + 40;
            Tn12(1:k3,2) = 4*(iEl3-1) + 43;
            Tn13(1:k1,1) = 4*(iEl-1) + 1;%
            Tn13(1:k1,2) = 4*iEl + 4;
            Tn14(1:k1,1) = 4*iEl;
            Tn14(1:k1,2) = 4*iEl + 1;
            Tn15(1:k1,1) = 4*(iEl-1) + 2;
            Tn15(1:k1,2) = 4*iEl + 3;
            Tn16(1:k1,1) = 4*iEl - 1;
            Tn16(1:k1,2) = 4*iEl + 2;
            Tn17(1,1:2)  = [4,2];
            Tn17(2,1:2)  = [1,3];
            Tn17(3,1:2)  = [73,75];
            Tn17(4,1:2)  = [76,74];
            Tn = [Tn1;Tn2;Tn3;Tn4;Tn5;Tn6;Tn7;Tn8;Tn9;Tn10;Tn11;Tn12;Tn13;Tn14;Tn15;Tn16;Tn17];
            obj.data.nodalconnec = Tn;
        end

    end

    methods (Static, Access = private)

        function nodes  = assemblyXandYandZsup(x,y12,y34,nodes,zsup)
            for i = 1:length(x)
                m = 4*(i-1);
                nodes(m+1,1) = x(i);
                nodes(m+1,2) = y12;
                nodes(m+2,1) = x(i);
                nodes(m+2,2) = y12;
                nodes(m+2,3) = zsup;
                nodes(m+3,1) = x(i);
                nodes(m+3,2) = y34;
                nodes(m+3,3) = zsup;
                nodes(m+4,1) = x(i);
                nodes(m+4,2) = y34;
            end
        end
        
    end
    
end