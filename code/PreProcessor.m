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

        function obj = compute(obj)
            obj.computeCoor();
%             obj.connectNodes();
        end
        
    end
    
    methods (Access = private)

        function obj = initFixedData(obj,fdatafile)
            run(fdatafile);
            obj.data = data;
            obj.dim  = dim;
        end
        
        function obj = computeCoor(obj)
            nNodes  = obj.dim.nnod;
            nDOF    = obj.dim.ni;
            nodes   = zeros(nNodes,nDOF);
            z_node  = obj.data.z_node;
            Nz_node = length(z_node);
            x       = obj.data.x;
            y12     = obj.data.y12;
            y34     = obj.data.y34;
            zsup    = obj.data.zsup;
            nodes   = obj.assemblyXandYandZsup(x,y12,y34,nodes,zsup);
            nodes   = obj.assemblyZinf(z_node,nNodes,nodes,Nz_node);
            obj.data.nodes = nodes;
        end

        function obj = connectNodes(obj)
                
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

        function nodes  = assemblyZinf(z_node,nNodes,nodes,Nz_node)
            for k = 1:Nz_node
                m = 4*(k-1);
                n = nNodes-3-m;
                nodes(m+1,3) = z_node(k);
                nodes(n,3)   = z_node(k);
                nodes(m+4,3) = z_node(k);
                nodes(n+3,3) = z_node(k);
            end
        end
        
    end
    
end