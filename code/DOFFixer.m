classdef DOFFixer < handle

    properties(SetAccess = private, GetAccess = public)
        fixedDisp
        fixedDOFs
        freeDOFs
    end

    properties(Access = private)
        dim
        fixnod
    end

    methods(Access = public)
        function obj = DOFFixer(cParams)
            obj.init(cParams);
        end

        function obj = fix(obj)
            obj.calculateFixedData();
            obj.calculateFreeDOFs();
        end
    end

    methods(Access = private)
        function init(obj, cParams)
            obj.dim    = cParams.dim;
            obj.fixnod = cParams.data.fixnod;
        end

        function calculateFixedData(obj)
            fnod = obj.fixnod;
            d = obj.dim;
            h = height(fnod);
            vr = zeros(h, 1);
            ur = zeros(h, 1);
            for j = 1:h
               val = fnod(j,3);
               converter = NodeDOFConverter(d, fnod, j);
               DOF = converter.convert();
               vr(j) = DOF;
               ur(j) = val;
            end
            obj.fixedDOFs = vr;
            obj.fixedDisp = ur;
        end

        function calculateFreeDOFs(obj)
            d = obj.dim;
            count = 1;
            for dof = 1:d.ndof
                if ~ismember(dof, obj.fixedDOFs)
                    vl(count) = dof;
                    count = count+1;
                end
            end
            obj.freeDOFs = vl;
        end     
    end
end