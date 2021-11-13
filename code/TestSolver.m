classdef TestSolver < handle
    
    properties(SetAccess = protected, GetAccess = public)
        solucio
    end
    
    properties(Access = protected)
        expectedResults
    end
    
    methods(Access = public)
        function obj = TestSolver()
        end
        
        function check(obj, results)
        end
    end
    
    methods (Static)

        function test = create(testType, dataFile)
            run(dataFile)
            switch testType
                case {'FEM'}
                    test = FEMTest();
                case {'STIFFNESS'}
                    test = StiffnessTest();
                otherwise
                    error('Invalid Test Type.')
            end
        end

    end
end