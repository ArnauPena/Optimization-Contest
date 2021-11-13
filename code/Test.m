classdef (Abstract) Test < handle
    
    properties(SetAccess = protected, GetAccess = public)
        
    end
    
    properties(Access = protected)
        fileName
    end
    
    methods(Access = public)
        function checkPassed(obj, fileName)
            obj.fileName = fileName;
            if obj.passed()
                obj.displayStatus(1, fileName)
            end
        end
    end
    
    methods(Access = protected)
         p = passed(obj)
    end
    
    methods (Static, Access = private)
        function displayStatus(passed, type)
            switch passed
                case {1}
                    fprintf('['); fprintf(type); fprintf('] ');
                    fprintf('Test '); cprintf('-comment', 'passed'); fprintf('!\n');
                case {0}
                    fprintf('['); fprintf(type); fprintf('] ');
                    fprintf('Test ') ; cprintf('-err', 'failed') ; fprintf('!\n') ;
            end
        end
    end
end