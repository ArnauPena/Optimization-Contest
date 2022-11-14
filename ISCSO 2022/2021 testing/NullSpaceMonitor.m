classdef NullSpaceMonitor
    
    properties
        monitor
    end
    
    methods (Access = public)
        function obj = NullSpaceMonitor()
            obj.monitor = figure();
        end

        function update(obj,iter,sections,cost,constraint,stepLength,evaluations)
            obj.updateObjectiveGraph(iter,cost);
            obj.updateStressConstraint(iter, constraint);
            obj.updateDispConstraint(iter, constraint);
            obj.updateSectionDistributionGraph(sections);
            obj.updateStepLength(iter,stepLength);
            obj.updateEvaluations(iter,evaluations);
            drawnow
        end
    end

    methods (Static, Access = private)

        function updateObjectiveGraph(iter,objective)
            subplot(2,3,1)
            plot(iter,objective, 'Color', '#0072BD')
            title('Objective function');
        end
        
        function updateSectionDistributionGraph(sections)
            subplot(2,3,4)
            bar(sections);
            ylim([1 37]);
            title('Bar sections');
        end

        function updateStressConstraint(iter, constraint)
            subplot(2,3,2)
            plot(iter,constraint(:,1),'Color', '#0072BD')
            title('Stress violation');
        end

        function updateDispConstraint(iter,constraint)
            subplot(2,3,3)
            plot(iter,constraint(:,2),'Color', '#D95319')
            title('Disp violation');
        end
        
        function updateStepLength(iter, t)
            subplot(2,3,5)
            plot(iter,t)
            title('Step length');
        end

        function updateEvaluations(iter,evaluations)
            subplot(2,3,6)
            plot(iter,evaluations)
            title('Function evaluations');
        end

    end

end
