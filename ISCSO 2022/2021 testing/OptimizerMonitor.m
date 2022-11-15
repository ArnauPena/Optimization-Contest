classdef OptimizerMonitor
    
    properties
        monitor
    end
    
    methods (Access = public)
        function obj = OptimizerMonitor()
            obj.monitor = figure();
        end

        function update(obj,iter,sections,cost,constraint,stepLength,evaluations,dualVar)
            obj.updateObjectiveGraph(iter,cost);
            obj.updateStressConstraint(iter, constraint);
            obj.updateDispConstraint(iter, constraint);
            obj.updateDualVariable(iter,dualVar);
            obj.updateSectionDistributionGraph(sections);
            obj.updateStepLength(iter,stepLength);
            obj.updateEvaluations(iter,evaluations);
            drawnow
        end
    end

    methods (Static, Access = private)

        function updateObjectiveGraph(iter,objective)
            subplot(2,4,1)
            plot(iter,objective, 'Color', '#0072BD')
            title('Objective function');
        end
        
        function updateSectionDistributionGraph(sections)
            subplot(2,4,4)
            bar(sections);
            ylim([1 37]);
            title('Bar sections');
        end

        function updateStressConstraint(iter, constraint)
            subplot(2,4,2)
            plot(iter,constraint(:,1),'Color', '#0072BD')
            title('Stress violation');
        end

        function updateDispConstraint(iter,constraint)
            subplot(2,4,3)
            plot(iter,constraint(:,2),'Color', '#D95319')
            title('Disp violation');
        end

        function updateDualVariable(iter, dualVar)
            subplot(2,4,5)
            plot(iter,dualVar(:,1),'Color', '#0072BD')
            title('Dual of stress the violation');
            subplot(2,4,6)
            plot(iter,dualVar(:,2),'Color', '#0072BD')
            title('Dual of the displacement violation');
        end

        
        function updateStepLength(iter, t)
            subplot(2,4,7)
            plot(iter,t)
            title('Step length');
        end

        function updateEvaluations(iter,evaluations)
            subplot(2,4,8)
            plot(iter,evaluations)
            title('Function evaluations');
        end

    end

end
