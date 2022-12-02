classdef AlternateDirectionsMonitor
    
    properties
        monitor
    end
    
    methods (Access = public)
        function obj = AlternateDirectionsMonitor()
            obj.monitor = figure();
        end

        function update(obj, iter, sections, objective, stressvio, dispvio,dc,dc_sig,dc_u,z)
            totalvio = stressvio + dispvio;
            obj.updateObjectiveGraph(iter,objective);
            obj.updateSectionDistributionGraph(sections);
            obj.updateVioGraph(iter, stressvio, dispvio);
            obj.updateTotalVioGraph(iter, totalvio);
            drawnow
        end
    end

    methods (Static, Access = private)

        function updateObjectiveGraph(iter,objective)
            subplot(2,2,1)
            plot(iter,objective, 'Color', '#0072BD')
            title('Objective function');
            hold off
        end
        
        function updateSectionDistributionGraph(sections)
            subplot(2,2,2)
            bar(sections);
            ylim([1 37]);
            title('Bar sections');
        end

        function updateVioGraph(iter, stressvio, dispvio)
            subplot(2,2,3)
            hold on
            plot(iter,stressvio,'Color', '#0072BD')
            plot(iter,dispvio,'Color', '#D95319')
            title('Stress violation');
            legend('stress', 'disp')
            hold off
        end
        
        function updateTotalVioGraph(iter, totalvio)
            subplot(2,2,4)
            plot(iter,totalvio)
            title('Total violation');
        end
    end

end

