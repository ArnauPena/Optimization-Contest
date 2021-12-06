function auglagr_monitoring(monitor, iter, sections, objective, stressvio, dispvio, totalvio, w0)
    createObjectiveGraph(iter,objective, w0);
    createSectionDistributionGraph(sections);
    createVioGraph(iter, stressvio, dispvio);
    createTotalVioGraph(iter, totalvio);
    drawnow
end

function createObjectiveGraph(iter,objective, w0)
    subplot(2,2,1)
    plot(iter,objective, 'Color', '#0072BD')
    w0plot = w0*ones(length(iter));
    plot(iter,w0plot,'--','Color', '#D95319')
    title('Objective function');
end

function createSectionDistributionGraph(sections)
    subplot(2,2,2)
    bar(sections);
    ylim([1 37]);
    title('Bar sections');
end

function createVioGraph(iter, stressvio, dispvio)
    subplot(2,2,3)
    hold on
    plot(iter,stressvio,'Color', '#0072BD')
    plot(iter,dispvio,'Color', '#D95319')
    title('Stress violation');
    legend('stress', 'disp')
    hold off
end

function createTotalVioGraph(iter, totalvio)
    subplot(2,2,4)
    plot(iter,totalvio)
    title('Total violation');
end

function f = freq(X)
    for i = [1:37]
        len = length(X(X==i));
        f(i) = len;
    end
end