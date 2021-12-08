function auglagr_monitoring(monitor, iter, sections, objective, ...
    stressvio, dispvio, lagrU, lagrS, tau, w0)

    totalvio = stressvio + dispvio;
    createObjectiveGraph(iter,objective, w0);
    createSectionDistributionGraph(sections);
    createVioGraph(iter, stressvio, dispvio);
    createTotalVioGraph(iter, totalvio);
    createLagrangianGraph(iter, lagrU, lagrS);
    createTauGraph(iter, tau)
    drawnow
end

function createObjectiveGraph(iter,objective, w0)
    subplot(2,3,1)
    hold on
    plot(iter,objective, 'Color', '#0072BD')
    w0plot = w0*ones(length(iter));
    plot(iter,w0plot,'--','Color', '#D95319')
    title('Objective function');
    hold off
end

function createSectionDistributionGraph(sections)
    subplot(2,3,2)
    bar(sections);
    ylim([1 37]);
    title('Bar sections');
end

function createLagrangianGraph(iter, lagrU, lagrS)
    subplot(2,3,3)
    plot(iter,lagrS,'Color', '#0072BD')
    hold on
    plot(iter,lagrU,'Color', '#D95319')
    title('Lagrangian multipliers');
    legend('stress', 'disp')
end

function createVioGraph(iter, stressvio, dispvio)
    subplot(2,3,4)
    hold on
    plot(iter,stressvio,'Color', '#0072BD')
    plot(iter,dispvio,'Color', '#D95319')
    title('Stress violation');
    legend('stress', 'disp')
    hold off
end

function createTotalVioGraph(iter, totalvio)
    subplot(2,3,5)
    plot(iter,totalvio)
    title('Total violation');
end

function createTauGraph(iter, tau)
    subplot(2,3,6)
    plot(iter,tau)
    title('Tau');
end

function f = freq(X)
    for i = [1:37]
        len = length(X(X==i));
        f(i) = len;
    end
end