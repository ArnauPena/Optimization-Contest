function plotMonitoring(cBest,gBest,iter,ffMin,ffConst,fmin)
    iterV = 1:iter;
    sections = round(gBest*26 + 1);
    sections = max(1,min(37,sections));
    cBest = ones(iter,2).*cBest; 
    subplot(1,2,1)
    hold on
    plot(iterV,gBest,'r',iterV,ffMin,'b')
    xlabel('Iteration')
    ylabel('Cost function')
    title('Cost function')
    hold off
    subplot(1,2,2)
    hold on
    plot(iterV,cBest(:,1),'r',iterV,cBest(:,2),':r',iterV,ffConst(:,1),'b',iterV,ffConst(:,2),':b')
    xlabel('Iteration')
    ylabel('Constraint')
    title('Constraint value')
    hold off
                subplot(2,2,2)
            bar(sections);
            ylim([1 37]);
            title('Bar sections');
    drawnow
end