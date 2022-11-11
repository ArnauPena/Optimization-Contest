% Main code

function [x,fff,rgBest] = computePSO()
    % PSO parameters
    ub      = 1;
    lb      = 0;
    xNumber = 345;
    wMax    = 0.09; % Intertia
    wMin    = 0.04; % Intertia
    k1      = 2.2;   % Acceleration factor
    k2      = 2.2;   % Acceleration factor
    n       = 20;  % population size
    maxIter = 5e3;
    maxRun  = 10;
    monitor = figure()
    for run = 1:maxRun
        for i = 1:n
            for j = 1:xNumber
                x0(i,j) = rand();
            end
        end
        x = x0;
        v = 0.01*x;
        for i = 1:n
            [f0(i,1),c1,c2] = computeCostFunction(x(i,:));
            const(i,1) = c1;
            const(i,2) = c2;
        end
        [fmin0,index0] = min(f0);
        pBest = x0;
        gBest = x0(index0,:);
        cBest(1,1) = const(index0,1);
        cBest(1,2) = const(index0,2);
        iter = 1;
        tol  = 1;
        while iter<=maxIter && tol > 1e-12
            w = wMax - (wMax - wMin)*iter/maxIter;
            for i = 1:n
                for j = 1:xNumber
                    v(i,j) = w*v(i,j) + k1*rand()*(pBest(i,j) - x(i,j))...
                        + k2*rand()*(gBest(1,j) - x(i,j));
                end
            end
            x = x + v;
            x = min(ub,max(lb,x));
            for i = 1:n
                [f(i,1),c1,c2] = computeCostFunction(x(i,:));
                if f(i,1) < f0(i,1)
                    pBest(i,:) = x(i,:);
                    f0(i,:)    = f(i,:);
                    const(i,1) = c1;
                    const(i,2) = c2;
                end
            end
            [fmin,index]      = min(f0);
            ffMin(1,iter)     = fmin;
            ffConst(iter,1:2) = const(index,:);
            ffIter(run)       = iter;
            if fmin < fmin0
                gBest      = pBest(index,:);
                cBest(1,1) = const(index,1);
                cBest(1,2) = const(index,2);
                fmin0      = fmin;
            end

            if iter > 100
                tol = abs(ffMin(iter - 100,run) - fmin0);
            end
%             disp('Iter: ');
%             disp(iter);
%             disp('Value: ');
%             disp(fmin0);
            iter = iter + 1;
            plotMonitoring(cBest,gBest,iter,ffMin,ffConst,fmin)
        end
        fval     = computeCostFunction(gBest);
        fff(run) = fval;
        rgBest(run,:) = gBest;
    end
end