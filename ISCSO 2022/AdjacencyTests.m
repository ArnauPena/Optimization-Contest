nBars = 336;
% First iteration
secs = 1    *ones(1,nBars);
[weightDefault, vioStressDefault, U] = ISCSO_2022(secs,0);

A = zeros(nBars,nBars);
for iBar = 1:nBars
    secsI = secs;
    secsI(iBar) = secsI(iBar) + 1;
    [dWi, dSi, Ui] = ISCSO_2022(secsI,0);
    for jBar = 1:nBars
        if (jBar>=iBar)
            secsJ  = secs;
            secsJ(jBar) = secsJ(jBar) + 1;
            [dWj, dSj, Uj] = ISCSO_2022(secsJ,0);
    
            if (iBar == jBar)
                secsIJ = secsI;
            else
                secsIJ = secsI;
                secsIJ(jBar) = secsIJ(jBar) + 1;
            end
    
            [dWij, dSij, Uij] = ISCSO_2022(secsIJ,0);
            incUi  = Ui-U;
            incUj  = Uj-U;
            incUij = Uij-U;
            sumUij = incUi + incUj;
            A(iBar,jBar) = (sumUij - incUij)/sumUij;
        else
            % ja transposarem al final
        end
    end
    
figure(1)
[x1, x2] = meshgrid(sort(1:nBars),sort(1:nBars));
scatter(x1(:), x2(:), [], A(:), 'filled');
drawnow
end
