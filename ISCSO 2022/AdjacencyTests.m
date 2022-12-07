function AdjacencyTests
nBars = 336;
% First iteration
secs = 2    *ones(1,nBars);
[weightDefault, vioStressDefault, U] = ISCSO_2022(secs,0);

A = zeros(nBars,nBars);
for iBar = 1:nBars
    secsI = secs;
    secsI(iBar) = secsI(iBar) + 1;
    [dWi, dSi, Ui] = ISCSO_2022(secsI,0);
    for jBar = 1:nBars
        secsJ =secsI;
        secsJ(jBar) = secsJ(jBar)-1;
        [dWj, dSj, Uj] = ISCSO_2022(secsJ,0);
        A(iBar,jBar) = (Ui - Uj)/Ui;  
    end
figure(1)
[x1, x2] = meshgrid(sort(1:nBars),sort(1:nBars));
scatter(x1(:), x2(:), [], A(:), 'filled');
drawnow
end

end
