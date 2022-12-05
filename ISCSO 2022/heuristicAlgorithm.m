%% Heuristic algorithm
clc; clear; close all;

nBars = 336;
Si = load('AvailableSections.mat','Si').Si;
barLengths = load('BarLengths2022.mat','barLengths').barLengths;

secs = 2*ones(1,nBars);
[weightDefault, vioStressDefault, vioDispDefault] = ISCSO_2022(secs,0);

s = secs;
comptador = 2;

h = figure();

% secs = load('MagicalHeuristics.mat').secs;
% s = secs;
% costVec = zeros(1,10000);
% stressVioVec = zeros(1,10000);
% dispVioVec = zeros(1,10000);
while true
    vioStresses = zeros(1,nBars);
    vioDispls   = zeros(1,nBars);
    
    for i = 1:nBars
        s(i) = secs(i) + 1;
        [weight, vioStress, vioDisp] = ISCSO_2022(s,0);
        vioStresses(i) = vioStress;
        vioDispls(i)   = vioDisp;
        s = secs;
    end
    
    sensitivityStrs =  vioStresses - vioStressDefault;
    sensitivityDisp =  vioDispls - vioDispDefault;

    [idx,col] = quinesBarresCanviem(vioStressDefault, vioDispDefault, sensitivityStrs, sensitivityDisp);

    % Agafem els cinc que pitjor ratio S/L i els incrementem
    secs(idx) = secs(idx) + 1;

    [weightDefault, vioStressDefault, vioDispDefault] = ISCSO_2022(secs,0);

    % Guardem cosetes
    iter = 1:comptador;
    costVec(comptador) = weightDefault;
    stressVioVec(comptador) = vioStressDefault;
    dispVioVec(comptador)   = vioDispDefault;

    plotCostsConstraints(iter, costVec, stressVioVec, dispVioVec, sensitivityStrs, sensitivityDisp, idx, col)
    comptador = comptador + 1;

end

function [idx,c] = quinesBarresCanviem(vioStress, vioDisp, sensitivityStrs, sensitivityDisp)
    barLengths = load('BarLengths2022.mat','barLengths').barLengths;
    perc = 0.25; % percentatge barres a canviar

    if (vioStress>= vioDisp)
        [~, idx] = sort(sensitivityStrs./barLengths*10^(-3), 'ascend');
        barresAcanviar = max(1,perc*length(find(sensitivityStrs<0)));
        idx = idx(1:barresAcanviar);
        idx(sensitivityStrs(idx)>=0) = [];
        c = ['og','or'];
    else
        [~, idx] = sort(sensitivityDisp./barLengths*10^(-3), 'ascend');
        barresAcanviar = max(1,perc*length(find(sensitivityDisp<0)));
        idx = idx(1:barresAcanviar);
        idx(sensitivityDisp(idx)>=0) = [];
        c = ['or','og'];
    end
    

end

function plotCostsConstraints(iter, weight, stressvio, dispvio, sStrs, sDisp, idx, col)

    barLengths = load('BarLengths2022.mat','barLengths').barLengths;

    subplot(2,3,1)
    plot(iter, weight)

    subplot(2,3,2)
    plot(iter, stressvio)

    subplot(2,3,3)
    plot(iter, dispvio)

    subplot(2,3,5)
    plot(barLengths, sStrs, 'ob', barLengths(idx), sStrs(idx), col(1:2))
    xlabel('approx length'), ylabel('stress sensitivity inf')
    grid minor

    subplot(2,3,6)
    plot(barLengths, sDisp, 'ob', barLengths(idx), sDisp(idx), col(3:4))
    xlabel('approx length'), ylabel('stress sensitivity inf')
    grid minor

    drawnow
end


