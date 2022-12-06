%% Heuristic algorithm
clc; clear; close all;

% Load known data
Si = load('AvailableSections.mat','Si').Si;
barLengths = load('BarLengths2022.mat','barLengths').barLengths;
nBars = 336;

% Parameters
lambda = 1000;
stepA = 3;


% First iteration
secs = 2*ones(1,nBars);
[weightDefault, vioStressDefault, vioDispDefault] = ISCSO_2022(secs,0);

% secs = min(37,max(1,secs));

s = secs;
comptador = 2;

h = figure();

% secs = load('MagicalHeuristics.mat').secs;
% secs = load('EndgameSections.mat').secs;
% s = secs;

%% Fem creixer arees
while true
    [sensitivityStrs, sensitivityDisp] = calculateSensitivities(secs, vioStressDefault, vioDispDefault);
    [newVals, lambda, col] = quinesBarresCanviem(vioStressDefault, vioDispDefault, sensitivityStrs, sensitivityDisp, secs, lambda, stepA);

    % Agafem els cinc que pitjor ratio S/L i els incrementem
    secs = secs - newVals;
    secs = min(37,max(1,secs));
    [weightDefault, vioStressDefault, vioDispDefault] = ISCSO_2022(secs,0);

%     rho = dJ + lambda*dc ;
%     lambda = lambda + rho*max(vioStressDefault,vioDispDefault);

    % Guardem cosetes
    iter = 1:comptador;
    costVec(comptador) = weightDefault;
    stressVioVec(comptador) = vioStressDefault;
    dispVioVec(comptador)   = vioDispDefault;
    lambdas(comptador)      = lambda;

    plotCostsConstraints(iter, secs, costVec, stressVioVec, dispVioVec, sensitivityStrs, sensitivityDisp, lambdas, col)
    comptador = comptador + 1;
    secs(secs==10) = 11;

end

%% Funcions


function dJ = costGradient(secs)
    barLengths = load('BarLengths2022.mat','barLengths').barLengths;
    Si = load('AvailableSections.mat').Si;
    dA = Si(secs+1) - Si(secs);
    dAdx = dA/(secs+1 - secs);
    dJ = 10^3*barLengths.*dAdx;
end

function [sens_Strs, sens_Disp] = calculateSensitivities(secs, vioStrsPre, vioDispPre)
    Si = load('AvailableSections.mat').Si;
    nBars = length(secs);
    s = secs;
    for i = 1:nBars
        s(i) = secs(i) + 1;
        [~, vioStress, vioDisp] = ISCSO_2022(s,0);
        dA = Si(s(i)) - Si(secs(i));
        dAdx = dA/(secs(i)+1 - secs(i));
        sens_Strs(i) = (vioStress - vioStrsPre)/(dA) * dAdx;
        sens_Disp(i) = (vioDisp - vioDispPre)/(dA) * dAdx;
        s = secs;
    end
end

function [newVals, lambda, col] = quinesBarresCanviem(vioStress, vioDisp, sensitivityStrs, sensitivityDisp, secs, lambda, stepA)
    rho = 1;
    dJ = costGradient(secs);
    if (vioStress>= vioDisp)
        dC = sensitivityStrs;
        c = vioStress;
        col = ['og','or'];
    else
        dC = sensitivityDisp;
        c = vioDisp;
        col = ['or','og'];
    end

    isdCNeg = dC<0;
    lambdaTrial = min(-dJ(isdCNeg)./dC(isdCNeg));

    lambda = max(lambdaTrial + 0.001, lambda+rho*c);
    dS = lambda*dC + dJ;
    newVals = -round(stepA*dS/min(dS));
    

end

function plotCostsConstraints(iter, secs, weight, stressvio, dispvio, sStrs, sDisp, lambdas, col)

    barLengths = load('BarLengths2022.mat','barLengths').barLengths;

    subplot(2,4,1)
    plot(iter, weight)
    title(num2str(weight(end)))

    subplot(2,4,2)
    plot(iter, stressvio)
    title(num2str(stressvio(end)))

    subplot(2,4,3)
    plot(iter, dispvio)
    title(num2str(dispvio(end)))

    subplot(2,4,4)
    plot(iter, lambdas)
    title(num2str(lambdas(end)))

    subplot(2,4,5)
    bar(secs)
    title(['336 barres canviades'])

    subplot(2,4,6)
    plot(barLengths, sStrs, col(1:2))
    xlabel('approx length'), ylabel('stress sensitivity inf')
    title([num2str(length(find(sStrs<0))),' sens negatives'])
    grid minor

    subplot(2,4,7)
    plot(barLengths, sDisp, col(3:4))
    xlabel('approx length'), ylabel('stress sensitivity inf')
    title([num2str(length(find(sDisp<0))),' sens negatives'])
    grid minor

    drawnow
end


