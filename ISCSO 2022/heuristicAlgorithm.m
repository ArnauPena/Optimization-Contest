%% Heuristic algorithm
clc; clear; close all;

nBars = 336;
Si = load('AvailableSections.mat','Si').Si;
barLengths = load('BarLengths2022.mat','barLengths').barLengths;

secs = 2*ones(1,nBars);
[weightDefault, vioStressDefault, vioDispDefault] = ISCSO_2022(secs,0);
[sensitivityStrs, sensitivityDisp] = calculateFullSensitivities(secs, vioStressDefault, vioDispDefault);
[idx,col] = quinesBarresCanviem(vioStressDefault, vioDispDefault, sensitivityStrs, sensitivityDisp);

s = secs;
comptador = 2;

h = figure();

% secs = load('MagicalHeuristics.mat').secs;
% secs = load('EndgameSections.mat').secs;
% s = secs;

%% Fem creixer arees
while true
    [sensitivityStrs, sensitivityDisp] = calculateFullSensitivities(secs, vioStressDefault, vioDispDefault);
%     [sensitivityStrs, sensitivityDisp] = calculateCertainSensitivities(secs, vioStressDefault, vioDispDefault,idx);
    [idx,col] = quinesBarresCanviem(vioStressDefault, vioDispDefault, sensitivityStrs, sensitivityDisp);

    % Agafem els cinc que pitjor ratio S/L i els incrementem
    secs(idx) = secs(idx) + 1;

    [weightDefault, vioStressDefault, vioDispDefault] = ISCSO_2022(secs,0);

    % Guardem cosetes
    iter = 1:comptador;
    costVec(comptador) = weightDefault;
    stressVioVec(comptador) = vioStressDefault;
    dispVioVec(comptador)   = vioDispDefault;

    plotCostsConstraints(iter, secs, costVec, stressVioVec, dispVioVec, sensitivityStrs, sensitivityDisp, idx, col)
    comptador = comptador + 1;
    secs(secs==10) = 11;

    if mod(comptador,5) == 0
        % Cada cinc iteracions recalculem
        idx = 1:336;
    end

end

%% Fem decreixer arees
% 
% while true
%     load('FinalHeuristicaWeight9081.mat')
% 
%     vioStresses = zeros(1,nBars);
%     vioDispls   = zeros(1,nBars);
%     
%     for i = 1:nBars
%         s(i) = secs(i) - 1;
%         [weight, vioStress, vioDisp] = ISCSO_2022(s,0);
%         vioStresses(i) = vioStress;
%         vioDispls(i)   = vioDisp;
%         s = secs;
%     end
% 
%     sensitivityStrs =  vioStresses - vioStressDefault;
%     sensitivityDisp =  vioDispls - vioDispDefault;
% 
%     [idx,col] = quinesBarresCanviem(vioStressDefault, vioDispDefault, sensitivityStrs, sensitivityDisp);
% 
%     % Agafem els cinc que pitjor ratio S/L i els incrementem
%     secs(idx) = secs(idx) + 1;
% 
%     [weightDefault, vioStressDefault, vioDispDefault] = ISCSO_2022(secs,0);
% 
%     % Guardem cosetes
%     iter = 1:comptador;
%     costVec(comptador) = weightDefault;
%     stressVioVec(comptador) = vioStressDefault;
%     dispVioVec(comptador)   = vioDispDefault;
% 
%     plotCostsConstraints(iter, secs, costVec, stressVioVec, dispVioVec, sensitivityStrs, sensitivityDisp, idx, col)
%     comptador = comptador + 1;
%     secs(secs==11) = 10;
% end


%% Funcions

function calculateSensitivities()
end

function [sens_Strs, sens_Disp] = calculateCertainSensitivities(secs, vioStrsPre, vioDispPre, idx)
    nBars = length(secs);
    vioStresses = zeros(1,nBars);
    vioDispls   = zeros(1,nBars);
    s = secs;
    for i = idx
        s(i) = secs(i) + 1;
        [~, vioStress, vioDisp] = ISCSO_2022(s,0);
        vioStresses(i) = vioStress;
        vioDispls(i)   = vioDisp;
        s = secs;
    end
    
    sens_Strs =  vioStresses - vioStrsPre;
    sens_Disp =  vioDispls   - vioDispPre;
end

function [sens_Strs, sens_Disp] = calculateFullSensitivities(secs, vioStrsPre, vioDispPre)
    nBars = length(secs);
    vioStresses = zeros(1,nBars);
    vioDispls   = zeros(1,nBars);
    s = secs;
    for i = 1:nBars
        s(i) = secs(i) + 1;
        [~, vioStress, vioDisp] = ISCSO_2022(s,0);
        vioStresses(i) = vioStress;
        vioDispls(i)   = vioDisp;
        s = secs;
    end
    
    sens_Strs =  vioStresses - vioStrsPre;
    sens_Disp =  vioDispls   - vioDispPre;
end

function [idx,c] = quinesBarresCanviem(vioStress, vioDisp, sensitivityStrs, sensitivityDisp)
    barLengths = load('BarLengths2022.mat','barLengths').barLengths;
    perc = 0.25; % percentatge barres a canviar
    perc2 = 0.35; % percentatge barres a canviar

%     if (vioStress > 100 && vioDisp> 100)
        if (vioStress>= vioDisp)
            [~, idx] = sort(sensitivityStrs./barLengths*10^(-3), 'ascend');
            barresAcanviar = round(max(1,perc*length(find(sensitivityStrs<0))));
            idx = idx(1:barresAcanviar);
            idx(sensitivityStrs(idx)>=0) = [];
            c = ['og','or'];
        else
            [~, idx] = sort(sensitivityDisp./barLengths*10^(-3), 'ascend');
            barresAcanviar = round(max(1,perc*length(find(sensitivityDisp<0))));
            idx = idx(1:barresAcanviar);
            idx(sensitivityDisp(idx)>=0) = [];
            c = ['or','og'];
        end
%     else
%         % We're in the endgame now
%         if (vioStress>= vioDisp)
%             [~, idx] = sort(sensitivityStrs./barLengths*10^(-3), 'ascend');
%             barresAcanviar = round(max(1,perc2*length(find(sensitivityDisp<0))));
%             idx = idx(1:barresAcanviar);
%             idx(sensitivityStrs(idx)>=0) = [];
%             c = ['og','or'];
%         else
%             [~, idx] = sort(sensitivityDisp./barLengths*10^(-3), 'ascend');
%             barresAcanviar = round(max(1,perc2*length(find(sensitivityDisp<0))));
%             idx = idx(1:barresAcanviar);
%             idx(sensitivityDisp(idx)>=0) = [];
%             c = ['or','og'];
%         end
%     end
    

end

function plotCostsConstraints(iter, secs, weight, stressvio, dispvio, sStrs, sDisp, idx, col)

    barLengths = load('BarLengths2022.mat','barLengths').barLengths;

    subplot(2,3,1)
    plot(iter, weight)
    title(num2str(weight(end)))

    subplot(2,3,2)
    plot(iter, stressvio)
    title(num2str(stressvio(end)))

    subplot(2,3,3)
    plot(iter, dispvio)
    title(num2str(dispvio(end)))

    subplot(2,3,4)
    bar(secs)
    title([num2str(length(idx)),' barres canviades'])

    subplot(2,3,5)
    plot(barLengths, sStrs, 'ob', barLengths(idx), sStrs(idx), col(1:2))
    xlabel('approx length'), ylabel('stress sensitivity inf')
    title([num2str(length(find(sStrs<0))),' sens negatives'])
    grid minor

    subplot(2,3,6)
    plot(barLengths, sDisp, 'ob', barLengths(idx), sDisp(idx), col(3:4))
    xlabel('approx length'), ylabel('stress sensitivity inf')
    title([num2str(length(find(sDisp<0))),' sens negatives'])
    grid minor

    drawnow
end


