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
h = figure();
costVect = 8235;
s = load("FinalHeuristicaWeight8235.mat","secs").secs;
it = 2;
while true
    [isDecreasable, sCost, sStress, sDisp] = checkIfDecreasable(s);
    if ~any(isDecreasable)
        break
    end
    s  = decreaseMostFavorableBar(s,isDecreasable,sCost);
    [cost,stress,disp] = ISCSO_2022(s,0);
    iter = 1:it;
    costVect(it) = cost;
    it = it + 1;
    plotDecreaseProcedure(sCost,sStress,sDisp,costVect,iter);
end
disp('Cannot descend more. Local minimum reached.')

%% Interior point method

load('minimumHeuristical7106.mat')
[~, sCost, sStress, sDisp] = checkIfDecreasable(s)
dJ = sCost;
dC = max(sStress,sDisp);
lambda = 100;
dL = dJ + lambda*dC;
[wPre, ~, ~] = ISCSO_2022(s,0);
[vals,idx] = sort(dL);
i = 0;
% seccions = s;
for b = 1:336
    seccions = s;
    seccions(idx(b)) = seccions(idx(b)) - 1;
    while true
        seccions(idx(end-i)) = seccions(idx(end-i)) + 1;
        [w,vS,vD] = ISCSO_2022(seccions, 0);
        if (vS ==0 && vD ==0 && w < wPre && i < 336)
            return
        end
        i = i+1;
    end
end



%% Plot

load('minimumHeuristical7106.mat')
% [isDec, sCost, sStrs, sDisp] = checkIfDecreasable(s);

figure()
subplot(3,1,1)
bar(s)
ylabel('Sections')

subplot(3,1,2)
bar(sCost)
ylabel('Cost')

% subplot(4,1,3)
% bar(sStress)
% ylabel('Stress')
% 
% subplot(4,1,4)
% bar(sDisp)
% ylabel('Displacement')

subplot(3,1,3)
bar(max(sDisp,sStrs))
ylabel('Max of both')
% subplot(4,1,4)
% bar(sCost + l*max(sDisp,sStress))

% subplot(4,1,4)
% bar(lInv)
% ylabel('\lambda')

%% Funcions

function plotDecreaseProcedure(sCost,sStress,sDisp,cost,iter)
    barLengths = load("BarLengths2022.mat","barLengths").barLengths;
    subplot(2,2,1)
    plot(iter, cost)
    title(num2str(cost(end)))

    subplot(2,2,2)
    plot(barLengths, sCost, 'ob');
    xlabel('approx length'), ylabel('cost sensitivity inf')
    grid minor

    subplot(2,2,3)
    plot(barLengths, sStress, 'ob');
    xlabel('approx length'), ylabel('stress sensitivity inf')
    title([num2str(length(find(sStress == 0))),' sens 0'])
    grid minor

    subplot(2,2,4)
    plot(barLengths, sDisp, 'ob')
    xlabel('approx length'), ylabel('disp sensitivity inf')
    title([num2str(length(find(sDisp == 0))),' sens 0'])
    grid minor

    drawnow
end

function [isDecreasable, sCost, sStress, sDisp] = checkIfDecreasable(s)
    sec = s;
    sCost   = zeros(length(s),1);
    sStress = sCost;
    sDisp   = sCost;
    [cost0,c10,c20] = ISCSO_2022(sec,0);
    for i = 1:length(s)
        sec(i) = s(i) - 1;
        sec(i) = max(1,sec(i));
        [cost,c1,c2] = ISCSO_2022(sec,0);
        if c1 == 0 && c2 == 0 && s(i) > 1
            isDecreasable(i) = true;
        else
            isDecreasable(i) = false;
        end
        sCost(i)   = cost - cost0;
        sStress(i) = c1 - c10;
        sDisp(i)   = c2 - c20;
        sec = s;
    end
end

function s = decreaseMostFavorableBar(s,isDecreasable,sCost)
    [~,idx] = sort(sCost,'ascend');
    for i = 1:length(idx)
        if isDecreasable(idx(i))
            s(idx(i)) = s(idx(i)) - 1;
            break
        end
    end
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


