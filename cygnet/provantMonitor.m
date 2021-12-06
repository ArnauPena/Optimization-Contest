clc; clear; close all;
x0 = [14	13	14	14	14	18	15	14	14	14	14	14	18	14	14	14	14, ...
    14	14	13	14	14	14	14	14	14	14	14	18	14	14	14	14	18, ...
    18	14	15	14	14	18	14	14	14	14	17	14	14	14	18	14	14, ...
    14	14	14	14	14	14	14	14	14	14	14	14	14	14	14	18	14, ...
    14	18	14	14	14	14	14	14	14	14	18	14	18	18	18	18	14, ...
    13	14	14	14	14	14	14	16	14	14	14	14	17	14	14	18	18, ...
    14	14	14	18	14	14	18	17	13	18	14	14	14	14	14	18	14, ...
    14	18	14	18	14	18	14	18	16	14	14	14	18	14	18	18	14, ...
    14	14	14	14	14	14	14	18	13	14	18	18	18	14	14	14	14, ...
    18	14	14	18	14	18	14	14	18	14	14	18	14	14	14	18	14, ...
    17	14	16	14	14	14	14	14	14	18	14	14	14	14	14	14	14, ...
    17	14	14	14	14	14	18	18	18	14	14	14	14	18	14	18	18, ...
    14	18	14	14	14	14	13	14	14	14	18	14	14	14	14	14	14, ...
    14	14	14	14	14	14	14	18	14	14	14	14	14	14	14	14	14, ...
    14	14	14	14	14	14	14	14	14	14	14	14	18	14	14	18	14, ...
    18	14	14	14	14	13	14	14	18	14	14	14	14	14	18	18	14, ...
    14	14	14	14	14	14	18	13	14	14	14	14	14	14	14	14	13, ...
    18	14	14	14	14	13	18	14	13	14	13	14	14	14	14	14	18, ...
    14	18	18	14	18	14	14	14	14	14	14	18	13	14	13	14	18, ...
    18	14	14	14	14	14	14	14	14	14	14	18	18	14	14	14	14, ...
    14	18	14	14	18];

[w0, v10, v20] = ISCSO_2021(x0,0);

hold on
monitor = createMonitor();
i = 1;

while w0 > 6000
% for iter = 1:600
    hold on;
    monitor = createMonitor();
    deltaSect = 1;
    [objective, totalvio] = modifySection(x0, deltaSect, monitor, w0);
    [~, index] = findMaxDelta(objective, totalvio, w0);
    [w0, x0] = updateIterationParameters(index, objective, x0, deltaSect);
    disp(w0)
    if isempty(w0)
        break
    end
end

% monitor = createMonitor();
% auglagr_monitoring(monitor, iter, x, objective, stressvio, dispvio, totalvio)


function monitor = createMonitor()
    monitor = figure(1);
    hold on;
end

function [objective, totalvio] = modifySection(x0, deltaSect, monitor, w0)
    for i = 1:1:length(x0)
%     for i = 1:1:5
        x = x0;
        x(i) = x0(i)-deltaSect;
        [w, v1, v2] = ISCSO_2021(x,0);
        iter(i)      = i;
        objective(i) = w;
        stressvio(i) = v1;
        dispvio(i)   = v2;
        totalvio(i)  = v1+v2;
        auglagr_monitoring(monitor, iter, x, objective, stressvio, dispvio, totalvio, w0)
%         hold on;
    end
end

function [M, index] = findMaxDelta(objective, totalvio, w0)
    okay = find(totalvio == 0);
    objective = objective(okay);
    delta = w0 - objective;
    [M, indexOkay] = max(delta);
    index = okay(indexOkay);
end

function [w0, x0] = updateIterationParameters(index, objective, x0, deltaSect)
    w0 = objective(index);
    x0(index) = x0(index) - deltaSect;
end

function grad = calculateObjectiveGradient(c, cpre, s, spre)
    grad = (c - cpre) / (s-spre);
end