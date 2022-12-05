%% Analyzing sensitivities
clc; clear; close all,

nBars = 336;
% secs = load('OptimalSections2021.mat','Sopt').Sopt;
secs = 2*ones(1,nBars);
[weightDefault, vioStressDefault, vioDispDefault] = ISCSO_2022(secs,1);

vioStressesInf = zeros(1,nBars);
vioDisplsInf   = zeros(1,nBars);
vioStressesSup = zeros(1,nBars);
vioDisplsSup   = zeros(1,nBars);

%% Find sensitivities
for i = 1:nBars
    secs(i) = secs(i) - 1;
    [~, vioStressInf, vioDispInf] = ISCSO_2022(secs,0);
    vioStressesInf(i) = vioStressInf;
    vioDisplsInf(i)   = vioDispInf;
    secs = resetSecs();
    secs(i) = secs(i) + 1;
    [weight, vioStressSup, vioDispSup] = ISCSO_2022(secs,0);
    vioStressesSup(i) = vioStressSup;
    vioDisplsSup(i)   = vioDispSup;
    secs = resetSecs();
end

sensitivityStrsInf = vioStressDefault - vioStressesInf;
sensitivityDispInf = vioDispDefault   - vioDisplsInf;
sensitivityStrsSup = vioStressDefault - vioStressesSup;
sensitivityDispSup = vioDispDefault   - vioDisplsSup;

% load('SensitivitiesStrsDisp.mat','sensitivityDisp')
% load('SensitivitiesStrsDisp.mat','sensitivityStrs')

%% Find approximate lengths
Si = load('AvailableSections.mat','Si').Si;
% barLengths = zeros(1,nBars);
% for i = 1:nBars
%     newSec = secs(i) - 1;
%     secs(i) = newSec;
%     [weight, ~, ~] = ISCSO_2022(secs,0);
%     barWeight = weight-weightDefault;
%     barArea   = Si(newSec);
%     barRho    = 1;
%     barLengths(i) = -barWeight/(barRho*barArea);
%     secs = resetSecs();
% end

barLengths = load('BarLengths2022.mat','barLengths').barLengths;
%% Figures
% figure()
%     subplot(1,2,1)
%     plot(barLengths, sensitivityStrs, 'o','MarkerEdgeColor', '#0072BD')
%     xlabel('approx length'), ylabel('stress sensitivity')
%     grid minor
%     
%     subplot(1,2,2)
%     plot(barLengths, sensitivityDisp, 'o','MarkerEdgeColor', '#D95319')
%     xlabel('approx length'), ylabel('displacement sensitivity')
%     grid minor

%% Plot normal
figure()
    subplot(1,2,1)
    hold on
    for i = 1:nBars
        plot(barLengths(i), -sensitivityStrsInf(i), 'o','MarkerEdgeColor', '#0072BD','MarkerSize',Si(secs(i))/1000*30)
    end
    hold off
    xlabel('approx length'), ylabel('stress sensitivity inf')
    grid minor

    subplot(1,2,2)
    hold on
    for i = 1:nBars
        plot(barLengths(i), -sensitivityDispInf(i), 'o','MarkerEdgeColor', '#D95319','MarkerSize',Si(secs(i))/1000*30)
    end
    hold off
    xlabel('approx length'), ylabel('displacement sensitivity inf')
    grid minor


figure()
    subplot(1,2,1)
    hold on
    for i = 1:nBars
        plot(barLengths(i), -sensitivityStrsSup(i), 'o','MarkerEdgeColor', '#0072BD','MarkerSize',Si(secs(i))/1000*30)
    end
    hold off
    xlabel('approx length'), ylabel('stress sensitivity sup')
    grid minor

    subplot(1,2,2)
    hold on
    for i = 1:nBars
        plot(barLengths(i), -sensitivityDispSup(i), 'o','MarkerEdgeColor', '#D95319','MarkerSize',Si(secs(i))/1000*30)
    end
    hold off
    xlabel('approx length'), ylabel('displacement sensitivity inf')
    grid minor

%% Plot de l'optim
% 
%     figure()
%         subplot(1,2,1)
%         hold on
%         for i = 1:nBars
%             plot(barLengths(i), -sensitivityStrsInf(i), 'o','MarkerEdgeColor', '#0072BD','MarkerSize',Si(secsOpt(i))/1000*5)
%         end
%         hold off
%         xlabel('approx length'), ylabel('stress sensitivity inf OPT')
%         grid minor
%     
%         subplot(1,2,2)
%         hold on
%         for i = 1:nBars
%             plot(barLengths(i), -sensitivityDispInf(i), 'o','MarkerEdgeColor', '#D95319','MarkerSize',Si(secsOpt(i))/1000*5)
%         end
%         hold off
%         xlabel('approx length'), ylabel('displacement sensitivity inf OPT')
%         grid minor
%     
%     
%     figure()
%         subplot(1,2,1)
%         hold on
%         for i = 1:nBars
%             plot(barLengths(i), -sensitivityStrsSup(i), 'o','MarkerEdgeColor', '#0072BD','MarkerSize',Si(secsOpt(i))/1000*5)
%         end
%         hold off
%         xlabel('approx length'), ylabel('stress sensitivity sup OPT')
%         grid minor
%     
%         subplot(1,2,2)
%         hold on
%         for i = 1:nBars
%             plot(barLengths(i), -sensitivityDispSup(i), 'o','MarkerEdgeColor', '#D95319','MarkerSize',Si(secsOpt(i))/1000*5)
%         end
%         hold off
%         xlabel('approx length'), ylabel('displacement sensitivity inf OPT')
%         grid minor

%% Frequencia acumulada en sensitivitat
% index = 1;
% freqac = length(0:0.01:42);
% for i = 0:0.01:42
%     freqac(index) = length(find(sensitivityStrsSup<=i))/nBars;
%     index = index +1;
% end
% figure()
% plot(0:0.01:42, freqac)

%% Functions
% Reset sections
function secs = resetSecs()
    secs = 2*ones(1,336);
%     secs = load('OptimalSections2021.mat','Sopt').Sopt;
end