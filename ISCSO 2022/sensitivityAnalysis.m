%% Analyzing sensitivities
clc; clear; close all,

nBars = 336;

secs = 1*ones(1,nBars);
[weightDefault, vioStressDefault, vioDispDefault] = ISCSO_2022(secs,1);

vioStresses = zeros(1,nBars);
vioDispls   = zeros(1,nBars);

%% Find sensitivities
% for i = 1:nBars
%     secs(i) = 2*secs(i);
%     [weight, vioStress, vioDisp] = ISCSO_2022(secs,0);
%     vioStresses(i) = vioStress;
%     vioDispls(i)   = vioDisp;
%     secs = resetSecs();
% end
% 
% sensitivityStrs = vioStressDefault - vioStresses;
% sensitivityDisp = vioDispDefault   - vioDispls;

load('SensitivitiesStrsDisp.mat','sensitivityDisp')
load('SensitivitiesStrsDisp.mat','sensitivityStrs')

%% Find approximate lengths
Si = load('AvailableSections.mat','Si').Si;
barLengths = zeros(1,nBars);
for i = 1:nBars
    secs(i) = 2*secs(i);
    [weight, ~, ~] = ISCSO_2022(secs,0);
    barWeight = weight-weightDefault;
    barArea   = Si(2);
    barRho    = 1;
    barLengths(i) = barWeight/(barRho*barArea);
    secs = resetSecs();
end

%% Figures
figure()
    subplot(1,2,1)
    plot(barLengths, sensitivityStrs, 'o','MarkerEdgeColor', '#0072BD')
    xlabel('approx length'), ylabel('stress sensitivity')
    grid minor
    
    subplot(1,2,2)
    plot(barLengths, sensitivityDisp, 'o','MarkerEdgeColor', '#D95319')
    xlabel('approx length'), ylabel('displacement sensitivity')
    grid minor

%% Functions
% Reset sections
function secs = resetSecs()
    secs = 1*ones(1,336);
end