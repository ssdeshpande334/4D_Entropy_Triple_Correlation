%%
clearvars
close all
clc

% SCRIPT FOR FIGURES 6-7

load batch1_variables_rat_cortical_cultures.mat
load batch2_variables_rat_cortical_cultures.mat
 
% two 4D variables, A and B, both of size 12x35x30x14
A = all_AT_batch1;
B = all_AT_batch2;

conc_AT = cat(3, A, B);
conc_AT(isinf(conc_AT)) = NaN;
avg_AT = mean(conc_AT, 3,'omitnan');

A = all_actual_batch1;
B = all_actual_batch2;

conc_actual = cat(3, A, B);
avg_actual = mean(conc_actual, 3,'omitnan');

A = all_entropy_4d_network_batch1;
B = all_entropy_4d_network_batch2;

conc_entropy_4d_network = cat(3, A, B);
rat_entropy_4d_network = mean(conc_entropy_4d_network,3,'omitnan');
%% compute the 4D entropy and plot it

figure;
hold on
set(gca,'fontsize',30)
boxplot(rat_entropy_4d_network,'Colors','g')
ylabel({'4D Network Entropy'})
xlabel('Days{\it in vitro}')
ylim([0 20])
title('Figure 6E')
%% compute the spike rate and plot it
% Assuming 'average_result' is your 4D variable of size 12x35x1x14
% Create 14 separate variables
for i = 1:14
    % Extract each slice along the fourth dimension
    slice_i = avg_actual(:, :, i);
    
    % Assign the slice to a separate variable (e.g., variable1, variable2, ...)
    variable_name = sprintf('actual_variable%d', i);
    assignin('base', variable_name, slice_i);
end

for i = 1:14
    % Extract each slice along the fourth dimension
    slice_i = avg_AT(:, :, i);
    
    % Assign the slice to a separate variable (e.g., variable1, variable2, ...)
    variable_name = sprintf('AT_variable%d', i);
    assignin('base', variable_name, slice_i);
end

%%
% plot motif class 0 actual contributions
figure;
hold on
set(gca,'fontsize',30)
boxplot((actual_variable1*numel(snippet_raster))./64,'Colors','r')
ylabel({'Spike rate (spikes/s)'})
xlabel('Days{\it in vitro}')
title('Figure 6B')

% plot motif class V A/T
figure;
hold on
set(gca,'fontsize',30)
boxplot(AT_variable6,'Colors','b')
yline(0,'k')
ylim([-1 1.5])
ylabel({'Motif Class V'})
xlabel('Days{\it in vitro}')
title('Figure 6D')

%%
ymax = 40;
figure; 
hold on
sgtitle('Figure 7 - top row')
subplot(1,7,1)
hold on;
A = AT_variable1;
A = A(:,~all(isnan(A)));
boxplot(A,'Colors','b')
ylim([-1 ymax])
med = median(A,1,'omitnan');
plot(med,'b','linewidth',4)


subplot(1,7,2)
hold on;
A = AT_variable2;
A = A(:,~all(isnan(A)));
A(isinf(A)) = NaN;
boxplot(A,'Colors','b')
yline(0,'k','linewidth',2)
ylim([-1 ymax])
med = median(A,1,'omitnan');
plot(med,'b','linewidth',4)

subplot(1,7,3)
hold on;
A = AT_variable3;
A(isinf(A)) = NaN;
A = A(:,~all(isnan(A)));
boxplot(A,'Colors','b')
yline(0,'k','linewidth',2)
ylim([-1 ymax])
med = median(A,1,'omitnan');
plot(med,'b','linewidth',4)

subplot(1,7,4)
hold on;
A = AT_variable4;
A = A(:,~all(isnan(A)));
A(isinf(A)) = NaN;
boxplot(A,'Colors','b')
yline(0,'k','linewidth',2)
ylim([-1 ymax])
med = median(A,1,'omitnan');
plot(med,'b','linewidth',4)

subplot(1,7,5)
hold on;
A = AT_variable5;
A = A(:,~all(isnan(A)));
A(isinf(A)) = NaN;
boxplot(A,'Colors','b')
yline(0,'k','linewidth',2)
ylim([-1 ymax])
med = median(A,1,'omitnan');
plot(med,'b','linewidth',4)

subplot(1,7,6)
hold on;
A = AT_variable6;
A = A(:,~all(isnan(A)));
A(isinf(A)) = NaN;
boxplot(A,'Colors','b')
yline(0,'k','linewidth',2)
ylim([-1 ymax])
med = median(A,1,'omitnan');
plot(med,'b','linewidth',4)

subplot(1,7,7)
hold on;
A = AT_variable7;
A = A(:,~all(isnan(A)));
A(:,1) = NaN;
A(:,2) = NaN;
A(isinf(A)) = NaN;
boxplot(A,'Colors','b')
yline(0,'k','linewidth',2)
ylim([-1 ymax])
med = median(A,1,'omitnan');
plot(med,'b','linewidth',4)

%%
ymax = 8;
figure;
hold on
sgtitle('Figure 7 - bottom row')
subplot(1,7,1)
hold on;
A = AT_variable8;
A = A(:,~all(isnan(A)));
A(isinf(A)) = NaN;
boxplot(A,'Colors','b')
yline(0,'k','linewidth',2)
ylim([-1 ymax])
A(:,1) = NaN;
A(:,2) = NaN;
med = median(A,1,'omitnan');
plot(med,'b','linewidth',4)


subplot(1,7,2)
hold on;
A = AT_variable9;
A = A(:,~all(isnan(A)));
A(isinf(A)) = NaN;
boxplot(A,'Colors','b')
yline(0,'k','linewidth',2)
ylim([-1 ymax])
med = median(A,1,'omitnan');
plot(med,'b','linewidth',4)

subplot(1,7,3)
hold on;
A = AT_variable10;
A = A(:,~all(isnan(A)));
A(isinf(A)) = NaN;
boxplot(A,'Colors','b')
yline(0,'k','linewidth',2)
ylim([-1 ymax])
med = median(A,1,'omitnan');
plot(med,'b','linewidth',4)

subplot(1,7,4)
hold on;
A = AT_variable11;
A = A(:,~all(isnan(A)));
A(isinf(A)) = NaN;
boxplot(A,'Colors','b')
yline(0,'k','linewidth',2)
ylim([-1 ymax])
med = median(A,1,'omitnan');
plot(med,'b','linewidth',4)

subplot(1,7,5)
hold on;
A = AT_variable12;
A = A(:,~all(isnan(A)));
A(isinf(A)) = NaN;
boxplot(A,'Colors','b')
yline(0,'k','linewidth',2)
ylim([-1 ymax])
med = median(A,1,'omitnan');
plot(med,'b','linewidth',4)

subplot(1,7,6)
hold on;
A = AT_variable13;
A = A(:,~all(isnan(A)));
A(isinf(A)) = NaN;
boxplot(A,'Colors','b')
yline(0,'k','linewidth',2)
ylim([-1 ymax])
med = median(A,1,'omitnan');
plot(med,'b','linewidth',4)

subplot(1,7,7)
hold on;
A = AT_variable14;
A = A(:,~all(isnan(A)));
A(isinf(A)) = NaN;
boxplot(A,'Colors','b')
yline(0,'k','linewidth',2)
ylim([-1 ymax])
med = median(A,1,'omitnan');
plot(med,'b','linewidth',4)