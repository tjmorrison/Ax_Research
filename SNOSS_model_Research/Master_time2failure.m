%Matlab script to open Dhiraj's figures from DEID processing and recreate a
%nice figure for the DEID report
%Author: Travis Morrison
%Last edits: Feb 2 2022
clear; close all;

%% Strip the figures for their data
%Complexity
open('C:\Users\tjmor\OneDrive\Research\Papers_writing\TARP_DEID_Report\Figures\Feb_3_2021\Complexity.fig');
c = findobj(gca,'Type','line');
complexity = get(c,'Xdata') ;
HS =get(c,'Ydata') ;

%density
open('C:\Users\tjmor\OneDrive\Research\Papers_writing\TARP_DEID_Report\Figures\Feb_3_2021\density.fig');
d = findobj(gca,'Type','line');
density = get(d,'Xdata');

%SDI
open('C:\Users\tjmor\OneDrive\Research\Papers_writing\TARP_DEID_Report\Figures\Feb_3_2021\SDI.fig');
S = findobj(gca,'Type','line');
SDI = get(S,'Xdata');

%tilt board
% open('C:\Users\tjmor\OneDrive\Research\Papers_writing\TARP_DEID_Report\Figures\Feb_3_2021\Feb_05_2021_tilt_board.fig');
% S = findobj(gca,'Type','line');
% percent_water = get(S,'Xdata');

close all;
%% Calcluate averages from data
delta_HS = max(HS)/length(HS);
avg_height = 50; %[mm]
N_averages = floor(max(HS)/avg_height);
N_index_to_avg = floor(avg_height./delta_HS);

start_i = 1;
end_i = N_index_to_avg;
for ii = 1:N_averages + 1
    if ii <= N_averages
        HS_avg(ii) = mean(HS(start_i:end_i));
        density_avg(ii) = mean(density(start_i:end_i));
        complexity_avg(ii) = mean(complexity(start_i:end_i));
        SDI_avg(ii) = mean(SDI(start_i:end_i));
        
        start_i = start_i + N_index_to_avg;
        end_i  = end_i + N_index_to_avg;
    else
        HS_avg(ii) = mean(HS(start_i:end));
        density_avg(ii) = mean(density(start_i:end));
        complexity_avg(ii) = mean(complexity(start_i:end));
        SDI_avg(ii) = mean(SDI(start_i:end));
    end
end

%

%% Make new figure
ft_size =12; 
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex'); 
set(0,'defaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',ft_size);% Plots

% Density figure
figure('Renderer', 'painters', 'Position', [10 10 1100 600])
ax1 = axes();
box(ax1);
barh(HS_avg, density_avg,'FaceColor',[0.5,0.5,0.5])
hold on
scatter(density, HS,40, 'Filled','k')
ylim([0 max(HS)])
xlim([min(density)-20 max(density)])
xticks([64 83 102 113 138 154])
xticklabels({'64','83','102','113','138','154'})
xlabel('$\rho$ (kg m$^{-3}$)')
ylabel('HS (mm)')
% handle second X-axis
ax2 = axes('Position', get(ax1,'Position'), ...
    'XAxisLocation','top', ...
    'Color','none', ...
    'XColor','k');
ax2.YAxis.Visible = 'off';
ax2.XLim = ax1.XLim;

xticks([64 83 102 113 138 154])
xticklabels({'F-','F','F+','4F','1F-','1F'})
ylim([0 max(HS)])
xlim([min(density)-20 max(density)])
grid on
%%
figure('Renderer', 'painters', 'Position', [10 10 400 600])
barh(HS_avg,complexity_avg,'FaceColor',[0.5,0.5,0.5])
hold on
scatter(complexity, HS,40, 'Filled','k')
ylim([0 max(HS)])
xlabel('Complexity')
grid on

%%
figure('Renderer', 'painters', 'Position', [10 10 400 600])
barh(HS_avg,SDI_avg,'FaceColor',[0.5,0.5,0.5])
hold on
scatter(SDI, HS,40, 'Filled','k')
xlabel('SDI')
ylim([0 max(HS)])
grid on

%% Validation figure of SDI and Complexity for shear strength
%Will use all storms for this data~ need data from Dhiraj

%Data harvested from Jameison Thesis in 1995 Fig 6.1
Precip_particles_x = [49.93324964165598, 49.88172304925008,54.89034204288887,69.13744484312497, ...
    74.98922625795149, 60.24910765310424, 80.23439915309018, 95.04595235195472, 105.08544045868035, ...
    100.13771652879399, 110.26503405439334];
Precip_particles_y = [0.17052023121387272,0.10693641618497107, 0.0375722543352599, 0.11849710982658923,...
    0.08959537572254295,0.1502890173410405,0.31213872832369915,0.08959537572254295,...
    0.2283236994219653, 0.3728323699421967,0.36994219653179217];

% Calculate Group 1 from Daniels curve from Jamieson thesis
density_range = [10:1:150];
density_ice = 1000;
Group_1 = 23.*(density_range./density_ice).^2;

Group_2 = 5.35.*(density./density_ice).^1.35;
Group_2_avg = 5.35.*(density_avg./density_ice).^1.35;

Feb_3_strength = SDI.*(density./density_ice).^(complexity);
Feb_3_strength_avg = SDI_avg.*(density_avg./density_ice).^(complexity_avg);
%% Figure: Jamieson vs Dhiraj vs height 
figure()
%barh(HS_avg, density_avg,'FaceColor',[0.5,0.5,0.5])
hold on
scatter(Feb_3_strength,HS, 40,'ko','filled')
scatter(Group_2,HS, 40,'ro','filled')
%plot([0 1],[152.4 152.4],'r-')
legend('SDI($\rho /rho_{ice}$)$^{complexity}$','5.35($\rho /rho_{ice}$)$^{1.35}$')
xlabel('Shear Strength (kPa)')
ylabel('HS (mm)')
grid on
%% Jamieson vs height
figure('Renderer', 'painters', 'Position', [10 10 400 600])
barh(HS_avg, Group_2_avg,'FaceColor',[0.5,0.5,0.5])
hold on
scatter(Group_2,HS, 40,'ko','filled')
%plot([0 1],[152.4 152.4],'r-')
%legend('$\sigma$ = 5.35($\rho /rho_{ice}$)$^{1.35}$')
xlabel('$\sigma$(kPa)')
ylabel('HS (mm)')
ylim([0 max(HS)])
grid on
%%
figure()
plot(Precip_particles_x, Precip_particles_y,'k+')
hold on
grid on
plot(density_range,Group_1,'k--')
scatter(density_avg, Feb_3_strength_avg, 40, HS_avg, 'Filled')
scatter(density, Group_2, 40, 'r', 'Filled')
legend('Jameison 1995, Precip Particles','DEID Method')


%% Time to failure

open('C:\Users\tjmor\OneDrive\Research\Papers_writing\TARP_DEID_Report\Figures\Feb_3_2021\Expectedtimetofailure.fig')
c = findobj(gca,'Type','line');
Time = get(c,'Xdata') ;
TimeToFailure =get(c,'Ydata') ;
%%
figure('Renderer', 'painters', 'Position', [10 10 1400 200])
semilogy(Time, TimeToFailure,'k-','LineWidth',2)
ylabel('Time to Failure (Hrs)')
xlabel('Time (Hrs)')
grid on
%%
open('C:\Users\tjmor\OneDrive\Research\Papers_writing\TARP_DEID_Report\Figures\Feb_3_2021\Shear_overbudden (1).fig')
c = findobj(gca,'Type','line');
tmp =get(c,'Ydata') ;
tmp2 = get(c,'Xdata') ;
Time = tmp2{1,:};
Strength = tmp{1,:};
overburden = tmp{2,:};
%%
figure('Renderer', 'painters', 'Position', [10 10 1400 200])
semilogy(Time, Strength./1000,'k-','LineWidth',2)
hold on
semilogy(Time, overburden./1000,'r-','LineWidth',2)
legend('Strength','Overburden')
ylabel('Stress (kPa)')
xlabel('Time (Hrs)')
grid on