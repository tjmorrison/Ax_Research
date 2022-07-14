

% Travis Morrison
% Date: Feb, 4 2022
%%
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

%% Calc. Dhiraj's values for strength
%Complexity
open('C:\Users\tjmor\OneDrive\Research\Papers_writing\TARP_DEID_Report\Figures\Feb_3_2021\Complexity.fig');
c = findobj(gca,'Type','line');
complexity = get(c,'Xdata') ;
HS =get(c,'Ydata') ;

%density
open('C:\Users\tjmor\OneDrive\Research\Papers_writing\TARP_DEID_Report\Figures\Feb_3_2021\density.fig');
d = findobj(gca,'Type','line');
density = get(d,'Xdata');

open('C:\Users\tjmor\OneDrive\Research\Papers_writing\TARP_DEID_Report\Figures\Feb_3_2021\SDI.fig');
S = findobj(gca,'Type','line');
SDI = get(S,'Xdata');


Feb_3_strength = SDI.*(density./density_ice).^(complexity);
plot(Feb_3_strength,HS)

%%
figure()
plot(Precip_particles_x, Precip_particles_y,'k+')
hold on
grid on
plot(density_range,Group_1,'k--')
scatter(density, Feb_3_strength, 40, HS, 'Filled')

legend('Jameison 1995, Precip Particles','DEID Method')

