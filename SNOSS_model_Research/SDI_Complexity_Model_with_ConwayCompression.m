%Matlab script to recompute the overburder and strength model for the DEID
%report. Raw DEID variables are provided by Dhiraj

%Author: Travis Morrison
%Last edits: Feb 17 2022
clear; close all;
%% Load storm of interest

start_time = datetime('12/14/2020 00:00:00'); %data collection start in MST
Data_path = 'E:\Travis\Data\DEID\Dec_14_2020_time_mass_density_swerate_snowheight_complx_SDI_sweheight.txt';
figure_path = 'E:\Travis\DEID_Ax_Research\SNOSS_model_Research\SDI_Complexity_Model\Dec_14_2020\';

% Read data
%time(sec),mass(kg),density(kg/m3),swe rate (mm/hr),snow height
%(mm),complexity,SDI,swe height (mm)
fileID = fopen(Data_path,'r');
formatSpec = '%f %f %f %f %f %f %f %f';
sizeA = [8 Inf];
data = fscanf(fileID,formatSpec,sizeA)';
fclose(fileID);

%pad with zeros for time 0 (for my sanity) - ie. time 0, no snow, no
%density.... this may be unneccary 
%data = [zeros(1,8);data]; 

% Get data real time (MST)
data_time = start_time + seconds(data(:,1));
%% Compute Stength and Overburden without compression
%Conway Stength and Overburden model

% Where do these coeffiecents from from?
%compute strength using exponential fitting
lambda = data(:,7)*1000; %[Pa] leading coeff in model 
beta = data(:,6); % Exponent in stength
rho_ice = 917; %[kg m^-3] density of ice

%Compute strength without compression model from layer 2, to top, so first
%layer with snow, since layer 1 was at time 0

%This is the 1D data of the snow strength
strength_wo_comp = lambda.*(data(1:end,3)./rho_ice).^beta; %Strength of each layer when it lands (no compression)
%preallocate to make in time and height, rows increase in time, col
%increase in height
strength_wo_compression = nan(length(data),length(data));

%Compute strength for all layers except the bottime (0) layer 
for zz = 1:length(data)
    strength_wo_compression(zz:end,zz) = strength_wo_comp(zz);
end

%compute overburden for a particular slope angle
slope_angle = 38; %degrees
delta_t = diff(data(:,1)); %delta_t in seconds
delta_t = [0;delta_t]; %pad the delta T to align with data size, 0 at time 0
swe = (data(:,4).*(delta_t./3600)); %divide swe rate by delta t in hours to get swe in mm
% Then divide by 1000 to get swe in m 
layer_weight = 9.8.*swe.*cosd(slope_angle).*sind(slope_angle); %this is the weight of each layer as a function of slope angle
layer_weight_zz = 9.8.*swe.*cosd(slope_angle).^2; %this is the weight neglecting the shear force used in the overburden model
%% Conway Model compression model
%Now program the model for every time step...
%coefficients for the kinetic term
B1 = 6.5e-7; %[Pa s] Dependent on density - functional form? khoshima paper
B2 = 19.3;
E = 67.3; %[kJ^-1 mol^-1]
R = 0.0083; %[kJ mol^-1 K^-1]
T = 273.15; %[K] will populate with data later on
sigma_m = 75; %[pa] stress for metamorphic processes

%Get time in seconds
seconds = data(:,1);
%% Conway density model forward Euler Method method

%Preallocate the density
density= nan(length(data),length(data));
%For all heights loop over the model
for zz = 1:length(data(:,3))
    %At each height, compute the density for all time, starting at hieght
    %zz. This gets smaller as we march in height 
    for tt = zz:length(data(:,3))
        %for the first time step allocate the density as the measured density 
        if tt == zz
            density(tt,zz) =data(zz,3);
        %for the remaining timesteps compute the densification due to the
        %overburden and metamorphic stresses
        else
            %compute the new overburden
            sigma_zz = sum(layer_weight_zz(zz:tt));
            %Compte eta at time t using the old density
            eta =  B1.*exp(B2.*(density(tt-1,zz)./rho_ice)).*exp(E./(R.*T));
            %compte the increase in density over time
            delta_density = (density(tt-1,zz)/eta).*(sigma_m + sigma_zz)*delta_t(tt); 
            %Add the increase in density to the density to the previous time
            %step
            density(tt,zz) = density(tt-1,zz) + delta_density;

        end
    end
end


%% Compute strength & overburden with compression
%Compute new strength with compression accounted for
for zz = 1:length(layer_weight)
    strength(:,zz) = lambda(zz).*(density(:,zz)./rho_ice).^beta(zz);
end

%Compute the overburden 
overburden = nan(length(data),length(data)); %note overburder will be the same for both cases wince we are not changing the height
for zz = 1:length(layer_weight)
    overburden(zz:end,zz) = cumsum(layer_weight(zz:end));
end

%Compute SI 
SI_wo_compression = strength_wo_compression./overburden;
SI = strength./overburden;

%% Time to failure model 
% need to add Netwon Raphtson root finding technique
% SInce diff works on the row (time) this works in a one liner
T2Fail_wo_compression= (SI_wo_compression(2:end,:) - 1)./abs((diff(SI_wo_compression)./diff(seconds)));
for zz =1:length(seconds)
    T2Fail(:,zz) = (SI(2:end,zz) - 1)./abs((diff(SI(:,zz))./diff(seconds)));
end
%% Calc change in height
%since we know how much each row increased in density, we can estimate the
%decrease in height over the storm
change_density = diff(density);
total_h_compression = (nansum(change_density(end,:)))^(1/3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make figures
ft_size =12; 
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex'); 
set(0,'defaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',ft_size);% Plots

%Set variables
height = [1:1:length(data)].*5; %[mm]

%% Density 

figure()
plot(data(2:end,3),height(1:end-1),'k-','LineWidth',3)
hold on
plot(density(end,:)',height(1:end),'k--','LineWidth',3)
legend('Observed','Modeled')
grid on
ylabel('HS (mm)')
xlabel('$\rho$ (kg m$^{-3}$)')
title('Observed $\rho$')
saveas(gcf,strcat(figure_path,'observed_density'),'png')

figure()
pcolor(data_time,height, (density)')
shading interp
colormap jet
colorbar
%caxis([0 10])
ylabel('HS (mm)')
xlabel('Time (MST)')
title('Modeled $\rho$')
hcb = colorbar;
hcb.Title;
hcb.Title.String = "\rho (kg m^{-3})";
saveas(gcf,strcat(figure_path,'modeled_density'),'png')

%Clabel....need

% Strength & Overburden
figure()
pcolor(data_time,height(2:end), strength(:,2:end)')
shading interp
colormap jet
colorbar
%caxis([30 60])
ylabel('HS (mm)')
xlabel('Time (MST)')
title('Strength')
hcb = colorbar;
hcb.Title;
hcb.Title.String = "\sigma_f (Pa)";
saveas(gcf,strcat(figure_path,'strength'),'png')

% Need to double check overburden
figure()
pcolor(data_time,height(2:end), overburden(:,2:end)')
shading interp
colormap jet
colorbar
%caxis([30 60])
ylabel('HS (mm)')
xlabel('Time (MST)')
title('Overburden')
hcb = colorbar;
hcb.Title;
hcb.Title.String = "\sigma_f (Pa)";
saveas(gcf,strcat(figure_path,'overburden'),'png')

% SI index
figure()
pcolor(data_time,height(2:end), SI(:,2:end)')
shading interp
colormap jet
colorbar
caxis([0 10])
ylabel('HS (mm)')
xlabel('Time (MST)')
title('SI')
saveas(gcf,strcat(figure_path,'SI'),'png')

figure()
semilogy(data_time, SI(:,1)','k-','LineWidth',2)
grid on
axis tight
ylabel('SI')
xlabel('Time (MST)')
title('SI Basal Layer')
saveas(gcf,strcat(figure_path,'SI_BasalLayer'),'png')


min_SI = nanmin(SI,[],1);
[min_SI,z_Index] = nanmin(min_SI,[],2)
min_SI_z = height(z_Index)

figure()
pcolor(data_time,height(2:end), SI_wo_compression(:,2:end)')
shading interp
colormap jet
colorbar
caxis([0 10])
ylabel('HS (mm)')
xlabel('Time (MST)')
title('SI w/o Compression')
saveas(gcf,strcat(figure_path,'SI_wo_compression'),'png')


min_SI_wo = nanmin(SI_wo_compression,[],1);
[min_SI_wo,z_Index_wo] = nanmin(min_SI_wo,[],2)
min_SIwo_z = height(z_Index_wo)

%% Time to failure
figure()
pcolor(data_time(2:end),height(2:end), abs(T2Fail(:,2:end))')
shading interp
colormap jet
colorbar

%set(c,'ColorScale','log')
caxis([1e3 10e3])
ylabel('HS (mm)')
xlabel('Time (MST)')
title('$t_f$')
hcb = colorbar;
hcb.Title;
hcb.Title.String = "(s)";
saveas(gcf,strcat(figure_path,'time2fail'),'png')

% Time to failure along basal layer
figure()
semilogy(data_time(2:end), (T2Fail(:,z_Index)),'k-','Linewidth', 2)
ylabel('$t_f$ (s)')
xlabel('Tims (MST)')
title(strcat('$t_f$($h = $', num2str(height(z_Index)), ' mm)'))
grid on
saveas(gcf,strcat(figure_path,'time2fail_weaklayer'),'png')

figure()
semilogy(data_time(2:end), (T2Fail(:,1)),'k-','Linewidth', 2)
ylabel('$t_f$ (s)')
xlabel('Tims (MST)')
title('$t_f$')
grid on
saveas(gcf,strcat(figure_path,'time2fail_Basallayer'),'png')

close all;
