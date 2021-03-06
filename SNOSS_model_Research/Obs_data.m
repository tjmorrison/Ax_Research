%Code to quickly look at UAC datalogger output
%written by Travis Morrison
%3/15/2021
%%
clear;close all;clc;
%%
obs_path = 'E:\Travis\DEID_Ax_Research\SNOSS_model_Research\12_1_2020_4_30_2021_avalanche_obs.csv';
fid = fopen(strcat(obs_path), 'rt');
b = textscan(fid, '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s', ...
      'Delimiter',',', 'CollectOutput',1, 'HeaderLines',1);
fclose(fid);

%%
bigdata = b{1};
t = datetime(bigdata(:,1),'InputFormat','MM/dd/yyyy');
index = find(~isnat(t)); %line up with weather data
t = t(index);
t1 = find(flipud(t) > datetime([2020, 12, 1, 00, 00, 00]), 1, 'first');
t2 = find(flipud(t) < datetime([2021, 4, 31, 00, 00, 00]), 1, 'last');
t = t(t1:t2);


Region = bigdata(t1:t2,2);
bigdata = bigdata(find(contains(Region,'Salt Lake')),:);
%fill empty entry with unknown
bigdata(find(strcmp(bigdata(:,4),'')),4) = {'Unknown'};
Avy_trigger = bigdata(:,4);
Natural_data = bigdata(find(contains(Avy_trigger,'Natural')),:);
Natural_data = timetable(datetime(Natural_data(:,1)),Natural_data(:,6),'VariableNames',{'Type'});
Skier_data = bigdata(find(contains(Avy_trigger,'Skier')),:);
Skier_data = timetable(datetime(Skier_data(:,1)),Skier_data(:,6),'VariableNames',{'Type'});
SnowBo_data = bigdata(find(contains(Avy_trigger,'Snowboarder')),:);
SnowBo_data = timetable(datetime(SnowBo_data(:,1)),SnowBo_data(:,6),'VariableNames',{'Type'});
Unknown_data = bigdata(find(contains(Avy_trigger,'Unknown')),:);
Unknown_data = timetable(datetime(Unknown_data(:,1)),Unknown_data(:,6),'VariableNames',{'Type'});
Snowmobiler_data = bigdata(find(contains(Avy_trigger,'Snowmobiler')),:);
Snowmobiler_data = timetable(datetime(Snowmobiler_data(:,1)),Snowmobiler_data(:,6),'VariableNames',{'Type'});
Explosive_data = bigdata(find(contains(Avy_trigger,'Explosive')),:);
Explosive_data = timetable(datetime(Explosive_data(:,1)),Explosive_data(:,6),'VariableNames',{'Type'});


%% Calculation for Table for DEID report
clc;
startTime = datetime('14-Dec-2020');
endTime = datetime('17-Dec-2020');

Natural_num = length(Natural_data.Time(startTime:endTime))
Unknown_num = length(Unknown_data.Time(startTime:endTime))
Explosive_num = length(Explosive_data.Time(startTime:endTime))
User_num = length(Skier_data.Time(startTime:endTime))+length(SnowBo_data.Time(startTime:endTime))+...
    length(Snowmobiler_data.Time(startTime:endTime))
total = Natural_num + Unknown_num + User_num + Explosive_num
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%
Skier_triggered_index = strcmp(Avy_trigger,'Skier');
Snowboarder_triggered_index = strcmp(Avy_trigger,'Snowboarder');
Snowmobiler_triggered_index = strcmp(Avy_trigger,'Snowmobiler');
User_triggered_index = Skier_triggered_index+Snowboarder_triggered_index+Snowmobiler_triggered_index;
Natural_triggered_index = strcmp(Avy_trigger,'Natural');
Unknown_triggered_index = strcmp(Avy_trigger,'Unknown');
Explosive_triggered_index = strcmp(Avy_trigger,'Explosive');

%%


User_Avy_Class.bed_surface.Facets = strcmp(bigdata(find(User_triggered_index),6),'Facets');
User_Avy_Class.bed_surface.DH = strcmp(bigdata(find(User_triggered_index),6),'Depth Hoar');
User_Avy_Class.bed_surface.NS_OS = strcmp(bigdata(find(User_triggered_index),6),'New Snow/Old Snow Interface');
User_Avy_Class.bed_surface.NS = strcmp(bigdata(find(User_triggered_index),6),'New Snow');

User_Trig = timetable(t(find(User_triggered_index)),User_Avy_Class.bed_surface.Facets,...
    User_Avy_Class.bed_surface.NS_OS,User_Avy_Class.bed_surface.DH,User_Avy_Class.bed_surface.NS,'VariableNames',{'Facets','NS_OS','DH','NS'});
User_Trig_daily_sum=retime(User_Trig,'daily',@sum);

Nat_Avy_Class.bed_surface.Facets = strcmp(bigdata(find(Natural_triggered_index),6),'Facets');
Nat_Avy_Class.bed_surface.DH = strcmp(bigdata(find(Natural_triggered_index),6),'Depth Hoar');
Nat_Avy_Class.bed_surface.NS_OS = strcmp(bigdata(find(Natural_triggered_index),6),'New Snow/Old Snow Interface');
Nat_Avy_Class.bed_surface.NS = strcmp(bigdata(find(Natural_triggered_index),6),'New Snow');

Nat_Trig = timetable(t(find(Natural_triggered_index)),Nat_Avy_Class.bed_surface.Facets,...
    Nat_Avy_Class.bed_surface.NS_OS,Nat_Avy_Class.bed_surface.DH,Nat_Avy_Class.bed_surface.NS,'VariableNames',{'Facets','NS_OS','DH','NS'});

Nat_Trig = timetable(t(find(Natural_triggered_index)),bigdata(find(Natural_triggered_index),:));

Nat_Trig_daily_sum=retime(Nat_Trig,'daily',@sum);

Exp_Avy_Class.bed_surface.Facets = strcmp(bigdata(find(Explosive_triggered_index),6),'Facets');
Exp_Avy_Class.bed_surface.DH = strcmp(bigdata(find(Explosive_triggered_index),6),'Depth Hoar');
Exp_Avy_Class.bed_surface.NS_OS = strcmp(bigdata(find(Explosive_triggered_index),6),'New Snow/Old Snow Interface');
Exp_Avy_Class.bed_surface.NS = strcmp(bigdata(find(Explosive_triggered_index),6),'New Snow');

Exp_Trig = timetable(t(find(Explosive_triggered_index)),Exp_Avy_Class.bed_surface.Facets,...
    Exp_Avy_Class.bed_surface.NS_OS,Exp_Avy_Class.bed_surface.DH,Exp_Avy_Class.bed_surface.NS, 'VariableNames',{'Facets','NS_OS','DH','NS'});
Exp_Trig_daily_sum=retime(Exp_Trig,'daily',@sum);


%Obs_TT = timetable(t,User_triggered_index,Natural_triggered_index);
%Obs_daily_sum=retime(Obs_TT,'daily',@sum);

%writetimetable(Obs_daily_sum,'DailyObs.csv')

%% plot Daily Avg of people, User trig avys, Fat
figure_path = '';
ft_size =12; 
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex'); 
set(0,'defaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',ft_size);% Plots


%% Figures

S = timerange(startTime,endTime);
%for daily new snow, powder day
day_prior_start = '12/12/2020 01:00:00'; 
day_prior_end = '2/15/2021 20:00:00';
Sprev = timerange(day_prior_start,day_prior_end);

hfig=figure;
a = 25; b= 15; 
hfig.Units='centimeters';hfig.Position=[10 10 a b];
set(hfig,'defaulttextinterpreter','latex','DefaultAxesFontSize',ft_size)
haxes=tight_subplot(1,1,[0.05 0.05],[2 1]/a,[1 1]/b);

axes(haxes(1))

for ii = 1:size(daily_sum.date,1)-1
   if isweekend(daily_sum.date(ii))
       b2 = area([daily_sum.date(ii)  daily_sum.date(ii+1)],[400 400],'FaceColor','k');
       b2.FaceAlpha=0.1;
       b2.EdgeAlpha = 0;
       hold on
       b2 = area([daily_sum.date(ii)  daily_sum.date(ii+1)],[-100 -100],'FaceColor','k');
       b2.FaceAlpha=0.1;
       b2.EdgeAlpha = 0;
       
   end
end

b1 = bar(daily_sum.date(S),daily_sum.uac_ppl(S),'b');
b1.FaceAlpha = 0.5;
ylabel('People [N]')
hold on
% b2 = bar(daily_sum.date(S),daily_sum.uac_beac(S),'r');
% b2.FaceAlpha=0.5;
b3 = bar(daily_sum.date(S),daily_sum.wba_ppl(S),'b');
b3.FaceAlpha=0.5;

p1 = plot([datetime(2021,2,6) datetime(2021,2,6)],[0 350],'r-','LineWidth',2);
plot([datetime(2021,1,30) datetime(2021,1,30)],[0 350],'r-','LineWidth',2)
plot([datetime(2021,1,8) datetime(2021,1,8)],[0 350],'r-','LineWidth',2)
ylim([0 340])

yyaxis right
p2 = plot(Obs_daily_sum.t(S),Obs_daily_sum.User_triggered_index(S),'k-o','LineWidth',1.5,'MarkerSize',5);
ylabel('User Triggered Avalaches [N]')
legend([b1 p1 p2],'People','Fatality','User Avy')
%xlabel('[MM/DD/YYYY]')
%title(strcat(plot_title,' Daily Avg'))
grid on
set(gca,'XMinorTick','on','YMinorTick','on')
axis tight
ax = gca;
ax.YColor = 'black';

set(hfig,'PaperUnits','centimeters');
set(hfig,'PaperSize',[a b]);
set(hfig,'PaperPosition', [0 0 a b]);
print(hfig,'-dpng',[figure_path,'Ppl_SkierTrig_Fat_WP_wba_uac_raw','.png']);
print(hfig,'-dpdf',[figure_path,'Ppl_SkierTrig_Fat_WPwba_uac_raw','.pdf']);
savefig(strcat(figure_path,'Ppl_SkierTrig_Fat_WP_wba_uac_raw'));
%% Weather figure
hfig=figure;
a = 25; b= 15; 
hfig.Units='centimeters';hfig.Position=[10 10 a b];
set(hfig,'defaulttextinterpreter','latex','DefaultAxesFontSize',ft_size)
haxes=tight_subplot(1,1,[0.05 0.05],[2 1]/a,[1 1]/b);

%weather figure
axes(haxes(1))
b3 = bar(daily_weather.Time(S)-hours(7),daily_weather.ISWR(S)./max(daily_weather.ISWR(S)),'b','LineWidth',1);
b3.FaceAlpha = 0.3;

ylabel('Blue Bird Index [ ]')
%xlabel('[MM/DD/YYYY]')
%title(strcat(plot_title,' Daily Avg'))
grid on
set(gca,'XMinorTick','on','YMinorTick','on')
axis tight
hold on
ylim([0 1])

yyaxis right
%calc new snow in cm
new_snow = (daily_weather.HS(S)-daily_weather.HS(Sprev)).*100;
new_snow(find(new_snow<0)) = 0;
b2 =bar(daily_weather.Time(S)-hours(7),new_snow,'k');
b2.FaceAlpha = 0.75;
ylim([0 80])
ylabel('Daily New Snow [cm]')
ax = gca;
ax.YColor = 'black';

p1 = plot([datetime(2021,2,6) datetime(2021,2,6)],[0 160],'r-','LineWidth',2);
plot([datetime(2021,1,30) datetime(2021,1,30)],[0 160],'r-','LineWidth',2)
plot([datetime(2021,1,8) datetime(2021,1,8)],[0 160],'r-','LineWidth',2)
legend([b3 b2 p1],'BBI','New Snow','Fatality')


set(hfig,'PaperUnits','centimeters');
set(hfig,'PaperSize',[a b]);
set(hfig,'PaperPosition', [0 0 a b]);
% print(hfig,'-dpng',[figure_path,'NewSnow_BBI','.png']);
% print(hfig,'-dpdf',[figure_path,'NewSnow_BBI','.pdf']);
% savefig(strcat(figure_path,'NewSnow_BBI'));
%% Avalache Size and type figure
close all
hfig=figure;
a = 30; b= 25; 
hfig.Units='centimeters';hfig.Position=[10 10 a b];
set(hfig,'defaulttextinterpreter','latex','DefaultAxesFontSize',ft_size)
haxes=tight_subplot(3,1,[0.02 0.02],[2 1]/a,[2 1.5]/b);
marker_size = 100;

%set axis 
st = stab_data.Time(S);


%avy timeline figure
axes(haxes(1))
for ii = 1:size(daily_sum.date,1)-1
   if isweekend(daily_sum.date(ii))
       b2 = area([daily_sum.date(ii)  daily_sum.date(ii+1)],[400 400],'FaceColor','k');
       b2.FaceAlpha=0.1;
       b2.EdgeAlpha = 0;
       hold on
       b2 = area([daily_sum.date(ii)  daily_sum.date(ii+1)],[-100 -100],'FaceColor','k');
       b2.FaceAlpha=0.1;
       b2.EdgeAlpha = 0;
       
   end
end

s1 = scatter(User_Trig_daily_sum.Time(find(User_Trig_daily_sum.Facets>0)),...
    User_Trig_daily_sum.Facets(find(User_Trig_daily_sum.Facets>0)),marker_size,'b','filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.5);
hold on
s2 = scatter(User_Trig_daily_sum.Time(find(User_Trig_daily_sum.NS_OS>0)),...
    User_Trig_daily_sum.NS_OS(find(User_Trig_daily_sum.NS_OS>0)),marker_size,'r','filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.5);
s3 = scatter(User_Trig_daily_sum.Time(find(User_Trig_daily_sum.DH>0)),...
    User_Trig_daily_sum.DH(find(User_Trig_daily_sum.DH>0)),marker_size,'k','filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.5);
s4 = scatter(User_Trig_daily_sum.Time(find(User_Trig_daily_sum.NS>0)),...
    User_Trig_daily_sum.NS(find(User_Trig_daily_sum.NS>0)),marker_size,'g','filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.5);

%naturals
s5 = scatter(Nat_Trig_daily_sum.Time(find(Nat_Trig_daily_sum.Facets>0)),...
    Nat_Trig_daily_sum.Facets(find(Nat_Trig_daily_sum.Facets>0)),marker_size,'b','filled','p','MarkerEdgeColor','k','MarkerFaceAlpha',0.5);
s6 = scatter(Nat_Trig_daily_sum.Time(find(Nat_Trig_daily_sum.NS_OS>0)),...
    Nat_Trig_daily_sum.NS_OS(find(Nat_Trig_daily_sum.NS_OS>0)),marker_size,'r','filled','p','MarkerEdgeColor','k','MarkerFaceAlpha',0.5);
s7 = scatter(Nat_Trig_daily_sum.Time(find(Nat_Trig_daily_sum.DH>0)),...
    Nat_Trig_daily_sum.DH(find(Nat_Trig_daily_sum.DH>0)),marker_size,'k','filled','p','MarkerEdgeColor','k','MarkerFaceAlpha',0.5);
s8 = scatter(Nat_Trig_daily_sum.Time(find(Nat_Trig_daily_sum.NS>0)),...
    Nat_Trig_daily_sum.NS(find(Nat_Trig_daily_sum.NS>0)),marker_size,'g','filled','p','MarkerEdgeColor','k','MarkerFaceAlpha',0.5);

%explosive work
% s9 = scatter(Exp_Trig_daily_sum.Time(find(Exp_Trig_daily_sum.Facets>0)),...
%     Exp_Trig_daily_sum.Facets(find(Exp_Trig_daily_sum.Facets>0)),marker_size,'b','filled','d','MarkerEdgeColor','k','MarkerFaceAlpha',0.5);
% s10 = scatter(Exp_Trig_daily_sum.Time(find(Exp_Trig_daily_sum.NS_OS>0)),...
%     Exp_Trig_daily_sum.NS_OS(find(Exp_Trig_daily_sum.NS_OS>0)),marker_size,'r','filled','d','MarkerEdgeColor','k','MarkerFaceAlpha',0.5);
% s11 = scatter(Exp_Trig_daily_sum.Time(find(Exp_Trig_daily_sum.DH>0)),...
%     Exp_Trig_daily_sum.DH(find(Exp_Trig_daily_sum.DH>0)),marker_size,'k','filled','d','MarkerEdgeColor','k','MarkerFaceAlpha',0.5);
% s12 = scatter(Exp_Trig_daily_sum.Time(find(Exp_Trig_daily_sum.NS>0)),...
%     Exp_Trig_daily_sum.NS(find(Exp_Trig_daily_sum.NS>0)),marker_size,'g','filled','d','MarkerEdgeColor','k','MarkerFaceAlpha',0.5);
%set(gca, 'YScale', 'log')

ylim([0 6])
grid on
ylabel('Number of Daily Avalaches [N]')

pfat = plot([datetime(2021,2,6) datetime(2021,2,6)],[0 160],'r-','LineWidth',2);
plot([datetime(2021,1,30) datetime(2021,1,30)],[0 160],'r-','LineWidth',2)
plot([datetime(2021,1,8) datetime(2021,1,8)],[0 160],'r-','LineWidth',2)

%dummy plot for symbols in legend
p_user_symbol = plot(datetime(2021,2,6), 200,'ko');
p_nat_symbol = plot(datetime(2021,2,6), 200,'kp');
%p_ext_symbol = plot(datetime(2021,2,6), 200,'kd');

p_facets = plot(datetime(2021,2,6), 200,'ks','MarkerFaceColor','b');
p_NSOS = plot(datetime(2021,2,6), 200,'ks','MarkerFaceColor','r');
p_DH_symbol = plot(datetime(2021,2,6), 200,'ks','MarkerFaceColor','k');
p_NS_symbol = plot(datetime(2021,2,6), 200,'ks','MarkerFaceColor','g');

DateString = {'12/19/2020';'02/20/2021'};
ticklocs = datetime(DateString); %R2014b or later
ticklocs.Format = 'dd-MMM-yyyy';
set(haxes([1]), 'XTick', ticklocs);    %R2016b or later
xlim([st(1) st(end)])

legend([p_user_symbol p_nat_symbol p_facets p_NSOS p_DH_symbol p_NS_symbol pfat],...
    'User Trig','Natural','Facets','New/Old Snow','Depth Hoar','New Snow','Fatality',...
'Location','northwest')

axes(haxes(3))
for ii = 1:size(daily_sum.date,1)-1
   if isweekend(daily_sum.date(ii))
       b2 = area([daily_sum.date(ii)  daily_sum.date(ii+1)],[400 400],'FaceColor','k');
       b2.FaceAlpha=0.1;
       b2.EdgeAlpha = 0;
       hold on
       b2 = area([daily_sum.date(ii)  daily_sum.date(ii+1)],[-100 -100],'FaceColor','k');
       b2.FaceAlpha=0.1;
       b2.EdgeAlpha = 0;
       
   end
end

stab_index_fail = find(stab_data.sk38<0.5);
stab_index_good = find(stab_data.sk38>=0.5);
scatter(stab_data.Time(stab_index_fail)-hours(7),stab_data.sk38(stab_index_fail),50,'ro','Filled')
hold on
scatter(stab_data.Time(stab_index_good)-hours(7),stab_data.sk38(stab_index_good),50,'ko','Filled')
ylabel('SK38 index [ ]')
grid on
xlim([st(1) st(end)])
ylim([0 2])
pfat = plot([datetime(2021,2,6) datetime(2021,2,6)],[0 160],'r-','LineWidth',2);
plot([datetime(2021,1,30) datetime(2021,1,30)],[0 160],'r-','LineWidth',2)
plot([datetime(2021,1,8) datetime(2021,1,8)],[0 160],'r-','LineWidth',2)
%yyaxis right
%stab_data.sk38(find(stab_data.sk38>2)) = nan;
%plot(stab_data.Time(S),stab_data.sk38(S),'b-')



axes(haxes(2))
for ii = 1:size(daily_sum.date,1)-1
   if isweekend(daily_sum.date(ii))
       b2 = area([daily_sum.date(ii)  daily_sum.date(ii+1)],[400 400],'FaceColor','k');
       b2.FaceAlpha=0.1;
       b2.EdgeAlpha = 0;
       hold on
       b2 = area([daily_sum.date(ii)  daily_sum.date(ii+1)],[-100 -100],'FaceColor','k');
       b2.FaceAlpha=0.1;
       b2.EdgeAlpha = 0;
       
   end
end

b3 = bar(daily_weather.Time(S)-hours(7),daily_weather.ISWR(S)./max(daily_weather.ISWR(S)),'b','LineWidth',1);
b3.FaceAlpha = 0.3;

ylabel('Blue Bird Index [ ]')
%xlabel('[MM/DD/YYYY]')
%title(strcat(plot_title,' Daily Avg'))
grid on
set(gca,'XMinorTick','on','YMinorTick','on')
axis tight
hold on
ylim([0 1])
xlim([st(1) st(end)])

yyaxis right
%calc new snow in cm
new_snow = (daily_weather.HS(S)-daily_weather.HS(Sprev)).*100;
new_snow(find(new_snow<0)) = 0;
b2 =bar(daily_weather.Time(S)-hours(7),new_snow,'k');
b2.FaceAlpha = 0.75;
ylim([0 80])
ylabel('Daily New Snow [cm]')
ax = gca;
ax.YColor = 'black';

p1 = plot([datetime(2021,2,6) datetime(2021,2,6)],[0 160],'r-','LineWidth',2);
plot([datetime(2021,1,30) datetime(2021,1,30)],[0 160],'r-','LineWidth',2)
plot([datetime(2021,1,8) datetime(2021,1,8)],[0 160],'r-','LineWidth',2)
legend([b3 b2 p1],'BBI','New Snow','Fatality')
xlim([st(1) st(end)])

set(haxes([1 2]),'XtickLabel',[])

set(hfig,'PaperUnits','centimeters');
set(hfig,'PaperSize',[a b]);
set(hfig,'PaperPosition', [0 0 a b]);
print(hfig,'-dpng',[figure_path,'AvyTimeline_weather','.png']);
print(hfig,'-dpdf',[figure_path,'AvyTimeline_weather','.pdf']);
savefig(strcat(figure_path,'AvyTimeline_weather'));
%%





