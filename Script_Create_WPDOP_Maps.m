%In this script we generate an ionspheric risk map
%of time defined by epoch_ini and epoch_end and then we estimate a weighted
% PDOP on the ground at time T
%this code was created by A. Koulouri, 15.11.2018
%Further details are given in the publication: Methodology to estimate ionospheric
%scintillation risk maps and their contribution to position dilution of precision on the ground 
%https://arxiv.org/abs/1911.08229?context=physics

close all;
clear all;
addpath(genpath(cd));



%load data
%Time interval of the collected data 
data_date_begin = '2014-11-01 00:00:00';
data_date_end   = '2014-11-04 24:00:00';
Split_date_begin = regexp(regexp(data_date_begin,' ','Split'),':','Split');
Split_date_end = regexp(regexp(data_date_end,' ','Split'),':','Split');

load_name = ['data_ar_Begin_',Split_date_begin{1}{1},'_',Split_date_begin{2}{1},'_',Split_date_begin{2}{2},'_',Split_date_begin{2}{3},'_End_',Split_date_end{1}{1},'_',Split_date_end{2}{1},'_',Split_date_end{2}{2},'_',Split_date_end{2}{3},'.mat'];
matfile = fullfile('data/', load_name);
load(matfile);


%Input data=[Svid_id,IPPs_lon,IPPs_lat,S4,S4t,St_id]
%IPPs(Ionospheric piece points at 350km) in (Lon,Lat)
%Svid_id (PRN of the satellite), St_id(code name of the ground receiver) (this information defines the satellite-ground station link) 
%S4 values for each IPP location
%S4t elapsed time in minutes (S4t(0)=0)
%The original downloaded data in given in data_orig_Begin_2014-11-01_00_00_00_End_2014-11-04_24_00_00

[Svid_id,IPPLon,IPPLat,S4,S4t,StationId]= deal(data_ar(:,1),data_ar(:,2),data_ar(:,3),data_ar(:,4),data_ar(:,5),data_ar(:,6));
load RefTime; %the actual time that the first data point was collected

%% Specify the interval to produce the static (average) ionospheric risk map
epoch_ini   = {'2014-11-02 23:00:00'};
epoch_end   = {'2014-11-03 03:00:00'};
           
% define the square region of the map
Res = 2;      % spatial resolution of the map in degrees
eps1 = 0.4;   % extensions of the map
LatLim_ion = [min(IPPLat)-eps1 max(IPPLat)+eps1];
LonLim_ion = [min(IPPLon)-eps1 max(IPPLon)+eps1];

% parameters for the ionospheric map
Dth = [1 5 10];    % scintillation duration threshold in min
Sth = [0.3 0.5 0.7];%+eps;  % scintillation intensity threshold 


%estimate the probability maps ->Ionospheric Risk Map!
t1= datenum(epoch_ini);
t2 = datenum(epoch_end);
date_begin = (t1- ReferenceTime)*1440; %estimate time elapsed from the starting point of data collection and convert into minutes
date_end   = (t2- ReferenceTime)*1440;

ind_ep_ini = find(abs(S4t-date_begin) == min(abs(S4t-date_begin)));%find(S4t <= date_begin & S4t >= 0.9*date_begin);
ind_ep_end = find(abs(S4t-date_end) == min(abs(S4t-date_end)));    %find(S4t >=date_end & S4t<=1.1*date_end);
Interval = [ind_ep_ini(1):ind_ep_end(end)];

data_for_risk_map = data_ar(Interval,:);%[IPPLon(Interval) IPPLat(Interval) Svid_ar(Interval) StationId(Interval) S4(Interval) S4t(Interval)];

%Estimate the ionospheric map given thresholds Sth(1) and Dth(1)
%This function returns the risk map and the Lon/Lat grid
IonosphericMap = IonosphericRiskMap(Res,data_for_risk_map,Sth(2),Dth(1),LatLim_ion,LonLim_ion);
           
        
%% Estimate the weighted PDOP MAPS

%define the time instant T to estimate the WPDOP ground map 
T = 0.650*(t2-t1)+t1; %you can select any minute between [t1-t2]  time interval

%%estimate the locations of the satellites at time T
%define the constellation (GPS or GNSS)
SatXYZ  = SatelliteLocations(T, ReferenceTime,data_for_risk_map,'GPS'); %CHECK THIS!

%Estimate  PDOP maps at time T
n = 30;
Gridsize = [n n]; %size of the ground map in degrees
[Receiver_loc,MaskGridLoc,Lon_grid_gr,Lat_grid_gr] = PlaceFictRec_SouthAmerica(n^2);%locations of the hypothetical receivers
[WPDOP,PDOP] = WPDOP_map(Receiver_loc,MaskGridLoc,SatXYZ,IonosphericMap);

%% Plot the ionspheric and ground maps
figure
set(gcf, 'Units','centimeters', 'Position',[5 5 30 10]);
subplot(1,4,1)%Ionospheric map
load coast_gr.mat
h = pcolor(IonosphericMap.grid.lon,IonosphericMap.grid.lat ,IonosphericMap.risk);
set(h, 'EdgeColor', [102 99 99]./255);
hold on;
plot(lon_cost,lat_cost,'color',[0.7 0.7 0.7],'linewidth',1.1)
% MagEqu = shaperead('figures/equador-brasil.shp');
% plot([MagEqu(1:end).X],[MagEqu(1:end).Y],'color',[1	80	32]./255,'marker','.')
axis square
axis xy
title('Ionospheric Scintillation Risks','FontSize',8)
cmap = hot(40);
cmap = cmap(end-4:-1:1,:);
cmap(2:4,:)=[];
colormap(cmap)
c1=colorbar('southoutside');
%c1.Label.String = '0<=Risk Values<=1';
x1=get(gca,'position');
x=get(c1,'Position');
x(4)=x(4)/2;
set(c1,'Position',x)
set(gca,'position',x1)



subplot(1,4,2)%ground map
load coast.mat
h=pcolor(Lon_grid_gr,Lat_grid_gr,WPDOP);
set(h, 'EdgeColor', [10 99 99]./255);
hold on;
plot(long,lat,'color',[0.1 0.1 0.1],'linewidth',1.8)
caxis([0, 10])

title('WPDOP')
axis square
axis xy

c1=colorbar('southoutside');
%c1.Label.String = '0<=Risk Values<=1';
x1=get(gca,'position');
x=get(c1,'Position');
x(4)=x(4)/2;
set(c1,'Position',x)
set(gca,'position',x1)



subplot(1,4,3)%ground map
load coast.mat
h=pcolor(Lon_grid_gr,Lat_grid_gr,PDOP);
set(h, 'EdgeColor', [10 99 99]./255);
hold on;
plot(long,lat,'color',[0.1 0.1 0.1],'linewidth',1.8)
title('PDOP') %white pixels designate no information
caxis([0, 10])
axis square
axis xy

c1=colorbar('southoutside');
%c1.Label.String = '0<=Risk Values<=1';
x1=get(gca,'position');
x=get(c1,'Position');
x(4)=x(4)/2;
set(c1,'Position',x)
set(gca,'position',x1)



subplot(1,4,4)%ground map
load coast.mat
h=pcolor(Lon_grid_gr,Lat_grid_gr,(WPDOP-PDOP)./WPDOP*100);
set(h, 'EdgeColor',[10 99 99]./255);
hold on;
plot(long,lat,'color',[0.1 0.1 0.1],'linewidth',1.8)
title('(WPDOP-PDOP)/WPDOP %') %white pixels designate no information
caxis([0, 100])
axis square
axis xy
c1=colorbar('southoutside');
%c2.Label.String = '%'
x1=get(gca,'position');
x=get(c1,'Position');
x(4)=x(4)/2;
set(c1,'Position',x)
set(gca,'position',x1)
