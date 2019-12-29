%In this script we generate ionspheric risk maps 
%following the publication: Methodology to estimate ionospheric
%scintillation risk maps and their contribution to position dilution of precision on the ground 
%https://arxiv.org/abs/1911.08229?context=physics
%this code was created by A. Koulouri, 01.09.2018
close all;
clear all;
addpath(genpath(cd));

%Ground stations (coordinates and codes)
load ('data/StationMAT.mat')


%load S4 data
%Time interval of the collected data 
data_date_begin = '2014-11-01 00:00:00';
data_date_end   = '2014-11-04 24:00:00';
Split_date_begin = regexp(regexp(data_date_begin,' ','Split'),':','Split');
Split_date_end = regexp(regexp(data_date_end,' ','Split'),':','Split');

load_name = ['data_ar_Begin_',Split_date_begin{1}{1},'_',Split_date_begin{2}{1},'_',Split_date_begin{2}{2},'_',Split_date_begin{2}{3},'_End_',Split_date_end{1}{1},'_',Split_date_end{2}{1},'_',Split_date_end{2}{2},'_',Split_date_end{2}{3},'.mat'];
matfile = fullfile('data/', load_name);
load(matfile);

[Svid_id,IPPLon,IPPLat,S4,S4t,StationId]= deal(data_ar(:,1),data_ar(:,2),data_ar(:,3),data_ar(:,4),data_ar(:,5),data_ar(:,6));

% Svid_ar = data_ar(:,1);
% IPPLat = data_ar(:,3);
% IPPLon = data_ar(:,2);
% S4 = data_ar(:,4);
% S4t = data_ar(:,5); % time interval Dt of the collected data S4, t(1) is 0 and t(end) show the length of the interval (time in minutes)
% StationId = data_ar(:,6);

load RefTime; %the time that the first data point was collected

%%  define the square region of the ionospheric map
Res = 2;      % spatial resolution of the map in degrees
eps1 = 0.4;   % extensions of the map
LatLim_ion = [min(IPPLat)-eps1 max(IPPLat)+eps1]; %limits of the map
LonLim_ion = [min(IPPLon)-eps1 max(IPPLon)+eps1];


%Specify a particular time interval to create static ionospheric maps
epoch_ini   = {'2014-11-01 23:00:00'};
epoch_end   = {'2014-11-02 03:00:00'};

% parameters
Dth = [1 5 10];    % scintillation duration threshold in min
Sth = [0.3 0.5 0.7];%+eps;  % scintillation intensity threshold 



%%estimate the probability maps ->Ionospheric Risk Map for different thresholds!
for i = 1:length(epoch_ini)
    t1= datenum(epoch_ini(i));
    t2 = datenum(epoch_end(i));
    date_begin = (t1- ReferenceTime)*1440; %estimate time elapsed from the starting point of data collection and convert into minutes
    date_end   = (t2- ReferenceTime)*1440;

    ind_ep_ini = find(abs(S4t-date_begin) == min(abs(S4t-date_begin)));%find(S4t <= date_begin & S4t >= 0.9*date_begin);
    ind_ep_end = find(abs(S4t-date_end) == min(abs(S4t-date_end)));    %find(S4t >=date_end & S4t<=1.1*date_end);
    Interval = [ind_ep_ini(1):ind_ep_end(end)];
    
    data_for_risk_map = data_ar(Interval,:);%[IPPLon(Interval) IPPLat(Interval) Svid_ar(Interval) StationId(Interval) S4(Interval) S4t(Interval)];
    for iterDur = 1:length(Dth)
        for iterTh = 1:length(Sth)
            
            
            IonosphericMap = IonosphericRiskMap(Res,data_for_risk_map,Sth(iterTh),Dth(iterDur),LatLim_ion,LonLim_ion);
                       

            figure
            load coast_gr.mat
            h = pcolor(IonosphericMap.grid.lon,IonosphericMap.grid.lat ,IonosphericMap.risk);
            set(h, 'EdgeColor', [102 99 99]./255);
            hold on;
            plot(lon_cost,lat_cost,'color',[0.7 0.7 0.7],'linewidth',1.1)
            axis square
            axis xy
            cmap = hot(40);
            cmap = cmap(end-4:-1:1,:);
            cmap(2:4,:)=[];
            colormap(cmap);
            colorbar

        
        end
        
        
        
    end

    
end
