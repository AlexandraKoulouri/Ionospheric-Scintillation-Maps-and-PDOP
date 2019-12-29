function SatXYZ = SatelliteLocations(TimeInstant, ReferenceTime,data,SAToption)

%Use the availabe IPP data and ground stations to estimate the locations of
%the satellites
%[IPP_lon IPP_lat Svid_GPS StationId_GPS] = deal(data(:,1),data(:,2),data(:,3),data(:,4));


[Svid,S4t]= deal(data(:,1),data(:,5));


epoch_ini_gr = (TimeInstant - ReferenceTime)*1440; %datenum(IPPs{1,5}) is the time that the first sample was estimated
epoch_end_gr   =  epoch_ini_gr+0.58;%30 seconds  


if (SAToption =='GPS')
    
    Ind_sat =find(Svid>700000 );
    
else 
    
    Ind_sat = find(Svid>700000 | Svid<6900); %GPS and GLONASS
end

S4t = S4t(Ind_sat);
ind_epoch_ini_gr = find(abs(S4t-epoch_ini_gr) == min(abs(S4t-epoch_ini_gr)));
ind_epoch_end_gr = find(abs(S4t-epoch_end_gr) == min(abs(S4t-epoch_end_gr)));
DataInterval = [ind_epoch_ini_gr(1):ind_epoch_end_gr(end)];
data_t = data(Ind_sat(DataInterval),:);
clear Svid;

[Svid,IPPLon,IPPLat,StationId]= deal(data_t(:,1),data_t(:,2),data_t(:,3),data_t(:,6));


Links = [Svid StationId];
%SatIds = unique(Svid);
SatIds = unique(Links(:,1));
%StIds = unique(Links(:,2));
SatXYZ = zeros(length(SatIds),3);



%load ground stations from South America
%Ground stations (coordinates and code)
load ('data/StationMAT.mat')


%find the satellite positions at this time interval (t1_gr to t2_gr)
for SatInd =  1:length(SatIds)

      AllStInd =  find(Links(:,1) == SatIds(SatInd)); %for this satellite find all the available ground stations

      PA = zeros(length(AllStInd),3);
      PB = zeros(length(AllStInd),3);
      for StIter = 1:length(AllStInd)
          StInd = find(Station_id == Links(AllStInd(StIter),2)); %find the longitude and latitude for each of the stations
          LatGround = Station_Lat(StInd);
          LonGround = Station_Lon(StInd); 
          [xgr,ygr,zgr] = Geodetic2Geocentric(LonGround, LatGround,0);
          PA(StIter,:) = [xgr ygr zgr];
          [xipp,yipp,zipp] = Geodetic2Geocentric(IPPLon(AllStInd(StIter)), IPPLat(AllStInd(StIter)),350);
          PB(StIter,:) = [xipp yipp zipp];


      end
      % %% Estimate satellite locations (not accurate approach but considering
      % that their trajectory is at 220.000km using the IPP information over
      % this period of time and the corresponding stations


      SatXYZ(SatInd,:) = lineIntersect3D(PA,PB); %Cartesian Coordinates of the satellites

end