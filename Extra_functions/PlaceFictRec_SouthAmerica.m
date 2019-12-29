   
%find the pixes which correspond to South America

function  [Receiver_loc,Mask,Lon_gr,Lat_gr] = PlaceFictRec_SouthAmerica(Receiver_number)

%insert the number of receivers
%Find the area of South America (mainly the land) and place hypothetical receivers


load coast.mat


LatLim_gr = [-40.5937   11.7574];
LonLim_gr = [-87.2475  -27.7632];

isouth = find(long>=-89 & long<=-30 &lat>=-35 &lat<=9);
lon_south=long(isouth);
lat_south=lat(isouth);



%Place receivers uniformly in the continent
LatRange = linspace(LatLim_gr(1),LatLim_gr(2),round(sqrt(Receiver_number)));
LonRange = linspace(LonLim_gr (1),LonLim_gr (2),round(sqrt(Receiver_number)));
[Lon_gr,Lat_gr]= meshgrid(LonRange,LatRange);

Lon_gr_vec = Lon_gr(:);
Lat_gr_vec = Lat_gr(:);

Mask = zeros(length(Lon_gr_vec),1);


Receiver_loc = [];
for i=1:length(Lon_gr(:))
   
    Slat= abs(lat_south-Lat_gr_vec(i));
    
    [Slat,ii]= sort(Slat);
    lat_south = lat_south(ii);
    lon_south = lon_south(ii);
   
   
    jj = find(Slat<1.5);
   
    if ~isempty(jj)
        ii_in = jj(1);
        ii_end= jj(2);


        for k=2:length(jj)
            Slon = abs(lon_south(ii_in)-lon_south(jj(k))); 
            if Slon>=7 
                ii_end = (jj(k));
                break;
            end
        end 


        Long_test = sort([lon_south(ii_in); lon_south(ii_end)]);
       % plot( Long_test,lat_south([ii_in,ii_end]),'ok')
        if (Lon_gr_vec(i)>=Long_test(1)-2 &&  Lon_gr_vec(i)<=Long_test(2)+10 && Lat_gr_vec(i)<=0 && Lat_gr_vec(i)>=-3)
           
            Mask(i)= 1;
           Receiver_loc=[Receiver_loc;Lon_gr_vec(i),Lat_gr_vec(i);];
        %  plot(Lon_gr_vec(i),Lat_gr_vec(i),'^','MarkerFaceColor',[.1 0.29 .03],'MarkerSize',3,'MarkerEdgeColor','none')
        elseif (Lon_gr_vec(i)>=Long_test(1)-2 &&  Lon_gr_vec(i)<=Long_test(2)+5 && Lat_gr_vec(i)<5)
      
            Mask(i)=1;
         % plot(Lon_gr_vec(i),Lat_gr_vec(i),'^','MarkerFaceColor',[.1 0.29 .03],'MarkerSize',3,'MarkerEdgeColor','none')
          Receiver_loc=[Receiver_loc;Lon_gr_vec(i),Lat_gr_vec(i);];
        elseif (Lon_gr_vec(i)>=Long_test(1)-1 &&  Lon_gr_vec(i)<=Long_test(2)+2)
             
            Mask(i)=1;
          %  plot(Lon_gr_vec(i),Lat_gr_vec(i),'^','MarkerFaceColor',[.1 0.29 .03],'MarkerSize',3,'MarkerEdgeColor','none')
             Receiver_loc=[Receiver_loc;Lon_gr_vec(i),Lat_gr_vec(i);];
        
         end
    end
    
end
% axis square
% title('Locations of Hypothetical Receivers')
