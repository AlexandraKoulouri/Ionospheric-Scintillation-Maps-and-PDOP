function  [x,y,z,lat,lon] = Project_Earth2IPP(Sat,Gr,R,m)
%inputs
%Sat = [x2,y2,z2] ->satellite Cartesian coordinates
%Gr = [x1,,y1,z1]->ground receiver Cartesian coordinates
%sphere (center of sphere m, radius h) ->ionospheric surface considered
%perfect sphere
 
%output 
% IPP [x,y,z] in Cartesian coordinates
% [lat,lon] also geodetic coordinates

%this code was created by A. Koulouri, 25.5.2018

%intersection line segment L_s-L_g with a sphere [(0,0,0),h]
%t takes values between (0,1)
x1 = Gr(:,1); y1 = Gr(:,2); z1 = Gr(:,3);
x2 = Sat(:,1); y2 = Sat(:,2); z2 = Sat(:,3);

%line representation
%x = x1 + (x2 - x1)*t = x1 + i*t
%y = y1 + (y2 - y1)*t = y1 + j*t
%z = z1 + (z2 - z1)*t = z1 + k*t

%and a sphere defined as:
%(x - m(1))^2 + (y - m(2))^2 + (z - m(3))^2 = R^2

%Substituting in gives the quadratic equation:
%a*t^2 + b*t + c = 0
%where:

x =NaN; y=x; z=x; lat=x;lon=x;

a = (x2 - x1).^2 + (y2 - y1).^2 + (z2 - z1).^2;
b = 2*(x2 - x1).*(x1 - m(1)) + 2.*(y2 - y1).*(y1 - m(2)) + 2.*(z2 - z1).*(z1 - m(3));
c = (x1-m(1)).^2 + (y1-m(2)).^2 + (z1-m(3)).^2 - R^2;

D = b.^2-4.*a.*c;
if (D>=0) %this condition is valid when there is intersection
    t= [-b+sqrt(D) -b-sqrt(D)]./(2*a);
    
    i = find(t>=0 & t<=1); %from the two option keep t in [0,1];

    if (~isempty(i))
        x = x1 + (x2 - x1).*t(i);
        y = y1 + (y2 - y1).*t(i);
        z = z1 + (z2 - z1).*t(i);
        % define ionhosheric spheroid (consider perfect sphere)
        s = oblateSpheroid; %
        s.SemimajorAxis = R;
        s.SemiminorAxis = R;
       wgs84 = wgs84Ellipsoid('kilometers');
       [lat,lon,~] = ecef2geodetic(wgs84,x,y,z);
      %  [lat,lon,h] = ecef2geodetic(s,x,y,z);
        
    end


end